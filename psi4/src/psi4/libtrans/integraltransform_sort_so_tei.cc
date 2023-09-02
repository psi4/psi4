/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "integraltransform_functors.h"
#include "mospace.h"
#include "integraltransform.h"

#include "psi4/libpsio/psio.hpp"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libmints/matrix.h"
#include "psi4/psifiles.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <cmath>
#include <numeric>

using namespace psi;

/**
 * @brief Computes Fock matrices, frozen core operators and other Fock-like quantities.  This shouldn't
 *        be needed because those quantities are computed during the SO integral presort.  However, in
 *        Brueckner CC and other orbital optimized methods, the Fock matrices need an update so this
 *        function may be called to extract them from the DPD buffers, allowing the IWL buffers to be
 *        nuked after the first integral transformation..
 * @param Hcore - the SO basis core Hamiltonian contribution to the Fock matrix.
 * @param Cmats - a list of the subset of C matrices for each Fock-like quantity, e.g. the nso*nfzc
 *                subset for the frozen core operator.
 * @return a list of the Fock matrices associated with each incoming C matrix.
 */
std::vector<SharedMatrix> IntegralTransform::compute_fock_like_matrices(SharedMatrix Hcore,
                                                                        std::vector<SharedMatrix> Cmats) {
    // This function is supposed to only be called after an initial presort, but we'll check to make sure.
    if (!alreadyPresorted_) presort_so_tei();

    // Some
    int nBuckets, thisBucketRows;
    size_t rowsPerBucket, rowsLeft, memFree;

    int nmats = Cmats.size();

    // Form the Density matrices associated with the C matrices, and allocate the F matrices.
    std::vector<SharedMatrix> Dmats, Fmats;
    for (size_t N = 0; N < nmats; ++N) {
        SharedMatrix Cmat = Cmats[N];
        if (Cmat->rowspi() != sopi_)
            throw PSIEXCEPTION(
                "Row dimension of C matrix is not equal to SOs per irrep in LibTrans::compute_fock_like_matrices()");
        auto Fmat = std::make_shared<Matrix>("F matrix", sopi_, sopi_);
        Fmats.push_back(Fmat);
        auto Dmat = std::make_shared<Matrix>("D matrix", sopi_, sopi_);
        Dmat->gemm(false, true, 1.0, Cmat, Cmat, 0.0);
        Dmats.push_back(Dmat);
    }

    psio_->open(PSIF_SO_PRESORT, PSIO_OPEN_OLD);

    // Grab control of DPD for now, but store the active number to restore it later
    int currentActiveDPD = psi::dpd_default;
    dpd_set_default(myDPDNum_);

    dpdbuf4 J;
    global_dpd_->buf4_init(&J, PSIF_SO_PRESORT, 0, DPD_ID("[n,n]"), DPD_ID("[n,n]"), DPD_ID("[n>=n]+"),
                           DPD_ID("[n>=n]+"), 0, "SO Ints (nn|nn)");

    for (int h = 0; h < nirreps_; h++) {
        if (J.params->coltot[h] && J.params->rowtot[h]) {
            memFree = static_cast<size_t>(dpd_memfree() - J.params->coltot[h]);
            rowsPerBucket = memFree / (J.params->coltot[h]);
            if (rowsPerBucket > J.params->rowtot[h]) rowsPerBucket = (size_t)J.params->rowtot[h];
            nBuckets = static_cast<int>(
                std::ceil(static_cast<double>(J.params->rowtot[h]) / static_cast<double>(rowsPerBucket)));
            rowsLeft = static_cast<size_t>(J.params->rowtot[h] % rowsPerBucket);
        } else {
            nBuckets = 0;
            rowsPerBucket = 0;
            rowsLeft = 0;
        }

        if (print_ > 1) {
            outfile->Printf("\th = %d; memfree         = %lu\n", h, memFree);
            outfile->Printf("\th = %d; rows_per_bucket = %lu\n", h, rowsPerBucket);
            outfile->Printf("\th = %d; rows_left       = %lu\n", h, rowsLeft);
            outfile->Printf("\th = %d; nbuckets        = %d\n", h, nBuckets);
        }

        // We assume that we're always building totally symmetric Fock matrices here.
        int sym = 0;

        for (int n = 0; n < nBuckets; n++) {
            if (nBuckets == 1)
                thisBucketRows = rowsPerBucket;
            else
                thisBucketRows = (n < nBuckets - 1) ? rowsPerBucket : rowsLeft;
            global_dpd_->buf4_mat_irrep_rd_block(&J, h, n * rowsPerBucket, thisBucketRows);
            for (int pq = 0; pq < thisBucketRows; pq++) {
                int pabs = J.params->roworb[h][pq][0];
                int qabs = J.params->roworb[h][pq][1];
                int psym = J.params->psym[pabs];
                int qsym = J.params->qsym[qabs];
                int prel = pabs - J.params->poff[psym];
                int qrel = qabs - J.params->qoff[qsym];
                for (int rs = 0; rs < J.params->coltot[h]; ++rs) {
                    int rabs = J.params->colorb[h][rs][0];
                    int sabs = J.params->colorb[h][rs][1];
                    int rsym = J.params->rsym[rabs];
                    int ssym = J.params->ssym[sabs];
                    int rrel = rabs - J.params->roff[rsym];
                    int srel = sabs - J.params->soff[ssym];
                    int pqsym = psym ^ qsym;
                    int rssym = rsym ^ ssym;
                    int qrsym = qsym ^ rsym;
                    int pssym = psym ^ ssym;
                    int prsym = psym ^ rsym;
                    int qssym = qsym ^ ssym;
                    double value = J.matrix[h][pq][rs];
                    for (int N = 0; N < nmats; ++N) {
                        SharedMatrix D = Dmats[N];
                        SharedMatrix F = Fmats[N];
                        if (pqsym == rssym && pqsym == sym) {
                            F->add(rsym, rrel, srel, D->get(psym, prel, qrel) * value);
                        }

                        if (qrsym == pssym && qrsym == sym) {
                            F->add(qsym, qrel, rrel, -0.5 * D->get(psym, prel, srel) * value);
                        }
                    } /* N */
                }     /* rs */
            }         /* pq */
        }
        global_dpd_->buf4_mat_irrep_close_block(&J, h, rowsPerBucket);
    } /* h */

    for (int N = 0; N < nmats; ++N) Fmats[N]->add(Hcore);

    global_dpd_->buf4_close(&J);

    psio_->close(PSIF_SO_PRESORT, keepDpdSoInts_);

    // Hand DPD control back to the user
    dpd_set_default(currentActiveDPD);

    return Fmats;
}

/**
 * Presort the two-electron integrals into DPD buffers to prepare them for
 * the transformation.  The frozen core operator is built simultaneously.
 * If this action has already been performed, it will just load the frozen
 * core operator from disk and return.
 */
void IntegralTransform::presort_so_tei() {
    check_initialized();

    if (alreadyPresorted_) {
        if (print_ > 5) outfile->Printf("\tSO integrals are already sorted, moving on...\n");
        return;
    }

    /// Initialize frozen core density
    auto aFzcD = std::vector<double>(nTriSo_); // Density matrix from core orbitals only
    std::vector<double> bFzcD, bFzcOp;
    if (transformationType_ != TransformationType::Restricted) {
        bFzcD = std::vector<double>(nTriSo_);
    }

    // Set frozen core density
    for (int h = 0, soOffset = 0; h < nirreps_; ++h) {
        auto pCa = Ca_->pointer(h);
        auto pCb = Cb_->pointer(h);
        for (int p = 0; p < sopi_[h]; ++p) {
            for (int q = 0; q <= p; ++q) {
                int pq = INDEX((p + soOffset), (q + soOffset));
                for (int i = 0; i < frzcpi_[h]; ++i) aFzcD[pq] += pCa[p][i] * pCa[q][i];
                if (transformationType_ != TransformationType::Restricted) {
                    for (int i = 0; i < frzcpi_[h]; ++i) bFzcD[pq] += pCb[p][i] * pCb[q][i];
                }
            }
        }
        soOffset += sopi_[h];
    }

    // Copy the OEI into the frozen core operator
    auto aoH = H_->to_lower_triangle();
    std::vector<double> aFzcOp(aoH, aoH + nTriSo_); // The frozen core operator. Not complete yet.
    if (transformationType_ != TransformationType::Restricted) {
        bFzcOp = std::vector<double>(aoH, aoH + nTriSo_);
    }

    /// Resume TEI sorting
    int currentActiveDPD = psi::dpd_default;
    dpd_set_default(myDPDNum_);

    if (print_) {
        outfile->Printf("\tPresorting SO-basis two-electron integrals.\n");
    }

    dpdfile4 I;
    psio_->open(PSIF_SO_PRESORT, PSIO_OPEN_NEW);
    global_dpd_->file4_init(&I, PSIF_SO_PRESORT, 0, DPD_ID("[n>=n]+"), DPD_ID("[n>=n]+"), "SO Ints (nn|nn)");

    size_t memoryd = memory_ / sizeof(double);

    int nump = 0, numq = 0;
    for (int h = 0; h < nirreps_; ++h) {
        nump += I.params->ppi[h];
        numq += I.params->qpi[h];
    }
    int **bucketMap = init_int_matrix(nump, numq);

    /* Room for one bucket to begin with */
    int **bucketOffset = (int **)malloc(sizeof(int *));
    bucketOffset[0] = init_int_array(nirreps_);
    int **bucketRowDim = (int **)malloc(sizeof(int *));
    bucketRowDim[0] = init_int_array(nirreps_);
    long int **bucketSize = (long int **)malloc(sizeof(long int *));
    bucketSize[0] = init_long_int_array(nirreps_);

    /* Figure out how many passes we need and where each p,q goes */
    int nBuckets = 1;
    size_t coreLeft = memoryd;
    psio_address next;
    for (int h = 0; h < nirreps_; ++h) {
        size_t rowLength = (size_t)I.params->coltot[h ^ (I.my_irrep)];
        for (int row = 0; row < I.params->rowtot[h]; ++row) {
            if (coreLeft >= rowLength) {
                coreLeft -= rowLength;
                bucketRowDim[nBuckets - 1][h]++;
                bucketSize[nBuckets - 1][h] += rowLength;
            } else {
                nBuckets++;
                coreLeft = memoryd - rowLength;
                /* Make room for another bucket */
                int **p;

                p = static_cast<int **>(realloc(static_cast<void *>(bucketOffset), nBuckets * sizeof(int *)));
                if (p == nullptr) {
                    throw PsiException("file_build: allocation error", __FILE__, __LINE__);
                } else {
                    bucketOffset = p;
                }
                bucketOffset[nBuckets - 1] = init_int_array(nirreps_);
                bucketOffset[nBuckets - 1][h] = row;

                p = static_cast<int **>(realloc(static_cast<void *>(bucketRowDim), nBuckets * sizeof(int *)));
                if (p == nullptr) {
                    throw PsiException("file_build: allocation error", __FILE__, __LINE__);
                } else {
                    bucketRowDim = p;
                }
                bucketRowDim[nBuckets - 1] = init_int_array(nirreps_);
                bucketRowDim[nBuckets - 1][h] = 1;

                long int **pp;
                pp = static_cast<long int **>(realloc(static_cast<void *>(bucketSize), nBuckets * sizeof(long int *)));
                if (pp == nullptr) {
                    throw PsiException("file_build: allocation error", __FILE__, __LINE__);
                } else {
                    bucketSize = pp;
                }
                bucketSize[nBuckets - 1] = init_long_int_array(nirreps_);
                bucketSize[nBuckets - 1][h] = rowLength;
            }
            int p = I.params->roworb[h][row][0];
            int q = I.params->roworb[h][row][1];
            bucketMap[p][q] = nBuckets - 1;
        }
    }

    if (print_) {
        outfile->Printf("\tSorting File: %s nbuckets = %d\n", I.label, nBuckets);
    }

    next = PSIO_ZERO;
    for (int n = 0; n < nBuckets; ++n) { /* nbuckets = number of passes */
        /* Prepare target matrix */
        for (int h = 0; h < nirreps_; h++) {
            I.matrix[h] = block_matrix(bucketRowDim[n][h], I.params->coltot[h]);
        }

        DPDFillerFunctor dpdfiller(&I, n, bucketMap, bucketOffset, false, true);
        NullFunctor null;
        IWL *iwl = new IWL(psio_.get(), soIntTEIFile_, tolerance_, 1, 1);
        // We need to feed the IWL integrals to construct the frozen core operator only once
        // If we're not on the first DPD bucket, skip it for efficiency.
        if (transformationType_ == TransformationType::Restricted) {
            FrozenCoreRestrictedFunctor frozencore(aFzcD.data(), aFzcOp.data());
            if (n)
                iwl_integrals(iwl, dpdfiller, null);
            else
                iwl_integrals(iwl, dpdfiller, frozencore);
        } else {
            FrozenCoreUnrestrictedFunctor frozencore(aFzcD.data(), bFzcD.data(), aFzcOp.data(), bFzcOp.data());
            if (n)
                iwl_integrals(iwl, dpdfiller, null);
            else
                iwl_integrals(iwl, dpdfiller, frozencore);
        }
        delete iwl;

        for (int h = 0; h < nirreps_; ++h) {
            if (bucketSize[n][h])
                psio_->write(I.filenum, I.label, (char *)I.matrix[h][0], bucketSize[n][h] * ((long int)sizeof(double)),
                             next, &next);
            free_block(I.matrix[h]);
        }
    } /* end loop over buckets/passes */

    /* Get rid of the input integral file */
    psio_->open(soIntTEIFile_, PSIO_OPEN_OLD);
    psio_->close(soIntTEIFile_, keepIwlSoInts_);

    free_int_matrix(bucketMap);

    for (int n = 0; n < nBuckets; ++n) {
        free(bucketOffset[n]);
        free(bucketRowDim[n]);
        free(bucketSize[n]);
    }
    free(bucketOffset);
    free(bucketRowDim);
    free(bucketSize);

    if (print_) outfile->Printf("\tConstructing frozen core operators\n");
    if (transformationType_ == TransformationType::Restricted) {
        // Compute frozen core energy
        size_t pq = 0;
        frozen_core_energy_ = 0.0;

        for (int p = 0; p < nso_; p++) {
            for (int q = 0; q <= p; q++, pq++) {
                double prefact = p == q ? 1.0 : 2.0;
                frozen_core_energy_ += prefact * aFzcD[pq] * (aoH[pq] + aFzcOp[pq]);
            }
        }

        auto aFzcMat = std::make_shared<Matrix>(PSIF_MO_FZC, sopi_, sopi_);
        aFzcMat->set(aFzcOp.data());
        aFzcMat->transform(Ca_);
        aFzcMat->save(psio_, PSIF_OEI, Matrix::SaveType::LowerTriangle);
    } else {
        // Compute frozen-core energy
        size_t pq = 0;
        frozen_core_energy_ = 0.0;
        for (int p = 0; p < nso_; p++) {
            for (int q = 0; q <= p; q++, pq++) {
                double prefact = p == q ? 0.5 : 1.0;
                frozen_core_energy_ += prefact * aFzcD[pq] * (aoH[pq] + aFzcOp[pq]);
                frozen_core_energy_ += prefact * bFzcD[pq] * (aoH[pq] + bFzcOp[pq]);
            }
        }

        auto aFzcMat = std::make_shared<Matrix>(PSIF_MO_A_FZC, sopi_, sopi_);
        aFzcMat->set(aFzcOp.data());
        aFzcMat->transform(Ca_);
        aFzcMat->save(psio_, PSIF_OEI, Matrix::SaveType::LowerTriangle);

        auto bFzcMat = std::make_shared<Matrix>(PSIF_MO_B_FZC, sopi_, sopi_);
        bFzcMat->set(bFzcOp.data());
        bFzcMat->transform(Cb_);
        bFzcMat->save(psio_, PSIF_OEI, Matrix::SaveType::LowerTriangle);
    }
    delete[] aoH;

    dpd_set_default(currentActiveDPD);

    alreadyPresorted_ = true;

    global_dpd_->file4_close(&I);
    psio_->close(PSIF_SO_PRESORT, 1);
}
