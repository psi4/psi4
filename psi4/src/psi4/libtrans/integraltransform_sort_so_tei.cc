/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libmints/matrix.h"
#include "psi4/psifiles.h"
#include "integraltransform_functors.h"
#include "mospace.h"
#define EXTERN
#include "psi4/libdpd/dpd.gbl"

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
std::vector<SharedMatrix>
IntegralTransform::compute_fock_like_matrices(SharedMatrix Hcore, std::vector<SharedMatrix> Cmats)
{
    // This function is supposed to only be called after an initial presort, but we'll check to make sure.
    if(!alreadyPresorted_)
        presort_so_tei();

    // Some
    int nBuckets, thisBucketRows;
    size_t rowsPerBucket, rowsLeft, memFree;

    int nmats = Cmats.size();

    // Form the Density matrices associated with the C matrices, and allocate the F matrices.
    std::vector<SharedMatrix> Dmats, Fmats;
    for (size_t N = 0; N < nmats; ++N) {
        SharedMatrix Cmat = Cmats[N];
        if(Cmat->rowspi() != sopi_)
            throw PSIEXCEPTION("Row dimension of C matrix is not equal to SOs per irrep in LibTrans::compute_fock_like_matrices()");
        SharedMatrix Fmat(new Matrix("F matrix", sopi_, sopi_));
        Fmats.push_back(Fmat);
        SharedMatrix Dmat(new Matrix("D matrix", sopi_, sopi_));
        Dmat->gemm(false, true, 1.0, Cmat, Cmat, 0.0);
        Dmats.push_back(Dmat);
    }

    psio_->open(PSIF_SO_PRESORT, PSIO_OPEN_OLD);

    // Grab control of DPD for now, but store the active number to restore it later
    int currentActiveDPD = psi::dpd_default;
    dpd_set_default(myDPDNum_);

    dpdbuf4 J;
    global_dpd_->buf4_init(&J, PSIF_SO_PRESORT, 0, DPD_ID("[n,n]"), DPD_ID("[n,n]"),
                           DPD_ID("[n>=n]+"), DPD_ID("[n>=n]+"), 0, "SO Ints (nn|nn)");

    for(int h=0; h < nirreps_; h++) {
        if(J.params->coltot[h] && J.params->rowtot[h]) {
            memFree = static_cast<size_t>(dpd_memfree() - J.params->coltot[h]);
            rowsPerBucket = memFree / (J.params->coltot[h]);
            if(rowsPerBucket > J.params->rowtot[h]) rowsPerBucket = (size_t) J.params->rowtot[h];
            nBuckets = static_cast<int>(ceil(static_cast<double>(J.params->rowtot[h])/
                                             static_cast<double>(rowsPerBucket)));
            rowsLeft = static_cast<size_t>(J.params->rowtot[h] % rowsPerBucket);
        }else{
            nBuckets = 0;
            rowsPerBucket = 0;
            rowsLeft = 0;
        }

        if(print_ > 1) {
            outfile->Printf( "\th = %d; memfree         = %lu\n", h, memFree);
            outfile->Printf( "\th = %d; rows_per_bucket = %lu\n", h, rowsPerBucket);
            outfile->Printf( "\th = %d; rows_left       = %lu\n", h, rowsLeft);
            outfile->Printf( "\th = %d; nbuckets        = %d\n", h, nBuckets);

        }

        // We assume that we're always building totally symmetric Fock matrices here.
        int sym = 0;

        for(int n=0; n < nBuckets; n++){
            if(nBuckets == 1)
                thisBucketRows = rowsPerBucket;
            else
                thisBucketRows = (n < nBuckets-1) ? rowsPerBucket : rowsLeft;
            global_dpd_->buf4_mat_irrep_rd_block(&J, h, n*rowsPerBucket, thisBucketRows);
            for(int pq=0; pq < thisBucketRows; pq++) {
                int pabs = J.params->roworb[h][pq][0];
                int qabs = J.params->roworb[h][pq][1];
                int psym = J.params->psym[pabs];
                int qsym = J.params->qsym[qabs];
                int prel = pabs - J.params->poff[psym];
                int qrel = qabs - J.params->qoff[qsym];
                for(int rs = 0; rs < J.params->coltot[h]; ++rs){
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
                    for(int N = 0; N < nmats; ++N){
                        SharedMatrix D = Dmats[N];
                        SharedMatrix F = Fmats[N];
                        if(pqsym == rssym && pqsym == sym){
                            F->add(rsym, rrel, srel, D->get(psym, prel, qrel) * value);
                        }

                        if(qrsym == pssym && qrsym == sym){
                            F->add(qsym, qrel, rrel, -0.5*D->get(psym, prel, srel) * value);
                        }
                    } /* N */
                } /* rs */
            } /* pq */
        }
        global_dpd_->buf4_mat_irrep_close_block(&J, h, rowsPerBucket);
    } /* h */

    for(int N = 0; N < nmats; ++N)
        Fmats[N]->add(Hcore);

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
void
IntegralTransform::presort_so_tei()
{
    check_initialized();

    if(alreadyPresorted_){
        if(print_>5)
            outfile->Printf( "\tSO integrals are already sorted, moving on...\n");
        return;
    }

    // Set aside some memory for the frozen core density and frozen core operator
    double *aFzcD  = init_array(nTriSo_);
    double *aFzcOp = init_array(nTriSo_);
    double *aD     = init_array(nTriSo_);
    double *aFock  = init_array(nTriSo_);
    double *aoH    = init_array(nTriSo_);
    double *bFzcD  = aFzcD;
    double *bFzcOp = aFzcOp;
    double *bD     = aD;
    double *bFock  = aFock;
    if(transformationType_ != Restricted){
        bFzcD  = init_array(nTriSo_);
        bFzcOp = init_array(nTriSo_);
        bD     = init_array(nTriSo_);
        bFock  = init_array(nTriSo_);
    }

    // Form the Density matrices
    for(int h = 0, soOffset = 0; h < nirreps_; ++h){
        double **pCa = Ca_->pointer(h);
        double **pCb = Cb_->pointer(h);
        for(int p = 0; p < sopi_[h]; ++p){
            for(int q = 0; q <= p; ++q){
                int pq = INDEX((p + soOffset), (q + soOffset));
                for(int i = 0; i < frzcpi_[h]; ++i)
                    aFzcD[pq] += pCa[p][i] * pCa[q][i];
                for(int i = 0; i < nalphapi_[h]; ++i)
                    aD[pq] += pCa[p][i] * pCa[q][i];
                if(transformationType_ != Restricted){
                    for(int i = 0; i < frzcpi_[h]; ++i)
                        bFzcD[pq] += pCb[p][i] * pCb[q][i];
                    for(int i = 0; i < nbetapi_[h]; ++i)
                        bD[pq] += pCb[p][i] * pCb[q][i];
                }
            }
        }
        soOffset += sopi_[h];
    }

    double *T = init_array(nTriSo_);
    if(print_>4) outfile->Printf( "The SO basis kinetic energy integrals\n");
    IWL::read_one(psio_.get(), PSIF_OEI, PSIF_SO_T,   T, nTriSo_, 0, print_ > 4, "outfile");
    if(print_>4) outfile->Printf( "The SO basis nuclear attraction integrals\n");
    IWL::read_one(psio_.get(), PSIF_OEI, PSIF_SO_V, aoH, nTriSo_, 0, print_ > 4, "outfile");

    for(int pq=0; pq < nTriSo_; ++pq){
        aoH[pq] += T[pq];
        aFzcOp[pq] = aoH[pq];
        aFock[pq]  = aoH[pq];
        if(transformationType_ != Restricted){
            bFock[pq]  = aoH[pq];
            bFzcOp[pq] = aoH[pq];
        }
    }
    free(T);

    int currentActiveDPD = psi::dpd_default;
    dpd_set_default(myDPDNum_);

    if(print_){
        outfile->Printf( "\tPresorting SO-basis two-electron integrals.\n");

    }

    dpdfile4 I;
    psio_->open(PSIF_SO_PRESORT, PSIO_OPEN_NEW);
    global_dpd_->file4_init(&I, PSIF_SO_PRESORT, 0, DPD_ID("[n>=n]+"), DPD_ID("[n>=n]+"), "SO Ints (nn|nn)");

    size_t memoryd = memory_ / sizeof(double);

    int nump = 0, numq = 0;
    for(int h=0; h < nirreps_; ++h){
        nump += I.params->ppi[h];
        numq += I.params->qpi[h];
    }
    int **bucketMap = init_int_matrix(nump, numq);

    /* Room for one bucket to begin with */
    int **bucketOffset = (int **) malloc(sizeof(int *));
    bucketOffset[0] = init_int_array(nirreps_);
    int **bucketRowDim = (int **) malloc(sizeof(int *));
    bucketRowDim[0] = init_int_array(nirreps_);
    long int **bucketSize = (long int **) malloc(sizeof(long int *));
    bucketSize[0] = init_long_int_array(nirreps_);

    /* Figure out how many passes we need and where each p,q goes */
    int nBuckets = 1;
    size_t coreLeft = memoryd;
    psio_address next;
    for(int h = 0; h < nirreps_; ++h){
        size_t rowLength = (size_t) I.params->coltot[h^(I.my_irrep)];
        for(int row=0; row < I.params->rowtot[h]; ++row) {
            if(coreLeft >= rowLength){
                coreLeft -= rowLength;
                bucketRowDim[nBuckets-1][h]++;
                bucketSize[nBuckets-1][h] += rowLength;
            } else {
                nBuckets++;
                coreLeft = memoryd - rowLength;
                /* Make room for another bucket */
		int **p;

		p = static_cast<int **>(realloc(static_cast<void *>(bucketOffset),
						nBuckets * sizeof(int *)));
		if(p == NULL) {
		  throw PsiException("file_build: allocation error", __FILE__, __LINE__);
		} else {
		  bucketOffset = p;
		}
                bucketOffset[nBuckets-1] = init_int_array(nirreps_);
                bucketOffset[nBuckets-1][h] = row;


		p = static_cast<int **>(realloc(static_cast<void *>(bucketRowDim),
						nBuckets * sizeof(int *)));
		if(p == NULL) {
		  throw PsiException("file_build: allocation error", __FILE__, __LINE__);
		} else {
		  bucketRowDim = p;
		}
		bucketRowDim[nBuckets-1] = init_int_array(nirreps_);
		bucketRowDim[nBuckets-1][h] = 1;


		long int **pp;
		pp = static_cast<long int **>(realloc(static_cast<void *>(bucketSize),
						nBuckets * sizeof(long int *)));
		if(pp == NULL) {
		  throw PsiException("file_build: allocation error", __FILE__, __LINE__);
		} else {
		  bucketSize = pp;
		}
                bucketSize[nBuckets-1] = init_long_int_array(nirreps_);
                bucketSize[nBuckets-1][h] = rowLength;
            }
            int p = I.params->roworb[h][row][0];
            int q = I.params->roworb[h][row][1];
            bucketMap[p][q] = nBuckets - 1;
        }
    }

    if(print_) {
        outfile->Printf( "\tSorting File: %s nbuckets = %d\n", I.label, nBuckets);

    }

    next = PSIO_ZERO;
    for(int n=0; n < nBuckets; ++n) { /* nbuckets = number of passes */
        /* Prepare target matrix */
        for(int h=0; h < nirreps_; h++) {
            I.matrix[h] = block_matrix(bucketRowDim[n][h], I.params->coltot[h]);
        }

        DPDFillerFunctor dpdfiller(&I,n,bucketMap,bucketOffset, false, true);
        NullFunctor null;
        IWL *iwl = new IWL(psio_.get(), soIntTEIFile_, tolerance_, 1, 1);
        // In the functors below, we only want to build the Fock matrix on the first pass
        if(transformationType_ == Restricted){
            FrozenCoreAndFockRestrictedFunctor fock(aD, aFzcD,aFock,aFzcOp);
            if(n)
                iwl_integrals(iwl, dpdfiller, null);
            else
                iwl_integrals(iwl, dpdfiller, fock);
        }else{
            FrozenCoreAndFockUnrestrictedFunctor fock(aD, bD, aFzcD, bFzcD,
                                                      aFock, bFock, aFzcOp, bFzcOp);
            if(n)
                iwl_integrals(iwl, dpdfiller, null);
            else
                iwl_integrals(iwl, dpdfiller, fock);
        }
        delete iwl;


        for(int h=0; h < nirreps_; ++h) {
            if(bucketSize[n][h])
                psio_->write(I.filenum, I.label, (char *) I.matrix[h][0],
                        bucketSize[n][h]*((long int) sizeof(double)), next, &next);
            free_block(I.matrix[h]);
        }
    } /* end loop over buckets/passes */

    /* Get rid of the input integral file */
    psio_->open(soIntTEIFile_, PSIO_OPEN_OLD);
    psio_->close(soIntTEIFile_, keepIwlSoInts_);

    free_int_matrix(bucketMap);

    for(int n=0; n < nBuckets; ++n) {
        free(bucketOffset[n]);
        free(bucketRowDim[n]);
        free(bucketSize[n]);
    }
    free(bucketOffset);
    free(bucketRowDim);
    free(bucketSize);

    double *moInts = init_array(nTriMo_);
    int *order = init_int_array(nmo_);
    // We want to keep Pitzer ordering, so this is just an identity mapping
    for(int n = 0; n < nmo_; ++n) order[n] = n;
    if(print_)
        outfile->Printf( "\tTransforming the one-electron integrals and constructing Fock matrices\n");
    if(transformationType_ == Restricted){

        // Compute frozen core energy
        size_t pq = 0;
        frozen_core_energy_ = 0.0;

        for(int p=0; p < nso_; p++) {
            for(int q=0; q <= p ; q++, pq++) {
                double prefact = p == q ? 1.0 : 2.0;
                frozen_core_energy_ += prefact * aFzcD[pq] * (aoH[pq] + aFzcOp[pq]);
            }
        }

        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCa = Ca_->pointer(h);
            trans_one(sopi_[h], mopi_[h], aoH, moInts, pCa, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            outfile->Printf( "The MO basis one-electron integrals\n");
            print_array(moInts, nmo_, "outfile");
        }
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_OEI, nTriMo_, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCa = Ca_->pointer(h);
            trans_one(sopi_[h], mopi_[h], aFzcOp, moInts, pCa, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            outfile->Printf( "The MO basis frozen core operator\n");
            print_array(moInts, nmo_, "outfile");
        }
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_FZC, nTriMo_, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCa = Ca_->pointer(h);
            trans_one(sopi_[h], mopi_[h], aFock, moInts, pCa, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            outfile->Printf( "The MO basis Fock operator\n");
            print_array(moInts, nmo_, "outfile");
        }

        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_FOCK, nTriMo_, moInts);
    }else{

        // Compute frozen-core energy
        size_t pq = 0;
        frozen_core_energy_ = 0.0;
        for(int p=0; p < nso_; p++) {
            for(int q=0; q <= p; q++, pq++) {
                double prefact = p == q ? 0.5 : 1.0;
                frozen_core_energy_ += prefact * aFzcD[pq] * (aoH[pq] + aFzcOp[pq]);
                frozen_core_energy_ += prefact * bFzcD[pq] * (aoH[pq] + bFzcOp[pq]);
            }
        }

        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCa = Ca_->pointer(h);
            trans_one(sopi_[h], mopi_[h], aoH, moInts, pCa, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            outfile->Printf( "The MO basis alpha one-electron integrals\n");
            print_array(moInts, nmo_, "outfile");
        }
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_A_OEI, nTriMo_, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCb = Cb_->pointer(h);
            trans_one(sopi_[h], mopi_[h], aoH, moInts, pCb, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            outfile->Printf( "The MO basis beta one-electron integrals\n");
            print_array(moInts, nmo_, "outfile");
        }
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_B_OEI, nTriMo_, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCa = Ca_->pointer(h);
            trans_one(sopi_[h], mopi_[h], aFzcOp, moInts, pCa, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            outfile->Printf( "The MO basis alpha frozen core operator\n");
            print_array(moInts, nmo_, "outfile");
        }
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_A_FZC, nTriMo_, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCb = Cb_->pointer(h);
            trans_one(sopi_[h], mopi_[h], bFzcOp, moInts, pCb, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            outfile->Printf( "The MO basis beta frozen core operator\n");
            print_array(moInts, nmo_, "outfile");
        }
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_B_FZC, nTriMo_, moInts);
        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCa = Ca_->pointer(h);
            trans_one(sopi_[h], mopi_[h], aFock, moInts, pCa, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            outfile->Printf( "The MO basis alpha Fock operator\n");
            print_array(moInts, nmo_, "outfile");
        }
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_A_FOCK, nTriMo_, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCb = Cb_->pointer(h);
            trans_one(sopi_[h], mopi_[h], bFock, moInts, pCb, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            outfile->Printf( "The MO basis beta Fock operator\n");
            print_array(moInts, nmo_, "outfile");
        }
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_B_FOCK, nTriMo_, moInts);
    }
    free(order);
    free(moInts);
    free(aFzcD);
    free(aFzcOp);
    free(aD);
    free(aoH);
    free(aFock);
    if(transformationType_ != Restricted){
        free(bFzcD);
        free(bFzcOp);
        free(bD);
        free(bFock);
    }

    dpd_set_default(currentActiveDPD);

    alreadyPresorted_ = true;

    global_dpd_->file4_close(&I);
    psio_->close(PSIF_SO_PRESORT, 1);
}
