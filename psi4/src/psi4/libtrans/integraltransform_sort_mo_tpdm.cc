/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#include "integraltransform.h"
#include "integraltransform_functors.h"
#include "mospace.h"

#include "psi4/libpsio/psio.hpp"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/psifiles.h"

using namespace psi;

/**
 * Presort the (restricted) MO TPDM into DPD buffers to prepare it
 * for the the transformation.
 */
void IntegralTransform::presort_mo_tpdm_restricted() {
    check_initialized();

    if (tpdmAlreadyPresorted_) {
        if (print_ > 5) outfile->Printf("\tMO TPDM already sorted, moving on...\n");
        return;
    }

    int currentActiveDPD = psi::dpd_default;
    dpd_set_default(myDPDNum_);

    if (print_) {
        outfile->Printf("\tPresorting MO-basis TPDM.\n");
    }

    dpdfile4 I;
    psio_->open(PSIF_TPDM_PRESORT, PSIO_OPEN_NEW);
    global_dpd_->file4_init(&I, PSIF_TPDM_PRESORT, 0, DPD_ID("[A>=A]+"), DPD_ID("[A>=A]+"), "MO TPDM (AA|AA)");

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
    //int **bucketSize = (int **)malloc(sizeof(int *));
    //bucketSize[0] = init_int_array(nirreps_);
    size_t **bucketSize = (size_t **) malloc(sizeof(size_t *));
    bucketSize[0] = init_size_t_array(nirreps_);


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

                //p = static_cast<int **>(realloc(static_cast<void *>(bucketSize), nBuckets * sizeof(int *)));
                size_t **P = static_cast<size_t **>(realloc(static_cast<void *>(bucketSize),
                                                    nBuckets * sizeof(size_t *)));

                if (P == nullptr) {
                    throw PsiException("file_build: allocation error", __FILE__, __LINE__);
                } else {
                    bucketSize = P;
                }
                //bucketSize[nBuckets - 1] = init_int_array(nirreps_);
                bucketSize[nBuckets - 1] = init_size_t_array(nirreps_);
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

        IWL *iwl = new IWL(psio_.get(), PSIF_MO_TPDM, tolerance_, 1, 0);
        DPDFillerFunctor dpdFiller(&I, n, bucketMap, bucketOffset, true, true);

        Label *lblptr = iwl->labels();
        Value *valptr = iwl->values();
        int lastbuf;
        /* Now run through the IWL buffers */
        do {
            iwl->fetch();
            lastbuf = iwl->last_buffer();
            for (int index = 0; index < iwl->buffer_count(); ++index) {
                int labelIndex = 4 * index;
                int p = aCorrToPitzer_[std::abs((int)lblptr[labelIndex++])];
                int q = aCorrToPitzer_[(int)lblptr[labelIndex++]];
                int r = aCorrToPitzer_[(int)lblptr[labelIndex++]];
                int s = aCorrToPitzer_[(int)lblptr[labelIndex++]];
                auto value = (double)valptr[index];
                // Check:
                //                outfile->Printf("\t%4d %4d %4d %4d = %20.10f\n", p, q, r, s, value);
                dpdFiller(p, q, r, s, value);
            }               /* end loop through current buffer */
        } while (!lastbuf); /* end loop over reading buffers */
        iwl->set_keep_flag(true);
        delete iwl;

        for (int h = 0; h < nirreps_; ++h) {
            if (bucketSize[n][h])
                psio_->write(I.filenum, I.label, (char *)I.matrix[h][0], bucketSize[n][h] * ((long int)sizeof(double)),
                             next, &next);
            free_block(I.matrix[h]);
        }
    } /* end loop over buckets/passes */

    /* Get rid of the input integral file */
    psio_->open(PSIF_MO_TPDM, PSIO_OPEN_OLD);
    psio_->close(PSIF_MO_TPDM, keepIwlMoTpdm_);

    free_int_matrix(bucketMap);

    for (int n = 0; n < nBuckets; ++n) {
        free(bucketOffset[n]);
        free(bucketRowDim[n]);
        free(bucketSize[n]);
    }
    free(bucketOffset);
    free(bucketRowDim);
    free(bucketSize);

    dpd_set_default(currentActiveDPD);

    tpdmAlreadyPresorted_ = true;

    global_dpd_->file4_close(&I);
    psio_->close(PSIF_TPDM_PRESORT, 1);
}

/**
 * Presort the (unrestricted) MO TPDMs into DPD buffers to prepare them
 * for the the transformation.
 */
void IntegralTransform::presort_mo_tpdm_unrestricted() {
    check_initialized();

    int currentActiveDPD = psi::dpd_default;
    dpd_set_default(myDPDNum_);

    if (print_) {
        outfile->Printf("\tPresorting MO-basis TPDMs.\n");
    }

    dpdfile4 I;
    psio_->open(PSIF_TPDM_PRESORT, PSIO_OPEN_NEW);
    global_dpd_->file4_init(&I, PSIF_TPDM_PRESORT, 0, DPD_ID("[A>=A]+"), DPD_ID("[A>=A]+"), "MO TPDM (AA|AA)");

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
    int **bucketSize = (int **)malloc(sizeof(int *));
    bucketSize[0] = init_int_array(nirreps_);

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

                p = static_cast<int **>(realloc(static_cast<void *>(bucketSize), nBuckets * sizeof(int *)));
                if (p == nullptr) {
                    throw PsiException("file_build: allocation error", __FILE__, __LINE__);
                } else {
                    bucketSize = p;
                }
                bucketSize[nBuckets - 1] = init_int_array(nirreps_);
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

    // The alpha - alpha spin case
    next = PSIO_ZERO;
    for (int n = 0; n < nBuckets; ++n) { /* nbuckets = number of passes */
        /* Prepare target matrix */
        for (int h = 0; h < nirreps_; h++) {
            I.matrix[h] = block_matrix(bucketRowDim[n][h], I.params->coltot[h]);
        }
        IWL *iwl = new IWL(psio_.get(), PSIF_MO_AA_TPDM, tolerance_, 1, 0);
        DPDFillerFunctor aaDpdFiller(&I, n, bucketMap, bucketOffset, true, true);

        Label *lblptr = iwl->labels();
        Value *valptr = iwl->values();
        int lastbuf;
        /* Now run through the IWL buffers */
        do {
            iwl->fetch();
            lastbuf = iwl->last_buffer();
            for (int index = 0; index < iwl->buffer_count(); ++index) {
                int labelIndex = 4 * index;
                int p = aCorrToPitzer_[std::abs((int)lblptr[labelIndex++])];
                int q = aCorrToPitzer_[(int)lblptr[labelIndex++]];
                int r = aCorrToPitzer_[(int)lblptr[labelIndex++]];
                int s = aCorrToPitzer_[(int)lblptr[labelIndex++]];
                auto value = (double)valptr[index];
                aaDpdFiller(p, q, r, s, value);
            }               /* end loop through current buffer */
        } while (!lastbuf); /* end loop over reading buffers */
        iwl->set_keep_flag(true);
        delete iwl;

        for (int h = 0; h < nirreps_; ++h) {
            if (bucketSize[n][h])
                psio_->write(I.filenum, I.label, (char *)I.matrix[h][0], bucketSize[n][h] * ((long int)sizeof(double)),
                             next, &next);
            free_block(I.matrix[h]);
        }
    } /* end loop over buckets/passes */

    /* Get rid of the input integral file */
    psio_->open(PSIF_MO_AA_TPDM, PSIO_OPEN_OLD);
    psio_->close(PSIF_MO_AA_TPDM, keepIwlMoTpdm_);

    // The alpha - beta spin case
    global_dpd_->file4_init(&I, PSIF_TPDM_PRESORT, 0, DPD_ID("[A>=A]+"), DPD_ID("[a>=a]+"), "MO TPDM (AA|aa)");
    if (print_) {
        outfile->Printf("\tSorting File: %s nbuckets = %d\n", I.label, nBuckets);
    }
    next = PSIO_ZERO;
    for (int n = 0; n < nBuckets; ++n) { /* nbuckets = number of passes */
        /* Prepare target matrix */
        for (int h = 0; h < nirreps_; h++) {
            I.matrix[h] = block_matrix(bucketRowDim[n][h], I.params->coltot[h]);
        }
        IWL *iwl = new IWL(psio_.get(), PSIF_MO_AB_TPDM, tolerance_, 1, 0);
        DPDFillerFunctor abDpdFiller(&I, n, bucketMap, bucketOffset, true, false);

        Label *lblptr = iwl->labels();
        Value *valptr = iwl->values();
        int lastbuf;
        /* Now run through the IWL buffers */
        do {
            iwl->fetch();
            lastbuf = iwl->last_buffer();
            for (int index = 0; index < iwl->buffer_count(); ++index) {
                int labelIndex = 4 * index;
                int p = aCorrToPitzer_[std::abs((int)lblptr[labelIndex++])];
                int q = aCorrToPitzer_[(int)lblptr[labelIndex++]];
                int r = bCorrToPitzer_[(int)lblptr[labelIndex++]];
                int s = bCorrToPitzer_[(int)lblptr[labelIndex++]];
                auto value = (double)valptr[index];
                // Check:
                //                outfile->Printf("\t%4d %4d %4d %4d = %20.10f\n", p, q, r, s, value);
                abDpdFiller(p, q, r, s, value);
            }               /* end loop through current buffer */
        } while (!lastbuf); /* end loop over reading buffers */
        iwl->set_keep_flag(true);
        delete iwl;

        for (int h = 0; h < nirreps_; ++h) {
            if (bucketSize[n][h])
                psio_->write(I.filenum, I.label, (char *)I.matrix[h][0], bucketSize[n][h] * ((long int)sizeof(double)),
                             next, &next);
            free_block(I.matrix[h]);
        }
    } /* end loop over buckets/passes */

    /* Get rid of the input integral file */
    psio_->open(PSIF_MO_AB_TPDM, PSIO_OPEN_OLD);
    psio_->close(PSIF_MO_AB_TPDM, keepIwlMoTpdm_);

    // The beta - beta spin case
    global_dpd_->file4_init(&I, PSIF_TPDM_PRESORT, 0, DPD_ID("[a>=a]+"), DPD_ID("[a>=a]+"), "MO TPDM (aa|aa)");
    if (print_) {
        outfile->Printf("\tSorting File: %s nbuckets = %d\n", I.label, nBuckets);
    }
    next = PSIO_ZERO;
    for (int n = 0; n < nBuckets; ++n) { /* nbuckets = number of passes */
        /* Prepare target matrix */
        for (int h = 0; h < nirreps_; h++) {
            I.matrix[h] = block_matrix(bucketRowDim[n][h], I.params->coltot[h]);
        }
        IWL *iwl = new IWL(psio_.get(), PSIF_MO_BB_TPDM, tolerance_, 1, 0);
        DPDFillerFunctor bbDpdFiller(&I, n, bucketMap, bucketOffset, true, true);

        Label *lblptr = iwl->labels();
        Value *valptr = iwl->values();
        int lastbuf;
        /* Now run through the IWL buffers */
        do {
            iwl->fetch();
            lastbuf = iwl->last_buffer();
            for (int index = 0; index < iwl->buffer_count(); ++index) {
                int labelIndex = 4 * index;
                int p = bCorrToPitzer_[std::abs((int)lblptr[labelIndex++])];
                int q = bCorrToPitzer_[(int)lblptr[labelIndex++]];
                int r = bCorrToPitzer_[(int)lblptr[labelIndex++]];
                int s = bCorrToPitzer_[(int)lblptr[labelIndex++]];
                auto value = (double)valptr[index];
                bbDpdFiller(p, q, r, s, value);
            }               /* end loop through current buffer */
        } while (!lastbuf); /* end loop over reading buffers */
        iwl->set_keep_flag(true);
        delete iwl;

        for (int h = 0; h < nirreps_; ++h) {
            if (bucketSize[n][h])
                psio_->write(I.filenum, I.label, (char *)I.matrix[h][0], bucketSize[n][h] * ((long int)sizeof(double)),
                             next, &next);
            free_block(I.matrix[h]);
        }
    } /* end loop over buckets/passes */

    /* Get rid of the input integral file */
    psio_->open(PSIF_MO_BB_TPDM, PSIO_OPEN_OLD);
    psio_->close(PSIF_MO_BB_TPDM, keepIwlMoTpdm_);

    free_int_matrix(bucketMap);

    for (int n = 0; n < nBuckets; ++n) {
        free(bucketOffset[n]);
        free(bucketRowDim[n]);
        free(bucketSize[n]);
    }
    free(bucketOffset);
    free(bucketRowDim);
    free(bucketSize);

    dpd_set_default(currentActiveDPD);

    tpdmAlreadyPresorted_ = true;

    global_dpd_->file4_close(&I);
    psio_->close(PSIF_TPDM_PRESORT, 1);
}
