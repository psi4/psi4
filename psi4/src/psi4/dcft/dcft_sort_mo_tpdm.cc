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

#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libtrans/integraltransform_functors.h"
#include "psi4/psifiles.h"
#include "psi4/libtrans/mospace.h"
#define EXTERN
#include "psi4/libdpd/dpd.gbl"
#include "dcft.h"
#include "defines.h"

namespace psi{ namespace dcft{

/** Presort TPDM (MO) for closed-shell odc-12
   ** In open-shell cases, TPDM (MO) is presorted by IntegralTransform::presort_mo_tpdm_restricted(),
   ** where AA (or BB) and AB contribution are treated differently.
   ** However, in closed-shell odc-12, AA, AB, and BB contribution are all saved in PSIF_MO_TPDM,
   ** where all values are treated as in AA matrix (which causes problems).
   ** Thus, we presort TPDM (MO) here and turn off IntegralTransform::presort_mo_tpdm_restricted().
   */
void DCFTSolver::presort_mo_tpdm_AB()
{
    int currentActiveDPD = psi::dpd_default;

    if(print_){
        outfile->Printf("\tPre-Presorting MO-basis TPDM: AB.\n\n");
    }

    dpdfile4 I;

    psio_->open(PSIF_TPDM_PRESORT, PSIO_OPEN_NEW);
    global_dpd_->file4_init(&I, PSIF_TPDM_PRESORT, 0, ID("[A>=A]+"), ID("[A>=A]+"), "MO TPDM (AA|AA)");

    size_t memoryd = Process::environment.get_memory() / sizeof(double);

    int nump = 0, numq = 0;
    for(int h = 0; h < nirrep_; ++h){
        nump += I.params->ppi[h];
        numq += I.params->qpi[h];
    }
    int **bucketMap = init_int_matrix(nump, numq);

    /* Room for one bucket to begin with */
    int **bucketOffset = (int **) malloc(sizeof(int *));
    bucketOffset[0] = init_int_array(nirrep_);
    int **bucketRowDim = (int **) malloc(sizeof(int *));
    bucketRowDim[0] = init_int_array(nirrep_);
    int **bucketSize = (int **) malloc(sizeof(int *));
    bucketSize[0] = init_int_array(nirrep_);

    /* Figure out how many passes we need and where each p,q goes */
    int nBuckets = 1;
    size_t coreLeft = memoryd;
    psio_address next;
    for(int h = 0; h < nirrep_; ++h){
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
                bucketOffset[nBuckets-1] = init_int_array(nirrep_);
                bucketOffset[nBuckets-1][h] = row;


		p = static_cast<int **>(realloc(static_cast<void *>(bucketRowDim),
						nBuckets * sizeof(int *)));
		if(p == NULL) {
		  throw PsiException("file_build: allocation error", __FILE__, __LINE__);
		} else {
		  bucketRowDim = p;
		}
		bucketRowDim[nBuckets-1] = init_int_array(nirrep_);
		bucketRowDim[nBuckets-1][h] = 1;


		p = static_cast<int **>(realloc(static_cast<void *>(bucketSize),
						nBuckets * sizeof(int *)));
		if(p == NULL) {
		  throw PsiException("file_build: allocation error", __FILE__, __LINE__);
		} else {
		  bucketSize = p;
		}
                bucketSize[nBuckets-1] = init_int_array(nirrep_);
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
        for(int h=0; h < nirrep_; h++) {
            I.matrix[h] = block_matrix(bucketRowDim[n][h], I.params->coltot[h]);
        }

        IWL *iwl = new IWL(psio_.get(), PSIF_MO_TPDM, 1.0E-16, 1, 0);
        DPDFillerFunctor dpdFiller(&I, n, bucketMap, bucketOffset, true, false);

        Label *lblptr = iwl->labels();
        Value *valptr = iwl->values();
        int lastbuf;
        /* Now run through the IWL buffers */
        do{
            iwl->fetch();
            lastbuf = iwl->last_buffer();
            for(int index = 0; index < iwl->buffer_count(); ++index){
                int labelIndex = 4*index;
                int p = _ints->alpha_corr_to_pitzer()[abs((int) lblptr[labelIndex++])];
                int q = _ints->alpha_corr_to_pitzer()[(int) lblptr[labelIndex++]];
                int r = _ints->alpha_corr_to_pitzer()[(int) lblptr[labelIndex++]];
                int s = _ints->alpha_corr_to_pitzer()[(int) lblptr[labelIndex++]];
                double value = (double) valptr[index];
                dpdFiller(p, q, r, s, value);
            } /* end loop through current buffer */
        } while(!lastbuf); /* end loop over reading buffers */
        iwl->set_keep_flag(1);
        delete iwl;

        for(int h=0; h < nirrep_; ++h) {
            if(bucketSize[n][h])
                psio_->write(I.filenum, I.label, (char *) I.matrix[h][0],
                bucketSize[n][h]*((long int) sizeof(double)), next, &next);
            free_block(I.matrix[h]);
        }
    } /* end loop over buckets/passes */

    /* Get rid of the input integral file */
    psio_->open(PSIF_MO_TPDM, PSIO_OPEN_OLD);
    psio_->close(PSIF_MO_TPDM, 1);

    free_int_matrix(bucketMap);

    for(int n=0; n < nBuckets; ++n) {
        free(bucketOffset[n]);
        free(bucketRowDim[n]);
        free(bucketSize[n]);
    }
    free(bucketOffset);
    free(bucketRowDim);
    free(bucketSize);

    dpd_set_default(currentActiveDPD);

    global_dpd_->file4_close(&I);
    psio_->close(PSIF_TPDM_PRESORT, 1);
}

void DCFTSolver::presort_mo_tpdm_AA()
{
    int currentActiveDPD = psi::dpd_default;

    if(print_){
        outfile->Printf("\tPre-Presorting MO-basis TPDM: AA and BB.\n\n");
    }

    dpdfile4 I;
    dpdbuf4 Ibuf, Itot;

    psio_->open(PSIF_TPDM_PRESORT, PSIO_OPEN_OLD);

    global_dpd_->buf4_init(&Ibuf, PSIF_TPDM_PRESORT, 0, ID("[A>=A]+"), ID("[A>=A]+"),
                           ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO TPDM (AA|AA)");
    global_dpd_->buf4_copy(&Ibuf, PSIF_TPDM_PRESORT, "MO TPDM (AA|AA) TEMP");
    global_dpd_->buf4_close(&Ibuf);

    global_dpd_->file4_init(&I, PSIF_TPDM_PRESORT, 0, ID("[A>=A]+"), ID("[A>=A]+"), "MO TPDM (AA|AA) TEMP");

    size_t memoryd = Process::environment.get_memory() / sizeof(double);

    int nump = 0, numq = 0;
    for(int h = 0; h < nirrep_; ++h){
        nump += I.params->ppi[h];
        numq += I.params->qpi[h];
    }
    int **bucketMap = init_int_matrix(nump, numq);

    /* Room for one bucket to begin with */
    int **bucketOffset = (int **) malloc(sizeof(int *));
    bucketOffset[0] = init_int_array(nirrep_);
    int **bucketRowDim = (int **) malloc(sizeof(int *));
    bucketRowDim[0] = init_int_array(nirrep_);
    int **bucketSize = (int **) malloc(sizeof(int *));
    bucketSize[0] = init_int_array(nirrep_);

    /* Figure out how many passes we need and where each p,q goes */
    int nBuckets = 1;
    size_t coreLeft = memoryd;
    psio_address next;
    for(int h = 0; h < nirrep_; ++h){
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
                bucketOffset[nBuckets-1] = init_int_array(nirrep_);
                bucketOffset[nBuckets-1][h] = row;


		p = static_cast<int **>(realloc(static_cast<void *>(bucketRowDim),
						nBuckets * sizeof(int *)));
		if(p == NULL) {
		  throw PsiException("file_build: allocation error", __FILE__, __LINE__);
		} else {
		  bucketRowDim = p;
		}
		bucketRowDim[nBuckets-1] = init_int_array(nirrep_);
		bucketRowDim[nBuckets-1][h] = 1;


		p = static_cast<int **>(realloc(static_cast<void *>(bucketSize),
						nBuckets * sizeof(int *)));
		if(p == NULL) {
		  throw PsiException("file_build: allocation error", __FILE__, __LINE__);
		} else {
		  bucketSize = p;
		}
                bucketSize[nBuckets-1] = init_int_array(nirrep_);
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

    outfile->Printf("\tnbuckets = %d\n\n", nBuckets);

    for(int n=0; n < nBuckets; ++n) { /* nbuckets = number of passes */
        /* Prepare target matrix */
        for(int h=0; h < nirrep_; h++) {
            I.matrix[h] = block_matrix(bucketRowDim[n][h], I.params->coltot[h]);
        }

        IWL *iwl = new IWL(psio_.get(), PSIF_MO_TPDM, 1.0E-16, 1, 0);
        DPDFillerFunctor dpdFiller(&I, n, bucketMap, bucketOffset, true, true);

        Label *lblptr = iwl->labels();
        Value *valptr = iwl->values();
        int lastbuf;
        /* Now run through the IWL buffers */
        do{
            iwl->fetch();
            lastbuf = iwl->last_buffer();
            for(int index = 0; index < iwl->buffer_count(); ++index){
                int labelIndex = 4*index;
                int p = _ints->alpha_corr_to_pitzer()[abs((int) lblptr[labelIndex++])];
                int q = _ints->alpha_corr_to_pitzer()[(int) lblptr[labelIndex++]];
                int r = _ints->alpha_corr_to_pitzer()[(int) lblptr[labelIndex++]];
                int s = _ints->alpha_corr_to_pitzer()[(int) lblptr[labelIndex++]];
                double value = (double) valptr[index];
                dpdFiller(p, q, r, s, value);

            } /* end loop through current buffer */
        } while(!lastbuf); /* end loop over reading buffers */
        iwl->set_keep_flag(1);
        delete iwl;

        for(int h=0; h < nirrep_; ++h) {
            if(bucketSize[n][h])
                psio_->write(I.filenum, I.label, (char *) I.matrix[h][0],
                bucketSize[n][h]*((long int) sizeof(double)), next, &next);
            free_block(I.matrix[h]);
        }
    } /* end loop over buckets/passes */


    /* Get rid of the input integral file */
    psio_->open(PSIF_MO_TPDM, PSIO_OPEN_OLD);
    psio_->close(PSIF_MO_TPDM, 1);

    free_int_matrix(bucketMap);

    for(int n=0; n < nBuckets; ++n) {
        free(bucketOffset[n]);
        free(bucketRowDim[n]);
        free(bucketSize[n]);
    }
    free(bucketOffset);
    free(bucketRowDim);
    free(bucketSize);

    dpd_set_default(currentActiveDPD);

    global_dpd_->file4_close(&I);

    global_dpd_->buf4_init(&Itot, PSIF_TPDM_PRESORT, 0, ID("[A>=A]+"), ID("[A>=A]+"),
                           ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO TPDM (AA|AA)");
    global_dpd_->buf4_init(&Ibuf, PSIF_TPDM_PRESORT, 0, ID("[A>=A]+"), ID("[A>=A]+"),
                           ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO TPDM (AA|AA) TEMP");
    global_dpd_->buf4_axpy(&Ibuf, &Itot, 1.0);
    global_dpd_->buf4_close(&Ibuf);
    global_dpd_->buf4_close(&Itot);

    psio_->close(PSIF_TPDM_PRESORT, 1);
}

}}
