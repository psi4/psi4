#include "integraltransform.h"
#include <libchkpt/chkpt.hpp>
#include <libpsio/psio.hpp>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libiwl/iwl.hpp>
#include "psifiles.h"
#include "mospace.h"
#define EXTERN
#include <libdpd/dpd.gbl>

using namespace psi;

/**
 * Presort the (restricted) MO TPDM into DPD buffers to prepare it
 * for the the transformation.
 */
void
IntegralTransform::presort_mo_tpdm_restricted()
{
    check_initialized();

    int nDocc = 0;
    for(int h = 0; h < nirreps_; ++h) nDocc += clsdpi_[h];

//    // The IWL TPDM comes in QT-ordered, but we're effectively using Pitzer order in the DPD structures
//    int *qtToPitzer = new int[_nmo];
//    int qtCount = 0;
//    // The frozen core orbitals
//    int pitzerOffset = 0;
//    for(int h = 0; h < _nirreps; ++h){
//        for(int pitzer = 0; pitzer < _frzcpi[h]; ++pitzer){
//            qtToPitzer[qtCount++] = pitzer + pitzerOffset;
//        }
//        pitzerOffset += _mopi[h];
//    }
//    // The doubly occupied orbitals
//    pitzerOffset = 0;
//    for(int h = 0; h < _nirreps; ++h){
//        pitzerOffset += _frzcpi[h];
//        for(int pitzer = 0; pitzer < _clsdpi[h] - _frzcpi[h]; ++pitzer){
//            qtToPitzer[qtCount++] = pitzer + pitzerOffset;
//        }
//        pitzerOffset += _mopi[h] - _frzcpi[h];
//    }
//    // The active virtuals
//    pitzerOffset = 0;
//    for(int h = 0; h < _nirreps; ++h){
//        pitzerOffset += _clsdpi[h];
//        for(int pitzer = 0; pitzer < _mopi[h] - _clsdpi[h] - _frzvpi[h]; ++pitzer){
//            qtToPitzer[qtCount++] = pitzer + pitzerOffset;
//        }
//        pitzerOffset += _mopi[h] - _clsdpi[h];
//    }

    if(tpdmAlreadyPresorted_){
        if(print_>5)
            fprintf(outfile, "\tMO TPDM already sorted, moving on...\n");
            return;
    }

    int currentActiveDPD = psi::dpd_default;
    dpd_set_default(myDPDNum_);

    if(print_){
        fprintf(outfile, "\tPresorting MO-basis TPDM.\n");
        fflush(outfile);
    }

    dpdfile4 I;
    psio_->open(PSIF_TPDM_PRESORT, PSIO_OPEN_NEW);
    dpd_file4_init(&I, PSIF_TPDM_PRESORT, 0, DPD_ID("[A>=A]+"), DPD_ID("[A>=A]+"), "MO TPDM (AA|AA)");

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
    int **bucketSize = (int **) malloc(sizeof(int *));
    bucketSize[0] = init_int_array(nirreps_);

    /* Figure out how many passes we need and where each p,q goes */
    int nBuckets = 1;
    size_t coreLeft = memoryd;
    psio_address next;
    for(int h = 0; h < nirreps_; ++h){
        size_t rowLength = (size_t) I.params->coltot[h^(I.my_irrep)];
        for(int row=0; row < I.params->rowtot[h]; ++row) {
            if((coreLeft - rowLength) >= 0){  // <-- This is aways true (unsigned - unsigned >= 0)
                coreLeft -= rowLength;
                bucketRowDim[nBuckets-1][h]++;
                bucketSize[nBuckets-1][h] += rowLength;
            }
            else {
                nBuckets++;
                coreLeft = memoryd - rowLength;
                /* Make room for another bucket */
                bucketOffset = (int **) realloc((void *) bucketOffset,
                                             nBuckets * sizeof(int *));
                bucketOffset[nBuckets-1] = init_int_array(nirreps_);
                bucketOffset[nBuckets-1][h] = row;

                bucketRowDim = (int **) realloc((void *) bucketRowDim,
                                             nBuckets * sizeof(int *));
                bucketRowDim[nBuckets-1] = init_int_array(nirreps_);
                bucketRowDim[nBuckets-1][h] = 1;

                bucketSize = (int **) realloc((void *) bucketSize,
                                                nBuckets * sizeof(int *));
                bucketSize[nBuckets-1] = init_int_array(nirreps_);
                bucketSize[nBuckets-1][h] = rowLength;
            }
            int p = I.params->roworb[h][row][0];
            int q = I.params->roworb[h][row][1];
            bucketMap[p][q] = nBuckets - 1;
        }
    }

    if(print_) {
        fprintf(outfile, "\tSorting File: %s nbuckets = %d\n", I.label, nBuckets);
        fflush(outfile);
    }

    next = PSIO_ZERO;
    for(int n=0; n < nBuckets; ++n) { /* nbuckets = number of passes */
        /* Prepare target matrix */
        for(int h=0; h < nirreps_; h++) {
            I.matrix[h] = block_matrix(bucketRowDim[n][h], I.params->coltot[h]);
        }
        IWL *iwl = new IWL(psio_.get(), PSIF_MO_TPDM, tolerance_, 1, 0);

        Label *lblptr = iwl->labels();
        Value *valptr = iwl->values();
        int lastbuf;
        /* Now run through the IWL buffers */
        do{
            iwl->fetch();
            lastbuf = iwl->last_buffer();
            for(int index = 0; index < iwl->buffer_count(); ++index){
                int labelIndex = 4*index;
                int p = aCorrToPitzer_[abs((int) lblptr[labelIndex++])];
                int q = aCorrToPitzer_[(int) lblptr[labelIndex++]];
                int r = aCorrToPitzer_[(int) lblptr[labelIndex++]];
                int s = aCorrToPitzer_[(int) lblptr[labelIndex++]];
                double value = (double) valptr[index];
//                fprintf(outfile, "%d %d %d %d -> %d %d %d %d %lf\n",p,q,r,s,qtToPitzer[p],qtToPitzer[q],qtToPitzer[r],qtToPitzer[s],value);
                idx_permute_presort(&I,n,bucketMap,bucketOffset,p,q,r,s,value, true);
            } /* end loop through current buffer */
        } while(!lastbuf); /* end loop over reading buffers */
        iwl->set_keep_flag(1);
        delete iwl;

        for(int h=0; h < nirreps_; ++h) {
            if(bucketSize[n][h])
                psio_->write(I.filenum, I.label, (char *) I.matrix[h][0],
                bucketSize[n][h]*((long int) sizeof(double)), next, &next);
            free_block(I.matrix[h]);
        }
    } /* end loop over buckets/passes */

    /* Get rid of the input integral file */
    psio_->open(PSIF_MO_TPDM, PSIO_OPEN_OLD);
    psio_->close(PSIF_MO_TPDM, keepIwlMoTpdm_);

//    delete [] qtToPitzer;
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

    tpdmAlreadyPresorted_ = true;

    dpd_file4_close(&I);
    psio_->close(PSIF_TPDM_PRESORT, 1);
}



/**
 * Presort the (unrestricted) MO TPDMs into DPD buffers to prepare them
 * for the the transformation.
 */
void
IntegralTransform::presort_mo_tpdm_unrestricted()
{
    check_initialized();

    int currentActiveDPD = psi::dpd_default;
    dpd_set_default(myDPDNum_);

    if(print_){
        fprintf(outfile, "\tPresorting MO-basis TPDMs.\n");
        fflush(outfile);
    }

    dpdfile4 I;
    psio_->open(PSIF_TPDM_PRESORT, PSIO_OPEN_NEW);
    dpd_file4_init(&I, PSIF_TPDM_PRESORT, 0, DPD_ID("[A>=A]+"), DPD_ID("[A>=A]+"), "MO TPDM (AA|AA)");

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
    int **bucketSize = (int **) malloc(sizeof(int *));
    bucketSize[0] = init_int_array(nirreps_);

    /* Figure out how many passes we need and where each p,q goes */
    int nBuckets = 1;
    size_t coreLeft = memoryd;
    psio_address next;
    for(int h = 0; h < nirreps_; ++h){
        size_t rowLength = (size_t) I.params->coltot[h^(I.my_irrep)];
        for(int row=0; row < I.params->rowtot[h]; ++row) {
            if((coreLeft - rowLength) >= 0){  // <-- This is always true (unsigned - unsigned >= 0)
                coreLeft -= rowLength;
                bucketRowDim[nBuckets-1][h]++;
                bucketSize[nBuckets-1][h] += rowLength;
            }
            else {
                nBuckets++;
                coreLeft = memoryd - rowLength;
                /* Make room for another bucket */
                bucketOffset = (int **) realloc((void *) bucketOffset,
                                             nBuckets * sizeof(int *));
                bucketOffset[nBuckets-1] = init_int_array(nirreps_);
                bucketOffset[nBuckets-1][h] = row;

                bucketRowDim = (int **) realloc((void *) bucketRowDim,
                                             nBuckets * sizeof(int *));
                bucketRowDim[nBuckets-1] = init_int_array(nirreps_);
                bucketRowDim[nBuckets-1][h] = 1;

                bucketSize = (int **) realloc((void *) bucketSize,
                                                nBuckets * sizeof(int *));
                bucketSize[nBuckets-1] = init_int_array(nirreps_);
                bucketSize[nBuckets-1][h] = rowLength;
            }
            int p = I.params->roworb[h][row][0];
            int q = I.params->roworb[h][row][1];
            bucketMap[p][q] = nBuckets - 1;
        }
    }

    if(print_) {
        fprintf(outfile, "\tSorting File: %s nbuckets = %d\n", I.label, nBuckets);
        fflush(outfile);
    }

    // The alpha - alpha spin case
    next = PSIO_ZERO;
    for(int n=0; n < nBuckets; ++n) { /* nbuckets = number of passes */
        /* Prepare target matrix */
        for(int h=0; h < nirreps_; h++) {
            I.matrix[h] = block_matrix(bucketRowDim[n][h], I.params->coltot[h]);
        }
        IWL *iwl = new IWL(psio_.get(), PSIF_MO_AA_TPDM, tolerance_, 1, 0);

        Label *lblptr = iwl->labels();
        Value *valptr = iwl->values();
        int lastbuf;
        /* Now run through the IWL buffers */
        do{
            iwl->fetch();
            lastbuf = iwl->last_buffer();
            for(int index = 0; index < iwl->buffer_count(); ++index){
                int labelIndex = 4*index;
                int p = aCorrToPitzer_[abs((int) lblptr[labelIndex++])];
                int q = aCorrToPitzer_[(int) lblptr[labelIndex++]];
                int r = aCorrToPitzer_[(int) lblptr[labelIndex++]];
                int s = aCorrToPitzer_[(int) lblptr[labelIndex++]];
                double value = (double) valptr[index];
                idx_permute_presort(&I,n,bucketMap,bucketOffset,p,q,r,s,value, true, true);
            } /* end loop through current buffer */
        } while(!lastbuf); /* end loop over reading buffers */
        iwl->set_keep_flag(1);
        delete iwl;

        for(int h=0; h < nirreps_; ++h) {
            if(bucketSize[n][h])
                psio_->write(I.filenum, I.label, (char *) I.matrix[h][0],
                bucketSize[n][h]*((long int) sizeof(double)), next, &next);
            free_block(I.matrix[h]);
        }
    } /* end loop over buckets/passes */

    /* Get rid of the input integral file */
    psio_->open(PSIF_MO_AA_TPDM, PSIO_OPEN_OLD);
    psio_->close(PSIF_MO_AA_TPDM, keepIwlMoTpdm_);

    // The alpha - beta spin case
    dpd_file4_init(&I, PSIF_TPDM_PRESORT, 0, DPD_ID("[A>=A]+"), DPD_ID("[a>=a]+"), "MO TPDM (AA|aa)");
    if(print_) {
        fprintf(outfile, "\tSorting File: %s nbuckets = %d\n", I.label, nBuckets);
        fflush(outfile);
    }
    next = PSIO_ZERO;
    for(int n=0; n < nBuckets; ++n) { /* nbuckets = number of passes */
        /* Prepare target matrix */
        for(int h=0; h < nirreps_; h++) {
            I.matrix[h] = block_matrix(bucketRowDim[n][h], I.params->coltot[h]);
        }
        IWL *iwl = new IWL(psio_.get(), PSIF_MO_AB_TPDM, tolerance_, 1, 0);

        Label *lblptr = iwl->labels();
        Value *valptr = iwl->values();
        int lastbuf;
        /* Now run through the IWL buffers */
        do{
            iwl->fetch();
            lastbuf = iwl->last_buffer();
            for(int index = 0; index < iwl->buffer_count(); ++index){
                int labelIndex = 4*index;
                int p = aCorrToPitzer_[abs((int) lblptr[labelIndex++])];
                int q = aCorrToPitzer_[(int) lblptr[labelIndex++]];
                int r = bCorrToPitzer_[(int) lblptr[labelIndex++]];
                int s = bCorrToPitzer_[(int) lblptr[labelIndex++]];
                double value = (double) valptr[index];
                idx_permute_presort(&I,n,bucketMap,bucketOffset,p,q,r,s,value, true, false);
            } /* end loop through current buffer */
        } while(!lastbuf); /* end loop over reading buffers */
        iwl->set_keep_flag(1);
        delete iwl;

        for(int h=0; h < nirreps_; ++h) {
            if(bucketSize[n][h])
                psio_->write(I.filenum, I.label, (char *) I.matrix[h][0],
                bucketSize[n][h]*((long int) sizeof(double)), next, &next);
            free_block(I.matrix[h]);
        }
    } /* end loop over buckets/passes */

    /* Get rid of the input integral file */
    psio_->open(PSIF_MO_AB_TPDM, PSIO_OPEN_OLD);
    psio_->close(PSIF_MO_AB_TPDM, keepIwlMoTpdm_);

    // The beta - beta spin case
    dpd_file4_init(&I, PSIF_TPDM_PRESORT, 0, DPD_ID("[a>=a]+"), DPD_ID("[a>=a]+"), "MO TPDM (aa|aa)");
    if(print_) {
        fprintf(outfile, "\tSorting File: %s nbuckets = %d\n", I.label, nBuckets);
        fflush(outfile);
    }
    next = PSIO_ZERO;
    for(int n=0; n < nBuckets; ++n) { /* nbuckets = number of passes */
        /* Prepare target matrix */
        for(int h=0; h < nirreps_; h++) {
            I.matrix[h] = block_matrix(bucketRowDim[n][h], I.params->coltot[h]);
        }
        IWL *iwl = new IWL(psio_.get(), PSIF_MO_BB_TPDM, tolerance_, 1, 0);

        Label *lblptr = iwl->labels();
        Value *valptr = iwl->values();
        int lastbuf;
        /* Now run through the IWL buffers */
        do{
            iwl->fetch();
            lastbuf = iwl->last_buffer();
            for(int index = 0; index < iwl->buffer_count(); ++index){
                int labelIndex = 4*index;
                int p = bCorrToPitzer_[abs((int) lblptr[labelIndex++])];
                int q = bCorrToPitzer_[(int) lblptr[labelIndex++]];
                int r = bCorrToPitzer_[(int) lblptr[labelIndex++]];
                int s = bCorrToPitzer_[(int) lblptr[labelIndex++]];
                double value = (double) valptr[index];
                idx_permute_presort(&I,n,bucketMap,bucketOffset,p,q,r,s,value, true, true);
            } /* end loop through current buffer */
        } while(!lastbuf); /* end loop over reading buffers */
        iwl->set_keep_flag(1);
        delete iwl;

        for(int h=0; h < nirreps_; ++h) {
            if(bucketSize[n][h])
                psio_->write(I.filenum, I.label, (char *) I.matrix[h][0],
                bucketSize[n][h]*((long int) sizeof(double)), next, &next);
            free_block(I.matrix[h]);
        }
    } /* end loop over buckets/passes */

    /* Get rid of the input integral file */
    psio_->open(PSIF_MO_BB_TPDM, PSIO_OPEN_OLD);
    psio_->close(PSIF_MO_BB_TPDM, keepIwlMoTpdm_);


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

    tpdmAlreadyPresorted_ = true;

    dpd_file4_close(&I);
    psio_->close(PSIF_TPDM_PRESORT, 1);
}
