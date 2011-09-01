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
 * Presort the MO TPDM into DPD buffers to prepare it for the
 * the transformation.
 */
void
IntegralTransform::presort_mo_tpdm_restricted()
{
    check_initialized();

    int nDocc = 0;
    for(int h = 0; h < _nirreps; ++h) nDocc += _clsdpi[h];

    // The IWL TPDM comes in QT-ordered, but we're effectively using Pitzer order in the DPD structures
    int *qtToPitzer = new int[_nmo];
    int qtCount = 0;
    // The frozen core orbitals
    int pitzerOffset = 0;
    for(int h = 0; h < _nirreps; ++h){
        for(int pitzer = 0; pitzer < _frzcpi[h]; ++pitzer){
            qtToPitzer[qtCount++] = pitzer + pitzerOffset;
        }
        pitzerOffset += _mopi[h];
    }
    // The doubly occupied orbitals
    pitzerOffset = 0;
    for(int h = 0; h < _nirreps; ++h){
        pitzerOffset += _frzcpi[h];
        for(int pitzer = 0; pitzer < _clsdpi[h] - _frzcpi[h]; ++pitzer){
            qtToPitzer[qtCount++] = pitzer + pitzerOffset;
        }
        pitzerOffset += _mopi[h] - _frzcpi[h];
    }
    // The active virtuals
    pitzerOffset = 0;
    for(int h = 0; h < _nirreps; ++h){
        pitzerOffset += _clsdpi[h];
        for(int pitzer = 0; pitzer < _mopi[h] - _clsdpi[h] - _frzvpi[h]; ++pitzer){
            qtToPitzer[qtCount++] = pitzer + pitzerOffset;
        }
        pitzerOffset += _mopi[h] - _clsdpi[h];
    }

    if(_tpdmAlreadyPresorted){
        if(_print>5)
            fprintf(outfile, "\tMO TPDM already sorted, moving on...\n");
            return;
    }

    int currentActiveDPD = psi::dpd_default;
    dpd_set_default(_myDPDNum);

    if(_print){
        fprintf(outfile, "\tPresorting MO-basis TPDM.\n");
        fflush(outfile);
    }

    dpdfile4 I;
    _psio->open(PSIF_TPDM_PRESORT, PSIO_OPEN_NEW);
    dpd_file4_init(&I, PSIF_TPDM_PRESORT, 0, DPD_ID("[A>=A]+"), DPD_ID("[A>=A]+"), "MO TPDM (AA|AA)");

    size_t memoryd = _memory / sizeof(double);

    int nump = 0, numq = 0;
    for(int h=0; h < _nirreps; ++h){
        nump += I.params->ppi[h];
        numq += I.params->qpi[h];
    }
    int **bucketMap = init_int_matrix(nump, numq);

    /* Room for one bucket to begin with */
    int **bucketOffset = (int **) malloc(sizeof(int *));
    bucketOffset[0] = init_int_array(_nirreps);
    int **bucketRowDim = (int **) malloc(sizeof(int *));
    bucketRowDim[0] = init_int_array(_nirreps);
    int **bucketSize = (int **) malloc(sizeof(int *));
    bucketSize[0] = init_int_array(_nirreps);

    /* Figure out how many passes we need and where each p,q goes */
    int nBuckets = 1;
    size_t coreLeft = memoryd;
    psio_address next;
    for(int h = 0; h < _nirreps; ++h){
        size_t rowLength = (size_t) I.params->coltot[h^(I.my_irrep)];
        for(int row=0; row < I.params->rowtot[h]; ++row) {
            if((coreLeft - rowLength) >= 0){
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
                bucketOffset[nBuckets-1] = init_int_array(_nirreps);
                bucketOffset[nBuckets-1][h] = row;

                bucketRowDim = (int **) realloc((void *) bucketRowDim,
                                             nBuckets * sizeof(int *));
                bucketRowDim[nBuckets-1] = init_int_array(_nirreps);
                bucketRowDim[nBuckets-1][h] = 1;

                bucketSize = (int **) realloc((void *) bucketSize,
                                                nBuckets * sizeof(int *));
                bucketSize[nBuckets-1] = init_int_array(_nirreps);
                bucketSize[nBuckets-1][h] = rowLength;
            }
            int p = I.params->roworb[h][row][0];
            int q = I.params->roworb[h][row][1];
            bucketMap[p][q] = nBuckets - 1;
        }
    }

    if(_print) {
        fprintf(outfile, "\tSorting File: %s nbuckets = %d\n", I.label, nBuckets);
        fflush(outfile);
    }

    next = PSIO_ZERO;
    for(int n=0; n < nBuckets; ++n) { /* nbuckets = number of passes */
        /* Prepare target matrix */
        for(int h=0; h < _nirreps; h++) {
            I.matrix[h] = block_matrix(bucketRowDim[n][h], I.params->coltot[h]);
        }
        IWL *iwl = new IWL(_psio.get(), PSIF_MO_TPDM, _tolerance, 1, 0);

        Label *lblptr = iwl->labels();
        Value *valptr = iwl->values();
        int lastbuf;
        /* Now run through the IWL buffers */
        do{
            iwl->fetch();
            lastbuf = iwl->last_buffer();
            for(int index = 0; index < iwl->buffer_count(); ++index){
                int labelIndex = 4*index;
                int p = qtToPitzer[abs((int) lblptr[labelIndex++])];
                int q = qtToPitzer[(int) lblptr[labelIndex++]];
                int r = qtToPitzer[(int) lblptr[labelIndex++]];
                int s = qtToPitzer[(int) lblptr[labelIndex++]];
                double value = (double) valptr[index];
//                fprintf(outfile, "%d %d %d %d -> %d %d %d %d %lf\n",p,q,r,s,qtToPitzer[p],qtToPitzer[q],qtToPitzer[r],qtToPitzer[s],value);
                idx_permute_presort(&I,n,bucketMap,bucketOffset,p,q,r,s,value, true);
            } /* end loop through current buffer */
        } while(!lastbuf); /* end loop over reading buffers */
        iwl->set_keep_flag(1);
        delete iwl;

//        // Add in the reference contribution
//        for(int qtI = 0; qtI < nDocc; ++qtI){
//            int i = qtToPitzer[qtI];
//            idx_permute_presort(&I,n,bucketMap,bucketOffset, i, i, i, i, 1.0, true);
//            for(int qtJ = 0; qtJ < qtI; ++qtJ){
//                int j = qtToPitzer[qtJ];
//                idx_permute_presort(&I,n,bucketMap,bucketOffset, i, j, j, i, -1.0, true);
//                idx_permute_presort(&I,n,bucketMap,bucketOffset, j, i, i, j, -1.0, true);
//                idx_permute_presort(&I,n,bucketMap,bucketOffset, i, i, j, j, 2.0, true);
//                idx_permute_presort(&I,n,bucketMap,bucketOffset, j, j, i, i, 2.0, true);
//            }
//        }

        for(int h=0; h < _nirreps; ++h) {
            if(bucketSize[n][h])
                _psio->write(I.filenum, I.label, (char *) I.matrix[h][0],
                bucketSize[n][h]*((long int) sizeof(double)), next, &next);
            free_block(I.matrix[h]);
        }
    } /* end loop over buckets/passes */

    /* Get rid of the input integral file */
    _psio->open(PSIF_MO_TPDM, PSIO_OPEN_OLD);
    _psio->close(PSIF_MO_TPDM, _keepIwlMoTpdm);

    delete [] qtToPitzer;
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

    _tpdmAlreadyPresorted = true;

    dpd_file4_close(&I);
    _psio->close(PSIF_TPDM_PRESORT, 1);
}
