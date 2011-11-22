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
 * Presort the two-electron integrals into DPD buffers to prepare them for
 * the transformation.  The frozen core operator is built simultaneously.
 * If this action has already been performed, it will just load the frozen
 * core operator from disk and return.
 */
void
IntegralTransform::presort_so_tei()
{
    check_initialized();

    if(_alreadyPresorted){
        if(_print>5)
            fprintf(outfile, "\tSO integrals are already sorted, moving on...\n");
            return;
    }

    // Set aside some memory for the frozen core density and frozen core operator
    double *aFzcD  = init_array(_nTriSo);
    double *aFzcOp = init_array(_nTriSo);
    double *aD     = init_array(_nTriSo);
    double *aFock  = init_array(_nTriSo);
    double *aoH    = init_array(_nTriSo);
    double *bFzcD  = aFzcD;
    double *bFzcOp = aFzcOp;
    double *bD     = aD;
    double *bFock  = aFock;
    if(_transformationType != Restricted){
        bFzcD  = init_array(_nTriSo);
        bFzcOp = init_array(_nTriSo);
        bD     = init_array(_nTriSo);
        bFock  = init_array(_nTriSo);
    }

    // Form the Density matrices
    for(int h = 0, soOffset = 0; h < _nirreps; ++h){
        for(int p = 0; p < _sopi[h]; ++p){
            for(int q = 0; q <= p; ++q){
                int pq = INDEX((p + soOffset), (q + soOffset));
                for(int i = 0; i < _frzcpi[h]; ++i)
                    aFzcD[pq] += _Ca[h][p][i] * _Ca[h][q][i];
                for(int i = 0; i < _clsdpi[h] + _openpi[h]; ++i)
                    aD[pq] += _Ca[h][p][i] * _Ca[h][q][i];
                if(_transformationType != Restricted){
                    for(int i = 0; i < _frzcpi[h]; ++i)
                        bFzcD[pq] += _Cb[h][p][i] * _Cb[h][q][i];
                    for(int i = 0; i < _clsdpi[h]; ++i)
                        bD[pq] += _Cb[h][p][i] * _Cb[h][q][i];
                }
            }
        }
        soOffset += _sopi[h];
    }

    double *T = init_array(_nTriSo);
    if(_print>4) fprintf(outfile, "The SO basis kinetic energy integrals\n");
    IWL::read_one(_psio.get(), PSIF_OEI, PSIF_SO_T,   T, _nTriSo, 0, _print > 4, outfile);
    if(_print>4) fprintf(outfile, "The SO basis nuclear attraction integrals\n");
    IWL::read_one(_psio.get(), PSIF_OEI, PSIF_SO_V, aoH, _nTriSo, 0, _print > 4, outfile);

    for(int pq=0; pq < _nTriSo; ++pq){
        aoH[pq] += T[pq];
        aFzcOp[pq] = aoH[pq];
        aFock[pq]  = aoH[pq];
        if(_transformationType != Restricted){
            bFock[pq]  = aoH[pq];
            bFzcOp[pq] = aoH[pq];
        }
    }
    free(T);

    int currentActiveDPD = psi::dpd_default;
    dpd_set_default(_myDPDNum);

    if(_print){
        fprintf(outfile, "\tPresorting SO-basis two-electron integrals.\n");
        fflush(outfile);
    }

    int soIntFile = PSIF_SO_TEI;

    dpdfile4 I;
    _psio->open(PSIF_SO_PRESORT, PSIO_OPEN_NEW);
    dpd_file4_init(&I, PSIF_SO_PRESORT, 0, 3, 3, "SO Ints (nn|nn)");

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

        IWL *iwl = new IWL(_psio.get(), soIntFile, _tolerance, 1, 1);

        Label *lblptr = iwl->labels();
        Value *valptr = iwl->values();
        int lastbuf   = iwl->last_buffer();

        for(int index = iwl->index(); index < iwl->buffer_count(); ++index){
            int labelIndex = 4*index;
            int p = abs((int) lblptr[labelIndex++]);
            int q = (int) lblptr[labelIndex++];
            int r = (int) lblptr[labelIndex++];
            int s = (int) lblptr[labelIndex++];
            double value = (double) valptr[index];
            idx_permute_presort(&I,n,bucketMap,bucketOffset,p,q,r,s,value);
            if(!n) /* build frozen-core operator and Fock matrix only on first pass*/
                build_fzc_and_fock(p, q, r, s, value, aFzcD, bFzcD,
                                   aFzcOp, bFzcOp, aD, bD, aFock, bFock);
        } /* end loop through current buffer */

        /* Now run through the rest of the buffers in the file */
        while(!lastbuf){
            iwl->fetch();
            lastbuf = iwl->last_buffer();
            for(int index = iwl->index(); index < iwl->buffer_count(); ++index){
                int labelIndex = 4*index;
                int p = abs((int) lblptr[labelIndex++]);
                int q = (int) lblptr[labelIndex++];
                int r = (int) lblptr[labelIndex++];
                int s = (int) lblptr[labelIndex++];
                double value = (double) valptr[index];
                idx_permute_presort(&I,n,bucketMap,bucketOffset,p,q,r,s,value);
                if(!n) /* build frozen-core operator and Fock matrix only on first pass*/
                    build_fzc_and_fock(p, q, r, s, value, aFzcD, bFzcD,
                                       aFzcOp, bFzcOp, aD, bD, aFock, bFock);
            } /* end loop through current buffer */
        } /* end loop over reading buffers */
        iwl->set_keep_flag(1);
        delete iwl;

        for(int h=0; h < _nirreps; ++h) {
            if(bucketSize[n][h])
                _psio->write(I.filenum, I.label, (char *) I.matrix[h][0],
                bucketSize[n][h]*((long int) sizeof(double)), next, &next);
            free_block(I.matrix[h]);
        }
    } /* end loop over buckets/passes */

    /* Get rid of the input integral file */
    _psio->open(soIntFile, PSIO_OPEN_OLD);
    _psio->close(soIntFile, _keepIwlSoInts);

    free_int_matrix(bucketMap);

    for(int n=0; n < nBuckets; ++n) {
        free(bucketOffset[n]);
        free(bucketRowDim[n]);
        free(bucketSize[n]);
    }
    free(bucketOffset);
    free(bucketRowDim);
    free(bucketSize);

    double *moInts = init_array(_nTriMo);
    int *order = init_int_array(_nmo);
    // We want to keep Pitzer ordering, so this is just an identity mapping
    for(int n = 0; n < _nmo; ++n) order[n] = n;
    if(_print)
        fprintf(outfile, "\tTransforming the one-electron integrals and constructing Fock matrices\n");
    if(_transformationType == Restricted){
        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_sopi[h], _mopi[h], aoH, moInts, _Ca[h], soOffset, &(order[moOffset]));
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The MO basis one-electron integrals\n");
            print_array(moInts, _nmo, outfile);
        }
        IWL::write_one(_psio.get(), PSIF_OEI, PSIF_MO_OEI, _nTriMo, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_sopi[h], _mopi[h], aFzcOp, moInts, _Ca[h], soOffset, &(order[moOffset]));
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The MO basis frozen core operator\n");
            print_array(moInts, _nmo, outfile);
        }
        IWL::write_one(_psio.get(), PSIF_OEI, PSIF_MO_FZC, _nTriMo, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_sopi[h], _mopi[h], aFock, moInts, _Ca[h], soOffset, &(order[moOffset]));
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The MO basis Fock operator\n");
            print_array(moInts, _nmo, outfile);
        }

        IWL::write_one(_psio.get(), PSIF_OEI, PSIF_MO_FOCK, _nTriMo, aFock);
    }else{
        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_sopi[h], _mopi[h], aoH, moInts, _Ca[h], soOffset, &(order[moOffset]));
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The MO basis alpha one-electron integrals\n");
            print_array(moInts, _nmo, outfile);
        }
        IWL::write_one(_psio.get(), PSIF_OEI, PSIF_MO_A_OEI, _nTriMo, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_sopi[h], _mopi[h], aoH, moInts, _Cb[h], soOffset, &(order[moOffset]));
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The MO basis beta one-electron integrals\n");
            print_array(moInts, _nmo, outfile);
        }
        IWL::write_one(_psio.get(), PSIF_OEI, PSIF_MO_B_OEI, _nTriMo, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_sopi[h], _mopi[h], aFzcOp, moInts, _Ca[h], soOffset, &(order[moOffset]));
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The MO basis alpha frozen core operator\n");
            print_array(moInts, _nmo, outfile);
        }
        IWL::write_one(_psio.get(), PSIF_OEI, PSIF_MO_A_FZC, _nTriMo, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_sopi[h], _mopi[h], bFzcOp, moInts, _Cb[h], soOffset, &(order[moOffset]));
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The MO basis beta frozen core operator\n");
            print_array(moInts, _nmo, outfile);
        }
        IWL::write_one(_psio.get(), PSIF_OEI, PSIF_MO_B_FZC, _nTriMo, moInts);
        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_sopi[h], _mopi[h], aFock, moInts, _Ca[h], soOffset, &(order[moOffset]));
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The MO basis alpha Fock operator\n");
            print_array(moInts, _nmo, outfile);
        }
        IWL::write_one(_psio.get(), PSIF_OEI, PSIF_MO_A_FOCK, _nTriMo, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_sopi[h], _mopi[h], bFock, moInts, _Cb[h], soOffset, &(order[moOffset]));
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The MO basis beta Fock operator\n");
            print_array(moInts, _nmo, outfile);
        }
        IWL::write_one(_psio.get(), PSIF_OEI, PSIF_MO_B_FOCK, _nTriMo, moInts);
    }
    free(order);
    free(moInts);
    free(aFzcD);
    free(aFzcOp);
    free(aD);
    free(aoH);
    free(aFock);
    if(_transformationType != Restricted){
        free(bFzcD);
        free(bFzcOp);
        free(bD);
        free(bFock);
    }

    dpd_set_default(currentActiveDPD);

    _alreadyPresorted = true;

    dpd_file4_close(&I);
    _psio->close(PSIF_SO_PRESORT, 1);
}

/**
 * Builds the frozen core operator and Fock matrix using the integral currently
 * in memory. N.B. all matrices are passed in lower triangular array form.
 *
 * @param p - the first index in the integral
 * @param q - the second index in the integral
 * @param r - the third index in the integral
 * @param s - the fourth index in the integral
 * @param value - the integral (pq|rs)
 * @param aFzcD - the alpha frozen core density
 * @param bFzcD - the beta frozen core density
 * @param aFzcOp - the alpha frozen core operator
 * @param bFzcOp - the beta frozen core operator
 * @param aD - the alpha density
 * @param bD - the beta density
 * @param aFock - the alpha Fock matrix
 * @param bFock - the beta Fock matrix
 */
void
IntegralTransform::build_fzc_and_fock(int p, int q, int r, int s, double value,
                  double *aFzcD, double *bFzcD, double *aFzcOp, double *bFzcOp,
                  double *aD, double *bD, double *aFock, double *bFock)
{
    int al[8], bl[8], cl[8], dl[8];
    int dum, found;

    if(_transformationType == Restricted) {
        int a = al[0] = p;
        int b = bl[0] = q;
        int c = cl[0] = r;
        int d = dl[0] = s;
        int ab = INDEX(a,b);
        int cd = INDEX(c,d);
        int bc = INDEX(b,c);
        int ad = INDEX(a,d);
        aFzcOp[cd] += 2.0 * aFzcD[ab] * value;
        aFock[cd]  += 2.0 * aD[ab] * value;
        if(b >= c){
            aFzcOp[bc] -= aFzcD[ad] * value;
            aFock[bc]  -= aD[ad] * value;
        }

        a = al[1] = q;
        b = bl[1] = p;
        c = cl[1] = r;
        d = dl[1] = s;
        if(!(a==al[0] && b==bl[0] && c==cl[0] && d==dl[0])) {
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                aFzcOp[cd] += 2.0 * aFzcD[ab] * value;
                aFock[cd]  += 2.0 * aD[ab] * value;
            }
            if(b >= c){
                aFzcOp[bc] -= aFzcD[ad] * value;
                aFock[bc]  -= aD[ad] * value;
            }
        }

        a = al[2] = p;
        b = bl[2] = q;
        c = cl[2] = s;
        d = dl[2] = r;
        for(dum=0,found=0; dum < 2 && !found; dum++)
            if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        if(!found){
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                aFzcOp[cd] += 2.0 * aFzcD[ab] * value;
                aFock[cd]  += 2.0 * aD[ab] * value;
            }
            if(b >= c){
                aFzcOp[bc] -= aFzcD[ad] * value;
                aFock[bc]  -= aD[ad] * value;
            }
        }

        a = al[3] = q;
        b = bl[3] = p;
        c = cl[3] = s;
        d = dl[3] = r;
        for(dum=0,found=0; dum < 3 && !found; dum++)
            if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        if(!found){
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                aFzcOp[cd] += 2.0 * aFzcD[ab] * value;
                aFock[cd]  += 2.0 * aD[ab] * value;
            }
            if(b >= c){
                aFzcOp[bc] -= aFzcD[ad] * value;
                aFock[bc]  -= aD[ad] * value;
            }
        }

        a = al[4] = r;
        b = bl[4] = s;
        c = cl[4] = p;
        d = dl[4] = q;
        for(dum=0,found=0; dum < 4 && !found; dum++)
            if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        if(!found){
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                aFzcOp[cd] += 2.0 * aFzcD[ab] * value;
                aFock[cd]  += 2.0 * aD[ab] * value;
            }
            if(b >= c){
                aFzcOp[bc] -= aFzcD[ad] * value;
                aFock[bc]  -= aD[ad] * value;
            }
        }

        a = al[5] = r;
        b = bl[5] = s;
        c = cl[5] = q;
        d = dl[5] = p;
        for(dum=0,found=0; dum < 5 && !found; dum++)
            if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        if(!found) {
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                aFzcOp[cd] += 2.0 * aFzcD[ab] * value;
                aFock[cd]  += 2.0 * aD[ab] * value;
            }
            if(b >= c){
                aFzcOp[bc] -= aFzcD[ad] * value;
                aFock[bc]  -= aD[ad] * value;
            }
        }

        a = al[6] = s;
        b = bl[6] = r;
        c = cl[6] = p;
        d = dl[6] = q;
        for(dum=0, found=0; dum < 6 && !found; ++dum)
            if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        if(!found){
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                aFzcOp[cd] += 2.0 * aFzcD[ab] * value;
                aFock[cd]  += 2.0 * aD[ab] * value;
            }
            if(b >= c){
                aFzcOp[bc] -= aFzcD[ad] * value;
                aFock[bc]  -= aD[ad] * value;
            }
        }

        a = al[7] = s;
        b = bl[7] = r;
        c = cl[7] = q;
        d = dl[7] = p;
        for(dum=0,found=0; dum < 7 && !found; dum++)
            if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        if(!found){
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                aFzcOp[cd] += 2.0 * aFzcD[ab] * value;
                aFock[cd]  += 2.0 * aD[ab] * value;
            }
            if(b >= c){
                aFzcOp[bc] -= aFzcD[ad] * value;
                aFock[bc]  -= aD[ad] * value;
            }
        }
    }else {
        /* Unrestricted */
        int a = al[0] = p;
        int b = bl[0] = q;
        int c = cl[0] = r;
        int d = dl[0] = s;
        int ab = INDEX(a,b);
        int cd = INDEX(c,d);
        int bc = INDEX(b,c);
        int ad = INDEX(a,d);
        aFzcOp[cd] += (aFzcD[ab] + bFzcD[ab]) * value;
        bFzcOp[cd] += (aFzcD[ab] + bFzcD[ab]) * value;
        aFock[cd]  += (aD[ab] + bD[ab]) * value;
        bFock[cd]  += (aD[ab] + bD[ab]) * value;
        if(b >= c) {
            aFzcOp[bc] -= aFzcD[ad] * value;
            aFock[bc]  -= aD[ad] * value;
            bFzcOp[bc] -= bFzcD[ad] * value;
            bFock[bc]  -= bD[ad] * value;
        }

        a = al[1] = q;
        b = bl[1] = p;
        c = cl[1] = r;
        d = dl[1] = s;
        if(!(a==al[0] && b==bl[0] && c==cl[0] && d==dl[0])){
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                aFzcOp[cd] += (aFzcD[ab] + bFzcD[ab]) * value;
                aFock[cd] += (aD[ab] + bD[ab]) * value;
                bFzcOp[cd] += (aFzcD[ab] + bFzcD[ab]) * value;
                bFock[cd] += (aD[ab] + bD[ab]) * value;
            }
            if(b >= c){
                aFzcOp[bc] -= aFzcD[ad] * value;
                aFock[bc]  -= aD[ad] * value;
                bFzcOp[bc] -= bFzcD[ad] * value;
                bFock[bc]  -= bD[ad] * value;
            }
        }

        a = al[2] = p;
        b = bl[2] = q;
        c = cl[2] = s;
        d = dl[2] = r;
        for(dum=0,found=0; dum < 2 && !found; dum++)
            if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        if(!found){
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                aFzcOp[cd] += (aFzcD[ab] + bFzcD[ab]) * value;
                aFock[cd] += (aD[ab] + bD[ab]) * value;
                bFzcOp[cd] += (aFzcD[ab] + bFzcD[ab]) * value;
                bFock[cd] += (aD[ab] + bD[ab]) * value;
            }
            if(b >= c){
                aFzcOp[bc] -= aFzcD[ad] * value;
                aFock[bc]  -= aD[ad] * value;
                bFzcOp[bc] -= bFzcD[ad] * value;
                bFock[bc]  -= bD[ad] * value;
            }
        }

        a = al[3] = q;
        b = bl[3] = p;
        c = cl[3] = s;
        d = dl[3] = r;
        for(dum=0,found=0; dum < 3 && !found; dum++)
            if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        if(!found){
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                aFzcOp[cd] += (aFzcD[ab] + bFzcD[ab]) * value;
                aFock[cd]  += (aD[ab] + bD[ab]) * value;
                bFzcOp[cd] += (aFzcD[ab] + bFzcD[ab]) * value;
                bFock[cd]  += (aD[ab] + bD[ab]) * value;
            }
            if(b >= c){
                aFzcOp[bc] -= aFzcD[ad] * value;
                aFock[bc]  -= aD[ad] * value;
                bFzcOp[bc] -= bFzcD[ad] * value;
                bFock[bc]  -= bD[ad] * value;
            }
        }

        a = al[4] = r;
        b = bl[4] = s;
        c = cl[4] = p;
        d = dl[4] = q;
        for(dum=0,found=0; dum < 4 && !found; dum++)
            if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        if(!found){
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                aFzcOp[cd] += (aFzcD[ab] + bFzcD[ab]) * value;
                aFock[cd]  += (aD[ab] + bD[ab]) * value;
                bFzcOp[cd] += (aFzcD[ab] + bFzcD[ab]) * value;
                bFock[cd]  += (aD[ab] + bD[ab]) * value;
            }
            if(b >= c){
                aFzcOp[bc] -= aFzcD[ad] * value;
                aFock[bc]  -= aD[ad] * value;
                bFzcOp[bc] -= bFzcD[ad] * value;
                bFock[bc]  -= bD[ad] * value;
            }
        }

        a = al[5] = r;
        b = bl[5] = s;
        c = cl[5] = q;
        d = dl[5] = p;
        for(dum=0,found=0; dum < 5 && !found; dum++)
            if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        if(!found){
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                aFzcOp[cd] += (aFzcD[ab] + bFzcD[ab]) * value;
                aFock[cd]  += (aD[ab] + bD[ab]) * value;
                bFzcOp[cd] += (aFzcD[ab] + bFzcD[ab]) * value;
                bFock[cd]  += (aD[ab] + bD[ab]) * value;
            }
            if(b >= c){
                aFzcOp[bc] -= aFzcD[ad] * value;
                aFock[bc]  -= aD[ad] * value;
                bFzcOp[bc] -= bFzcD[ad] * value;
                bFock[bc]  -= bD[ad] * value;
            }
        }

        a = al[6] = s;
        b = bl[6] = r;
        c = cl[6] = p;
        d = dl[6] = q;
        for(dum=0,found=0; dum < 6 && !found; dum++)
            if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        if(!found) {
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                aFzcOp[cd] += (aFzcD[ab] + bFzcD[ab]) * value;
                aFock[cd]  += (aD[ab] + bD[ab]) * value;
                bFzcOp[cd] += (aFzcD[ab] + bFzcD[ab]) * value;
                bFock[cd]  += (aD[ab] + bD[ab]) * value;
            }
            if(b >= c){
                aFzcOp[bc] -= aFzcD[ad] * value;
                aFock[bc]  -= aD[ad] * value;
                bFzcOp[bc] -= bFzcD[ad] * value;
                bFock[bc]  -= bD[ad] * value;
            }
        }

        a = al[7] = s;
        b = bl[7] = r;
        c = cl[7] = q;
        d = dl[7] = p;
        for(dum=0,found=0; dum < 7 && !found; dum++)
            if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        if(!found){
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d){
                aFzcOp[cd] += (aFzcD[ab] + bFzcD[ab]) * value;
                aFock[cd]  += (aD[ab] + bD[ab]) * value;
                bFzcOp[cd] += (aFzcD[ab] + bFzcD[ab]) * value;
                bFock[cd]  += (aD[ab] + bD[ab]) * value;
            }
            if(b >= c){
                aFzcOp[bc] -= aFzcD[ad] * value;
                aFock[bc]  -= aD[ad] * value;
                bFzcOp[bc] -= bFzcD[ad] * value;
                bFock[bc]  -= bD[ad] * value;
            }
        }
    }
}

