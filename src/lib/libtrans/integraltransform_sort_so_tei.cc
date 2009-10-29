#include "integraltransform.h"
#include <libchkpt/chkpt.hpp>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "psifiles.h"
#include "mospace.h"
#include "spaceinfo.h"

namespace psi{ namespace libtrans{


/**
 * Presort the two-electron integrals into DPD buffers to prepare them for
 * the transformation.  The frozen core operator is built simultaneously.
 * If this action has already been performed, it will just load the frozen
 * core operator from disk and return.
 */
void
IntegralTransform::presort_so_tei()
{
    static bool alreadyPresorted = false;

    if(alreadyPresorted){
        if(_print>5)
            fprintf(outfile, "SO integrals are already sorted, moving on...\n");
            return;
    }

    // Set aside some memory for the frozen core density and frozen core operator
    _aFzcD  = init_array(_nTriSo);
    _aFzcOp = init_array(_nTriSo);
    if(_transformationType != Restricted){
        _bFzcD  = init_array(_nTriSo);
        _bFzcOp = init_array(_nTriSo);
    }else{
        _bFzcD  = _aFzcD;
        _bFzcOp = _aFzcOp;
    }
//TODO form frozen core density


    timer_on("presort");
    if(_print){
        fprintf(outfile, "\n\tPresorting SO-basis two-electron integrals.\n");
        fflush(outfile);
    }
    
    int soIntFile = PSIF_SO_TEI;
    
    dpdfile4 I;
    psio_open(PSIF_SO_PRESORT, 0);
    dpd_file4_init(&I, PSIF_SO_PRESORT, 0, 3, 3, "SO Ints (nn,nn)");

    size_t memoryd = _memory / sizeof(double);

    int nump = 0;
    int numq = 0;
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

        struct iwlbuf InBuf;
        iwl_buf_init(&InBuf, soIntFile, _tolerance, 1, 1);

        Label *lblptr = InBuf.labels;
        Value *valptr = InBuf.values;
        int lastbuf   = InBuf.lastbuf;

        for(int idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
            int p = abs((int) lblptr[idx++]);
            int q = (int) lblptr[idx++];
            int r = (int) lblptr[idx++];
            int s = (int) lblptr[idx++];
            double value = (double) valptr[InBuf.idx];
            idx_permute_presort(&I,n,bucketMap,bucketOffset,p,q,r,s,value);
            if(_nfzc && !n) /* build frozen-core operator only on first pass*/
                frozen_core(p,q,r,s,value);
        } /* end loop through current buffer */

        /* Now run through the rest of the buffers in the file */
        while(!lastbuf){
            iwl_buf_fetch(&InBuf);
            lastbuf = InBuf.lastbuf;
            for(int idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
                int p = abs((int) lblptr[idx++]);
                int q = (int) lblptr[idx++];
                int r = (int) lblptr[idx++];
                int s = (int) lblptr[idx++];
                double value = (double) valptr[InBuf.idx];
                idx_permute_presort(&I, n, bucketMap, bucketOffset, p, q, r, s, value);
                if(_nfzc && !n) /* build frozen-core operator only on first pass */
                    frozen_core(p,q,r,s,value);
            } /* end loop through current buffer */
        } /* end loop over reading buffers */

        iwl_buf_close(&InBuf, 1); /* close buffer for next pass */

        for(int h=0; h < _nirreps; ++h) {
            if(bucketSize[n][h])
                psio_write(I.filenum, I.label, (char *) I.matrix[h][0],
                bucketSize[n][h]*((long int) sizeof(double)), next, &next);
            free_block(I.matrix[h]);
        }
    } /* end loop over buckets/passes */

    /* Get rid of the input integral file */
    psio_open(soIntFile, PSIO_OPEN_OLD);
    psio_close(soIntFile, !_deleteIwlSoTei);

    free_int_matrix(bucketMap);

    for(int n=0; n < nBuckets; ++n) {
        free(bucketOffset[n]);
        free(bucketRowDim[n]);
        free(bucketSize[n]);
    }
    free(bucketOffset);
    free(bucketRowDim);
    free(bucketSize);

    alreadyPresorted = true;
    
    dpd_file4_close(&I);
    psio_close(PSIF_SO_PRESORT, 1);
    timer_off("presort");
}

/**
 * Builds the frozen core operator using the integral currently in memory
 * @param p - the first index in the integral
 * @param q - the second index in the integral
 * @param r - the third index in the integral
 * @param s - the fourth index in the integral
 * @param value - the integral (pq|rs)
 */
void
IntegralTransform::frozen_core(int p, int q, int r, int s, double value)
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
        _aFzcOp[cd] += 2.0 * _aFzcD[ab] * value;
        if(b >= c) _aFzcOp[bc] -= _aFzcD[ad] * value;

        a = al[1] = q;
        b = bl[1] = p;
        c = cl[1] = r;
        d = dl[1] = s;
        if(!(a==al[0] && b==bl[0] && c==cl[0] && d==dl[0])) {
            ab = INDEX(a,b);
            cd = INDEX(c,d);
            bc = INDEX(b,c);
            ad = INDEX(a,d);
            if(c >= d) _aFzcOp[cd] += 2.0 * _aFzcD[ab] * value;
            if(b >= c) _aFzcOp[bc] -= _aFzcD[ad] * value;
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
            if(c >= d) _aFzcOp[cd] += 2.0 * _aFzcD[ab] * value;
            if(b >= c) _aFzcOp[bc] -= _aFzcD[ad] * value;
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
            if(c >= d) _aFzcOp[cd] += 2.0 * _aFzcD[ab] * value;
            if(b >= c) _aFzcOp[bc] -= _aFzcD[ad] * value;
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
            if(c >= d) _aFzcOp[cd] += 2.0 * _aFzcD[ab] * value;
            if(b >= c) _aFzcOp[bc] -= _aFzcD[ad] * value;
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
            if(c >= d) _aFzcOp[cd] += 2.0 * _aFzcD[ab] * value;
            if(b >= c) _aFzcOp[bc] -= _aFzcD[ad] * value;
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
            if(c >= d) _aFzcOp[cd] += 2.0 * _aFzcD[ab] * value;
            if(b >= c) _aFzcOp[bc] -= _aFzcD[ad] * value;
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
            if(c >= d) _aFzcOp[cd] += 2.0 * _aFzcD[ab] * value;
            if(b >= c) _aFzcOp[bc] -= _aFzcD[ad] * value;
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
        _aFzcOp[cd] += (_aFzcD[ab] + _bFzcD[ab]) * value;
        _bFzcOp[cd] += (_aFzcD[ab] + _bFzcD[ab]) * value;
        if(b >= c) {
            _aFzcOp[bc] -= _aFzcD[ad] * value;
            _bFzcOp[bc] -= _bFzcD[ad] * value;
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
                _aFzcOp[cd] += (_aFzcD[ab] + _bFzcD[ab]) * value;
                _bFzcOp[cd] += (_aFzcD[ab] + _bFzcD[ab]) * value;
            }
            if(b >= c){
                _aFzcOp[bc] -= _aFzcD[ad] * value;
                _bFzcOp[bc] -= _bFzcD[ad] * value;
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
                _aFzcOp[cd] += (_aFzcD[ab] + _bFzcD[ab]) * value;
                _bFzcOp[cd] += (_aFzcD[ab] + _bFzcD[ab]) * value;
            }
            if(b >= c){
                _aFzcOp[bc] -= _aFzcD[ad] * value;
                _bFzcOp[bc] -= _bFzcD[ad] * value;
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
                _aFzcOp[cd] += (_aFzcD[ab] + _bFzcD[ab]) * value;
                _bFzcOp[cd] += (_aFzcD[ab] + _bFzcD[ab]) * value;
            }
            if(b >= c){
                _aFzcOp[bc] -= _aFzcD[ad] * value;
                _bFzcOp[bc] -= _bFzcD[ad] * value;
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
                _aFzcOp[cd] += (_aFzcD[ab] + _bFzcD[ab]) * value;
                _bFzcOp[cd] += (_aFzcD[ab] + _bFzcD[ab]) * value;
            }
            if(b >= c){
                _aFzcOp[bc] -= _aFzcD[ad] * value;
                _bFzcOp[bc] -= _bFzcD[ad] * value;
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
                _aFzcOp[cd] += (_aFzcD[ab] + _bFzcD[ab]) * value;
                _bFzcOp[cd] += (_aFzcD[ab] + _bFzcD[ab]) * value;
            }
            if(b >= c){
                _aFzcOp[bc] -= _aFzcD[ad] * value;
                _bFzcOp[bc] -= _bFzcD[ad] * value;
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
                _aFzcOp[cd] += (_aFzcD[ab] + _bFzcD[ab]) * value;
                _bFzcOp[cd] += (_aFzcD[ab] + _bFzcD[ab]) * value;
            }
            if(b >= c){
                _aFzcOp[bc] -= _aFzcD[ad] * value;
                _bFzcOp[bc] -= _bFzcD[ad] * value;
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
                _aFzcOp[cd] += (_aFzcD[ab] + _bFzcD[ab]) * value;
                _bFzcOp[cd] += (_aFzcD[ab] + _bFzcD[ab]) * value;
            }
            if(b >= c) {
                _aFzcOp[bc] -= _aFzcD[ad] * value;
                _bFzcOp[bc] -= _bFzcD[ad] * value;
            }
        }
    }
}

void
IntegralTransform::idx_permute_presort(dpdfile4 *File, int &thisBucket, int **&bucketMap,
                                       int **&bucketOffset, int &p, int &q, int &r, int &s,
                                       double &value)
{
    dpdparams4 *Params = File->params;
    int perm_pq = Params->perm_pq;
    int perm_rs = Params->perm_rs;

    /* Get the orbital symmetries */
    int p_sym = Params->psym[p];
    int q_sym = Params->qsym[q];
    int r_sym = Params->rsym[r];
    int s_sym = Params->ssym[s];
    int pq_sym = p_sym^q_sym;
    int rs_sym = r_sym^s_sym;

    /* The allowed (Mulliken) permutations are very simple in this case */
    if(bucketMap[p][q] == thisBucket) {
        /* Get the row and column indices and assign the value */
        int pq = Params->rowidx[p][q];
        int rs = Params->colidx[r][s];
        if((pq >= Params->rowtot[pq_sym]) || (rs >= Params->coltot[rs_sym]))
            idx_error("MP Params_make: pq, rs", p,q,r,s,pq,rs,pq_sym,rs_sym);

        int offset = bucketOffset[thisBucket][pq_sym];
        File->matrix[pq_sym][pq-offset][rs] = value;
    }

    if(bucketMap[r][s] == thisBucket) {

        int rs = Params->rowidx[r][s];
        int pq = Params->colidx[p][q];
        if((rs >= Params->rowtot[rs_sym])||(pq >= Params->coltot[pq_sym]))
            idx_error("MP Params_make: rs, pq", p,q,r,s,rs,pq,rs_sym,pq_sym);

        int offset = bucketOffset[thisBucket][rs_sym];
        File->matrix[rs_sym][rs-offset][pq] = value;
    }
}


void
IntegralTransform::idx_error(const char *message, int p, int q, int r, int s,
                             int pq, int rs, int pq_sym, int rs_sym)
{

    fprintf(outfile, "\n\tDPD Parameter Error in %s\n", message);
    fprintf(outfile,"\t-------------------------------------------------\n");
    fprintf(outfile,"\t    p      q      r      s  [   pq]  [   rs] pq_symm rs_symm\n");
    fprintf(outfile,"\t%5d  %5d  %5d  %5d  [%5d]  [%5d]   %1d   %1d\n", p,q,r,s,
      pq,rs,pq_sym,rs_sym);
    throw PsiException("DPD idx failure.", __FILE__, __LINE__);
}

}} // Namespaces