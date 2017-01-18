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

/*! \file
    \ingroup DPD
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libqt/qt.h"
#include "dpd.h"

namespace psi {



/* dpd_contract244(): Contracts a two-index quantity with a
** four-index quantity to compute the contribution to another
** four-index quantity, beta * Z = alpha * X * Y.
**
** Arguments:
**   dpdfile2 *X: A pointer to the two-index file.
**   dpdbuf4 *Y: A pointer to the four-index buffer.
**   dpdbuf4 *Z: A pointer to the target four-index buffer.
**   int sum_X: Indicates which index on X is to be summed (takes a value of
**              0 or 1).
**   int sum_Y: Indicates which index on Y is to be summed (takes a value of
**              0, 1, 2, or 3).
**   int Ztrans: A boolean which indicates whether the indices in Z
**                must be transposed in order to match the external
**                index ordering in the product X * Y.
**   double alpha: A prefactor for the product alpha * X * Y.
**   double beta: A prefactor for the target beta * Z.
*/

#define DPD_BIGNUM 2147483647 /* the four-byte signed int limit */

int DPD::contract244(dpdfile2 *X, dpdbuf4 *Y, dpdbuf4 *Z, int sum_X, int sum_Y,
                     int Ztrans, double alpha, double beta)
{
    int h, h0, Hx, hybuf, hzbuf, Hy, Hz, nirreps, GX, GY, GZ, bra_y;
    int rking=0, *yrow, *ycol, symlink;
    int Xtrans, Ytrans=0;
    int incore;
    int rowx, rowz, colx, colz;
    int pq, Gr, GsY, GsZ, Gs, GrZ, GrY;
    int ncols, nrows, nlinks;
    long int core, memoryd, core_total, rowtot, coltot, maxrows, Z_core;
    int *numlinks, *numrows, *numcols;
    dpdtrans4 Yt, Zt;
    double ***Ymat, ***Zmat;
#ifdef DPD_DEBUG
    int *xrow, *xcol, *zrow, *zcol;
#endif

    nirreps = Y->params->nirreps;
    GX = X->my_irrep;
    GY = Y->file.my_irrep;
    GZ = Z->file.my_irrep;

    memoryd = dpd_main.memory;
    incore = 1; /* default */

    file2_mat_init(X);
    file2_mat_rd(X);

    if(sum_X == 0) { Xtrans = 1; numlinks = X->params->rowtot; symlink = 0; }
    else if(sum_X == 1) { Xtrans = 0; numlinks = X->params->coltot; symlink = GX; }
    else { outfile->Printf( "Junk X index %d\n", sum_X); exit(PSI_RETURN_FAILURE); }

    if((sum_Y == 1) || (sum_Y == 2)) trans4_init(&Yt, Y);

    if(Ztrans) trans4_init(&Zt, Z);

    /*  if(fabs(beta) > 0.0) dpd_buf4_scm(Z, beta); */
    buf4_scm(Z, beta);

#ifdef DPD_DEBUG
    if(Xtrans) { xrow = X->params->coltot; xcol = X->params->rowtot; }
    else { xrow = X->params->rowtot; xcol = X->params->coltot; }
#endif

    for(hzbuf=0; hzbuf < nirreps; hzbuf++) {

        incore = 1; /* default */

        if (sum_Y < 2) {
            if (Ztrans) hybuf = hzbuf^GY; else hybuf = hzbuf^GZ^GY;
        }
        else {
            if (Ztrans) hybuf = hzbuf; else hybuf = hzbuf^GZ;
        }

        /* Compute the core requirements for the straight contraction */
        core_total = 0;
        /** Y terms **/
        coltot = Y->params->coltot[hybuf^GY];
        if(coltot) {
            maxrows = DPD_BIGNUM/coltot;
            if(maxrows < 1) {
                outfile->Printf( "\nLIBDPD Error: cannot compute even the number of rows in contract244.\n");
                dpd_error("contract244", "outfile");
            }
        }
        else maxrows = DPD_BIGNUM;
        rowtot = Y->params->rowtot[hybuf];
        for(; rowtot > maxrows; rowtot -= maxrows) {
            if(core_total > (core_total + maxrows*coltot)) incore = 0;
            else core_total += maxrows * coltot;
        }
        if(core_total > (core_total + rowtot*coltot)) incore = 0;
        core_total += rowtot * coltot;

        if(sum_Y == 1 || sum_Y == 2) core_total *= 2;  /* we need room to transpose the Y buffer */

        /** Z terms **/
        coltot = Z->params->coltot[hzbuf^GZ];
        if(coltot) {
            maxrows = DPD_BIGNUM/coltot;
            if(maxrows < 1) {
                outfile->Printf( "\nLIBDPD Error: cannot compute even the number of rows in contract244.\n");
                dpd_error("contract244", "outfile");
            }
        }
        else maxrows = DPD_BIGNUM;
        rowtot = Z->params->rowtot[hzbuf];
        Z_core = maxrows * coltot;
        if(Ztrans) Z_core *= 2;
        for(; rowtot > maxrows; rowtot -= maxrows) {
            if(core_total > (core_total + Z_core)) incore = 0;
            else core_total += Z_core;
        }
        Z_core = rowtot * coltot;
        if(Ztrans) Z_core *= 2;
        if(core_total > (core_total + Z_core)) incore = 0;
        core_total += Z_core;

        if(core_total > memoryd) incore = 0;

        /* Force incore for all but a "normal" 244 contraction for now */
        if(!Ztrans || sum_Y == 0 || sum_Y == 1 || sum_Y == 3) incore = 1;

        if(incore) {
            /*       dpd_buf4_scm(Z, beta); */
            buf4_mat_irrep_init(Z, hzbuf);
            if(fabs(beta) > 0.0) buf4_mat_irrep_rd(Z, hzbuf);
            if(Ztrans) {
                trans4_mat_irrep_init(&Zt, hzbuf);
                trans4_mat_irrep_rd(&Zt, hzbuf);
                buf4_mat_irrep_close(Z, hzbuf);
                trans4_mat_irrep_shift13(&Zt, hzbuf);
                numrows = Zt.shift.rowtot[hzbuf];
                numcols = Zt.shift.coltot[hzbuf];
                Zmat = Zt.shift.matrix[hzbuf];
#ifdef DPD_DEBUG
                zrow = Zt.shift.rowtot[hzbuf];
                zcol = Zt.shift.coltot[hzbuf];
#endif
            }
            else {
                buf4_mat_irrep_shift13(Z, hzbuf);
                numrows = Z->shift.rowtot[hzbuf];
                numcols = Z->shift.coltot[hzbuf];
                Zmat = Z->shift.matrix[hzbuf];
#ifdef DPD_DEBUG
                zrow = Z->shift.rowtot[hzbuf];
                zcol = Z->shift.coltot[hzbuf];
#endif
            }

            if (sum_Y < 2) {
                if (Ztrans) hybuf = hzbuf^GY; else hybuf = hzbuf^GZ^GY;
            }
            else {
                if (Ztrans) hybuf = hzbuf; else hybuf = hzbuf^GZ;
            }

            if(sum_Y == 0) {
                buf4_mat_irrep_init(Y, hybuf);
                buf4_mat_irrep_rd(Y, hybuf);
                buf4_mat_irrep_shift13(Y, hybuf);
                Ymat = Y->shift.matrix[hybuf];
                Ytrans = 0;
#ifdef DPD_DEBUG
                yrow = Y->shift.rowtot[hybuf];
                ycol = Y->shift.coltot[hybuf];
#endif
            }
            else if(sum_Y == 1) {
                buf4_mat_irrep_init(Y, hybuf);
                buf4_mat_irrep_rd(Y, hybuf);
                trans4_mat_irrep_init(&Yt, hybuf);
                trans4_mat_irrep_rd(&Yt, hybuf);
                buf4_mat_irrep_close(Y, hybuf);
                trans4_mat_irrep_shift31(&Yt, hybuf);
                rking = 1;
                Ytrans = 1;
                Ymat = Yt.shift.matrix[hybuf];
#ifdef DPD_DEBUG
                yrow = Yt.shift.coltot[hybuf];
                ycol = Yt.shift.rowtot[hybuf];
#endif
            }
            else if(sum_Y == 2) {
                buf4_mat_irrep_init(Y, hybuf);
                buf4_mat_irrep_rd(Y, hybuf);
                trans4_mat_irrep_init(&Yt, hybuf);
                trans4_mat_irrep_rd(&Yt, hybuf);
                buf4_mat_irrep_close(Y, hybuf);
                trans4_mat_irrep_shift13(&Yt, hybuf);
                Ymat = Yt.shift.matrix[hybuf];
                Ytrans = 0;
#ifdef DPD_DEBUG
                yrow = Yt.shift.rowtot[hybuf];
                ycol = Yt.shift.coltot[hybuf];
#endif
            }
            else if(sum_Y == 3) {
                buf4_mat_irrep_init(Y, hybuf);
                buf4_mat_irrep_rd(Y, hybuf);
                buf4_mat_irrep_shift31(Y, hybuf);
                rking = 1;
                Ytrans = 1;
                Ymat = Y->shift.matrix[hybuf];
#ifdef DPD_DEBUG
                yrow = Y->shift.coltot[hybuf];
                ycol = Y->shift.rowtot[hybuf];
#endif
            }

            if(rking)
                for(Hz=0; Hz < nirreps; Hz++) {
                    if      (!Xtrans && !Ytrans) {Hx=Hz;       Hy = Hz^GX; }
                    else if (!Xtrans &&  Ytrans) {Hx=Hz;       Hy = Hz^GX^GY; }
                    else if ( Xtrans && !Ytrans) {Hx=Hz^GX;    Hy = Hz^GX; }
                    else if ( Xtrans &&  Ytrans) {Hx=Hz^GX;    Hy = Hz^GX^GY; }
#ifdef DPD_DEBUG
                    if((xrow[Hz] != zrow[Hz]) || (ycol[Hz] != zcol[Hz]) ||
                            (xcol[Hz] != yrow[Hz])) {
                        outfile->Printf( "** Alignment error in contract244 **\n");
                        outfile->Printf( "** Irrep: %d; Subirrep: %d **\n",hzbuf,Hz);
                        dpd_error("dpd_contract244", "outfile");
                    }
#endif
                    newmm_rking(X->matrix[Hx],Xtrans, Ymat[Hy], Ytrans,
                                Zmat[Hz], numrows[Hz], numlinks[Hx^symlink],
                            numcols[Hz], alpha, 1.0);
                }
            else
                for(Hz=0; Hz < nirreps; Hz++) {
                    if      (!Xtrans && !Ytrans ) {Hx=Hz;       Hy = Hz^GX; }
                    else if (!Xtrans &&  Ytrans ) {Hx=Hz;       Hy = Hz^GX^GY; }
                    else if ( Xtrans && !Ytrans ) {Hx=Hz^GX;    Hy = Hz^GX; }
                    else if ( Xtrans &&  Ytrans ) {Hx=Hz^GX;    Hy = Hz^GX^GY; }

#ifdef DPD_DEBUG
                    if((xrow[Hz] != zrow[Hz]) || (ycol[Hz] != zcol[Hz]) ||
                            (xcol[Hz] != yrow[Hz])) {
                        outfile->Printf( "** Alignment error in contract244 **\n");
                        outfile->Printf( "** Irrep: %d; Subirrep: %d **\n",hzbuf,Hz);
                        dpd_error("dpd_contract244", "outfile");
                    }
#endif
                    /* outfile->Printf("Hz %d, Hx %d, Hy %d, numrows %d, numlinks %d, numcols %d\n",
         Hz, Hx, Hy, numrows[Hz],numlinks[Hx],numcols[Hz]); */

                    if(numrows[Hz] && numcols[Hz] && numlinks[Hx^symlink]) {
                        if(!Xtrans && !Ytrans) {
                            C_DGEMM('n','n',numrows[Hz],numcols[Hz],numlinks[Hx^symlink],
                                    alpha, &(X->matrix[Hx][0][0]),numlinks[Hz^symlink],
                                    &(Ymat[Hy][0][0]),numcols[Hz],1.0,
                                    &(Zmat[Hz][0][0]),numcols[Hz]);
                        }
                        else if(Xtrans && !Ytrans) {
                            C_DGEMM('t','n',numrows[Hz],numcols[Hz],numlinks[Hx^symlink],
                                    alpha, &(X->matrix[Hx][0][0]),numrows[Hz],
                                    &(Ymat[Hy][0][0]),numcols[Hz],1.0,
                                    &(Zmat[Hz][0][0]),numcols[Hz]);
                        }
                        else if(!Xtrans && Ytrans) {
                            C_DGEMM('n','t',numrows[Hz],numcols[Hz],numlinks[Hx^symlink],
                                    alpha, &(X->matrix[Hx][0][0]),numlinks[Hx^symlink],
                                    &(Ymat[Hy][0][0]),numlinks[Hx^symlink],1.0,
                                    &(Zmat[Hz][0][0]),numcols[Hz]);
                        }
                        else {
                            C_DGEMM('t','t',numrows[Hz],numcols[Hz],numlinks[Hx^symlink],
                                    alpha, &(X->matrix[Hx][0][0]),numrows[Hz],
                                    &(Ymat[Hy][0][0]),numlinks[Hx^symlink],1.0,
                                    &(Zmat[Hz][0][0]),numcols[Hz]);
                        }
                    }

                    /*
        newmm(X->matrix[Hx], Xtrans, Ymat[Hy], Ytrans,
        Zmat[Hz], numrows[Hz], numlinks[Hx^symlink],
        numcols[Hz], alpha, 1.0);
      */
                }

            if(sum_Y == 0) buf4_mat_irrep_close(Y, hybuf);
            else if(sum_Y == 1) trans4_mat_irrep_close(&Yt, hybuf);
            else if(sum_Y == 2) trans4_mat_irrep_close(&Yt, hybuf);
            else if(sum_Y == 3) buf4_mat_irrep_close(Y, hybuf);

            if(Ztrans) {
                buf4_mat_irrep_init(Z, hzbuf);
                trans4_mat_irrep_wrt(&Zt, hzbuf);
                trans4_mat_irrep_close(&Zt, hzbuf);
            }

            buf4_mat_irrep_wrt(Z, hzbuf);
            buf4_mat_irrep_close(Z, hzbuf);

        }  /* end if(incore) */
        else { /* out-of-core for "normal" 244 contractions */
            /* Prepare the input buffer for the X factor and the target*/

#ifdef DPD_DEBUG
            outfile->Printf( "\t244 out-of-core: %d\n", hybuf);
#endif
            buf4_mat_irrep_row_init(Y, hybuf);
            buf4_mat_irrep_row_init(Z, hzbuf);

            /* Loop over rows of the Y factor and the target */
            for(pq=0; pq < Z->params->rowtot[hzbuf]; pq++) {

                buf4_mat_irrep_row_zero(Y, hybuf, pq);
                buf4_mat_irrep_row_rd(Y, hybuf, pq);

                buf4_mat_irrep_row_zero(Z, hzbuf, pq);

                if(fabs(beta) > 0.0)
                    buf4_mat_irrep_row_rd(Z, hzbuf, pq);

                for(Gs=0; Gs < nirreps; Gs++) {
                    GrY = Gs^hybuf^GY;
                    GrZ = Gs^hzbuf^GZ;

                    nrows = Z->params->rpi[GrZ];
                    ncols = Z->params->spi[Gs];
                    nlinks = Y->params->rpi[GrY];

                    rowx = Y->params->rpi[GrY];
                    colx = Y->params->spi[Gs];
                    rowz = Z->params->rpi[GrZ];
                    colz = Z->params->spi[Gs];

                    if(nrows && ncols && nlinks) {
                        C_DGEMM(Xtrans?'t':'n','n',nrows,ncols,nlinks,alpha,
                                &(X->matrix[Xtrans?GrY:GrZ][0][0]),Xtrans?nrows:nlinks,
                                &(Y->matrix[hybuf][0][Y->col_offset[hybuf][GrY]]),ncols,1.0,
                                &(Z->matrix[hzbuf][0][Z->col_offset[hzbuf][GrZ]]),ncols);
                    }
                }
                buf4_mat_irrep_row_wrt(Z, hzbuf, pq);
            }

            buf4_mat_irrep_row_close(Y, hybuf);
            buf4_mat_irrep_row_close(Z, hzbuf);
        }
    }

    if(((sum_Y == 1) || (sum_Y == 2)) && incore) trans4_close(&Yt);

    if(Ztrans && incore) trans4_close(&Zt);

    file2_mat_close(X);

    return 0;
}

}
