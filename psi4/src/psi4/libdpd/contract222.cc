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
#include "psi4/psi4-dec.h"
namespace psi {

int DPD::contract222(dpdfile2 *X, dpdfile2 *Y, dpdfile2 *Z, int target_X,
                     int target_Y, double alpha, double beta)
{
    int h, nirreps, Xtrans, Ytrans, *numlinks;
    int GX, GY, GZ;
    int Hx, Hy, Hz;
    int symlink;

#ifdef DPD_DEBUG
    int *xrow, *xcol, *yrow, *ycol, *zrow, *zcol;
#endif

    nirreps = X->params->nirreps;
    GX = X->my_irrep;
    GY = Y->my_irrep;
    GZ = Z->my_irrep;

    file2_mat_init(X);
    file2_mat_rd(X);
    file2_mat_init(Y);
    file2_mat_rd(Y);
    file2_mat_init(Z);
    if(fabs(beta) > 0.0) file2_mat_rd(Z);

    if(target_X == 0) { Xtrans = 0; numlinks = X->params->coltot; symlink = GX; }
    else if(target_X == 1) { Xtrans = 1; numlinks = X->params->rowtot; symlink = 0; }
    else {
        outfile->Printf( "Junk X index %d in contract222\n", target_X);
        exit(PSI_RETURN_FAILURE);
    }
    if(target_Y == 0) Ytrans = 1;
    else if(target_Y == 1) Ytrans = 0;
    else {
        outfile->Printf( "Junk Y index %d in contract222\n", target_Y);
        exit(PSI_RETURN_FAILURE);
    }

#ifdef DPD_DEBUG
    if(Xtrans) { xrow = X->params->coltot; xcol = X->params->rowtot; }
    else { xrow = X->params->rowtot; xcol = X->params->coltot; }

    if(Ytrans) { yrow = Y->params->coltot; ycol = Y->params->rowtot; }
    else { yrow = Y->params->rowtot; ycol = Y->params->coltot; }

    zrow = Z->params->rowtot; zcol = Z->params->coltot;

    if((zrow != xrow) || (zcol != ycol) || (xcol != yrow)) {
        outfile->Printf( "** Alignment error in contract222 **\n");
        dpd_error("dpd_contract222", "outfile");
    }
#endif

    /* loop over row irreps of X */
    for(Hx=0; Hx < nirreps; Hx++) {
        if      ((!Xtrans)&&(!Ytrans)) {Hy = Hx^GX;    Hz = Hx;    }
        else if ((!Xtrans)&&( Ytrans)) {Hy = Hx^GX^GY; Hz = Hx;    }
        else if (( Xtrans)&&(!Ytrans)) {Hy = Hx;       Hz = Hx^GX; }
        else /*(( Xtrans)&&( Ytrans))*/{Hy = Hx^GY;    Hz = Hx^GX; }

        if(Z->params->rowtot[Hz] &&
                Z->params->coltot[Hz^GZ] &&
                numlinks[Hx^symlink]) {
            C_DGEMM(Xtrans?'t':'n',Ytrans?'t':'n',Z->params->rowtot[Hz],
                    Z->params->coltot[Hz^GZ],numlinks[Hx^symlink],alpha,
                    &(X->matrix[Hx][0][0]),X->params->coltot[Hx^GX],
                    &(Y->matrix[Hy][0][0]),Y->params->coltot[Hy^GY],beta,
                    &(Z->matrix[Hz][0][0]),Z->params->coltot[Hz^GZ]);
        }
        /*
    newmm(X->matrix[Hx], Xtrans, Y->matrix[Hy], Ytrans, Z->matrix[Hz],
        Z->params->rowtot[Hz], numlinks[Hx^symlink], Z->params->coltot[Hz^GZ],
        alpha, beta);
    */
    }

    file2_mat_wrt(Z);
    file2_mat_close(X);
    file2_mat_close(Y);
    file2_mat_close(Z);

    return 0;
}

}
