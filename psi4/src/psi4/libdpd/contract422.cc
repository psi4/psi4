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
#include <cmath>
#include "psi4/libqt/qt.h"
#include "dpd.h"

namespace psi {

/* dpd_contract422(): Contracts a four-index dpd with a two-index
** dpd to give another two-index dpd.  Both indices of the
** two-index factor must be summed, and both must be ket indices in
** the four-index buffer.
**
** Arguments:
**   dpdbuf4 *X: A pointer to the four-index buffer.
**   dpdfile2 *Y: A pointer to the two-index factor file.
**   dpdfile2 *Z: A pointer to the two-index target file.
**   int trans_Y: A boolean to indicate whether the indices in Y are
**                transposed relative to those in the ket of X.
**   int trans_Z: A boolean to indicate whether the indices in Z are
**                transposed relative to those in the bra of X.
**   double alpha: A prefactor for the product alpha * X * Y.
**   double beta: A prefactor for the target beta * Z.
*/

int DPD::contract422(dpdbuf4 *X, dpdfile2 *Y, dpdfile2 *Z, int trans_Y,
                     int trans_Z, double alpha, double beta)
{
    int nirreps, GX, GY, GZ, hxbuf;
    int row,p,q,r,s, psym, qsym, Gr, Gs, P, Q, R, S, col;
    double **TMP;
    double value;
#ifdef DPD_DEBUG
    int *yrow, *ycol, *zrow, *zcol;
#endif

    nirreps = X->params->nirreps;
    GX = X->file.my_irrep;
    GY = Y->my_irrep;
    GZ = Z->my_irrep;

    file2_mat_init(Y);
    file2_mat_rd(Y);
    file2_mat_init(Z);
    if(fabs(beta) > 0.0) file2_mat_rd(Z);

#ifdef DPD_DEBUG
    if(trans_Z) { zrow = Z->params->coltot; zcol = Z->params->rowtot; }
    else { zrow = Z->params->rowtot; zcol = Z->params->coltot; }

    if(trans_Y) { yrow = Y->params->coltot; ycol = Y->params->rowtot; }
    else { yrow = Y->params->rowtot; ycol = Y->params->coltot; }

    if((zrow != X->params->ppi) || (zcol != X->params->qpi) ||
            (yrow != X->params->rpi) || (ycol != X->params->spi)) {
        outfile->Printf( "** Alignment error in contract422 **\n");
        dpd_error("dpd_contract422", "outfile");
    }
#endif

    /* read in block of X whose row irrep is same as target Gpq=GZ */
    hxbuf = GZ;
    buf4_mat_irrep_init(X, hxbuf);
    buf4_mat_irrep_rd(X, hxbuf);

    for(row=0; row < X->params->rowtot[hxbuf]; row++) {
        p = X->params->roworb[hxbuf][row][0];
        psym = X->params->psym[p];
        P = p - X->params->poff[psym];
        q = X->params->roworb[hxbuf][row][1];
        qsym = X->params->qsym[q];
        Q = q - X->params->qoff[qsym];

        value = 0.0;

        for(Gr=0; Gr < nirreps; Gr++) {
            Gs = Gr^GY;

            if(X->params->rpi[Gr] && X->params->spi[Gs]) {
                if(trans_Y) {
                    TMP = dpd_block_matrix(X->params->spi[Gs],X->params->rpi[Gr]);
                }
                else {
                    TMP = dpd_block_matrix(X->params->rpi[Gr],X->params->spi[Gs]);
                }
            }

            for(r=0; r < X->params->rpi[Gr]; r++) {
                R = X->params->roff[Gr] + r;
                for(s=0; s < X->params->spi[Gs]; s++) {
                    S = X->params->soff[Gs] + s;

                    col = X->params->colidx[R][S];

                    if(trans_Y)
                        TMP[s][r] = X->matrix[GZ][row][col];
                    else
                        TMP[r][s] = X->matrix[GZ][row][col];
                }
            }

            if(trans_Y) {
                value += dot_block(TMP, Y->matrix[Gs], X->params->spi[Gs],
                                   X->params->rpi[Gr], alpha);
            } else {
                value += dot_block(TMP, Y->matrix[Gr], X->params->rpi[Gr],
                                   X->params->spi[Gs], alpha);
            }

            if(X->params->rpi[Gr] && X->params->spi[Gs]) {
                if(trans_Y) {
                    free_dpd_block(TMP, X->params->spi[Gs], X->params->rpi[Gr]);
                } else {
                    free_dpd_block(TMP,X->params->rpi[Gr],X->params->spi[Gs]);
                }
            }
        }

        if(trans_Z)
            Z->matrix[qsym][Q][P] = beta*Z->matrix[qsym][Q][P] + value;
        else
            Z->matrix[psym][P][Q] = beta*Z->matrix[psym][P][Q] + value;

    }

    buf4_mat_irrep_close(X, GZ);

    file2_mat_close(Y);
    file2_mat_wrt(Z);
    file2_mat_close(Z);

    return 0;
}

}
