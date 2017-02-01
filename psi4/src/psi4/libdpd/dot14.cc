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
#include "psi4/libqt/qt.h"
#include "dpd.h"

namespace psi {

/* the dot functions have not been tested for transposed cases */

int DPD::dot14(dpdfile2 *T, dpdbuf4 *I, dpdfile2 *Z,
               int transt, int transz, double alpha, double beta)
{
    int h, Gp, Gq, Gr, Gs, GT, GI, GZ, Tblock, Zblock;
    int p, q, r, s;
    int P, Q, R, S;
    int row, col;
    int nirreps;
    double **X;
    double value;

    nirreps = T->params->nirreps;
    GT = T->my_irrep;
    GI = I->file.my_irrep;
    GZ = Z->my_irrep;

    /* Get the two-index quantities from disk */
    file2_mat_init(T);
    file2_mat_rd(T);
    file2_scm(Z, beta);
    file2_mat_init(Z);
    file2_mat_rd(Z);

#ifdef DPD_TIMER
    timer_on("dot14");
#endif

    /* loop over irreps of buffer I; h = Gpq */
    for(h=0; h < nirreps; h++) {

        buf4_mat_irrep_init(I, h);
        buf4_mat_irrep_rd(I, h);

        /* Loop over irreps of the target */
        for(Gq=0; Gq < nirreps; Gq++) {
            /*Gr = Gq;  Gp = Gs = h^Gq; */
            Gp = h^Gq; Gr = Gq^GZ; Gs = h^Gq^GT;
            if (!transt) Tblock = Gp; else Tblock = Gs;
            if (!transz) Zblock = Gq; else Zblock = Gr;

            /* Allocate space for the X buffer */
            if(T->params->ppi[Gp] && T->params->qpi[Gs])
                X = dpd_block_matrix(T->params->ppi[Gp],T->params->qpi[Gs]);

            /* Loop over orbitals of the target */
            for(q=0; q < Z->params->ppi[Gq]; q++) {
                Q = Z->params->poff[Gq] + q;
                for(r=0; r < Z->params->qpi[Gr]; r++) {
                    R = Z->params->qoff[Gr] + r;

                    /* Loop over orbitals of the two-index term */
                    for(p=0; p < T->params->ppi[Gp]; p++) {
                        P = T->params->poff[Gp] + p;
                        for(s=0; s < T->params->qpi[Gs]; s++) {
                            S = T->params->qoff[Gs] + s;

                            /* Calculate row and column indices in I */
                            if(!transt && !transz) {
                                row = I->params->rowidx[P][Q];
                                col = I->params->colidx[R][S];
                            }
                            else if(transt && !transz) {
                                row = I->params->rowidx[S][Q];
                                col = I->params->colidx[R][P];
                            }
                            else if(!transt && transz) {
                                row = I->params->rowidx[P][R];
                                col = I->params->colidx[Q][S];
                            }
                            else if(transt && transz) {
                                row = I->params->rowidx[S][R];
                                col = I->params->colidx[Q][P];
                            }

                            /* Build the X buffer */
                            X[p][s] = I->matrix[h][row][col];

                        }
                    }

                    value = dot_block(T->matrix[Tblock], X, T->params->ppi[Gp],
                                      T->params->qpi[Gs], alpha);

                    Z->matrix[Zblock][q][r] += value;
                }
            }
            if(T->params->ppi[Gp] && T->params->qpi[Gs])
                free_dpd_block(X,T->params->ppi[Gp],T->params->qpi[Gs]);
        }
        buf4_mat_irrep_close(I, h);
    }

#ifdef DPD_TIMER
    timer_off("dot14");
#endif

    /* Close the two-index quantities */
    file2_mat_close(T);
    file2_mat_wrt(Z);
    file2_mat_close(Z);

    return 0;
}

}
