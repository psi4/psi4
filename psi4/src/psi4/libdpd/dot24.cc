/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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

/*! \file
    \ingroup DPD
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libqt/qt.h"
#include "dpd.h"

namespace psi {

/* non-symmetric states have not been tested for transposed cases */

int DPD::dot24(dpdfile2 *T, dpdbuf4 *I, dpdfile2 *Z, int transt, int transz, double alpha, double beta) {
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
    timer_on("dot24");
#endif

    /* loop over row irreps of buf4, h = Gpq = Grs^GI (symm: Gpq = Grs) */
    for (h = 0; h < nirreps; h++) {
        buf4_mat_irrep_init(I, h);
        buf4_mat_irrep_rd(I, h);

        /* Loop over row irreps of the target Z, GZ = Gpr */
        for (Gp = 0; Gp < nirreps; Gp++) {
            /* Gr = Gp;
       Gq = Gs = h^Gp;
     */
            Gq = h ^ Gp;
            Gr = Gp ^ GZ;
            Gs = h ^ Gp ^ GT;
            if (!transt)
                Tblock = Gq;
            else
                Tblock = Gs;
            if (!transz)
                Zblock = Gp;
            else
                Zblock = Gr;

            /* Allocate space for the X buffer */
            if (T->params->ppi[Gq] && T->params->qpi[Gs]) X = dpd_block_matrix(T->params->ppi[Gq], T->params->qpi[Gs]);

            /* Loop over orbitals of the target */
            for (p = 0; p < Z->params->ppi[Gp]; p++) {
                P = Z->params->poff[Gp] + p;
                for (r = 0; r < Z->params->qpi[Gr]; r++) {
                    R = Z->params->qoff[Gr] + r;

                    /* Loop over orbitals of the two-index term */
                    for (q = 0; q < T->params->ppi[Gq]; q++) {
                        Q = T->params->poff[Gq] + q;
                        for (s = 0; s < T->params->qpi[Gs]; s++) {
                            S = T->params->qoff[Gs] + s;

                            /* Calculate row and column indices in I */
                            if (!transt && !transz) {
                                row = I->params->rowidx[P][Q];
                                col = I->params->colidx[R][S];
                            } else if (transt && !transz) {
                                row = I->params->rowidx[P][S];
                                col = I->params->colidx[R][Q];
                            } else if (!transt && transz) {
                                row = I->params->rowidx[R][Q];
                                col = I->params->colidx[P][S];
                            } else if (transt && transz) {
                                row = I->params->rowidx[R][S];
                                col = I->params->colidx[P][Q];
                            }

                            /* Build the X buffer */
                            X[q][s] = I->matrix[h][row][col];
                        }
                    }

                    value = dot_block(T->matrix[Tblock], X, T->params->ppi[Gq], T->params->qpi[Gs], alpha);

                    Z->matrix[Zblock][p][r] += value;
                }
            }
            if (T->params->ppi[Gq] && T->params->qpi[Gs]) free_dpd_block(X, T->params->ppi[Gq], T->params->qpi[Gs]);
        }
        buf4_mat_irrep_close(I, h);
    }

#ifdef DPD_TIMER
    timer_off("dot24");
#endif

    /* Close the two-index quantities */
    file2_mat_close(T);
    file2_mat_wrt(Z);
    file2_mat_close(Z);

    return 0;
}

}  // namespace psi
