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

/* dpd_trace42_13(): Take a "trace" of the specified indices of a buf4 and put the
** result into a file2.
**
** In this version of the trace42 functions, the summation indices are 1 and 3:
**
**         B(q,s) = alpha * \Sum_{p==r} A(pq,rs) + beta * B(q,s)
**
** Arguments:
**   dpdbuf4 *A: Pointer to the input four-index buffer.
**   dpdfile2 *B: Pointer to the target two-index file.
**   int transb: Boolean to indicate whether B should be transposed to match the
**               target index ordering of A.
**   double alpha: Prefactor for A.
**   double beta: Prefactor for B.
*/

int DPD::trace42_13(dpdbuf4 *A, dpdfile2 *B, int transb, double alpha, double beta)
{
    int h, Gp, Gq, Gr, Gs;
    int p, q, r, s;
    int P, Q, R, S;
    int pq, rs;
    int nirreps;

    nirreps = A->params->nirreps;

    file2_scm(B, beta);
    file2_mat_init(B);
    file2_mat_rd(B);

#ifdef DPD_TIMER
    timer_on("trace42");
#endif

    /* Read all of A into core */
    for(h=0; h < nirreps; h++) {
        buf4_mat_irrep_init(A, h);
        buf4_mat_irrep_rd(A, h);
    }

    for(h=0; h < nirreps; h++) {

        for(Gp=0; Gp < nirreps; Gp++) {
            Gq = Gp ^ h;
            Gr = Gp;
            Gs = Gq;

            /* Loop over target indices */
            for(q=0; q < A->params->qpi[Gq]; q++) {
                Q = A->params->qoff[Gq] + q;

                for(s=0; s < A->params->spi[Gs]; s++) {
                    S = A->params->soff[Gs] + s;

                    /* Loop over elements for which P==R */
                    for(p=0; p < A->params->ppi[Gp]; p++) {
                        P = A->params->poff[Gp] + p;
                        R = P;

                        pq = A->params->rowidx[P][Q];
                        rs = A->params->colidx[R][S];

                        if(!transb)  B->matrix[Gq][q][s] += alpha * A->matrix[h][pq][rs];
                        else  B->matrix[Gq][s][q] += alpha * A->matrix[h][pq][rs];

                    }
                }
            }

        }

    }

    for(h=0; h < nirreps; h++) {
        buf4_mat_irrep_close(A, h);
    }

#ifdef DPD_TIMER
    timer_off("trace42");
#endif

    /* Close the two-index quantities */
    file2_mat_wrt(B);
    file2_mat_close(B);

    return 0;
}

}
