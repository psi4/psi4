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
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "dpd.h"

namespace psi {

int DPD::trans4_mat_irrep_shift31(dpdtrans4 *Trans, int buf_block)
{
    int h, h_rsp, h_pqr, pq, Gr, Gs, r, nirreps, cnt, all_buf_irrep;
    int rowtot, coltot;
    int *count;
    int *blocklen, *rowoff;
    double *data;
    long int pqcol;

#ifdef DPD_TIMER
    timer_on("shift");
#endif
    if(Trans->shift.shift_type) {
        outfile->Printf( "\n\tShift is already on! %d\n",
                Trans->shift.shift_type);
        exit(PSI_RETURN_FAILURE);
    }
    else Trans->shift.shift_type = 31;

    nirreps = Trans->buf.params->nirreps;
    all_buf_irrep = Trans->buf.file.my_irrep;
    rowtot = Trans->buf.params->coltot[buf_block^all_buf_irrep];
    coltot = Trans->buf.params->rowtot[buf_block];
    if (rowtot == 0 || coltot == 0) data = 0;
    else data = Trans->matrix[buf_block][0];

    /* Calculate row and column dimensions of each new sub-block */
    for(h_rsp=0; h_rsp < nirreps; h_rsp++) {
        Trans->shift.coltot[buf_block][h_rsp] = Trans->buf.params->qpi[h_rsp^all_buf_irrep];
        Trans->shift.rowtot[buf_block][h_rsp] = rowtot
                * Trans->buf.params->ppi[h_rsp^all_buf_irrep^buf_block];
    }

    /* Malloc the pointers to the rows for the shifted access matrix */
    Trans->shift.matrix[buf_block] = (double ***) malloc(nirreps*sizeof(double **));
    for(h=0; h < nirreps; h++)
        Trans->shift.matrix[buf_block][h] =
                ((!Trans->shift.rowtot[buf_block][h]) ? NULL :
                                                        (double **) malloc(Trans->shift.rowtot[buf_block][h] * sizeof(double *)));

    /* Calculate the row offsets */
    blocklen = init_int_array(nirreps);
    for(h_rsp=0; h_rsp < nirreps; h_rsp++)
        blocklen[h_rsp] = Trans->buf.params->ppi[h_rsp^all_buf_irrep^buf_block] *
                Trans->buf.params->qpi[h_rsp^all_buf_irrep];

    rowoff = init_int_array(nirreps);
    cnt = 0;
    for (h=0;h<nirreps;++h) {  /* loop over Gp */
        h_rsp = h^buf_block^all_buf_irrep;
        rowoff[h_rsp] = cnt;
        cnt += blocklen[h_rsp];
    }

    /* The row counter for each sub-block */
    count = init_int_array(nirreps);

    /* Loop over rows of original DPD matrix */
    for(pq=0; pq < Trans->buf.params->coltot[buf_block^all_buf_irrep]; pq++) {

        pqcol = ((long) pq) * ((long) coltot);

        /* Loop over irreps of s */
        for(h_pqr=0; h_pqr < nirreps; h_pqr++) { /* loop over rsp of original dpd*/
            Gs = h_pqr^all_buf_irrep;            /* Gq of original dpd */
            Gr = h_pqr^buf_block^all_buf_irrep;  /* Gp of original dpd */

            /* Loop over orbitals in Gr */
            for(r=0; (r < Trans->buf.params->ppi[Gr]) &&
                Trans->buf.params->qpi[Gs]; r++) {

                /* Re-assign the row pointer */
                Trans->shift.matrix[buf_block][h_pqr][count[h_pqr]] =
                        &(data[pqcol + rowoff[h_pqr] +
                        (r * Trans->buf.params->qpi[Gs])]);
                count[h_pqr]++;

            }
        }
    }

    free(count); free(rowoff); free (blocklen);

#ifdef DPD_TIMER
    timer_off("shift");
#endif

    return 0;
}

}
