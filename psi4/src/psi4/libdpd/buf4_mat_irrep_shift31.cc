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
#include "psi4/psi4-dec.h"
namespace psi {

int DPD::buf4_mat_irrep_shift31(dpdbuf4 *Buf, int buf_block)
{
    int h, pq, Gr, Gs, r, nirreps, all_buf_irrep, h_pqr;
    int rowtot, coltot,cnt;
    int *count;
    int *blocklen, *rowoff;
    double *data;
    long int pqcol;

#ifdef DPD_TIMER
    timer_on("shift");
#endif

    all_buf_irrep = Buf->file.my_irrep;

    if(Buf->shift.shift_type) {
        outfile->Printf( "\n\tShift is already on! %d\n",
                Buf->shift.shift_type);
        exit(PSI_RETURN_FAILURE);
    }
    else Buf->shift.shift_type = 31;

    nirreps = Buf->params->nirreps;
    rowtot = Buf->params->rowtot[buf_block];
    coltot = Buf->params->coltot[buf_block^all_buf_irrep];
    if (rowtot == 0 || coltot == 0) data = 0;
    else data = Buf->matrix[buf_block][0];

    /* Calculate row and column dimensions of each new sub-block */
    /* loop over h_pqr */
    for(h_pqr=0; h_pqr < nirreps; h_pqr++) {
        Buf->shift.rowtot[buf_block][h_pqr] = rowtot * Buf->params->rpi[h_pqr^buf_block];
        Buf->shift.coltot[buf_block][h_pqr] = Buf->params->spi[h_pqr^all_buf_irrep];
    }

    /* Malloc the pointers to the rows for the shifted access matrix */
    Buf->shift.matrix[buf_block] = (double ***) malloc(nirreps*sizeof(double **));
    for(h_pqr=0; h_pqr < nirreps; h_pqr++)
        Buf->shift.matrix[buf_block][h_pqr] =
                ((!Buf->shift.rowtot[buf_block][h_pqr]) ? NULL :
                                                          (double **) malloc(Buf->shift.rowtot[buf_block][h_pqr] * sizeof(double *)));

    /* Calculate the row offsets */
    blocklen = init_int_array(nirreps);
    for(h_pqr=0; h_pqr < nirreps; h_pqr++)
        blocklen[h_pqr] = Buf->params->rpi[h_pqr^buf_block] *
                Buf->params->spi[h_pqr^all_buf_irrep];

    rowoff = init_int_array(nirreps);
    cnt = 0;
    for (h=0;h<nirreps;++h) {  /* loop over Gr */
        h_pqr = buf_block^h;
        rowoff[h_pqr] = cnt;
        cnt += blocklen[h_pqr];
    }

    /* The row counter for each sub-block */
    count = init_int_array(nirreps);

    /* Loop over rows of original DPD matrix */
    for(pq=0; pq < Buf->params->rowtot[buf_block]; pq++) {
        pqcol = ((long) pq) * ((long) coltot);

        /* Loop over irreps of pqr */
        for(h_pqr=0; h_pqr < nirreps; h_pqr++) {
            Gr = h_pqr^buf_block;
            Gs = h_pqr^all_buf_irrep;

            /* Loop over orbitals in Gr */
            for(r=0; (r < Buf->params->rpi[Gr]) && Buf->params->spi[Gs]; r++) {

                /* Re-assign the row pointer */
                Buf->shift.matrix[buf_block][h_pqr][count[h_pqr]] =
                        &(data[pqcol + rowoff[h_pqr] + (r * Buf->params->spi[Gs])]);

                count[h_pqr]++;

            }
        }
    }

    free(count); free(rowoff); free(blocklen);

#ifdef DPD_TIMER
    timer_off("shift");
#endif

    return 0;
}

}
