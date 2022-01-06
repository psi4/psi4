/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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
#include <cstdlib>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "dpd.h"

namespace psi {

/* buf_block is the block of buffer to shift (also the row-irrep) */

int DPD::trans4_mat_irrep_shift13(dpdtrans4 *Trans, int buf_block) {
    int h, i, nirreps, all_buf_irrep;
    int *count;
    int *dataoff;
    int rowtot, coltot;
    double *data;

#ifdef DPD_TIMER
    timer_on("shift");
#endif

    all_buf_irrep = Trans->buf.file.my_irrep;

    if (Trans->shift.shift_type) {
        outfile->Printf("\n\tShift is already on! %d\n", Trans->shift.shift_type);
        exit(PSI_RETURN_FAILURE);
    } else
        Trans->shift.shift_type = 13;

    nirreps = Trans->buf.params->nirreps;
    rowtot = Trans->buf.params->coltot[buf_block ^ all_buf_irrep];
    coltot = Trans->buf.params->rowtot[buf_block];
    if (rowtot == 0 || coltot == 0)
        data = nullptr;
    else
        data = Trans->matrix[buf_block][0];

    /* Calculate row and column dimensions of each new sub-block */
    for (h = 0; h < nirreps; h++) {
        Trans->shift.rowtot[buf_block][h] = Trans->buf.params->rpi[h];
        Trans->shift.coltot[buf_block][h] = coltot * Trans->buf.params->spi[h ^ buf_block ^ all_buf_irrep];
    }

    /* Malloc the pointers to the rows for the shifted access matrix */
    Trans->shift.matrix[buf_block] = (double ***)malloc(nirreps * sizeof(double **));
    for (h = 0; h < nirreps; h++)
        Trans->shift.matrix[buf_block][h] =
            ((!Trans->shift.rowtot[buf_block][h])
                 ? nullptr
                 : (double **)malloc(Trans->shift.rowtot[buf_block][h] * sizeof(double *)));

    /* Calculate the data offset */
    dataoff = init_int_array(nirreps);
    dataoff[0] = 0;
    for (h = 1; h < nirreps; h++)
        dataoff[h] = dataoff[h - 1] +
                     ((long)Trans->shift.rowtot[buf_block][h - 1]) * ((long)Trans->shift.coltot[buf_block][h - 1]);

    /* The row counter for each sub-block */
    count = init_int_array(nirreps);

    /* Loop over irreps of isolated index */
    for (h = 0; h < nirreps; h++) {
        for (i = 0; (i < Trans->shift.rowtot[buf_block][h]) && Trans->shift.coltot[buf_block][h]; i++, count[h]++) {
            Trans->shift.matrix[buf_block][h][count[h]] =
                &(data[dataoff[h] + ((long)Trans->shift.coltot[buf_block][h]) * ((long)i)]);
        }
    }

    free(count);
    free(dataoff);

#ifdef DPD_TIMER
    timer_off("shift");
#endif

    return 0;
}

}  // namespace psi
