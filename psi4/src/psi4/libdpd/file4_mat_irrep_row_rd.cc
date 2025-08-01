/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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
#include "psi4/libpsio/psio.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libqt/qt.h"
#include "dpd.h"

namespace psi {

int DPD::file4_mat_irrep_row_rd(dpdfile4 *File, int irrep, int row) {
    int coltot, my_irrep, seek_block;
    psio_address row_ptr, next_address;

    if (File->incore) return 0; /* We already have this data in core */

#ifdef DPD_TIMER
    timer_on("f4_rowrd");
#endif

    my_irrep = File->my_irrep;

    row_ptr = File->lfiles[irrep];
    coltot = File->params->coltot[irrep ^ my_irrep];

    /* Advance file pointer to current row --- careful about overflows! */
    if (coltot) {
        seek_block = DPD_BIGNUM / (coltot * sizeof(double)); /* no. of rows for which we can compute the address */
        if (seek_block < 1) {
            outfile->Printf("\nLIBDPD Error: each row of %s is too long to compute an address.\n", File->label);
            dpd_error("dpd_file4_mat_irrep_row_rd", "outfile");
        }
        for (; row > seek_block; row -= seek_block)
            row_ptr = psio_get_address(row_ptr, sizeof(double) * seek_block * coltot);
        row_ptr = psio_get_address(row_ptr, sizeof(double) * row * coltot);
    }

    if (coltot)
        psio_read(File->filenum, File->label, (char *)File->matrix[irrep][0], coltot * sizeof(double), row_ptr,
                  &next_address);

#ifdef DPD_TIMER
    timer_off("f4_rowrd");
#endif

    return 0;
}

}  // namespace psi
