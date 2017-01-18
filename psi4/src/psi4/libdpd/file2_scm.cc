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
#include "dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.h"

namespace psi {

int DPD::file2_scm(dpdfile2 *InFile, double alpha)
{
    int h, nirreps, new_file2, my_irrep;
    int row, col, length;
    double *X;

    nirreps = InFile->params->nirreps;
    my_irrep = InFile->my_irrep;
    file2_mat_init(InFile);

    /* Look first for the TOC entry on disk */
    if(psio_tocscan(InFile->filenum, InFile->label) == NULL)
        new_file2 = 1;
    else new_file2 = 0;

    if(!new_file2) file2_mat_rd(InFile);

    for(h=0; h < nirreps; h++) {

        length = InFile->params->rowtot[h] * InFile->params->coltot[h^my_irrep];
        if(length) {
            X = &(InFile->matrix[h][0][0]);
            C_DSCAL(length, alpha, X, 1);
        }
    }

    file2_mat_wrt(InFile);
    file2_mat_close(InFile);

    return 0;
}

}
