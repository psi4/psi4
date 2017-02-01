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
#include "psi4/libciomr/libciomr.h"
#include "dpd.h"

namespace psi {

/* dpd_file4_mat_irrep_close(): Releases memory for a matrix for a
** single irrep of a dpd four-index file.
**
** Arguments:
**   dpdfile4 *File: A pointer to the input dpdfile.
**   int irrep: The irrep number to be freed.
*/

int DPD::file4_mat_irrep_close(dpdfile4 *File, int irrep)
{
    int my_irrep, rowtot, coltot;
    long int size;

    my_irrep = File->my_irrep;

    rowtot = File->params->rowtot[irrep];
    coltot = File->params->coltot[irrep^my_irrep];
    size = ((long) rowtot) * ((long) coltot);

    if(File->incore) return 0;  /* We need to keep the memory */

    if(size) free_dpd_block(File->matrix[irrep], rowtot, coltot);

    return 0;
}

}
