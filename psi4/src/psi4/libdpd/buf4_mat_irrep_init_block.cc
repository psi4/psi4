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
#include "psi4/libciomr/libciomr.h"
#include "dpd.h"

namespace psi {

/* dpd_buf4_mat_irrep_init_block(): Allocates and initializes memory
** for a subblock of a matrix for a single irrep of a dpd four-index
** buffer.
**
** Arguments:
**   dpdbuf4 *Buf: A pointer to the input dpdbuf.
**   int irrep: The irrep number to be prepared.
**   int num_pq: The number of rows needed.
**
*/

int DPD::buf4_mat_irrep_init_block(dpdbuf4 *Buf, int irrep, int num_pq) {
    int all_buf_irrep;
    all_buf_irrep = Buf->file.my_irrep;

    Buf->matrix[irrep] = dpd_block_matrix(num_pq, Buf->params->coltot[irrep ^ all_buf_irrep]);

    return 0;
}

}  // namespace psi
