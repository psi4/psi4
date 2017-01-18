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
#include "psi4/libqt/qt.h"
#include "dpd.h"

namespace psi {

/* dpd_buf4_mat_irrep_init(): Allocates and initializes memory for a
** matrix for a single irrep of a dpd four-index buffer.
**
** Arguments:
**   dpdbuf4 *Buf: A pointer to the input dpdbuf.
**   int irrep: The irrep number to be prepared.
*/

int DPD::buf4_mat_irrep_init(dpdbuf4 *Buf, int irrep)
{
    int rowtot, coltot, all_buf_irrep;
    long int size;

    all_buf_irrep = Buf->file.my_irrep;
    rowtot = Buf->params->rowtot[irrep];
    coltot = Buf->params->coltot[irrep^all_buf_irrep];

    size = ((long) rowtot) * ((long) coltot);

#ifdef DPD_TIMER
    timer_on("buf4_init");
#endif

    if(size) {

        /* If the file member is already in cache and its ordering is the
       same as the parent buffer, don't malloc() memory, just assign
       the pointer */
        if(Buf->file.incore && !(Buf->anti) &&
                (Buf->params->pqnum == Buf->file.params->pqnum) &&
                (Buf->params->rsnum == Buf->file.params->rsnum))
            Buf->matrix[irrep] = Buf->file.matrix[irrep];
        else {
            Buf->matrix[irrep] = dpd_block_matrix(rowtot,coltot);
        }

    }

#ifdef DPD_TIMER
    timer_off("buf4_init");
#endif

    return 0;

}

}
