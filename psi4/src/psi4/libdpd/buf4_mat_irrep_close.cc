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
#include "dpd.h"

namespace psi {

/* dpd_buf4_mat_irrep_close(): Releases memory for a matrix for a
** single irrep of a dpd four-index buffer.
**
** Arguments:
**   dpdbuf4 *Buf: A pointer to the input dpdbuf.
**   int irrep: The irrep number to be freed.
**
** Note that shift information is freed here as well.
*/

int DPD::buf4_mat_irrep_close(dpdbuf4 *Buf, int irrep)
{
    int h, nirreps, rowtot, coltot, my_irrep;
    long int size;

    my_irrep = Buf->file.my_irrep;
    rowtot = Buf->params->rowtot[irrep];
    coltot = Buf->params->coltot[irrep^my_irrep];

    size = ((long) rowtot) * ((long) coltot);

    nirreps = Buf->params->nirreps;

    /* Free the shift structure for this irrep if used */
    if(Buf->shift.shift_type) {
        for(h=0; h < nirreps; h++)
            if(Buf->shift.rowtot[irrep][h])
                free(Buf->shift.matrix[irrep][h]);
        free(Buf->shift.matrix[irrep]);
        Buf->shift.shift_type = 0;
    }

    if(size) {
        /* If the file member is already in cache and its ordering is the
         same as the buffer, then we just copied the pointer in
         buf4_mat_irrep_init(); don't free! */
        //      if(Buf->file.incore && !(Buf->anti) &&
        //          (Buf->params->pqnum == Buf->file.params->pqnum) &&
        //          (Buf->params->rsnum == Buf->file.params->rsnum))
        //          0;
        //      else
        //          dpd_free_block(Buf->matrix[irrep], rowtot, coltot);
        if(!(Buf->file.incore && !(Buf->anti) &&
             (Buf->params->pqnum == Buf->file.params->pqnum) &&
             (Buf->params->rsnum == Buf->file.params->rsnum)))
            free_dpd_block(Buf->matrix[irrep], rowtot, coltot);
    }

    return 0;
}

}
