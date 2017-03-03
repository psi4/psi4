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

namespace psi {

/* dpd_buf4_symm(): Symmetrizes a dpdbuf4 by taking,
  ** I'(pq,rs) = 1/2 [I(pq,rs) + I(rs,pq)].  Users should keep in mind
  ** that the original buffer will be overwritten when this function is
  ** called.  Also note that this routine will NOT check to see if the
  ** row and column dimensions of the input buffer are identical, which
  ** is necessary for this to work.
  **
  ** Arguments:
  **   dpdbuf4 *Buf: A pointer to the dpdbuf4 to be symmetrized.
  */

int DPD::buf4_symm(dpdbuf4 *Buf)
{
    int h, row, col, all_buf_irrep;
    double value;

    all_buf_irrep = Buf->file.my_irrep;

    for(h=0; h < Buf->params->nirreps; h++) {
        buf4_mat_irrep_init(Buf, h);
        buf4_mat_irrep_rd(Buf, h);

        for(row=0; row < Buf->params->rowtot[h]; row++)
            for(col=0; col < Buf->params->coltot[h^all_buf_irrep]; col++) {
                value = 0.5*(Buf->matrix[h][row][col]+Buf->matrix[h][col][row]);
                Buf->matrix[h][row][col] = Buf->matrix[h][col][row] = value;
            }

        buf4_mat_irrep_wrt(Buf, h);
        buf4_mat_irrep_close(Buf, h);
    }

    return 0;
}


}
