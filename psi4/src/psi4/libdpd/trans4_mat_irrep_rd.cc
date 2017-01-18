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
#include "psi4/libqt/qt.h"
#include "dpd.h"

namespace psi {

int DPD::trans4_mat_irrep_rd(dpdtrans4 *Trans, int irrep)
{
    int pq, rs, all_buf_irrep;
    int rows, cols;
    dpdbuf4 *Buf;
    double *A, *B;

    Buf = &(Trans->buf);
    all_buf_irrep = Buf->file.my_irrep;

#ifdef DPD_TIMER
    timer_on("trans4_rw");
#endif

    /* Loop over rows of transpose */
    /*
    for(pq=0; pq < Trans->buf.params->coltot[irrep]; pq++) {
    for(rs=0; rs < Trans->buf.params->rowtot[irrep]; rs++) {
    Trans->matrix[irrep][pq][rs] = Buf->matrix[irrep][rs][pq];
    }
    }
  */
    /* Transpose using BLAS DCOPY */
    rows = Buf->params->rowtot[irrep];
    cols = Buf->params->coltot[irrep^all_buf_irrep];

    if(rows && cols) {
        for(rs=0; rs < cols; rs++) {
            A = &(Buf->matrix[irrep][0][rs]);
            B = &(Trans->matrix[irrep][rs][0]);
            C_DCOPY(rows, A, cols, B, 1);
        }
    }

#ifdef DPD_TIMER
    timer_off("trans4_rw");
#endif

    return 0;
}

}
