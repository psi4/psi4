/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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
#include "psi4/libqt/qt.h"
#include "dpd.h"

namespace psi {

int DPD::trans4_mat_irrep_wrt(dpdtrans4 *Trans, int irrep) {
    int pq, rs, all_buf_irrep;
    dpdbuf4 *Buf;

    Buf = &(Trans->buf);
    all_buf_irrep = Buf->file.my_irrep;

#ifdef DPD_TIMER
    timer_on("trans4_rw");
#endif

    /* Loop over rows of transpose */
    for (pq = 0; pq < Trans->buf.params->coltot[irrep ^ all_buf_irrep]; pq++) {
        for (rs = 0; rs < Trans->buf.params->rowtot[irrep]; rs++) {
            Buf->matrix[irrep][rs][pq] = Trans->matrix[irrep][pq][rs];
        }
    }

#ifdef DPD_TIMER
    timer_off("trans4_rw");
#endif

    return 0;
}

}  // namespace psi
