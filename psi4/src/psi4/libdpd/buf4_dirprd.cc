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
#include "psi4/libqt/qt.h"
#include "dpd.h"

namespace psi {

/* dpd_buf4_dirprd(): Computes the direct product between two dpd four-index
** buffers.
**
** Arguments:
**   dpdbuf4 *BufA, *BufB: Pointers to the dpd four-index buffers.
**  The results is written to FileB.
*/

int DPD::buf4_dirprd(dpdbuf4 *BufA, dpdbuf4 *BufB) {
    int h, nirreps, my_irrep;

    nirreps = BufA->params->nirreps;
    my_irrep = BufA->file.my_irrep;

    for (h = 0; h < nirreps; h++) {
        buf4_mat_irrep_init(BufA, h);
        buf4_mat_irrep_init(BufB, h);
        buf4_mat_irrep_rd(BufA, h);
        buf4_mat_irrep_rd(BufB, h);

        dirprd_block(BufA->matrix[h], BufB->matrix[h], BufA->params->rowtot[h], BufA->params->coltot[h ^ my_irrep]);

        buf4_mat_irrep_wrt(BufB, h);
        buf4_mat_irrep_close(BufA, h);
        buf4_mat_irrep_close(BufB, h);
    }

    return 0;
}

}  // namespace psi
