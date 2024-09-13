/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#include "dpd.h"
#include "psi4/libqt/qt.h"

#include <cstring>

namespace psi {

int dpdbuf4::zero() {

    for (int h = 0; h < params->nirreps; ++h) {
        global_dpd_->buf4_mat_irrep_init(this, h);
        global_dpd_->buf4_mat_irrep_rd(this, h);
        size_t size = params->rowtot[h] * params->coltot[h ^ file.my_irrep] * sizeof(double);
        if (size) {
            memset(&(matrix[h][0][0]), 0, size);
        }
        global_dpd_->buf4_mat_irrep_wrt(this, h);
        global_dpd_->buf4_mat_irrep_close(this, h);
    }

    return 0;
}

}
