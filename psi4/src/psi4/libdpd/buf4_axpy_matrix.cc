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

#include "dpd.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libqt/qt.h"

namespace psi {

int dpdbuf4::axpy_matrix(const Matrix& MatX, double alpha) {
    if (params->nirreps != MatX.nirrep()) {
        throw PSIEXCEPTION("dpdbuf4 and Matrix have different numbers of irreps");
    }
    if (file.my_irrep != MatX.symmetry()) {
        throw PSIEXCEPTION("dpdbuf4 and Matrix have different symmetries");
    }
    auto symmetry = MatX.symmetry();

    for (int h = 0; h < params->nirreps; ++h) {
        global_dpd_->buf4_mat_irrep_init(this, h);
        global_dpd_->buf4_mat_irrep_rd(this, h);
        size_t sizeX = MatX.rowspi()[h] * MatX.colspi()[h ^ symmetry];
        size_t sizeY = params->rowtot[h] * params->coltot[h ^ symmetry];
        if (sizeX != sizeY) {
            throw PSIEXCEPTION("dpdbuf4 and Matrix have different size");
        }
        if (sizeX) {
            auto X = MatX.pointer(h)[0]; 
            auto Y = &(matrix[h][0][0]);
            C_DAXPY(sizeX, alpha, X, 1, Y, 1);
        }
        global_dpd_->buf4_mat_irrep_wrt(this, h);
        global_dpd_->buf4_mat_irrep_close(this, h);
    }

    return 0;
}

}
