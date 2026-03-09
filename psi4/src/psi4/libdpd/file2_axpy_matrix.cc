/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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
#include "dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/matrix.h"

namespace psi {

/* file2_axpy(): Evaluates the standard operation a * X + Y -> Y for
 ** dpdfile2's.
 **
 ** Arguments:
 **   dpdfile2 *FileA: A pointer to the leftmost dpdfile2.
 **   dpdfile2 *FileB: A pointer to the rightmost (and target) dpdfile2.
 **   double alpha: The scalar prefactor in the multiplication.
 **   int transA: A boolean indicating that we should use the transpose of
 **               FileA
 */

int dpdfile2::axpy_matrix(const Matrix& MatX, double alpha) {
    if (params->nirreps != MatX.nirrep()) {
        throw PSIEXCEPTION("dpdbuf4 and Matrix have different numbers of irreps");
    }
    if (my_irrep != MatX.symmetry()) {
        throw PSIEXCEPTION("dpdbuf4 and Matrix have different symmetries");
    }

    global_dpd_->file2_mat_init(this);
    global_dpd_->file2_mat_rd(this);

    for (int h = 0; h < params->nirreps; h++) {
        size_t sizeX = MatX.rowspi()[h] * MatX.colspi()[h ^ my_irrep];
        size_t sizeY = params->rowtot[h] * params->coltot[h ^ my_irrep];
        if (sizeX != sizeY) {
            throw PSIEXCEPTION("dpdfile2 and Matrix have different size");
        }
        if (sizeX) {
            auto X = MatX.pointer(h)[0]; 
            auto Y = &(matrix[h][0][0]);
            C_DAXPY(sizeX, alpha, X, 1, Y, 1);
        }
    }

    global_dpd_->file2_mat_wrt(this);
    global_dpd_->file2_mat_close(this);

    return 0;
}

}  // namespace psi
