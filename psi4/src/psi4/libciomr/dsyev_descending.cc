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

/*!
** \file
** \brief Diagnoalize a symmetrix square matrix
** \ingroup CIOMR
*/

#include "psi4/libqt/qt.h"
#include "libciomr.h"

namespace psi {
/*!
** DSYEV_descending(): diagonalize a symmetric square matrix ('array') using LAPACK DSYEV, with results reversed
** Please note that descending-order results are obtained by running an ascending-order diagonalizer and reordering
** the results. This may have a small performance penalty for large matrices, especially if eigenvectors are required.
**
** \param n      = number of rows (and columns)
** \param array  = matrix to diagonalize (2D row major array)
** \param e_vals = array to hold eigenvalues (returned in descending order)
** \param e_vecs = (optional) matrix of eigenvectors (2D row major array, one column for each eigvector)
**
** \ingroup CIOMR
*/
[[nodiscard]] int DSYEV_descending(const int N, const double* const* const array, double* e_vals,
                                   double* const* const e_vecs /* = nullptr*/) {
    const auto info = DSYEV_ascending(N, array, e_vals, e_vecs);
    // Reverse the order of eigenvalues
    for (int64_t i = 0; i < N / 2; i++) {
        std::swap(e_vals[i], e_vals[N - i - 1]);
    }
    // Reverse the order of columns of the row-major 2D eigenvector array
    for (int64_t i = 0; i < N; i++) {
        for (int64_t j = 0; j < N / 2; j++) {
            std::swap(e_vecs[i][j], e_vecs[i][N - j - 1]);
        }
    }
    return info;
}
}  // namespace psi
