/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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
#include <vector>

namespace psi {
/*!
** DSYEV_ascending(): diagonalize a symmetric square matrix ('array') using LAPACK DSYEV
**
** \param n      = number of rows (and columns)
** \param array  = matrix to diagonalize (2D row major array)
** \param e_vals = array to hold eigenvalues (returned in ascending order)
** \param e_vecs = (optional) matrix of eigenvectors (2D row major array, one column for each eigvector)
**
** \ingroup CIOMR
*/
[[nodiscard]] int DSYEV_ascending(const int N, const double* const* const array, double* e_vals,
                                  double* const* const e_vecs /* = nullptr*/) {
    // We need to make a copy of the matrix before diagonalization, because LAPACK overwrites it.
    // LAPACK also needs the mtx to be flattened to a 1D array, so a copy is inevitable.
    // The new 1D array will correspond to a column-major array, suitable for LAPACK.
    std::vector<double> tmp_matrix(N * N);
    for (int i = 0, ij = 0; i < N; i++) {
        for (int j = 0; j < N; j++, ij++) {
            tmp_matrix[ij] = array[j][i];
        }
    }
    // LAPACK also needs some extra memory to store temporaries in
    // TODO: query C_DSYEV for optimal workspace size
    const int64_t workspace_size = 3 * N;
    std::vector<double> tmp_work(workspace_size);
    const char jobtype = (e_vecs != nullptr) ? 'V' : 'N';
    const auto info = C_DSYEV(jobtype, 'U', N, tmp_matrix.data(), N, e_vals, tmp_work.data(), workspace_size);
    if ((info == 0) && (e_vecs != nullptr)) {
        // tmp_matrix has now been overwritten with the eigenvecs as the columns, flattened as column-major
        // Copy them to the columns of a row-major 2D array
        for (int j = 0, ij = 0; j < N; j++) {
            for (int i = 0; i < N; i++, ij++) {
                e_vecs[i][j] = tmp_matrix[ij];
            }
        }
    }
    return info;
}
}  // namespace psi
