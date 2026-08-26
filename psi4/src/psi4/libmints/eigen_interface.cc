/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#include "psi4/libmints/eigen_interface.h"

namespace psi {
namespace linalg {

Eigen::Map<Eigen::MatrixXd> eigen_map(Matrix& matrix) {
    if (matrix.nirrep() != 1) {
        std::string message = "linalg::eigen_map() called, but matrix only has one irrep! ";
        message += "Use linalg::eigen_maps() instead.";
        throw PSIEXCEPTION(message);
    }

    return Eigen::Map<Eigen::MatrixXd>(matrix.get_pointer(), matrix.nrow(), matrix.ncol());
}

std::vector<Eigen::Map<Eigen::MatrixXd>> eigen_maps(Matrix& matrix) {
    std::vector<Eigen::Map<Eigen::MatrixXd>> maps;
    maps.reserve(matrix.nirrep());

    for (int h = 0; h != matrix.nirrep(); ++h) {
        maps.emplace_back(matrix.get_pointer(h), matrix.rowdim(h), matrix.coldim(h));
    }

    return maps;
}

Matrix matrix_from_eigen(const Eigen::MatrixXd& eigen_mat, const std::string& name) {
    Matrix new_matrix(name, static_cast<int>(eigen_mat.rows()), static_cast<int>(eigen_mat.cols()) );
    new_matrix.copy_from(eigen_mat.data());
    return new_matrix;
}

}  // namespace linalg
}  // namespace psi
