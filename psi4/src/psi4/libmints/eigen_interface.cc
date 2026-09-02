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

Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> generate_permutation_to_cca(const BasisSet& basis) {
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> permutation_matrix(basis.nbf());

    bool needs_reorder = basis.has_puream();
#if psi4_SHGSHELL_ORDERING == LIBINT_SHGSHELL_ORDERING_STANDARD
    needs_reorder = false;
#elif psi4_SHGSHELL_ORDERING != LIBINT_SHGSHELL_ORDERING_GAUSSIAN
    #error "unknown value of macro psi4_SHGSHELL_ORDERING"
#endif

    // maps an index to the am of the associated basis fn in
    // CCA order
    std::vector<int> cca_integral_order(2*basis.max_am() + 1, 0);

    for (size_t l = 1, idx = 1; l <= basis.max_am(); idx += 2, ++l) {
        cca_integral_order[idx] = l;
        cca_integral_order[idx + 1] = -l;
    }

    for (int ish = 0, ibf = 0; ish != basis.nshell(); ++ish) {
        auto& sh = basis.shell(ish);
        auto am = sh.am();

        auto ibf_base = ibf;
        for (int ishbf = 0; ishbf != 2*am + 1; ++ishbf, ++ibf) {
            auto tgt = needs_reorder ? ibf_base + cca_integral_order[ishbf] + am : ibf;
            permutation_matrix.indices()[ibf] = tgt;
        }
    }

    return permutation_matrix;

}

}  // namespace linalg
}  // namespace psi
