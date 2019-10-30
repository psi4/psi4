/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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

#include <algorithm>
#include <set>
#include "integraltransform.h"
#include "psi4/libmints/matrix.h"

using namespace psi;

/* libmints's grad_two_center_computer expects libtrans to send  G_mn = 0.5 C_pm C_qn (G_pq + G_qp).
 * How we do this differs depending on the block. For blocks with the same label (e.g., "O", "O"), we transform
 * 0.5 (G_pq + G_qp). For blocks with different labels (e.g., "O", "V"), we pick one block and transform (G_pq + G_qp).
 * We account for the 0.5 by not transforming the block with the reversed labels (e.g., "V", "O").
 * The main advantage of transforming blocks over transforming the entire matrix is that we no longer
 * force the calling code to construct the entire matrix. They normally have blocks instead. */
SharedMatrix IntegralTransform::backtransform_two_index(std::map<std::array<char, 2>, SharedMatrix> opdm_blocks,
                                                        std::string name) {
    auto backtransformed_matrix = std::make_shared<Matrix>(name, nirreps_, sopi_, sopi_);
    // We want to skip backtransformation of any blocks whose transpose we already backtransformed, because we
    // backtransformed those two blocks together.
    std::set<std::array<char, 2>> labels_to_skip;
    for (const auto& kv : opdm_blocks) {
        auto submatrix_label = kv.first;
        if (labels_to_skip.find(submatrix_label) != labels_to_skip.end()) continue;
        auto label_first = submatrix_label[0];
        auto C_first = MOCoefficients_[label_first];
        auto label_second = submatrix_label[1];
        auto C_second = MOCoefficients_[label_second];
        auto opdm_submatrix = kv.second->clone();
        auto reverse_label = submatrix_label;
        std::reverse(reverse_label.begin(), reverse_label.end());
        labels_to_skip.insert(reverse_label);
        if (submatrix_label == reverse_label) {
            opdm_submatrix->hermitivitize();
        } else {
            auto transpose_pair = opdm_blocks.find(reverse_label);
            if (transpose_pair != opdm_blocks.end()) {
                opdm_submatrix->add(transpose_pair->second->transpose());
            } else {
                // Assume G Hermitian. G_qp = G_pq => G_pq + G_qp = 2 G_pq
                opdm_submatrix->scale(2);
            }
        }
        backtransformed_matrix->add(linalg::triplet(C_first, opdm_submatrix, C_second, false, false, true));
    }
    return backtransformed_matrix;
}
