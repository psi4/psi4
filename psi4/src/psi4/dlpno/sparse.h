/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
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

#ifndef PSI4_SRC_DLPNO_SPARSE_H_
#define PSI4_SRC_DLPNO_SPARSE_H_

#include "psi4/libmints/matrix.h"

#include <vector>

typedef std::vector<std::vector<int>> SparseMap;

namespace psi{

std::vector<int> merge_lists(const std::vector<int> &l1, const std::vector<int> &l2);
std::vector<int> contract_lists(const std::vector<int> &y, const std::vector<std::vector<int>> &A_to_y);
std::vector<int> block_list(const std::vector<int> &x_list, const std::vector<int> &x_to_y_map);

SparseMap invert_map(const SparseMap &x_to_y, int ny);
SparseMap chain_maps(const SparseMap &x_to_y, const SparseMap &y_to_z);
SparseMap extend_maps(const SparseMap &i_to_y, const std::vector<std::pair<int,int>> &ij_to_i_j);

SharedMatrix submatrix_rows(SharedMatrix mat, const std::vector<int> &row_inds);
SharedMatrix submatrix_cols(SharedMatrix mat, const std::vector<int> &col_inds);
SharedMatrix submatrix_rows_and_cols(SharedMatrix mat, const std::vector<int> &row_inds, const std::vector<int> &col_inds);

} // namespace psi

#endif // PSI4_SRC_DLPNO_SPARSE_H_
