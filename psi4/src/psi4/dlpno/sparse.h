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

#ifndef PSI4_SRC_DLPNO_SPARSE_H_
#define PSI4_SRC_DLPNO_SPARSE_H_

#include "psi4/libmints/matrix.h"

#include <vector>

typedef std::vector<std::vector<int>> SparseMap;

namespace psi{

/* Args: sorted lists l1 and l2
 * Returns: sorted union of l1 and l2
 */
std::vector<int> merge_lists(const std::vector<int> &l1, const std::vector<int> &l2);

/* Args: sorted list of y, sparse map from A to another list of y (assume sorted, each possible value of y appears exactly once in entire map)
 * Returns: the union of lists in A_to_y where at least one element is in y
 */
std::vector<int> contract_lists(const std::vector<int> &y, const std::vector<std::vector<int>> &A_to_y);

/* Args: x is a list of values (sorted), y is a map from values of x to values of y
 * Returns: a list of y values
 *
 * Multiple values in x may map to the same value in y (i.e. x is a list of bf, y is atoms)
 */
std::vector<int> block_list(const std::vector<int> &x_list, const std::vector<int> &x_to_y_map);

/* Args: SparseMap from x to y, maximum possible y value
 * Returns: SparseMap from y to x
 */
SparseMap invert_map(const SparseMap &x_to_y, int ny);

/* Args: SparseMap from x to y, SparseMap from y to z
 * Returns: SparseMap from x to z
 */
SparseMap chain_maps(const SparseMap &x_to_y, const SparseMap &y_to_z);

/* Args: SparseMap from x to y, list of pairs of type x
 * Returns: extended SparseMap from x to y
 */
SparseMap extend_maps(const SparseMap &i_to_y, const std::vector<std::pair<int,int>> &ij_to_i_j);

/* Args: Matrix, list of row indices */
SharedMatrix submatrix_rows(const Matrix &mat, const std::vector<int> &row_inds);

/* Args: Matrix, list of column indices */
SharedMatrix submatrix_cols(const Matrix &mat, const std::vector<int> &col_inds);

/* Args: Matrix, list of row and column indices */
SharedMatrix submatrix_rows_and_cols(const Matrix &mat, const std::vector<int> &row_inds, const std::vector<int> &col_inds);

} // namespace psi

#endif // PSI4_SRC_DLPNO_SPARSE_H_
