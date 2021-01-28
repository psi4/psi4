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

#include "sparse.h"

#include <tuple>
#include <algorithm>

namespace psi {

/* Args: sorted lists l1 and l2
 * Returns: sorted union of l1 and l2
 */
std::vector<int> merge_lists(const std::vector<int> &l1, const std::vector<int> &l2) {

    std::vector<int> l12;

    int i1 = 0, i2 = 0;
    while(i1 < l1.size() || i2 < l2.size()) {
        if(i1 == l1.size()) {
            l12.push_back(l2[i2]);
            ++i2;
        } else if(i2 == l2.size()) {
            l12.push_back(l1[i1]);
            ++i1;
        } else if(l1[i1] == l2[i2]) {
            l12.push_back(l1[i1]);
            ++i1;
            ++i2;
        } else if(l1[i1] < l2[i2]) {
            l12.push_back(l1[i1]);
            ++i1;
        } else {
            l12.push_back(l2[i2]);
            ++i2;
        }
    }

    return l12;

}

/* Args: sorted list of y, sparse map from A to y (assume sorted)
 * Returns: sorted list yA
 *
 * For all a, every value in A_to_y[a] is included in yA if at least one is present in y
 */
std::vector<int> contract_lists(const std::vector<int> &y, const std::vector<std::vector<int>> &A_to_y) {

    // TODO: runtime is proportional to A_to_y size (system size, O(N))
    // could maybe reduce to &y size (domain size, O(1)), probably doesn't matter
    std::vector<int> yA;

    for(int a = 0, y_ind = 0; a < A_to_y.size(); ++a) {

        bool is_a = false;
        for(auto y_val : A_to_y[a]) {
            if (y_ind < y.size() && y[y_ind] == y_val) {
                y_ind++;
                is_a = true;
            }
        }

        if(is_a) {
            for(auto y_val : A_to_y[a]) {
                yA.push_back(y_val);
            }
        }

    }

    return yA;

}

/* Args: x is a list of values (sorted), y is a map from values of x to values of y
 * Returns: a list of y values
 *
 * Multiple values in x may map to the same value in y (i.e. x is a list of bf, y is atoms)
 */
std::vector<int> block_list(const std::vector<int> &x_list, const std::vector<int> &x_to_y_map) {

    std::vector<int> y_list;

    for(int x_val : x_list) {
        int y_val = x_to_y_map[x_val];
        if(y_list.size() == 0) {
            y_list.push_back(y_val);
        } else if(y_list[y_list.size() - 1] != y_val) {
            y_list.push_back(y_val);
        }
    }

    return y_list;

}

/* Args: SparseMap from x to y, maximum possible y value
 * Returns: SparseMap from y to x
 */
std::vector<std::vector<int>> invert_map(const std::vector<std::vector<int>> &x_to_y, int ny) {

    int nx = x_to_y.size();
    std::vector<std::vector<int>> y_to_x(ny);

    for(int x = 0; x < nx; x++) {
        for(auto y : x_to_y[x]) {
            y_to_x[y].push_back(x);
        }
    }

    return y_to_x;

}

/* Args: SparseMap from x to y, SparseMap from y to z
 * Returns: SparseMap from x to z
 */
std::vector<std::vector<int>> chain_maps(const std::vector<std::vector<int>> &x_to_y, const std::vector<std::vector<int>> &y_to_z) {

    int nx = x_to_y.size();
    std::vector<std::vector<int>> x_to_z(nx);

    for(int x = 0; x < nx; x++) {
        for(auto y : x_to_y[x]) {
            for(auto z : y_to_z[y]) {
                x_to_z[x].push_back(z);
            }
        }
        std::sort(x_to_z[x].begin(), x_to_z[x].end());
        x_to_z[x].erase(std::unique(x_to_z[x].begin(), x_to_z[x].end()), x_to_z[x].end());
    }

    return x_to_z;

}

/* Args: SparseMap from x to y, list of pairs of type x
 * Returns: extended SparseMap from x to y
 */
std::vector<std::vector<int>> extend_maps(const std::vector<std::vector<int>> &x_to_y, const std::vector<std::pair<int,int>> &xpairs) {

    int nx = x_to_y.size();
    std::vector<std::vector<int>> xext_to_y(nx);

    for(auto xpair : xpairs) {
        size_t x1, x2;
        std::tie(x1,x2) = xpair;
        xext_to_y[x1] = merge_lists(xext_to_y[x1], x_to_y[x2]);
    }

    return xext_to_y;

}

SharedMatrix submatrix_rows(SharedMatrix mat, const std::vector<int> &row_inds) {

    SharedMatrix mat_new = std::make_shared<Matrix>(mat->name(), row_inds.size(), mat->colspi(0));
    for(int r_new = 0; r_new < row_inds.size(); r_new++) {
        int r_old = row_inds[r_new];
        for(int c = 0; c < mat->colspi(0); c++) {
            mat_new->set(r_new, c, mat->get(r_old, c));
        }
    }
    return mat_new;
}

SharedMatrix submatrix_cols(SharedMatrix mat, const std::vector<int> &col_inds) {

    SharedMatrix mat_new = std::make_shared<Matrix>(mat->name(), mat->rowspi(0), col_inds.size());
    for(int r = 0; r < mat->rowspi(0); r++) {
        for(int c_new = 0; c_new < col_inds.size(); c_new++) {
            int c_old = col_inds[c_new];
            mat_new->set(r, c_new, mat->get(r, c_old));
        }
    }
    return mat_new;
}

SharedMatrix submatrix_rows_and_cols(SharedMatrix mat, const std::vector<int> &row_inds, const std::vector<int> &col_inds) {

    SharedMatrix mat_new = std::make_shared<Matrix>(mat->name(), row_inds.size(), col_inds.size());
    for(int r_new = 0; r_new < row_inds.size(); r_new++) {
        int r_old = row_inds[r_new];
        for(int c_new = 0; c_new < col_inds.size(); c_new++) {
            int c_old = col_inds[c_new];
            mat_new->set(r_new, c_new, mat->get(r_old, c_old));
        }
    }
    return mat_new;
}

} // namespace psi

