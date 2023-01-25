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

#include "sparse.h"

#include "psi4/libqt/qt.h"

#include <tuple>
#include <algorithm>

namespace psi {

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

SharedMatrix submatrix_rows(const Matrix &mat, const std::vector<int> &row_inds) {

    SharedMatrix mat_new = std::make_shared<Matrix>(mat.name(), row_inds.size(), mat.colspi(0));
    double** mat_newp = mat_new->pointer();
    double** matp = mat.pointer();
    for(int r_new = 0; r_new < row_inds.size(); r_new++) {
        int r_old = row_inds[r_new];
        ::memcpy(&mat_newp[r_new][0], &matp[r_old][0], sizeof(double) * mat.colspi(0));
    }
    return mat_new;
}

SharedMatrix submatrix_cols(const Matrix &mat, const std::vector<int> &col_inds) {

    SharedMatrix mat_new = std::make_shared<Matrix>(mat.name(), mat.rowspi(0), col_inds.size());
    double** mat_newp = mat_new->pointer();
    double** matp = mat.pointer();
    for(int c_new = 0; c_new < col_inds.size(); c_new++) {
        int c_old = col_inds[c_new];
        C_DCOPY(mat.rowspi(0), &matp[0][c_old], mat.colspi(0), &mat_newp[0][c_new], col_inds.size());
    }
    return mat_new;
}

SharedMatrix submatrix_rows_and_cols(const Matrix &mat, const std::vector<int> &row_inds, const std::vector<int> &col_inds) {

    SharedMatrix mat_new = std::make_shared<Matrix>(mat.name(), row_inds.size(), col_inds.size());
    for(int r_new = 0; r_new < row_inds.size(); r_new++) {
        int r_old = row_inds[r_new];
        for(int c_new = 0; c_new < col_inds.size(); c_new++) {
            int c_old = col_inds[c_new];
            mat_new->set(r_new, c_new, mat.get(r_old, c_old));
        }
    }
    return mat_new;
}

} // namespace psi

