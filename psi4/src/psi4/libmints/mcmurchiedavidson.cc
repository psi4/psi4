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
#include "mcmurchiedavidson.h"

namespace mdintegrals {

inline int cart_dim(int L) { return (L + 1) * (L + 2) / 2; }

inline int cumulative_cart_dim(int L) { return ((L + 1) * (L + 2) * (L + 3)) / 6; }

std::vector<std::array<int, 4>> generate_am_components_cca(int am) {
    std::vector<std::array<int, 4>> ret(cart_dim(am));
    int index = 0;
    for (int l = am; l > -1; --l) {
        for (int n = 0; n < am - l + 1; ++n) {
            int m = am - l - n;
            ret[index] = {{l, m, n, index}};
            index++;
        }
    }
    return ret;
}

void fill_E_matrix(int maxam1, int maxam2, const Point& P, const Point& A, const Point& B, double a, double b,
                   std::vector<double>& Ex, std::vector<double>& Ey, std::vector<double>& Ez) {
    // computes the Hermite Gaussian expansion coefficients E_t^{ij} (eq 9.5.1)
    // equation numbers from Molecular Electronic-Structure Theory (10.1002/9781119019572)

    // make sure buffers are zeroed out
    std::fill(Ex.begin(), Ex.end(), 0.0);
    std::fill(Ey.begin(), Ey.end(), 0.0);
    std::fill(Ez.begin(), Ez.end(), 0.0);
    int dim1 = maxam1 + 1;
    int dim2 = maxam2 + 1;
    int dim3 = maxam1 + maxam2 + 2;

    // eq 9.2.11: total exponent
    double p = a + b;
    // eq 9.2.12: reduced exponent
    double mu = a * b / p;
    // eq 9.2.14: relative coordinate
    Point AB = point_diff(A, B);
    Point PA = point_diff(P, A);
    Point PB = point_diff(P, B);
    // eq 9.2.15: pre-exponential factor
    Ex[0] = std::exp(-1.0 * mu * AB[0] * AB[0]);
    Ey[0] = std::exp(-1.0 * mu * AB[1] * AB[1]);
    Ez[0] = std::exp(-1.0 * mu * AB[2] * AB[2]);

    // one over 2p
    double oo2p = 1.0 / (2.0 * p);
    // compute the upward recursion (eqs 9.5.6, 9.5.7)
    for (int i = 0; i < dim1; ++i) {
        for (int j = 0; j < dim2; ++j) {
            // special case: t = 0
            int idx_t0 = address_3d(i, j, 0, dim2, dim3);  // E_0^{ij}
            if (i > 0) {
                int idxi_t0 = address_3d(i - 1, j, 0, dim2, dim3);  // E_0^{i-1,j}
                int idxi_t1 = address_3d(i - 1, j, 1, dim2, dim3);  // E_1^{i-1,j}
                // terms 2 and 3 in eq 9.5.6
                Ex[idx_t0] += PA[0] * Ex[idxi_t0] + Ex[idxi_t1];
                Ey[idx_t0] += PA[1] * Ey[idxi_t0] + Ey[idxi_t1];
                Ez[idx_t0] += PA[2] * Ez[idxi_t0] + Ez[idxi_t1];
            } else if (j > 0) {
                int idxj_t0 = address_3d(i, j - 1, 0, dim2, dim3);  // E_0^{i,j-1}
                int idxj_t1 = address_3d(i, j - 1, 1, dim2, dim3);  // E_1^{i,j-1}
                // terms 2 and 3 in eq 9.5.7
                Ex[idx_t0] += PB[0] * Ex[idxj_t0] + Ex[idxj_t1];
                Ey[idx_t0] += PB[1] * Ey[idxj_t0] + Ey[idxj_t1];
                Ez[idx_t0] += PB[2] * Ez[idxj_t0] + Ez[idxj_t1];
            }
            // t > 0
            for (int t = 1; t < i + j + 1; ++t) {
                int idx = address_3d(i, j, t, dim2, dim3);  // E_t^{ij}
                if (i > 0) {
                    int idxi_t = address_3d(i - 1, j, t, dim2, dim3);       // E_t^{i-1,j}
                    int idxi_tm = address_3d(i - 1, j, t - 1, dim2, dim3);  // E_{t-1}^{i-1,j}
                    int idxi_tp = address_3d(i - 1, j, t + 1, dim2, dim3);  // E_{t+1}^{i-1,j}
                    // eq 9.5.6
                    Ex[idx] += PA[0] * Ex[idxi_t] + oo2p * Ex[idxi_tm] + (t + 1) * Ex[idxi_tp];
                    Ey[idx] += PA[1] * Ey[idxi_t] + oo2p * Ey[idxi_tm] + (t + 1) * Ey[idxi_tp];
                    Ez[idx] += PA[2] * Ez[idxi_t] + oo2p * Ez[idxi_tm] + (t + 1) * Ez[idxi_tp];
                } else if (j > 0) {
                    int idxj_t = address_3d(i, j - 1, t, dim2, dim3);       // E_t^{i,j-1}
                    int idxj_tm = address_3d(i, j - 1, t - 1, dim2, dim3);  // E_{t-1}^{i,j-1}
                    int idxj_tp = address_3d(i, j - 1, t + 1, dim2, dim3);  // E_{t+1}^{i,j-1}
                    // eq 9.5.7
                    Ex[idx] += PB[0] * Ex[idxj_t] + oo2p * Ex[idxj_tm] + (t + 1) * Ex[idxj_tp];
                    Ey[idx] += PB[1] * Ey[idxj_t] + oo2p * Ey[idxj_tm] + (t + 1) * Ey[idxj_tp];
                    Ez[idx] += PB[2] * Ez[idxj_t] + oo2p * Ez[idxj_tm] + (t + 1) * Ez[idxj_tp];
                }
            }
        }
    }
}

}  // namespace mdintegrals