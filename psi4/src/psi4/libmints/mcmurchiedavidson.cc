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
#include <cmath>

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

inline double norm(const Point& A) { return std::sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]); }

void fill_E_matrix(int maxam1, int maxam2, const Point& P, const Point& A, const Point& B, double a, double b,
                   std::vector<double>& Ex, std::vector<double>& Ey, std::vector<double>& Ez) {
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

    double oo2p = 1.0 / (2.0 * p);
    for (int i = 0; i < dim1; ++i) {
        for (int j = 0; j < dim2; ++j) {
            // t = 0 case
            int idx0 = address_3d(i, j, 0, dim2, dim3);
            if (i > 0) {
                int idxi0 = address_3d(i - 1, j, 0, dim2, dim3);
                int idxitp0 = address_3d(i - 1, j, 1, dim2, dim3);
                Ex[idx0] += PA[0] * Ex[idxi0] + Ex[idxitp0];
                Ey[idx0] += PA[1] * Ey[idxi0] + Ey[idxitp0];
                Ez[idx0] += PA[2] * Ez[idxi0] + Ez[idxitp0];
            } else if (j > 0) {
                int idxj0 = address_3d(i, j - 1, 0, dim2, dim3);
                int idxjtp0 = address_3d(i, j - 1, 1, dim2, dim3);
                Ex[idx0] += PB[0] * Ex[idxj0] + Ex[idxjtp0];
                Ey[idx0] += PB[1] * Ey[idxj0] + Ey[idxjtp0];
                Ez[idx0] += PB[2] * Ez[idxj0] + Ez[idxjtp0];
            }
            for (int t = 1; t < i + j + 1; ++t) {
                int idx = address_3d(i, j, t, dim2, dim3);
                if (i > 0) {
                    int idxi = address_3d(i - 1, j, t, dim2, dim3);
                    int idxit = address_3d(i - 1, j, t - 1, dim2, dim3);
                    int idxitp = address_3d(i - 1, j, t + 1, dim2, dim3);
                    Ex[idx] += PA[0] * Ex[idxi] + oo2p * Ex[idxit] + (t + 1) * Ex[idxitp];
                    Ey[idx] += PA[1] * Ey[idxi] + oo2p * Ey[idxit] + (t + 1) * Ey[idxitp];
                    Ez[idx] += PA[2] * Ez[idxi] + oo2p * Ez[idxit] + (t + 1) * Ez[idxitp];
                } else if (j > 0) {
                    int idxj = address_3d(i, j - 1, t, dim2, dim3);
                    int idxjt = address_3d(i, j - 1, t - 1, dim2, dim3);
                    int idxjtp = address_3d(i, j - 1, t + 1, dim2, dim3);
                    Ex[idx] += PB[0] * Ex[idxj] + oo2p * Ex[idxjt] + (t + 1) * Ex[idxjtp];
                    Ey[idx] += PB[1] * Ey[idxj] + oo2p * Ey[idxjt] + (t + 1) * Ey[idxjtp];
                    Ez[idx] += PB[2] * Ez[idxj] + oo2p * Ez[idxjt] + (t + 1) * Ez[idxjtp];
                }
            }
        }
    }
}

}  // namespace mdintegrals