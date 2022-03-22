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
#include "psi4/libmints/mcmurchiedavidson.h"

#include <libint2/boys.h>

namespace mdintegrals {

inline int cart_dim(int L) { return (L + 1) * (L + 2) / 2; }

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

void fill_M_matrix(int maxam, int maxpow, const Point& PC, double a, double b, std::vector<double>& Mx,
                   std::vector<double>& My, std::vector<double>& Mz) {
    // Generate multipole intermediates using eqs 9.5.31 to 9.5.36
    // from Molecular Electronic-Structure Theory (10.1002/9781119019572)

    // zero out buffers
    std::fill(Mx.begin(), Mx.end(), 0.0);
    std::fill(My.begin(), My.end(), 0.0);
    std::fill(Mz.begin(), Mz.end(), 0.0);

    int dim1 = std::max(maxam, maxpow) + 2;

    double p = a + b;
    // one over 2p
    double oo2p = 1.0 / (2.0 * p);
    double sqrtpip = std::sqrt(M_PI / p);
    // eq 9.5.32: M_t^0 = delta_{0t} \sqrt{\pi/p}
    Mx[0] = sqrtpip;
    My[0] = sqrtpip;
    Mz[0] = sqrtpip;
    for (int e = 1; e <= maxpow; ++e) {
        // t = 0 case
        int idx0 = e * dim1; // M_0^e
        int idx0_em = (e - 1) * dim1; // M_0^{e-1}
        // last two terms of eq 9.5.36
        Mx[idx0] += PC[0] * Mx[idx0_em] + oo2p * Mx[idx0_em + 1];
        My[idx0] += PC[1] * My[idx0_em] + oo2p * My[idx0_em + 1];
        Mz[idx0] += PC[2] * Mz[idx0_em] + oo2p * Mz[idx0_em + 1];
        // t > 0 case
        int upper_t = std::min(e + 1, std::max(maxam, maxpow) + 1);
        for (int t = 1; t < upper_t; ++t) {
            int idx = e * dim1 + t; // target index M_t^e
            int idx_em = (e - 1) * dim1 + t; // M_t^{e-1}
            int idx_em_tm = (e - 1) * dim1 + (t - 1); // M_{t-1}^{e-1}
            int idx_em_tp = (e - 1) * dim1 + (t + 1); // M_{t+1}^{e-1}
            // eq 9.5.36
            Mx[idx] += t * Mx[idx_em_tm] + PC[0] * Mx[idx_em] + oo2p * Mx[idx_em_tp];
            My[idx] += t * My[idx_em_tm] + PC[1] * My[idx_em] + oo2p * My[idx_em_tp];
            Mz[idx] += t * Mz[idx_em_tm] + PC[2] * Mz[idx_em] + oo2p * Mz[idx_em_tp];
        }
    }
}

void fill_R_matrix(int maxam, double p, const Point& P, const Point& C, std::vector<double>& R,
                   std::shared_ptr<const libint2::FmEval_Chebyshev7<double>> fm_eval) {
    // zero out buffer
    std::fill(R.begin(), R.end(), 0.0);
    // eq 9.9.13
    auto PC = point_diff(P, C);
    auto RPC = point_norm(PC);
    double T = p * RPC * RPC;

    std::vector<double> fmvals(maxam + 1);

    // evaluate Boys function
    fm_eval->eval(fmvals.data(), T, maxam);

    int dim1 = maxam + 1;
    int dim2 = dim1 * dim1 * dim1;

    // NOTE: avoiding std::pow(-2.0 * p, n)
    double fac = 1.0;
    double mult = -2.0 * p;
    for (int n = 0; n < dim1; ++n) {
        // eq 9.9.14
        R[n * dim2] = fac * fmvals[n];
        fac *= mult;
    }
    for (int v = 0; v < dim1; ++v) {
        for (int u = 0; u < dim1 - v; ++u) {
            for (int t = 0; t < dim1 - v - u; ++t) {
                for (int n = 0; n < maxam; ++n) {
                    double val = 0.0;
                    int noffset = (n + 1) * dim2;
                    if (t == 0 && u == 0 && v == 0) {
                        continue;
                    } else if (t == 0 && u == 0) {
                        // eq 9.9.20
                        if (v > 1) {
                            val += (v - 1) * R[noffset + address_3d(t, u, v - 2, dim1, dim1)];
                        }
                        val += PC[2] * R[noffset + address_3d(t, u, v - 1, dim1, dim1)];
                    } else if (t == 0) {
                        // eq 9.9.19
                        if (u > 1) {
                            val += (u - 1) * R[noffset + address_3d(t, u - 2, v, dim1, dim1)];
                        }
                        val += PC[1] * R[noffset + address_3d(t, u - 1, v, dim1, dim1)];
                    } else {
                        // eq 9.9.18
                        if (t > 1) {
                            val += (t - 1) * R[noffset + address_3d(t - 2, u, v, dim1, dim1)];
                        }
                        val += PC[0] * R[noffset + address_3d(t - 1, u, v, dim1, dim1)];
                    }
                    R[n * dim2 + address_3d(t, u, v, dim1, dim1)] = val;
                }
            }
        }
    }
}

}  // namespace mdintegrals