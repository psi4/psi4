/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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
#pragma once

#include <array>
#include <vector>
#include <cmath>
#include <memory>

namespace libint2 {
template <typename Real>
class FmEval_Chebyshev7;
}

namespace mdintegrals {

using Point = std::array<double, 3>;
inline Point point_diff(const Point& A, const Point& B) { return {A[0] - B[0], A[1] - B[1], A[2] - B[2]}; }
inline double point_norm(const Point& A) { return std::sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]); }

void fill_E_matrix(int maxam1, int maxam2, const Point& P, const Point& A, const Point& B, double a, double b,
                   std::vector<double>& Ex, std::vector<double>& Ey, std::vector<double>& Ez);
void fill_M_matrix(int maxam, int maxpow, const Point& PC, double a, double b, std::vector<double>& Mx,
                   std::vector<double>& My, std::vector<double>& Mz);
void fill_R_matrix(int maxam, double p, const Point& P, const Point& C, std::vector<double>& R,
                   std::shared_ptr<const libint2::FmEval_Chebyshev7<double>> fm_eval);

std::vector<std::array<int, 3>> generate_am_components_cca(int am);

inline int cumulative_cart_dim(int L) { return ((L + 1) * (L + 2) * (L + 3)) / 6; }

inline int address_3d(int i, int j, int k, int dim2, int dim3) { return k + dim3 * (j + dim2 * i); }

class MDHelper {
   protected:
    int maxam1_;
    int maxam2_;
    std::vector<double> Ex;
    std::vector<double> Ey;
    std::vector<double> Ez;
    std::vector<std::vector<std::array<int, 3>>> am_comps_;

   public:
    MDHelper(int maxam1, int maxam2) : maxam1_(maxam1), maxam2_(maxam2) {
        int max_am = std::max(maxam1_, maxam2_);
        // CCA-ordered Cartesian components (lx, ly, lz) for all required angular momenta
        am_comps_ = std::vector<std::vector<std::array<int, 3>>>(max_am + 1);
        for (int am = 0; am < max_am + 1; ++am) {
            am_comps_[am] = generate_am_components_cca(am);
        }
        // maximum dimensions of the E matrix
        int edim1 = maxam1_ + 1;
        int edim2 = maxam2_ + 2;
        int edim3 = edim1 + edim2;
        int esize = edim1 * edim2 * edim3;
        // pre-allocate E-matrix
        Ex = std::vector<double>(esize);
        Ey = std::vector<double>(esize);
        Ez = std::vector<double>(esize);
    };
};

}  // namespace mdintegrals