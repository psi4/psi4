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

#include "psi4/libmints/batchedpotential.h"

#include <algorithm>
#include <cmath>

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mcmurchiedavidson.h"
#include "psi4/libmints/snk/snk_kernels.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsi4util/process.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <libint2/boys.h>
#include <libint2/shell.h>

namespace psi {

using namespace mdintegrals;

namespace {

/// Number of doubles the Hermite integrals of one chunk of points may take: chunks are
/// as many points as fit, between 8 and 64.
constexpr size_t R_budget = 65536;

/// Primitive pairs whose Gaussian product prefactor is below this are skipped
constexpr double prim_threshold = 1.0e-15;

/// asymptotic_T[L]: the argument above which F_0 .. F_L follow the upward recursion from the asymptotic
/// F_0 = sqrt(pi / (4T)) to a relative error below 1e-16 (computed in 50-digit arithmetic against the series,
/// plus a margin of one). For higher L, Libint2's FmEval_Chebyshev7::T_crit = 117 holds.
constexpr double asymptotic_T[] = {36, 40, 44, 47, 49, 52, 54, 57, 59, 61, 64,
                                   66, 68, 70, 72, 74, 76, 78, 80, 82, 84};

/*! The Boys function F_L(T) and exp(-T) for 0 <= T <= asymptotic_T[nmax], L <= nmax, as in the preprint of DiCenso,
 *  Laqua and Head-Gordon (arXiv:2610.21075, sec. III B 1): third-order Chebyshev interpolation on intervals of width
 *  1/512, for many arguments at a time. Measured against the series in extended precision, the relative error is
 *  below 4e-15 for F_0 .. F_8, and 1.1e-14 (absolute 5e-15) for exp(-T). The tables are built once, in long
 *  double, on first use.
 */
class BoysChebyshevTable {
   public:
    static constexpr int nmax = 8;
    static constexpr double inv_delta = 512.0;
    static constexpr double Tmax = 59.0;  // asymptotic_T[nmax]
    static constexpr int nint = static_cast<int>(Tmax * inv_delta) + 1;

    static const BoysChebyshevTable& instance() {
        static const BoysChebyshevTable table;
        return table;
    }

    /// F_n(T[k]) into F[n * nk + k], n <= L <= nmax, for the nk arguments T, by interpolating F_L and exp(-T)
    /// and the downward recursion F_n = (2T F_{n+1} + exp(-T)) / (2n + 1); arguments beyond Tmax are clamped,
    /// and the caller discards their values. e is nk doubles of work space.
    void eval(int L, size_t nk, const double* T, double* F, double* e) const {
        const double* cF = &coef_[static_cast<size_t>(L) * nint * 4];
        const double* cE = &coef_[static_cast<size_t>(nmax + 1) * nint * 4];
        double* FL = F + L * nk;
#pragma omp simd
        for (size_t k = 0; k < nk; ++k) {
            const double x = std::min(T[k], Tmax) * inv_delta;
            const int i = std::min(static_cast<int>(x), nint - 1);
            const double t = 2.0 * (x - i) - 1.0;
            const double* c = cF + 4 * i;
            const double* d = cE + 4 * i;
            FL[k] = c[0] + t * (c[1] + t * (c[2] + t * c[3]));
            e[k] = d[0] + t * (d[1] + t * (d[2] + t * d[3]));
        }
        for (int n = L - 1; n >= 0; --n) {
            const double inv = 1.0 / (2 * n + 1);
            const double* up = F + (n + 1) * nk;
            double* cur = F + n * nk;
#pragma omp simd
            for (size_t k = 0; k < nk; ++k) cur[k] = (2.0 * T[k] * up[k] + e[k]) * inv;
        }
    }

   private:
    /// [n][interval][4] power coefficients in t in [-1, 1] of F_n, n <= nmax, then of exp(-T)
    std::vector<double> coef_;

    BoysChebyshevTable() {
        coef_.resize(static_cast<size_t>(nmax + 2) * nint * 4);
        // Chebyshev nodes of degree 3, and the inverse of the Vandermonde matrix in t at them
        long double tn[4], V[4][4], Vinv[4][4];
        for (int j = 0; j < 4; ++j) tn[j] = std::cos((2 * j + 1) * 3.14159265358979323846264338327950288L / 8);
        for (int j = 0; j < 4; ++j)
            for (int m = 0; m < 4; ++m) V[j][m] = m == 0 ? 1.0L : V[j][m - 1] * tn[j];
        invert4(V, Vinv);
        for (int i = 0; i < nint; ++i) {
            long double vals[nmax + 2][4];
            for (int j = 0; j < 4; ++j) {
                const long double x = (i + 0.5L + 0.5L * tn[j]) / inv_delta;
                // F_nmax by the series, all terms positive, then the downward recursion
                const long double ex = std::exp(-x);
                long double term = 1.0L / (2 * nmax + 1), sum = term;
                for (int k = 1; k < 10000; ++k) {
                    term *= 2 * x / (2 * nmax + 2 * k + 1);
                    sum += term;
                    if (term < sum * 1.0e-21L) break;
                }
                long double f = ex * sum;
                vals[nmax][j] = f;
                for (int n = nmax - 1; n >= 0; --n) vals[n][j] = f = (2 * x * f + ex) / (2 * n + 1);
                vals[nmax + 1][j] = ex;
            }
            for (int n = 0; n <= nmax + 1; ++n)
                for (int m = 0; m < 4; ++m) {
                    long double c = 0;
                    for (int j = 0; j < 4; ++j) c += Vinv[m][j] * vals[n][j];
                    coef_[(static_cast<size_t>(n) * nint + i) * 4 + m] = static_cast<double>(c);
                }
        }
    }

    static void invert4(const long double (&M)[4][4], long double (&Minv)[4][4]) {
        long double a[4][8];
        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 8; ++c) a[r][c] = c < 4 ? M[r][c] : (c - 4 == r ? 1.0L : 0.0L);
        for (int c = 0; c < 4; ++c) {
            int piv = c;
            for (int r = c + 1; r < 4; ++r)
                if (std::abs(a[r][c]) > std::abs(a[piv][c])) piv = r;
            for (int k = 0; k < 8; ++k) std::swap(a[c][k], a[piv][k]);
            const long double d = a[c][c];
            for (int k = 0; k < 8; ++k) a[c][k] /= d;
            for (int r = 0; r < 4; ++r)
                if (r != c) {
                    const long double f = a[r][c];
                    for (int k = 0; k < 8; ++k) a[r][k] -= f * a[c][k];
                }
        }
        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 4; ++c) Minv[r][c] = a[r][c + 4];
    }
};

/// Points per chunk for total angular momentum L.
size_t chunk_size(int L) {
    const size_t dim1 = L + 1;
    const size_t per_point = dim1 * dim1 * dim1 * dim1;
    const size_t fit = R_budget / per_point;
    return std::max<size_t>(8, std::min<size_t>(64, fit - fit % 8));
}

}  // namespace

BatchedPotentialInt::BatchedPotentialInt(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2)
    : bs1_(bs1), bs2_(bs2) {
    const int am1 = bs1_->max_am(), am2 = bs2_->max_am();
    max_L_ = am1 + am2;
    if (max_L_ >= 64) throw PSIEXCEPTION("BatchedPotentialInt: shell pairs of total angular momentum 64 or more.");
    fm_eval_ = libint2::FmEval_Chebyshev7<double>::instance(max_L_);
    for (int l = 0; l <= std::max(am1, am2); ++l) {
        sph_trans_.emplace_back(l);
        comps_.push_back(generate_am_components_cca(l));
    }

    const size_t esize = static_cast<size_t>(am1 + 1) * (am2 + 2) * (am1 + am2 + 3);
    Ex_.resize(esize);
    Ey_.resize(esize);
    Ez_.resize(esize);
}

BatchedPotentialInt::~BatchedPotentialInt() {}

void BatchedPotentialInt::boys(int L, size_t nk, const double* T, double* F, double* inv_T) {
    // Above T_crit, F_0 .. F_L follow the asymptotic upward recursion F_0 = sqrt(pi/4) sqrt(1/T),
    // F_n = F_{n-1} (n - 1/2) / T, which omits exp(-T) / (2T), to double precision (relative error below 1e-16);
    // the omitted term weighs more at higher orders. Evaluate it, as Libint2 does, for all the arguments at once
    // (at smaller T the values are meaningless, T = 0 included, and are overwritten); below T_crit,
    // BoysChebyshevTable for L <= BoysChebyshevTable::nmax, and Libint2's FmEval_Chebyshev7 argument by argument
    // for higher L.
    const double T_crit = L < static_cast<int>(sizeof(asymptotic_T) / sizeof(double)) ? asymptotic_T[L] : 117.0;
    for (size_t k = 0; k < nk; ++k) inv_T[k] = 1 / T[k];
    for (size_t k = 0; k < nk; ++k) F[k] = 0.88622692545275801365 * std::sqrt(inv_T[k]);
    for (int n = 1; n <= L; ++n) {
        const double ihalf = double(n) - double(0.5);
        const double* prev = &F[(n - 1) * nk];
        double* cur = &F[n * nk];
        for (size_t k = 0; k < nk; ++k) cur[k] = prev[k] * ihalf * inv_T[k];
    }
    double Tmin = T[0];
    for (size_t k = 1; k < nk; ++k) Tmin = std::min(Tmin, T[k]);
    if (Tmin > T_crit) return;
    if (L <= BoysChebyshevTable::nmax) {
        // the preprint's regime (3), for mixed batches as well: interpolate for all the arguments, and keep the
        // interpolated values where T <= T_crit
        table_F_.resize((L + 2) * nk);
        double* G = table_F_.data();
        BoysChebyshevTable::instance().eval(L, nk, T, G, G + (L + 1) * nk);
        for (int n = 0; n <= L; ++n) {
            double* f = F + n * nk;
            const double* g = G + n * nk;
#pragma omp simd
            for (size_t k = 0; k < nk; ++k) f[k] = T[k] <= T_crit ? g[k] : f[k];
        }
        return;
    }
    double Fm[64];
    for (size_t k = 0; k < nk; ++k) {
        if (T[k] > T_crit) continue;
        fm_eval_->eval(Fm, T[k], L);
        for (int n = 0; n <= L; ++n) F[n * nk + k] = Fm[n];
    }
}

void BatchedPotentialInt::compute_cart_vrr(const libint2::Shell& sh1, const libint2::Shell& sh2, size_t npoints,
                                           const double* x, const double* y, const double* z) {
    // Cartesian integrals by the generated Obara-Saika VRR or HRR kernels of snk_integral_codegen (see
    // snk/README.md), for the significant primitive pairs, over chunks of the points
    const int am1 = sh1.contr[0].l, am2 = sh2.contr[0].l;
    const int L = am1 + am2;
    const size_t ncart12 = comps_[am1].size() * comps_[am2].size();
    const double A[3] = {sh1.O[0], sh1.O[1], sh1.O[2]};
    const double B[3] = {sh2.O[0], sh2.O[1], sh2.O[2]};
    const bool same_center = A[0] == B[0] && A[1] == B[1] && A[2] == B[2];
    const double AB2 = (A[0] - B[0]) * (A[0] - B[0]) + (A[1] - B[1]) * (A[1] - B[1]) + (A[2] - B[2]) * (A[2] - B[2]);

    // the significant primitive pairs: exponent p, center P, P - A, and seed prefactor (2 pi / p) U_p
    prim_p_.clear();
    prim_P_.clear();
    prim_PA_.clear();
    prim_U_.clear();
    for (size_t p1 = 0; p1 < sh1.nprim(); ++p1) {
        const double a = sh1.alpha[p1], ca = sh1.contr[0].coeff[p1];
        for (size_t p2 = 0; p2 < sh2.nprim(); ++p2) {
            const double b = sh2.alpha[p2], cb = sh2.contr[0].coeff[p2];
            const double p = a + b;
            const double prefac = 2.0 * M_PI * ca * cb / p;
            const double U = prefac * std::exp(-a * b / p * AB2);
            if (std::abs(U) < prim_threshold) continue;
            prim_p_.push_back(p);
            for (int i = 0; i < 3; ++i) {
                const double Pi = (a * A[i] + b * B[i]) / p;
                prim_P_.push_back(Pi);
                prim_PA_.push_back(Pi - A[i]);
            }
            prim_U_.push_back(U);
        }
    }
    const int nprim = static_cast<int>(prim_p_.size());

    // The kernel, by the rules of the preprint (sec. III B): HRR for more than two primitive pairs; for two, if
    // am(A) + am(B) > 3 (and am(B) > 1 on the same center); for one, if am(A) + am(B) > 6; VRR otherwise. HRR
    // kernels exist for am(A), am(B) >= 1, in either order; VRR kernels for am(A) >= am(B).
    const bool hrr_exists = am1 >= 1 && am2 >= 1;
    const bool use_hrr = hrr_exists && (nprim > 2 || (nprim == 2 && L > 3 && (!same_center || am2 > 1)) ||
                                        (nprim == 1 && L > 6));
    if (!use_hrr && am1 < am2) {
        // VRR for the swapped pair, [c2][c1][k], transposed
        compute_cart_vrr(sh2, sh1, npoints, x, y, z);
        const size_t ncart1 = comps_[am1].size(), ncart2 = comps_[am2].size();
        vrr_out_.assign(cart_.begin(), cart_.end());
        for (size_t c1 = 0; c1 < ncart1; ++c1)
            for (size_t c2 = 0; c2 < ncart2; ++c2)
                std::copy(&vrr_out_[(c2 * ncart1 + c1) * npoints], &vrr_out_[(c2 * ncart1 + c1) * npoints] + npoints,
                          &cart_[(c1 * ncart2 + c2) * npoints]);
        return;
    }

    cart_.assign(ncart12 * npoints, 0.0);
    if (nprim == 0) return;

    // P - A by direction, as the kernels take them
    std::vector<double>& PAx = vrr_PA_;
    PAx.resize(3 * nprim);
    double *PAy = PAx.data() + nprim, *PAz = PAy + nprim;
    for (int q = 0; q < nprim; ++q) {
        PAx[q] = prim_PA_[3 * q];
        PAy[q] = prim_PA_[3 * q + 1];
        PAz[q] = prim_PA_[3 * q + 2];
    }

    constexpr size_t vrr_chunk = 128;
    const size_t nseed = static_cast<size_t>(L + 1);
    for (size_t k0 = 0; k0 < npoints; k0 += vrr_chunk) {
        const size_t nk = std::min(vrr_chunk, npoints - k0);
        // per primitive pair: P - G for the points, T, and the seeds (2 pi / p) U_p F_N(T)
        vrr_PG_.resize(3 * nprim * nk);
        vrr_seed_.resize(nprim * nseed * nk);
        vrr_T_.resize(2 * nk);
        if (boys_.size() < nseed * nk) boys_.resize(nseed * nk);
        double* PGx = vrr_PG_.data();
        double* PGy = PGx + nprim * nk;
        double* PGz = PGy + nprim * nk;
        for (int q = 0; q < nprim; ++q) {
            const double p = prim_p_[q];
            const double* P = &prim_P_[3 * q];
            double* gx = PGx + q * nk;
            double* gy = PGy + q * nk;
            double* gz = PGz + q * nk;
            double* T = vrr_T_.data();
            for (size_t k = 0; k < nk; ++k) {
                gx[k] = P[0] - x[k0 + k];
                gy[k] = P[1] - y[k0 + k];
                gz[k] = P[2] - z[k0 + k];
                T[k] = p * (gx[k] * gx[k] + gy[k] * gy[k] + gz[k] * gz[k]);
            }
            boys(L, nk, T, boys_.data(), T + nk);
            double* seed = &vrr_seed_[q * nseed * nk];
            const double U = prim_U_[q];
            for (size_t i = 0; i < nseed * nk; ++i) seed[i] = U * boys_[i];
        }

        vrr_out_.resize(ncart12 * nk);
        const int gb = static_cast<int>(nk);
        if (use_hrr) {
            hrr_base_.resize(static_cast<size_t>(snk::hrr_base_size(am1, am2, same_center)) * nk);
            if (same_center)
                snk::hrr_same_center[am1][am2](vrr_out_.data(), hrr_base_.data(), vrr_seed_.data(), prim_p_.data(), PGx,
                                               PGy, PGz, nprim, gb);
            else
                snk::hrr[am1][am2](vrr_out_.data(), hrr_base_.data(), vrr_seed_.data(), prim_p_.data(), PGx, PGy, PGz,
                                   PAx.data(), PAy, PAz, A[0] - B[0], A[1] - B[1], A[2] - B[2], nprim, gb);
        } else if (same_center) {
            snk::vrr_same_center[am1][am2](vrr_out_.data(), vrr_seed_.data(), prim_p_.data(), PGx, PGy, PGz, nprim, gb);
        } else {
            snk::vrr[am1][am2](vrr_out_.data(), vrr_seed_.data(), prim_p_.data(), PGx, PGy, PGz, PAx.data(), PAy, PAz,
                               A[0] - B[0], A[1] - B[1], A[2] - B[2], nprim, gb);
        }
        for (size_t c = 0; c < ncart12; ++c)
            std::copy(&vrr_out_[c * nk], &vrr_out_[c * nk] + nk, &cart_[c * npoints + k0]);
    }
}

void BatchedPotentialInt::compute_cart_md(const libint2::Shell& sh1, const libint2::Shell& sh2, size_t npoints,
                                          const double* x, const double* y, const double* z) {
    // Cartesian integrals by the McMurchie-Davidson scheme, for shells of any angular momentum
    const int am1 = sh1.contr[0].l, am2 = sh2.contr[0].l;
    const int L = am1 + am2;
    const auto& comps1 = comps_[am1];
    const auto& comps2 = comps_[am2];
    const size_t ncart1 = comps1.size(), ncart2 = comps2.size();

    cart_.assign(ncart1 * ncart2 * npoints, 0.0);

    const size_t nchunk_max = chunk_size(L);
    const int dim1 = L + 1;
    const size_t dim2 = static_cast<size_t>(dim1) * dim1 * dim1;  // (t, u, v) box
    if (R_.size() < dim1 * dim2 * nchunk_max) R_.resize(dim1 * dim2 * nchunk_max);
    if (boys_.size() < (dim1 + 2) * nchunk_max) boys_.resize((dim1 + 2) * nchunk_max);

    // E coefficient layout of fill_E_matrix for these angular momenta
    const int edim2 = am2 + 1, edim3 = am1 + am2 + 2;

    const Point A{sh1.O[0], sh1.O[1], sh1.O[2]};
    const Point B{sh2.O[0], sh2.O[1], sh2.O[2]};

    for (size_t p1 = 0; p1 < sh1.nprim(); ++p1) {
        const double a = sh1.alpha[p1], ca = sh1.contr[0].coeff[p1];
        for (size_t p2 = 0; p2 < sh2.nprim(); ++p2) {
            const double b = sh2.alpha[p2], cb = sh2.contr[0].coeff[p2];
            const double p = a + b;
            const Point P{(a * A[0] + b * B[0]) / p, (a * A[1] + b * B[1]) / p, (a * A[2] + b * B[2]) / p};
            const double prefac = 2.0 * M_PI * ca * cb / p;

            // Skip primitive pairs whose Gaussian product is negligible: their integrals are bounded by
            // |prefac| exp(-ab/p |A - B|^2) times polynomial factors of order one
            const double AB2 = (A[0] - B[0]) * (A[0] - B[0]) + (A[1] - B[1]) * (A[1] - B[1]) + (A[2] - B[2]) * (A[2] - B[2]);
            if (std::abs(prefac) * std::exp(-a * b / p * AB2) < prim_threshold) continue;

            // Hermite expansion coefficients of the primitive pair: independent of the points
            fill_E_matrix(am1, am2, P, A, B, a, b, Ex_, Ey_, Ez_);

            for (size_t k0 = 0; k0 < npoints; k0 += nchunk_max) {
                const size_t nk = std::min(nchunk_max, npoints - k0);
                double* R = R_.data();
                auto Rat = [&](int n, int t, int u, int v) -> double* {
                    return R + ((n * dim2) + address_3d(t, u, v, dim1, dim1)) * nk;
                };

                // Boys function F_n(T_k), T_k = p |P - C_k|^2, into boys_[n][k]
                double* T = &boys_[(L + 1) * nk];
                for (size_t k = 0; k < nk; ++k) {
                    const double PCx = P[0] - x[k0 + k], PCy = P[1] - y[k0 + k], PCz = P[2] - z[k0 + k];
                    T[k] = p * (PCx * PCx + PCy * PCy + PCz * PCz);
                }
                boys(L, nk, T, boys_.data(), T + nk);
                // R^n_000 = (-2p)^n F_n (eq 9.9.14)
                double fac = 1.0;
                for (int n = 0; n <= L; ++n) {
                    double* out = Rat(n, 0, 0, 0);
                    const double* F = &boys_[n * nk];
                    for (size_t k = 0; k < nk; ++k) out[k] = fac * F[k];
                    fac *= -2.0 * p;
                }

                // Hermite Coulomb integrals R^n_tuv, t + u + v <= L - n, by the recursions of
                // fill_R_matrix (eqs 9.9.18-20), with the points innermost
                const double* Cx = x + k0;
                const double* Cy = y + k0;
                const double* Cz = z + k0;
                for (int v = 1; v <= L; ++v) {
                    for (int n = 0; n + v <= L; ++n) {
                        double* out = Rat(n, 0, 0, v);
                        const double* r1 = Rat(n + 1, 0, 0, v - 1);
                        const double* r2 = v > 1 ? Rat(n + 1, 0, 0, v - 2) : nullptr;
                        for (size_t k = 0; k < nk; ++k) {
                            double val = (P[2] - Cz[k]) * r1[k];
                            if (r2) val += (v - 1) * r2[k];
                            out[k] = val;
                        }
                    }
                }
                for (int v = 0; v <= L; ++v) {
                    for (int u = 1; u + v <= L; ++u) {
                        for (int n = 0; n + u + v <= L; ++n) {
                            double* out = Rat(n, 0, u, v);
                            const double* r1 = Rat(n + 1, 0, u - 1, v);
                            const double* r2 = u > 1 ? Rat(n + 1, 0, u - 2, v) : nullptr;
                            for (size_t k = 0; k < nk; ++k) {
                                double val = (P[1] - Cy[k]) * r1[k];
                                if (r2) val += (u - 1) * r2[k];
                                out[k] = val;
                            }
                        }
                    }
                }
                for (int v = 0; v <= L; ++v) {
                    for (int u = 0; u + v <= L; ++u) {
                        for (int t = 1; t + u + v <= L; ++t) {
                            for (int n = 0; n + t + u + v <= L; ++n) {
                                double* out = Rat(n, t, u, v);
                                const double* r1 = Rat(n + 1, t - 1, u, v);
                                const double* r2 = t > 1 ? Rat(n + 1, t - 2, u, v) : nullptr;
                                for (size_t k = 0; k < nk; ++k) {
                                    double val = (P[0] - Cx[k]) * r1[k];
                                    if (r2) val += (t - 1) * r2[k];
                                    out[k] = val;
                                }
                            }
                        }
                    }
                }

                // (a|C|b) = 2 pi / p sum_tuv E^ab_t E^ab_u E^ab_v R^0_tuv (eq 9.9.32)
                for (size_t c1 = 0; c1 < ncart1; ++c1) {
                    const int lx1 = comps1[c1][0], ly1 = comps1[c1][1], lz1 = comps1[c1][2];
                    for (size_t c2 = 0; c2 < ncart2; ++c2) {
                        const int lx2 = comps2[c2][0], ly2 = comps2[c2][1], lz2 = comps2[c2][2];
                        const double* Ex = &Ex_[edim3 * (lx2 + edim2 * lx1)];
                        const double* Ey = &Ey_[edim3 * (ly2 + edim2 * ly1)];
                        const double* Ez = &Ez_[edim3 * (lz2 + edim2 * lz1)];
                        double* dst = &cart_[(c1 * ncart2 + c2) * npoints + k0];
                        for (int t = 0; t <= lx1 + lx2; ++t) {
                            for (int u = 0; u <= ly1 + ly2; ++u) {
                                const double Etu = prefac * Ex[t] * Ey[u];
                                for (int v = 0; v <= lz1 + lz2; ++v) {
                                    const double E = Etu * Ez[v];
                                    const double* src = Rat(0, t, u, v);
                                    for (size_t k = 0; k < nk; ++k) dst[k] += E * src[k];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void BatchedPotentialInt::compute_shell(int s1, int s2, size_t npoints, const double* x, const double* y,
                                      const double* z) {
    const libint2::Shell& sh1 = bs1_->l2_shell(s1);
    const libint2::Shell& sh2 = bs2_->l2_shell(s2);
    const int am1 = sh1.contr[0].l, am2 = sh2.contr[0].l;
    const size_t ncart1 = comps_[am1].size(), ncart2 = comps_[am2].size();

    // => Cartesian integrals, [c1][c2][k] <= //

    if (am1 <= snk::max_am && am2 <= snk::max_am) {
        compute_cart_vrr(sh1, sh2, npoints, x, y, z);
    } else {
        compute_cart_md(sh1, sh2, npoints, x, y, z);
    }

    // => Spherical transforms, one index at a time <= //

    const bool pure1 = sh1.contr[0].pure && am1 > 0, pure2 = sh2.contr[0].pure && am2 > 0;
    const size_t n1 = pure1 ? 2 * am1 + 1 : ncart1, n2 = pure2 ? 2 * am2 + 1 : ncart2;

    const std::vector<double>* in = &cart_;
    if (pure2) {
        half_.assign(ncart1 * n2 * npoints, 0.0);
        const SphericalTransform& st = sph_trans_[am2];
        for (size_t c1 = 0; c1 < ncart1; ++c1)
            for (int i = 0; i < st.n(); ++i) {
                const double coef = st.coef(i);
                const double* src = &cart_[(c1 * ncart2 + st.cartindex(i)) * npoints];
                double* dst = &half_[(c1 * n2 + st.pureindex(i)) * npoints];
                for (size_t k = 0; k < npoints; ++k) dst[k] += coef * src[k];
            }
        in = &half_;
    }
    if (pure1) {
        out_.assign(n1 * n2 * npoints, 0.0);
        const SphericalTransform& st = sph_trans_[am1];
        for (int i = 0; i < st.n(); ++i) {
            const double coef = st.coef(i);
            for (size_t i2 = 0; i2 < n2; ++i2) {
                const double* src = &(*in)[(st.cartindex(i) * n2 + i2) * npoints];
                double* dst = &out_[(st.pureindex(i) * n2 + i2) * npoints];
                for (size_t k = 0; k < npoints; ++k) dst[k] += coef * src[k];
            }
        }
    } else {
        out_.assign(in->begin(), in->begin() + n1 * n2 * npoints);
    }
}

namespace {

/// Points per BatchedPotentialInt call in the contraction drivers, to bound the size of its buffers
constexpr size_t driver_chunk = 512;

/// Shell pairs of bs1 x bs2: with s2 <= s1 if they are the same basis set.
std::vector<std::pair<int, int>> shell_pairs(const BasisSet& bs1, const BasisSet& bs2, bool same) {
    std::vector<std::pair<int, int>> pairs;
    for (int s1 = 0; s1 < bs1.nshell(); ++s1)
        for (int s2 = 0; s2 < (same ? s1 + 1 : bs2.nshell()); ++s2) pairs.emplace_back(s1, s2);
    return pairs;
}

/// The points as three coordinate arrays.
void split_points(const Matrix& points, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z) {
    if (points.nirrep() != 1 || points.coldim() != 3)
        throw PSIEXCEPTION("BatchedPotentialInt: points must be an npoints x 3 matrix.");
    const size_t npoints = points.rowdim();
    x.resize(npoints);
    y.resize(npoints);
    z.resize(npoints);
    for (size_t k = 0; k < npoints; ++k) {
        x[k] = points.get(k, 0);
        y[k] = points.get(k, 1);
        z[k] = points.get(k, 2);
    }
}

int driver_threads() {
#ifdef _OPENMP
    return Process::environment.get_n_threads();
#else
    return 1;
#endif
}

}  // namespace

void batched_potential_contract_density(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2, const Matrix& D,
                                        const Matrix& points, double* V) {
    if (D.nirrep() != 1 || D.rowdim() != bs1->nbf() || D.coldim() != bs2->nbf())
        throw PSIEXCEPTION("batched_potential_contract_density: D must be nbf(bs1) x nbf(bs2).");
    std::vector<double> x, y, z;
    split_points(points, x, y, z);
    const size_t npoints = x.size();
    const bool same = bs1 == bs2;
    const auto pairs = shell_pairs(*bs1, *bs2, same);
    double** Dp = D.pointer();

    const int nthreads = driver_threads();
    std::vector<std::vector<double>> VT(nthreads, std::vector<double>(npoints, 0.0));

#pragma omp parallel num_threads(nthreads)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        BatchedPotentialInt ints(bs1, bs2);
        double* Vt = VT[rank].data();
#pragma omp for schedule(dynamic)
        for (size_t ip = 0; ip < pairs.size(); ++ip) {
            const int P = pairs[ip].first, Q = pairs[ip].second;
            const int n1 = bs1->shell(P).nfunction(), n2 = bs2->shell(Q).nfunction();
            const int o1 = bs1->shell_to_basis_function(P), o2 = bs2->shell_to_basis_function(Q);
            const bool both = same && P != Q;  // (Q, P) is not in the list: it contributes the same integrals
            for (size_t k0 = 0; k0 < npoints; k0 += driver_chunk) {
                const size_t nk = std::min(driver_chunk, npoints - k0);
                ints.compute_shell(P, Q, nk, &x[k0], &y[k0], &z[k0]);
                const double* A = ints.buffer();
                for (int i1 = 0; i1 < n1; ++i1) {
                    for (int i2 = 0; i2 < n2; ++i2, A += nk) {
                        double d = Dp[o1 + i1][o2 + i2];
                        if (both) d += Dp[o2 + i2][o1 + i1];
                        for (size_t k = 0; k < nk; ++k) Vt[k0 + k] += d * A[k];
                    }
                }
            }
        }
    }

    for (int t = 0; t < nthreads; ++t)
        for (size_t k = 0; k < npoints; ++k) V[k] += VT[t][k];
}

void batched_potential_contract_charges(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2, const double* q,
                                        const Matrix& points, Matrix& F) {
    if (F.nirrep() != 1 || F.rowdim() != bs1->nbf() || F.coldim() != bs2->nbf())
        throw PSIEXCEPTION("batched_potential_contract_charges: F must be nbf(bs1) x nbf(bs2).");
    std::vector<double> x, y, z;
    split_points(points, x, y, z);
    const size_t npoints = x.size();
    const bool same = bs1 == bs2;
    const auto pairs = shell_pairs(*bs1, *bs2, same);
    double** Fp = F.pointer();

    const int nthreads = driver_threads();
#pragma omp parallel num_threads(nthreads)
    {
        BatchedPotentialInt ints(bs1, bs2);
        std::vector<double> block;
        // each shell pair writes only its own blocks of F
#pragma omp for schedule(dynamic)
        for (size_t ip = 0; ip < pairs.size(); ++ip) {
            const int P = pairs[ip].first, Q = pairs[ip].second;
            const int n1 = bs1->shell(P).nfunction(), n2 = bs2->shell(Q).nfunction();
            const int o1 = bs1->shell_to_basis_function(P), o2 = bs2->shell_to_basis_function(Q);
            block.assign(static_cast<size_t>(n1) * n2, 0.0);
            for (size_t k0 = 0; k0 < npoints; k0 += driver_chunk) {
                const size_t nk = std::min(driver_chunk, npoints - k0);
                ints.compute_shell(P, Q, nk, &x[k0], &y[k0], &z[k0]);
                const double* A = ints.buffer();
                for (int i12 = 0; i12 < n1 * n2; ++i12, A += nk) {
                    double f = 0.0;
                    for (size_t k = 0; k < nk; ++k) f += q[k0 + k] * A[k];
                    block[i12] += f;
                }
            }
            for (int i1 = 0; i1 < n1; ++i1)
                for (int i2 = 0; i2 < n2; ++i2) {
                    Fp[o1 + i1][o2 + i2] += block[i1 * n2 + i2];
                    if (same && P != Q) Fp[o2 + i2][o1 + i1] += block[i1 * n2 + i2];
                }
        }
    }
}

}  // namespace psi
