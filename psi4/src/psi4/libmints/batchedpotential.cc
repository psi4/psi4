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

void BatchedPotentialInt::compute_shell(int s1, int s2, size_t npoints, const double* x, const double* y,
                                      const double* z) {
    const libint2::Shell& sh1 = bs1_->l2_shell(s1);
    const libint2::Shell& sh2 = bs2_->l2_shell(s2);
    const int am1 = sh1.contr[0].l, am2 = sh2.contr[0].l;
    const int L = am1 + am2;
    const auto& comps1 = comps_[am1];
    const auto& comps2 = comps_[am2];
    const size_t ncart1 = comps1.size(), ncart2 = comps2.size();

    // => Cartesian integrals, [c1][c2][k] <= //

    cart_.assign(ncart1 * ncart2 * npoints, 0.0);

    const size_t nchunk_max = chunk_size(L);
    const int dim1 = L + 1;
    const size_t dim2 = static_cast<size_t>(dim1) * dim1 * dim1;  // (t, u, v) box
    if (R_.size() < dim1 * dim2 * nchunk_max) R_.resize(dim1 * dim2 * nchunk_max);
    std::vector<double> Fm(dim1);

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

            // Hermite expansion coefficients of the primitive pair: independent of the points
            fill_E_matrix(am1, am2, P, A, B, a, b, Ex_, Ey_, Ez_);

            for (size_t k0 = 0; k0 < npoints; k0 += nchunk_max) {
                const size_t nk = std::min(nchunk_max, npoints - k0);
                double* R = R_.data();
                auto Rat = [&](int n, int t, int u, int v) -> double* {
                    return R + ((n * dim2) + address_3d(t, u, v, dim1, dim1)) * nk;
                };

                // Boys function and R^n_000 = (-2p)^n F_n(p |P - C|^2) (eq 9.9.14)
                for (size_t k = 0; k < nk; ++k) {
                    const double PCx = P[0] - x[k0 + k], PCy = P[1] - y[k0 + k], PCz = P[2] - z[k0 + k];
                    fm_eval_->eval(Fm.data(), p * (PCx * PCx + PCy * PCy + PCz * PCz), L);
                    double fac = 1.0;
                    for (int n = 0; n <= L; ++n) {
                        Rat(n, 0, 0, 0)[k] = fac * Fm[n];
                        fac *= -2.0 * p;
                    }
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
