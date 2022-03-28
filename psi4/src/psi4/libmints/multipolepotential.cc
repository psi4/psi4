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

#include "psi4/libmints/multipolepotential.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"

#include <libint2/shell.h>
#include <libint2/boys.h>

using namespace psi;
using namespace mdintegrals;

MultipolePotentialInt::MultipolePotentialInt(std::vector<SphericalTransform>& spherical_transforms,
                                             std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2, int order,
                                             int deriv)
    : OneBodyAOInt(spherical_transforms, bs1, bs2, deriv), MDHelper(bs1->max_am(), bs2->max_am()), order_(order) {
    int maxnao1 = INT_NCART(maxam1_);
    int maxnao2 = INT_NCART(maxam2_);

    if (deriv > 0) {
        throw FeatureNotImplemented("LibMints", "MultipolePotentialInts called with deriv > 0", __FILE__, __LINE__);
    }

    // pre-allocate R matrix
    int am = maxam1_ + maxam2_;
    int rdim1 = am + order_ + 1;
    int rdim2 = rdim1 * rdim1 * rdim1;
    R = std::vector<double>(rdim1 * rdim2);

    // set up Boys function evaluator
    fm_eval_ = libint2::FmEval_Chebyshev7<double>::instance(am + order_);

    comps_der_ = std::vector<std::vector<std::array<int, 4>>>(order_ + 1);
    for (int d = 0; d < order_ + 1; ++d) {
        comps_der_[d] = generate_am_components_cca(d);
    }
    int nchunks = cumulative_cart_dim(order_);
    buffer_ = new double[nchunks * maxnao1 * maxnao2];
    set_chunks(nchunks);
    buffers_.resize(nchunk_);
}

MultipolePotentialInt::~MultipolePotentialInt() { delete[] buffer_; }

void MultipolePotentialInt::compute_pair(const libint2::Shell& s1, const libint2::Shell& s2) {
    double Cx = origin_[0];
    double Cy = origin_[1];
    double Cz = origin_[2];
    Point C{Cx, Cy, Cz};

    int am1 = s1.contr[0].l;
    int am2 = s2.contr[0].l;
    int am = am1 + am2;

    const auto& comps_am1 = am_comps_[am1];
    const auto& comps_am2 = am_comps_[am2];

    auto A = s1.O;
    auto B = s2.O;
    int nprim1 = s1.nprim();
    int nprim2 = s2.nprim();

    // output buffer dimensions
    int dim1 = INT_NCART(am1);
    int dim2 = INT_NCART(am2);
    int size = dim1 * dim2;
    memset(buffer_, 0, nchunk_ * size * sizeof(double));

    // pre-allocate R-matrix
    int r_am = am + order_;
    int rdim1 = r_am + 1;
    int rdim2 = rdim1 * rdim1 * rdim1;
    std::vector<double> Rmat(rdim1 * rdim2);

    // dimensions of the E matrix
    int edim2 = am2 + 1;
    int edim3 = am1 + am2 + 2;

    int ao12 = 0;
    for (int p1 = 0; p1 < nprim1; ++p1) {
        double a = s1.alpha[p1];
        double ca = s1.contr[0].coeff[p1];
        for (int p2 = 0; p2 < nprim2; ++p2) {
            double b = s2.alpha[p2];
            double cb = s2.contr[0].coeff[p2];
            double p = a + b;
            Point P{(a * A[0] + b * B[0]) / p, (a * A[1] + b * B[1]) / p, (a * A[2] + b * B[2]) / p};
            double prefac = 2.0 * M_PI * ca * cb / p;
            fill_E_matrix(am1, am2, P, A, B, a, b, Ex, Ey, Ez);
            fill_R_matrix(r_am, p, P, C, R, fm_eval_);
            int der_count = 0;
            double sign_prefac = prefac;
            for (int der = 0; der < order_ + 1; ++der) {
                const auto& comps = comps_der_[der];
                for (const auto& [ex, ey, ez, index0] : comps) {
                    ao12 = 0;
                    for (const auto& [l1, m1, n1, index1] : comps_am1) {
                        for (const auto& [l2, m2, n2, index2] : comps_am2) {
                            double val = 0.0;
                            int maxt = l1 + l2;
                            int maxu = m1 + m2;
                            int maxv = n1 + n2;
                            const double* ex_p = &Ex.data()[edim3 * (l2 + edim2 * l1)];
                            const double* ey_p = &Ey.data()[edim3 * (m2 + edim2 * m1)];
                            const double* ez_p = &Ez.data()[edim3 * (n2 + edim2 * n1)];
                            for (int t = 0; t <= maxt; ++t) {
                                for (int u = 0; u <= maxu; ++u) {
                                    for (int v = 0; v <= maxv; ++v) {
                                        // eq 9.9.32
                                        val += ex_p[t] * ey_p[u] * ez_p[v] * R[address_3d(t + ex, u + ey, v + ez, rdim1, rdim1)];
                                    }
                                }
                            }
                            buffer_[ao12 + size * der_count] += sign_prefac * val;
                            ++ao12;
                        }
                    }
                    der_count++;
                }
                sign_prefac *= -1.0;
            }
        }
    }
    pure_transform(s1, s2, cumulative_cart_dim(order_));
    for (int chunk = 0; chunk < cumulative_cart_dim(order_); ++chunk) {
        buffers_[chunk] = buffer_ + chunk * s1.size() * s2.size();
    }
}