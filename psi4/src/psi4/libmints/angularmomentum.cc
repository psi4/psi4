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

#include "psi4/libmints/angularmomentum.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"

#include <memory>
#include <stdexcept>

#include <libint2/shell.h>

using namespace psi;
using namespace mdintegrals;

AngularMomentumInt::AngularMomentumInt(std::vector<SphericalTransform>& spherical_transforms,
                                       std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2, int nderiv)
    : OneBodyAOInt(spherical_transforms, bs1, bs2, nderiv), MDHelper(bs1->max_am(), bs2->max_am()) {
    int maxnao1 = INT_NCART(maxam1_);
    int maxnao2 = INT_NCART(maxam2_);

    if (deriv_ == 0) {
        // one chunk for each Cartesian component
        buffer_ = new double[3 * maxnao1 * maxnao2];
        set_chunks(3);
    } else {
        throw PSIEXCEPTION("AngularMomentumInt does not provide derivatives.");
    }
    buffers_.resize(nchunk_);
}

AngularMomentumInt::~AngularMomentumInt() { delete[] buffer_; }

void AngularMomentumInt::compute_pair(const libint2::Shell& s1, const libint2::Shell& s2) {
    // Computes the angular momentum integrals for a pair of shells using eq 9.3.33
    // equation numbers from Molecular Electronic-Structure Theory (10.1002/9781119019572)
    int am1 = s1.contr[0].l;
    int am2 = s2.contr[0].l;

    const auto& comps_am1 = am_comps_[am1];
    const auto& comps_am2 = am_comps_[am2];

    auto A = s1.O;
    auto B = s2.O;
    // TODO: settable origin? (not supported previously)
    Point C = {0.0, 0.0, 0.0};
    int nprim1 = s1.nprim();
    int nprim2 = s2.nprim();

    int dim1 = INT_NCART(am1);
    int dim2 = INT_NCART(am2);
    int size = dim1 * dim2;
    memset(buffer_, 0, 3 * size * sizeof(double));

    // dimensions of the E matrix
    int edim2 = am2 + 2;
    int edim3 = (am1 + 1) + edim2;

    int ao12 = 0;
    for (int p1 = 0; p1 < nprim1; ++p1) {
        double a = s1.alpha[p1];
        double ca = s1.contr[0].coeff[p1];
        for (int p2 = 0; p2 < nprim2; ++p2) {
            double b = s2.alpha[p2];
            double cb = s2.contr[0].coeff[p2];
            double p = a + b;
            Point P{(a * A[0] + b * B[0]) / p, (a * A[1] + b * B[1]) / p, (a * A[2] + b * B[2]) / p};
            auto PC = point_diff(P, C);
            double prefac = ca * cb * std::pow(M_PI / p, 1.5);
            fill_E_matrix(am1, am2 + 1, P, A, B, a, b, Ex, Ey, Ez);
            ao12 = 0;
            for (const auto& [l1, m1, n1] : comps_am1) {
                for (const auto& [l2, m2, n2] : comps_am2) {
                    int xidx = address_3d(l1, l2, 0, edim2, edim3);
                    int yidx = address_3d(m1, m2, 0, edim2, edim3);
                    int zidx = address_3d(n1, n2, 0, edim2, edim3);
                    double Sx0 = Ex[xidx];
                    double Sy0 = Ey[yidx];
                    double Sz0 = Ez[zidx];
                    // Multipole integral, eq 9.5.43
                    double Sx1 = Ex[xidx + 1] + PC[0] * Ex[xidx];
                    double Sy1 = Ey[yidx + 1] + PC[1] * Ey[yidx];
                    double Sz1 = Ez[zidx + 1] + PC[2] * Ez[zidx];
                    // Use eq 9.3.30 to define derivative w.r.t. ket center
                    double Dx1 = -2.0 * b * Ex[address_3d(l1, l2 + 1, 0, edim2, edim3)];
                    double Dy1 = -2.0 * b * Ey[address_3d(m1, m2 + 1, 0, edim2, edim3)];
                    double Dz1 = -2.0 * b * Ez[address_3d(n1, n2 + 1, 0, edim2, edim3)];
                    if (l2 > 0) {
                        Dx1 += l2 * Ex[address_3d(l1, l2 - 1, 0, edim2, edim3)];
                    }
                    if (m2 > 0) {
                        Dy1 += m2 * Ey[address_3d(m1, m2 - 1, 0, edim2, edim3)];
                    }
                    if (n2 > 0) {
                        Dz1 += n2 * Ez[address_3d(n1, n2 - 1, 0, edim2, edim3)];
                    }
                    // Lx
                    buffer_[ao12 + size * 0] += (Sx0 * Dy1 * Sz1 - Sx0 * Sy1 * Dz1) * prefac;
                    // Ly
                    buffer_[ao12 + size * 1] += (Sx1 * Sy0 * Dz1 - Dx1 * Sy0 * Sz1) * prefac;
                    // Lz
                    buffer_[ao12 + size * 2] += (Dx1 * Sy1 * Sz0 - Sx1 * Dy1 * Sz0) * prefac;
                    ++ao12;
                }
            }
        }
    }
    pure_transform(s1, s2, 3);
    buffers_[0] = buffer_ + 0 * s1.size() * s2.size();
    buffers_[1] = buffer_ + 1 * s1.size() * s2.size();
    buffers_[2] = buffer_ + 2 * s1.size() * s2.size();
}