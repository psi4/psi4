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

#include "psi4/libmints/multipoles.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/vector.h"

#include <libint2/engine.h>

using namespace psi;
using namespace mdintegrals;

uint64_t binomial(int n, int c1);  // From solidharmonics.cc

MultipoleInt::MultipoleInt(std::vector<SphericalTransform>& spherical_transforms, std::shared_ptr<BasisSet> bs1,
                           std::shared_ptr<BasisSet> bs2, int order, int nderiv)
    : OneBodyAOInt(spherical_transforms, bs1, bs2, nderiv), order_(order), MDHelper(bs1->max_am(), bs2->max_am()) {
    int maxnao1 = INT_NCART(maxam1_);
    int maxnao2 = INT_NCART(maxam2_);

    // The number of multipole components to compute. N.B. we don't compute the 0th one (overlap)
    int n_mult = cumulative_cart_dim(order_) - 1;

    // pre-allocate M-matrix
    int am = maxam1_ + maxam2_;
    int mdim0 = order_ + 1;
    int mdim1 = std::max(am, order_) + 2;
    int msize = mdim0 * mdim1;
    Mx = std::vector<double>(msize);
    My = std::vector<double>(msize);
    Mz = std::vector<double>(msize);

    // pre-allocate S-matrix
    int sdim0 = maxam1_ + 1;
    int sdim1 = maxam2_ + 1;
    int sdim2 = order_ + 1;
    int ssize = sdim0 * sdim1 * sdim2;
    Sx = std::vector<double>(ssize);
    Sy = std::vector<double>(ssize);
    Sz = std::vector<double>(ssize);

    // CCA-ordered Cartesian components for the multipoles
    comps_mul_ = std::vector<std::vector<std::array<int, 4>>>(order_ + 1);
    for (int d = 0; d < order_ + 1; ++d) {
        comps_mul_[d] = generate_am_components_cca(d);
    }

    // Increase buffer size to handle all Cartesian components
    if (deriv_ == 0) {
        buffer_ = new double[n_mult * maxnao1 * maxnao2];
        set_chunks(n_mult);
    } else {
        throw PSIEXCEPTION("Derivatives are NYI for arbitrary-order multipoles");
    }
    buffers_.resize(nchunk_);
}

MultipoleInt::~MultipoleInt() { delete[] buffer_; }

SharedVector MultipoleInt::nuclear_contribution(std::shared_ptr<Molecule> mol, int order, const Vector3& origin) {
    int ntot = cumulative_cart_dim(order) - 1;
    auto sret = std::make_shared<Vector>(ntot);
    double* ret = sret->pointer();

    int address = 0;
    for (int l = 1; l <= order; ++l) {
        for (int ii = 0; ii <= l; ii++) {
            int lx = l - ii;
            for (int lz = 0; lz <= ii; lz++) {
                int ly = ii - lz;
                for (int atom = 0; atom < mol->natom(); ++atom) {
                    Vector3 geom = mol->xyz(atom) - origin;
                    ret[address] += mol->Z(atom) * pow(geom[0], lx) * pow(geom[1], ly) * pow(geom[2], lz);
                }
                ++address;
            }
        }
    }

    return sret;
}

void MultipoleInt::compute_pair(const libint2::Shell& s1, const libint2::Shell& s2) {
    int am1 = s1.contr[0].l;
    int am2 = s2.contr[0].l;
    int am = am1 + am2;
    int nprim1 = s1.nprim();
    int nprim2 = s2.nprim();

    // zero out buffer
    memset(buffer_, 0, nchunk_ * INT_NCART(am1) * INT_NCART(am2) * sizeof(double));

    const auto& comps_am1 = am_comps_[am1];
    const auto& comps_am2 = am_comps_[am2];

    auto A = s1.O;
    auto B = s2.O;
    const Point C = {origin_[0], origin_[1], origin_[2]};

    int dim1 = INT_NCART(am1);
    int dim2 = INT_NCART(am2);
    // The number of bf components in each shell pair
    int size = dim1 * dim2;
    
    // dimensions of M and S matrix
    int mdim1 = std::max(am, order_) + 2;
    int sdim1 = am2 + 1;
    int sdim2 = order_ + 1;

    // dimensions of the E matrix
    int edim1 = am2 + 1;
    int edim2 = (am1 + 1) + edim1;

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
            fill_E_matrix(am1, am2, P, A, B, a, b, Ex, Ey, Ez);
            fill_M_matrix(am, order_, PC, a, b, Mx, My, Mz);

            std::fill(Sx.begin(), Sx.end(), 0);
            std::fill(Sy.begin(), Sy.end(), 0);
            std::fill(Sz.begin(), Sz.end(), 0);
            for (int i = 0; i <= am1; ++i) {
                for (int j = 0; j <= am2; ++j) {
                    for (int e = 0; e <= order_; ++e) {
                        int uppert = std::min(i + j, e);
                        int idx = address_3d(i, j, e, sdim1, sdim2);
                        for (int t = 0; t <= uppert; ++t) {
                            int idxt = address_3d(i, j, t, edim1, edim2);
                            // eq 9.5.39
                            Sx[idx] += Ex[idxt] * Mx[e * mdim1 + t];
                            Sy[idx] += Ey[idxt] * My[e * mdim1 + t];
                            Sz[idx] += Ez[idxt] * Mz[e * mdim1 + t];
                        }
                    }
                }
            }
            int m_count = 0;
            for (int mul = 1; mul < order_ + 1; ++mul) {
                const auto& comps_mul = comps_mul_[mul];
                for (const auto& [ex, ey, ez, index0] : comps_mul) {
                    ao12 = 0;
                    for (const auto& [l1, m1, n1, index1] : comps_am1) {
                        for (const auto& [l2, m2, n2, index2] : comps_am2) {
                            // -1.0 for consistency with dipole/quadrupole implementation
                            buffer_[ao12 + size * m_count] += -1.0 *
                                (ca * cb * Sx[address_3d(l1, l2, ex, sdim1, sdim2)] *
                                 Sy[address_3d(m1, m2, ey, sdim1, sdim2)] * Sz[address_3d(n1, n2, ez, sdim1, sdim2)]);
                            ao12++;
                        }
                    }
                    m_count++;
                }
            }
        }
    }
    pure_transform(s1, s2, nchunk_);
    for (int chunk = 0; chunk < nchunk_; ++chunk) {
        buffers_[chunk] = buffer_ + chunk * s1.size() * s2.size();
    }
}
