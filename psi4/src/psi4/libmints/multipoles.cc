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

#include "psi4/libmints/multipoles.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/vector.h"

#include <libint2/shell.h>

using namespace psi;
using namespace mdintegrals;


MultipoleInt::MultipoleInt(std::vector<SphericalTransform>& spherical_transforms, std::shared_ptr<BasisSet> bs1,
                           std::shared_ptr<BasisSet> bs2, int order, int nderiv)
    : OneBodyAOInt(spherical_transforms, bs1, bs2, nderiv),
      order_(order),
      MDHelper(bs1->max_am() + nderiv, bs2->max_am() + nderiv) {
    int maxnao1 = INT_NCART(maxam1_);
    int maxnao2 = INT_NCART(maxam2_);

    if (order_ == 0) {
        throw PSIEXCEPTION(
            "Zeroth order multipoles are not accessible via MultipoleInt."
            " Use OverlapInt instead.");
    }

    // The number of multipole components to compute. N.B. we don't compute the 0th one (overlap)
    int n_mult = cumulative_cart_dim(order_) - 1;

    // Increase buffer size to handle all Cartesian components
    if (deriv_ == 0) {
        buffer_ = new double[n_mult * maxnao1 * maxnao2];
        set_chunks(n_mult);
    } else if (deriv_ == 1) {
        buffer_ = new double[6 * n_mult * maxnao1 * maxnao2];
        set_chunks(6 * n_mult);
    } else {
        throw PSIEXCEPTION("Only first derivatives are available for arbitrary-order multipoles.");
    }
    buffers_.resize(nchunk_);

    // pre-allocate M-matrix
    int am = maxam1_ + maxam2_;
    int mdim0 = order_ + 1;
    int mdim1 = std::max(am, order_) + 2;
    int msize = mdim0 * mdim1;
    Mx = std::vector<double>(msize);
    My = std::vector<double>(msize);
    Mz = std::vector<double>(msize);

    // pre-allocate S-matrix
    int sdim0 = maxam1_ + 1 + deriv_;
    int sdim1 = maxam2_ + 1 + deriv_;
    int sdim2 = order_ + 1;
    int ssize = sdim0 * sdim1 * sdim2;
    Sx = std::vector<double>(ssize);
    Sy = std::vector<double>(ssize);
    Sz = std::vector<double>(ssize);

    // CCA-ordered Cartesian components for the multipoles
    comps_mul_ = std::vector<std::vector<std::array<int, 3>>>(order_ + 1);
    for (int d = 0; d < order_ + 1; ++d) {
        comps_mul_[d] = generate_am_components_cca(d);
    }
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
            // compute S matrix according to eq 9.5.39
            for (int i = 0; i <= am1; ++i) {
                for (int j = 0; j <= am2; ++j) {
                    for (int e = 0; e <= order_; ++e) {
                        int uppert = std::min(i + j, e);
                        int idx = address_3d(i, j, e, sdim1, sdim2);  // S_{ij}^e
                        for (int t = 0; t <= uppert; ++t) {
                            int idxt = address_3d(i, j, t, edim1, edim2);  // E_t^{ij}
                            int idxm = e * mdim1 + t;                      // M_t^e
                            // eq 9.5.39
                            Sx[idx] += Ex[idxt] * Mx[idxm];
                            Sy[idx] += Ey[idxt] * My[idxm];
                            Sz[idx] += Ez[idxt] * Mz[idxm];
                        }
                    }
                }
            }
            // -1.0 for consistency with dipole/quadrupole implementation
            double prefac = -1.0 * ca * cb;
            int m_count = 0;
            for (int mul = 1; mul < order_ + 1; ++mul) {
                const auto& comps_mul = comps_mul_[mul];
                for (const auto& [ex, ey, ez] : comps_mul) {
                    ao12 = 0;
                    for (const auto& [l1, m1, n1] : comps_am1) {
                        for (const auto& [l2, m2, n2] : comps_am2) {
                            // multiply separable x, y, and z components (eq 9.3.12)
                            buffer_[ao12 + size * m_count] += prefac * Sx[address_3d(l1, l2, ex, sdim1, sdim2)] *
                                                              Sy[address_3d(m1, m2, ey, sdim1, sdim2)] *
                                                              Sz[address_3d(n1, n2, ez, sdim1, sdim2)];
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

void MultipoleInt::compute_pair_deriv1(const libint2::Shell& s1, const libint2::Shell& s2) {
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
    int sdim1 = am2 + 2;
    int sdim2 = order_ + 1;

    // dimensions of the E matrix
    int edim1 = am2 + 2;
    int edim2 = (am1 + 2) + edim1;

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
            // increase E matrix angular momentum for the derivative
            fill_E_matrix(am1 + 1, am2 + 1, P, A, B, a, b, Ex, Ey, Ez);
            fill_M_matrix(am, order_, PC, a, b, Mx, My, Mz);

            std::fill(Sx.begin(), Sx.end(), 0);
            std::fill(Sy.begin(), Sy.end(), 0);
            std::fill(Sz.begin(), Sz.end(), 0);
            // compute S matrix according to eq 9.5.39,
            // here with increased angular momentum for the derivatives
            for (int i = 0; i <= am1 + 1; ++i) {
                for (int j = 0; j <= am2 + 1; ++j) {
                    for (int e = 0; e <= order_; ++e) {
                        int uppert = std::min(i + j, e);
                        int idx = address_3d(i, j, e, sdim1, sdim2);  // S_{ij}^e
                        for (int t = 0; t <= uppert; ++t) {
                            int idxt = address_3d(i, j, t, edim1, edim2);  // E_t^{ij}
                            int idxm = e * mdim1 + t;                      // M_t^e
                            // eq 9.5.39
                            Sx[idx] += Ex[idxt] * Mx[idxm];
                            Sy[idx] += Ey[idxt] * My[idxm];
                            Sz[idx] += Ez[idxt] * Mz[idxm];
                        }
                    }
                }
            }
            double prefac = ca * cb;
            int m_count = 0;
            // loop over multipole order
            for (int mul = 1; mul < order_ + 1; ++mul) {
                const auto& comps_mul = comps_mul_[mul];
                // loop over components of the given multipole, i.e.,
                // x, y, z for dipole (mul = 1); xx, xy, xz, yy, ...
                // for quadrupole (mul = 2) and so on
                for (const auto& [ex, ey, ez] : comps_mul) {
                    ao12 = 0;
                    // loop over primitive angular momentum components
                    for (const auto& [l1, m1, n1] : comps_am1) {
                        for (const auto& [l2, m2, n2] : comps_am2) {
                            // get the 'non-differentiated' contributions to the
                            // multipole integral
                            double sx = Sx[address_3d(l1, l2, ex, sdim1, sdim2)];
                            double sy = Sy[address_3d(m1, m2, ey, sdim1, sdim2)];
                            double sz = Sz[address_3d(n1, n2, ez, sdim1, sdim2)];
                            // form derivatives (eq 9.3.30): since the Cartesian components of
                            // the multipole integral are separable, apply the product rule
                            // Ax (derivative wrt bra center x coordinate)
                            double DAx = -2.0 * a * Sx[address_3d(l1 + 1, l2, ex, sdim1, sdim2)];
                            if (l1) DAx += l1 * Sx[address_3d(l1 - 1, l2, ex, sdim1, sdim2)];
                            buffer_[ao12 + size * (m_count + 0)] += prefac * DAx * sy * sz;
                            // Ay
                            double DAy = -2.0 * a * Sy[address_3d(m1 + 1, m2, ey, sdim1, sdim2)];
                            if (m1) DAy += m1 * Sy[address_3d(m1 - 1, m2, ey, sdim1, sdim2)];
                            buffer_[ao12 + size * (m_count + 1)] += prefac * sx * DAy * sz;
                            // Az
                            double DAz = -2.0 * a * Sz[address_3d(n1 + 1, n2, ez, sdim1, sdim2)];
                            if (n1) DAz += n1 * Sz[address_3d(n1 - 1, n2, ez, sdim1, sdim2)];
                            buffer_[ao12 + size * (m_count + 2)] += prefac * sx * sy * DAz;
                            // Bx (derivative wrt ket center x coordinate)
                            double DBx = -2.0 * b * Sx[address_3d(l1, l2 + 1, ex, sdim1, sdim2)];
                            if (l2) DBx += l2 * Sx[address_3d(l1, l2 - 1, ex, sdim1, sdim2)];
                            buffer_[ao12 + size * (m_count + 3)] += prefac * DBx * sy * sz;
                            // By
                            double DBy = -2.0 * b * Sy[address_3d(m1, m2 + 1, ey, sdim1, sdim2)];
                            if (m2) DBy += m2 * Sy[address_3d(m1, m2 - 1, ey, sdim1, sdim2)];
                            buffer_[ao12 + size * (m_count + 4)] += prefac * sx * DBy * sz;
                            // Bz
                            double DBz = -2.0 * b * Sz[address_3d(n1, n2 + 1, ez, sdim1, sdim2)];
                            if (n2) DBz += n2 * Sz[address_3d(n1, n2 - 1, ez, sdim1, sdim2)];
                            buffer_[ao12 + size * (m_count + 5)] += prefac * sx * sy * DBz;
                            ao12++;
                        }
                    }
                    // 6 derivative components for each multipole components (Ax to Bz above)
                    m_count += 6;
                }
            }
        }
    }
    pure_transform(s1, s2, nchunk_);
    for (int chunk = 0; chunk < nchunk_; ++chunk) {
        buffers_[chunk] = buffer_ + chunk * s1.size() * s2.size();
    }
}
