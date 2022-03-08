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
#include <algorithm>
#include "psi4/libmints/dipole.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libciomr/libciomr.h"

#include <libint2/engine.h>

using namespace psi;

// Initialize overlap_recur_ to +1 basis set angular momentum, +1 on each center is sufficient
// to compute the dipole derivatives
DipoleInt::DipoleInt(std::vector<SphericalTransform> &spherical_transforms, std::shared_ptr<BasisSet> bs1,
                     std::shared_ptr<BasisSet> bs2, int nderiv)
    : OneBodyAOInt(spherical_transforms, bs1, bs2, nderiv), overlap_recur_(bs1->max_am() + 1, bs2->max_am() + 1) {
    // ACS delete these when L2 can provide the correct derivatives
    int maxam1 = bs1_->max_am();
    int maxam2 = bs2_->max_am();
    int maxnao1 = (maxam1 + 1) * (maxam1 + 2) / 2;
    int maxnao2 = (maxam2 + 1) * (maxam2 + 2) / 2;

    // // Increase buffer size to handle x, y, and z components
    // if (deriv_ == 0) {
    //     buffer_ = new double[3 * maxnao1 * maxnao2];
    //     set_chunks(3);
    // } else if (deriv_ == 1) {
    //     natom_ = bs1_->molecule()->natom();
    //     buffer_ = new double[6 * 3 * maxnao1 * maxnao2];
    //     set_chunks(6 * 3);
    // }

    int max_am = std::max(basis1()->max_am(), basis2()->max_am());
    int max_nprim = std::max(basis1()->max_nprimitive(), basis2()->max_nprimitive());

    if (nderiv == 0) {
        set_chunks(3);

        engine0_ =
            std::make_unique<libint2::Engine>(libint2::Operator::emultipole1, max_nprim, max_am, 0);
    } else if (nderiv == 1) {
        // We set chunk count for normalize_am and pure_transform
        set_chunks(18);

        engine0_ =
            std::make_unique<libint2::Engine>(libint2::Operator::emultipole1, max_nprim, max_am, 0);
        engine1_ =
            std::make_unique<libint2::Engine>(libint2::Operator::emultipole1, max_nprim, max_am, 1);
        buffer_ = new double[6 * 3 * maxnao1 * maxnao2];
    }

    //buffer_ = nullptr;
    buffers_.resize(nchunk_);
}

DipoleInt::~DipoleInt() { delete[] buffer_; }

SharedVector DipoleInt::nuclear_contribution(std::shared_ptr<Molecule> mol, const Vector3 &origin) {
    auto sret = std::make_shared<Vector>(3);
    double *ret = sret->pointer();

    for (int i = 0; i < mol->natom(); ++i) {
        Vector3 geom = mol->xyz(i) - origin;
        ret[0] += mol->Z(i) * geom[0];
        ret[1] += mol->Z(i) * geom[1];
        ret[2] += mol->Z(i) * geom[2];
    }

    return sret;
}

SharedMatrix DipoleInt::nuclear_gradient_contribution(std::shared_ptr<Molecule> mol) {
    auto sret = std::make_shared<Matrix>("Nuclear dipole derivative (3Nx3)", 3 * mol->natom(), 3);
    double **ret = sret->pointer();

    for (int i = 0; i < mol->natom(); ++i) {
        ret[3 * i + 0][0] = mol->Z(i);
        ret[3 * i + 1][1] = mol->Z(i);
        ret[3 * i + 2][2] = mol->Z(i);
    }

    return sret;
}

void DipoleInt::compute_pair(const libint2::Shell &s1, const libint2::Shell &s2) {
    engine0_->compute(s1, s2);

    size_t nints = s1.size() * s2.size();
    // Libint gives us the overlap, mu_x, mu_y, mu_z in the buffers.
    // We don't care about the overlap here so we just skip over it.
    for (int chunk = 1; chunk < 4; chunk++) {
        double * ptr = const_cast<double*>(engine0_->results()[chunk]);
        std::transform(ptr, ptr + nints, ptr, [](double val) -> double { return -val; });
        buffers_[chunk - 1] = engine0_->results()[chunk];
    }
}

// This function should replace the function below it when we can link against a libint2 that provides dipole derivatives https://github.com/evaleev/libint/issues/236
//void DipoleInt::compute_pair_deriv1(const libint2::Shell &s1, const libint2::Shell &s2) {
//    engine1_->compute(s1, s2);
//
//    size_t nints = s1.size() * s2.size();
//    for (int i = 0; i < 6; i) {
//        double * ptr = const_cast<double*>(engine1_->results()[4 * i + 1]);
//        std::transform(ptr, ptr + nints, ptr, [](double val) -> double { return -val; });
//        ptr = const_cast<double*>(engine1_->results()[4 * i + 1]);
//        std::transform(ptr, ptr + nints, ptr, [](double val) -> double { return -val; });
//        ptr = const_cast<double*>(engine1_->results()[4 * i + 1]);
//        std::transform(ptr, ptr + nints, ptr, [](double val) -> double { return -val; });
//        buffers_[3 * i + 0] = engine1_->results()[4 * i + 1];
//        buffers_[3 * i + 1] = engine1_->results()[4 * i + 2];
//        buffers_[3 * i + 2] = engine1_->results()[4 * i + 3];
//    }
//}

// The engine only supports segmented basis sets
void DipoleInt::compute_pair_deriv1(const libint2::Shell &s1, const libint2::Shell &s2) {
    int ao12;
    int am1 = s1.contr[0].l;
    int am2 = s2.contr[0].l;
    int nprim1 = s1.nprim();
    int nprim2 = s2.nprim();
    size_t length = INT_NCART(am1) * INT_NCART(am2);
    double A[3], B[3];

    A[0] = s1.O[0];
    A[1] = s1.O[1];
    A[2] = s1.O[2];
    B[0] = s2.O[0];
    B[1] = s2.O[1];
    B[2] = s2.O[2];

    size_t xaxdisp = 0;
    size_t xaydisp = 1 * length;
    size_t xazdisp = 2 * length;
    size_t xbxdisp = 3 * length;
    size_t xbydisp = 4 * length;
    size_t xbzdisp = 5 * length;
    size_t yaxdisp = 6 * length;
    size_t yaydisp = 7 * length;
    size_t yazdisp = 8 * length;
    size_t ybxdisp = 9 * length;
    size_t ybydisp = 10 * length;
    size_t ybzdisp = 11 * length;
    size_t zaxdisp = 12 * length;
    size_t zaydisp = 13 * length;
    size_t zazdisp = 14 * length;
    size_t zbxdisp = 15 * length;
    size_t zbydisp = 16 * length;
    size_t zbzdisp = 17 * length;

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    // Zero out the buffer
    memset(buffer_, 0, 6 * 3 * length * sizeof(double));

    double **x = overlap_recur_.x();
    double **y = overlap_recur_.y();
    double **z = overlap_recur_.z();
    double v1, v2, v3, v4;  // temporary value storage

    for (int p1 = 0; p1 < nprim1; ++p1) {
        double a1 = s1.alpha[p1];
        double c1 = s1.contr[0].coeff[p1];
        for (int p2 = 0; p2 < nprim2; ++p2) {
            double a2 = s2.alpha[p2];
            double c2 = s2.contr[0].coeff[p2];
            double gamma = a1 + a2;
            double oog = 1.0 / gamma;

            double PA[3], PB[3];
            double P[3];

            P[0] = (a1 * A[0] + a2 * B[0]) * oog;
            P[1] = (a1 * A[1] + a2 * B[1]) * oog;
            P[2] = (a1 * A[2] + a2 * B[2]) * oog;
            PA[0] = P[0] - A[0];
            PA[1] = P[1] - A[1];
            PA[2] = P[2] - A[2];
            PB[0] = P[0] - B[0];
            PB[1] = P[1] - B[1];
            PB[2] = P[2] - B[2];

            double over_pf = exp(-a1 * a2 * AB2 * oog) * sqrt(M_PI * oog) * M_PI * oog * c1 * c2;

            // Do recursion, this is sufficient information to compute dipole derivatives
            overlap_recur_.compute(PA, PB, gamma, am1 + 1, am2 + 1);

            ao12 = 0;
            for (int ii = 0; ii <= am1; ii++) {
                int l1 = am1 - ii;
                for (int jj = 0; jj <= ii; jj++) {
                    int m1 = ii - jj;
                    int n1 = jj;
                    /*--- create all am components of sj ---*/
                    for (int kk = 0; kk <= am2; kk++) {
                        int l2 = am2 - kk;
                        for (int ll = 0; ll <= kk; ll++) {
                            int m2 = kk - ll;
                            int n2 = ll;

                            // mu-x derivatives
                            v1 = v2 = v3 = v4 = 0.0;

                            //
                            // A derivatives with mu-x
                            //

                            // (a+1_x|mx|b+1_x)
                            v1 = x[l1 + 1][l2 + 1] * y[m1][m2] * z[n1][n2];
                            // (a+1_x|mx|b)
                            v2 = x[l1 + 1][l2] * y[m1][m2] * z[n1][n2];
                            if (l1) {
                                // (a-1_x|mx|b+1_x)
                                v3 = x[l1 - 1][l2 + 1] * y[m1][m2] * z[n1][n2];
                                // (a-1_x|mx|b)
                                v4 = x[l1 - 1][l2] * y[m1][m2] * z[n1][n2];
                            }
                            buffer_[ao12 + xaxdisp] -= (2.0 * a1 * (v1 + B[0] * v2) - l1 * (v3 + B[0] * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_y|mx|b+1_x)
                            v1 = x[l1][l2 + 1] * y[m1 + 1][m2] * z[n1][n2];
                            // (a+1_y|mx|b)
                            v2 = x[l1][l2] * y[m1 + 1][m2] * z[n1][n2];
                            if (m1) {
                                // (a-1_y|mx|b+1_x)
                                v3 = x[l1][l2 + 1] * y[m1 - 1][m2] * z[n1][n2];
                                // (a-1_y|mx|b)
                                v4 = x[l1][l2] * y[m1 - 1][m2] * z[n1][n2];
                            }
                            buffer_[ao12 + xaydisp] -= (2.0 * a1 * (v1 + B[0] * v2) - m1 * (v3 + B[0] * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_z|mx|b+1_x)
                            v1 = x[l1][l2 + 1] * y[m1][m2] * z[n1 + 1][n2];
                            // (a+1_z|mx|b)
                            v2 = x[l1][l2] * y[m1][m2] * z[n1 + 1][n2];
                            if (n1) {
                                // (a-1_z|mx|b+1_x)
                                v3 = x[l1][l2 + 1] * y[m1][m2] * z[n1 - 1][n2];
                                // (a-1_z|mx|b)
                                v4 = x[l1][l2] * y[m1][m2] * z[n1 - 1][n2];
                            }
                            buffer_[ao12 + xazdisp] -= (2.0 * a1 * (v1 + B[0] * v2) - n1 * (v3 + B[0] * v4)) * over_pf;

                            //
                            // B derivatives with mu-x
                            //

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_x|mx|b+1_x)
                            v1 = x[l1 + 1][l2 + 1] * y[m1][m2] * z[n1][n2];
                            // (a|mx|b+1_x)
                            v2 = x[l1][l2 + 1] * y[m1][m2] * z[n1][n2];
                            if (l2) {
                                // (a+1_x|mx|b-1_x)
                                v3 = x[l1 + 1][l2 - 1] * y[m1][m2] * z[n1][n2];
                                // (a|mx|b-1_x)
                                v4 = x[l1][l2 - 1] * y[m1][m2] * z[n1][n2];
                            }
                            buffer_[ao12 + xbxdisp] -= (2.0 * a2 * (v1 + A[0] * v2) - l2 * (v3 + A[0] * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_x|mx|b+1_y)
                            v1 = x[l1 + 1][l2] * y[m1][m2 + 1] * z[n1][n2];
                            // (a|mx|b+1_y)
                            v2 = x[l1][l2] * y[m1][m2 + 1] * z[n1][n2];
                            if (m2) {
                                // (a+1_x|mx|b-1_y)
                                v3 = x[l1 + 1][l2] * y[m1][m2 - 1] * z[n1][n2];
                                // (a|mx|b-1_y)
                                v4 = x[l1][l2] * y[m1][m2 - 1] * z[n1][n2];
                            }
                            buffer_[ao12 + xbydisp] -= (2.0 * a2 * (v1 + A[0] * v2) - m2 * (v3 + A[0] * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_x|mx|b+1_z)
                            v1 = x[l1 + 1][l2] * y[m1][m2] * z[n1][n2 + 1];
                            // (a|mx|b+1_z)
                            v2 = x[l1][l2] * y[m1][m2] * z[n1][n2 + 1];
                            if (n2) {
                                // (a+1_x|mx|b-1_z)
                                v3 = x[l1 + 1][l2] * y[m1][m2] * z[n1][n2 - 1];
                                // (a|mx|b-1_z)
                                v4 = x[l1][l2] * y[m1][m2] * z[n1][n2 - 1];
                            }
                            buffer_[ao12 + xbzdisp] -= (2.0 * a2 * (v1 + A[0] * v2) - n2 * (v3 + A[0] * v4)) * over_pf;

                            //
                            // A derivatives with mu-y
                            //

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_x|my|b+1_y)
                            v1 = x[l1 + 1][l2] * y[m1][m2 + 1] * z[n1][n2];
                            // (a+1_x|my|b)
                            v2 = x[l1 + 1][l2] * y[m1][m2] * z[n1][n2];
                            if (l1) {
                                // (a-1_x|my|b+1_y)
                                v3 = x[l1 - 1][l2] * y[m1][m2 + 1] * z[n1][n2];
                                // (a-1_x|my|b)
                                v4 = x[l1 - 1][l2] * y[m1][m2] * z[n1][n2];
                            }
                            buffer_[ao12 + yaxdisp] -= (2.0 * a1 * (v1 + B[1] * v2) - l1 * (v3 + B[1] * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_y|my|b+1_y)
                            v1 = x[l1][l2] * y[m1 + 1][m2 + 1] * z[n1][n2];
                            // (a+1_y|my|b)
                            v2 = x[l1][l2] * y[m1 + 1][m2] * z[n1][n2];
                            if (m1) {
                                // (a-1_y|my|b+1_y)
                                v3 = x[l1][l2] * y[m1 - 1][m2 + 1] * z[n1][n2];
                                // (a-1_y|my|b)
                                v4 = x[l1][l2] * y[m1 - 1][m2] * z[n1][n2];
                            }
                            buffer_[ao12 + yaydisp] -= (2.0 * a1 * (v1 + B[1] * v2) - m1 * (v3 + B[1] * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_z|my|b+1_y)
                            v1 = x[l1][l2] * y[m1][m2 + 1] * z[n1 + 1][n2];
                            // (a+1_z|my|b)
                            v2 = x[l1][l2] * y[m1][m2] * z[n1 + 1][n2];
                            if (n1) {
                                // (a-1_z|my|b+1_y)
                                v3 = x[l1][l2] * y[m1][m2 + 1] * z[n1 - 1][n2];
                                // (a-1_z|my|b)
                                v4 = x[l1][l2] * y[m1][m2] * z[n1 - 1][n2];
                            }
                            buffer_[ao12 + yazdisp] -= (2.0 * a1 * (v1 + B[1] * v2) - n1 * (v3 + B[1] * v4)) * over_pf;

                            //
                            // B derivatives with mu-y
                            //

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_y|my|b+1_x)
                            v1 = x[l1][l2 + 1] * y[m1 + 1][m2] * z[n1][n2];
                            // (a|my|b+1_x)
                            v2 = x[l1][l2 + 1] * y[m1][m2] * z[n1][n2];
                            if (l2) {
                                // (a+1_y|my|b-1_x)
                                v3 = x[l1][l2 - 1] * y[m1 + 1][m2] * z[n1][n2];
                                // (a|my|b-1_x)
                                v4 = x[l1][l2 - 1] * y[m1][m2] * z[n1][n2];
                            }
                            buffer_[ao12 + ybxdisp] -= (2.0 * a2 * (v1 + A[1] * v2) - l2 * (v3 + A[1] * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_y|my|b+1_y)
                            v1 = x[l1][l2] * y[m1 + 1][m2 + 1] * z[n1][n2];
                            // (a|my|b+1_y)
                            v2 = x[l1][l2] * y[m1][m2 + 1] * z[n1][n2];
                            if (m2) {
                                // (a+1_y|my|b-1_y)
                                v3 = x[l1][l2] * y[m1 + 1][m2 - 1] * z[n1][n2];
                                // (a|my|b-1_y)
                                v4 = x[l1][l2] * y[m1][m2 - 1] * z[n1][n2];
                            }
                            buffer_[ao12 + ybydisp] -= (2.0 * a2 * (v1 + A[1] * v2) - m2 * (v3 + A[1] * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_y|my|b+1_z)
                            v1 = x[l1][l2] * y[m1 + 1][m2] * z[n1][n2 + 1];
                            // (a|my|b+1_z)
                            v2 = x[l1][l2] * y[m1][m2] * z[n1][n2 + 1];
                            if (n2) {
                                // (a+1_y|my|b-1_z)
                                v3 = x[l1][l2] * y[m1 + 1][m2] * z[n1][n2 - 1];
                                // (a|my|b-1_z)
                                v4 = x[l1][l2] * y[m1][m2] * z[n1][n2 - 1];
                            }
                            buffer_[ao12 + ybzdisp] -= (2.0 * a2 * (v1 + A[1] * v2) - n2 * (v3 + A[1] * v4)) * over_pf;

                            //
                            // A derivatives with mu-z
                            //

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_x|mz|b+1_z)
                            v1 = x[l1 + 1][l2] * y[m1][m2] * z[n1][n2 + 1];
                            // (a+1_x|mz|b)
                            v2 = x[l1 + 1][l2] * y[m1][m2] * z[n1][n2];
                            if (l1) {
                                // (a-1_x|mz|b+1_z)
                                v3 = x[l1 - 1][l2] * y[m1][m2] * z[n1][n2 + 1];
                                // (a-1_x|mz|b)
                                v4 = x[l1 - 1][l2] * y[m1][m2] * z[n1][n2];
                            }
                            buffer_[ao12 + zaxdisp] -= (2.0 * a1 * (v1 + B[2] * v2) - l1 * (v3 + B[2] * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_y|mz|b+1_z)
                            v1 = x[l1][l2] * y[m1 + 1][m2] * z[n1][n2 + 1];
                            // (a+1_y|mz|b)
                            v2 = x[l1][l2] * y[m1 + 1][m2] * z[n1][n2];
                            if (m1) {
                                // (a-1_y|mz|b+1_z)
                                v3 = x[l1][l2] * y[m1 - 1][m2] * z[n1][n2 + 1];
                                // (a-1_y|mz|b)
                                v4 = x[l1][l2] * y[m1 - 1][m2] * z[n1][n2];
                            }
                            buffer_[ao12 + zaydisp] -= (2.0 * a1 * (v1 + B[2] * v2) - m1 * (v3 + B[2] * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_z|mz|b+1_z)
                            v1 = x[l1][l2] * y[m1][m2] * z[n1 + 1][n2 + 1];
                            // (a+1_z|mz|b)
                            v2 = x[l1][l2] * y[m1][m2] * z[n1 + 1][n2];
                            if (n1) {
                                // (a-1_z|mz|b+1_z)
                                v3 = x[l1][l2] * y[m1][m2] * z[n1 - 1][n2 + 1];
                                // (a-1_z|mz|b)
                                v4 = x[l1][l2] * y[m1][m2] * z[n1 - 1][n2];
                            }
                            buffer_[ao12 + zazdisp] -= (2.0 * a1 * (v1 + B[2] * v2) - n1 * (v3 + B[2] * v4)) * over_pf;

                            //
                            // B derivates with mu-z
                            //

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_z|mz|b+1_x)
                            v1 = x[l1][l2 + 1] * y[m1][m2] * z[n1 + 1][n2];
                            // (a|mz|b+1_x)
                            v2 = x[l1][l2 + 1] * y[m1][m2] * z[n1][n2];
                            if (l2) {
                                // (a+1_z|mz|b-1_x)
                                v3 = x[l1][l2 - 1] * y[m1][m2] * z[n1 + 1][n2];
                                // (a|mz|b-1_x)
                                v4 = x[l1][l2 - 1] * y[m1][m2] * z[n1][n2];
                            }
                            buffer_[ao12 + zbxdisp] -= (2.0 * a2 * (v1 + A[2] * v2) - l2 * (v3 + A[2] * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_z|mz|b+1_y)
                            v1 = x[l1][l2] * y[m1][m2 + 1] * z[n1 + 1][n2];
                            // (a|mz|b+1_y)
                            v2 = x[l1][l2] * y[m1][m2 + 1] * z[n1][n2];
                            if (m2) {
                                // (a+1_z|mz|b-1_y)
                                v3 = x[l1][l2] * y[m1][m2 - 1] * z[n1 + 1][n2];
                                // (a|mz|b-1_y)
                                v4 = x[l1][l2] * y[m1][m2 - 1] * z[n1][n2];
                            }
                            buffer_[ao12 + zbydisp] -= (2.0 * a2 * (v1 + A[2] * v2) - m2 * (v3 + A[2] * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_z|mz|b+1_z)
                            v1 = x[l1][l2] * y[m1][m2] * z[n1 + 1][n2 + 1];
                            // (a|mz|b+1_z)
                            v2 = x[l1][l2] * y[m1][m2] * z[n1][n2 + 1];
                            if (n2) {
                                // (a+1_z|mz|b-1_z)
                                v3 = x[l1][l2] * y[m1][m2] * z[n1 + 1][n2 - 1];
                                // (a|mz|b-1_z)
                                v4 = x[l1][l2] * y[m1][m2] * z[n1][n2 - 1];
                            }
                            buffer_[ao12 + zbzdisp] -= (2.0 * a2 * (v1 + A[2] * v2) - n2 * (v3 + A[2] * v4)) * over_pf;

                            ao12++;
                        }
                    }
                }
            }
        }
    }

    pure_transform(s1, s2, 18);
    for (int i = 0; i < 18; ++i) {
        buffers_[i] = buffer_ + i * s1.size() * s2.size();
    }
}
