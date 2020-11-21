/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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

;
using namespace psi;

#define have_moment(k, max_k) (k <= max_k);

namespace {
int number_of_chunks(int max_k) {
    int nc = 0;
    for (size_t i = 0; i <= max_k; i++) {
        nc += (i + 1) * (i + 2) / 2;
    }
    return nc;
}
}  // namespace

MultipolePotentialInt::MultipolePotentialInt(std::vector<SphericalTransform> &spherical_transforms,
                                             std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2, int max_k,
                                             int deriv)
    : OneBodyAOInt(spherical_transforms, bs1, bs2, deriv),
      max_k_(max_k) {

    if (deriv == 0) {
        mvi_recur_ = new ObaraSaikaTwoCenterMultipolePotentialRecursion(bs1->max_am(), bs2->max_am(), max_k);
    } else if (deriv == 1) {
        mvi_recur_ = new ObaraSaikaTwoCenterMultipolePotentialRecursion(bs1->max_am() + 1, bs2->max_am() + 1, max_k);
    } else {
        throw FeatureNotImplemented("LibMints", "MultipolePotentialInts called with deriv > 1", __FILE__, __LINE__);
     }
    int maxam1 = bs1_->max_am();
    int maxam2 = bs2_->max_am();

    int maxnao1 = INT_NCART(maxam1);
    int maxnao2 = INT_NCART(maxam2);

    if (max_k > 3) {
        throw FeatureNotImplemented("LibMints", "MultipolePotentialInts called with max_k > 3", __FILE__, __LINE__);
    }

    int nchunks = number_of_chunks(max_k_);
    if (deriv == 1) {
        set_chunks(3 * natom_ * nchunks);
        maxnao1 *= 3 * natom_;
    } else {
        set_chunks(nchunks);
    }

    buffer_ = new double[nchunks * maxnao1 * maxnao2];
}

MultipolePotentialInt::~MultipolePotentialInt() { delete[] buffer_; }

void MultipolePotentialInt::compute_pair(const GaussianShell &s1, const GaussianShell &s2) {
    int ao12;
    int am1 = s1.am();
    int am2 = s2.am();
    int nprim1 = s1.nprimitive();
    int nprim2 = s2.nprimitive();
    double A[3], B[3];
    A[0] = s1.center()[0];
    A[1] = s1.center()[1];
    A[2] = s1.center()[2];
    B[0] = s2.center()[0];
    B[1] = s2.center()[1];
    B[2] = s2.center()[2];

    int izm = 1;
    int iym = am1 + 1;
    int ixm = iym * iym;
    int jzm = 1;
    int jym = am2 + 1;
    int jxm = jym * jym;

    // Not sure if these are needed.
    int size = INT_NCART(am1) * INT_NCART(am2);

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    size_t n_elements = static_cast<size_t>(number_of_chunks(max_k_) * size);
    memset(buffer_, 0, n_elements * sizeof(double));

    bool have_dipole = have_moment(1, max_k_);
    bool have_quadrupole = have_moment(2, max_k_);
    bool have_octupole = have_moment(3, max_k_);

    // Charge
    double ***q = mvi_recur_->q();
    // Dipole
    double ***x = mvi_recur_->x();
    double ***y = mvi_recur_->y();
    double ***z = mvi_recur_->z();
    // Quadrupole
    double ***xx = mvi_recur_->xx();
    double ***yy = mvi_recur_->yy();
    double ***zz = mvi_recur_->zz();
    double ***xy = mvi_recur_->xy();
    double ***xz = mvi_recur_->xz();
    double ***yz = mvi_recur_->yz();

    // Octupole
    double ***xxx = mvi_recur_->xxx();
    double ***yyy = mvi_recur_->yyy();
    double ***zzz = mvi_recur_->zzz();
    double ***xxy = mvi_recur_->xxy();
    double ***xxz = mvi_recur_->xxz();
    double ***xyy = mvi_recur_->xyy();
    double ***yyz = mvi_recur_->yyz();
    double ***xzz = mvi_recur_->xzz();
    double ***yzz = mvi_recur_->yzz();
    double ***xyz = mvi_recur_->xyz();

    double Cx = origin_[0];
    double Cy = origin_[1];
    double Cz = origin_[2];

    for (int p1 = 0; p1 < nprim1; ++p1) {
        double a1 = s1.exp(p1);
        double c1 = s1.coef(p1);
        for (int p2 = 0; p2 < nprim2; ++p2) {
            double a2 = s2.exp(p2);
            double c2 = s2.coef(p2);
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
            double PC[3];

            PC[0] = P[0] - Cx;
            PC[1] = P[1] - Cy;
            PC[2] = P[2] - Cz;

            // Get recursive
            mvi_recur_->compute(PA, PB, PC, gamma, am1, am2);

            // Gather contributions.
            ao12 = 0;
            for (int ii = 0; ii <= am1; ++ii) {
                int l1 = am1 - ii;
                for (int jj = 0; jj <= ii; ++jj) {
                    int m1 = ii - jj;
                    int n1 = jj;

                    for (int kk = 0; kk <= am2; ++kk) {
                        int l2 = am2 - kk;
                        for (int ll = 0; ll <= kk; ++ll) {
                            int m2 = kk - ll;
                            int n2 = ll;

                            // Compute location in the recursion
                            int iind = l1 * ixm + m1 * iym + n1 * izm;
                            int jind = l2 * jxm + m2 * jym + n2 * jzm;

                            buffer_[ao12 + 0 * size] += q[iind][jind][0] * over_pf;
                            if (have_dipole) {
                                buffer_[ao12 + 1 * size] += x[iind][jind][0] * over_pf;
                                buffer_[ao12 + 2 * size] += y[iind][jind][0] * over_pf;
                                buffer_[ao12 + 3 * size] += z[iind][jind][0] * over_pf;
                            }
                            if (have_quadrupole) {
                                buffer_[ao12 + 4 * size] += xx[iind][jind][0] * over_pf;
                                buffer_[ao12 + 5 * size] += xy[iind][jind][0] * over_pf;
                                buffer_[ao12 + 6 * size] += xz[iind][jind][0] * over_pf;
                                buffer_[ao12 + 7 * size] += yy[iind][jind][0] * over_pf;
                                buffer_[ao12 + 8 * size] += yz[iind][jind][0] * over_pf;
                                buffer_[ao12 + 9 * size] += zz[iind][jind][0] * over_pf;
                            }
                            if (have_octupole) {
                                buffer_[ao12 + 10 * size] += xxx[iind][jind][0] * over_pf;
                                buffer_[ao12 + 11 * size] += xxy[iind][jind][0] * over_pf;
                                buffer_[ao12 + 12 * size] += xxz[iind][jind][0] * over_pf;
                                buffer_[ao12 + 13 * size] += xyy[iind][jind][0] * over_pf;
                                buffer_[ao12 + 14 * size] += xyz[iind][jind][0] * over_pf;
                                buffer_[ao12 + 15 * size] += xzz[iind][jind][0] * over_pf;
                                buffer_[ao12 + 16 * size] += yyy[iind][jind][0] * over_pf;
                                buffer_[ao12 + 17 * size] += yyz[iind][jind][0] * over_pf;
                                buffer_[ao12 + 18 * size] += yzz[iind][jind][0] * over_pf;
                                buffer_[ao12 + 19 * size] += zzz[iind][jind][0] * over_pf;
                            }
                            ao12++;
                        }
                    }
                }
            }
        }
    }
}

void MultipolePotentialInt::compute_pair_deriv1(const GaussianShell &s1, const GaussianShell &s2) {
    int ao12;
    int am1 = s1.am();
    int am2 = s2.am();
    int nprim1 = s1.nprimitive();
    int nprim2 = s2.nprimitive();
    const int ncenteri = s1.ncenter();
    const int ncenterj = s2.ncenter();
    double A[3], B[3];
    A[0] = s1.center()[0];
    A[1] = s1.center()[1];
    A[2] = s1.center()[2];
    B[0] = s2.center()[0];
    B[1] = s2.center()[1];
    B[2] = s2.center()[2];

    int izm1 = 1;
    int iym1 = am1 + 1 + 1;
    int ixm1 = iym1 * iym1;
    int jzm1 = 1;
    int jym1 = am2 + 1 + 1;
    int jxm1 = jym1 * jym1;

    // Not sure if these are needed.
    const size_t size = s1.ncartesian() * s2.ncartesian();
    const int center_i = ncenteri * 3 * size;
    const int center_j = ncenterj * 3 * size;

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    size_t n_elements = static_cast<size_t>(number_of_chunks(max_k_) * size);
    memset(buffer_, 0, 3 * natom_ * n_elements * sizeof(double));

    bool have_dipole = have_moment(1, max_k_);
    bool have_quadrupole = have_moment(2, max_k_);
    bool have_octupole = have_moment(3, max_k_);

    // Charge
    double ***q = mvi_recur_->q();
    // Dipole
    double ***x = mvi_recur_->x();
    double ***y = mvi_recur_->y();
    double ***z = mvi_recur_->z();
    // Quadrupole
    double ***xx = mvi_recur_->xx();
    double ***yy = mvi_recur_->yy();
    double ***zz = mvi_recur_->zz();
    double ***xy = mvi_recur_->xy();
    double ***xz = mvi_recur_->xz();
    double ***yz = mvi_recur_->yz();

    // Octupole
    double ***xxx = mvi_recur_->xxx();
    double ***yyy = mvi_recur_->yyy();
    double ***zzz = mvi_recur_->zzz();
    double ***xxy = mvi_recur_->xxy();
    double ***xxz = mvi_recur_->xxz();
    double ***xyy = mvi_recur_->xyy();
    double ***yyz = mvi_recur_->yyz();
    double ***xzz = mvi_recur_->xzz();
    double ***yzz = mvi_recur_->yzz();
    double ***xyz = mvi_recur_->xyz();

    double Cx = origin_[0];
    double Cy = origin_[1];
    double Cz = origin_[2];

    for (int p1 = 0; p1 < nprim1; ++p1) {
        double a1 = s1.exp(p1);
        double c1 = s1.coef(p1);
        for (int p2 = 0; p2 < nprim2; ++p2) {
            double a2 = s2.exp(p2);
            double c2 = s2.coef(p2);
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
            double PC[3];

            PC[0] = P[0] - Cx;
            PC[1] = P[1] - Cy;
            PC[2] = P[2] - Cz;

            // Get recursive
            mvi_recur_->compute(PA, PB, PC, gamma, am1 + 1, am2 + 1);

            // Gather contributions.
            ao12 = 0;
            for (int ii = 0; ii <= am1; ++ii) {
                int l1 = am1 - ii;
                for (int jj = 0; jj <= ii; ++jj) {
                    int m1 = ii - jj;
                    int n1 = jj;

                    for (int kk = 0; kk <= am2; ++kk) {
                        int l2 = am2 - kk;
                        for (int ll = 0; ll <= kk; ++ll) {
                            int m2 = kk - ll;
                            int n2 = ll;

                            // Compute location in the recursion
                            int iind = l1 * ixm1 + m1 * iym1 + n1 * izm1;
                            int jind = l2 * jxm1 + m2 * jym1 + n2 * jzm1;
                            
                            double Z = moments_[0];
                            double pfac = over_pf * Z;

                            // x
                            double temp = 2.0 * a1 * q[iind + ixm1][jind][0];
                            if (l1) temp -= l1 * q[iind - ixm1][jind][0];
                            buffer_[center_i + (0 * size) + ao12] -= temp * pfac;

                            temp = 2.0 * a2 * q[iind][jind + jxm1][0];
                            if (l2) temp -= l2 * q[iind][jind - jxm1][0];
                            buffer_[center_j + (0 * size) + ao12] -= temp * pfac;

                            // y
                            temp = 2.0 * a1 * q[iind + iym1][jind][0];
                            if (m1) temp -= m1 * q[iind - iym1][jind][0];
                            buffer_[center_i + (1 * size) + ao12] -= temp * pfac;

                            temp = 2.0 * a2 * q[iind][jind + jym1][0];
                            if (m2) temp -= m2 * q[iind][jind - jym1][0];
                            buffer_[center_j + (1 * size) + ao12] -= temp * pfac;

                            // z
                            temp = 2.0 * a1 * q[iind + izm1][jind][0];
                            if (n1) temp -= n1 * q[iind - izm1][jind][0];
                            buffer_[center_i + (2 * size) + ao12] -= temp * pfac;

                            temp = 2.0 * a2 * q[iind][jind + jzm1][0];
                            if (n2) temp -= n2 * q[iind][jind - jzm1][0];
                            buffer_[center_j + (2 * size) + ao12] -= temp * pfac;

                            int off;
                            if (have_dipole) {
                                off = 1;
                                const std::vector<double ***> dips{x, y, z};
                                for (int cc = 0; cc < 3; ++cc) {
                                    double*** dip = dips[cc];
                                    // x
                                    temp = 2.0 * a1 * dip[iind + ixm1][jind][0];
                                    if (l1) temp -= l1 * dip[iind - ixm1][jind][0];
                                    buffer_[center_i + (0 * size) + ao12] -= temp * over_pf * moments_[cc + off];

                                    temp = 2.0 * a2 * dip[iind][jind + jxm1][0];
                                    if (l2) temp -= l2 * dip[iind][jind - jxm1][0];
                                    buffer_[center_j + (0 * size) + ao12] -= temp * over_pf * moments_[cc + off];

                                    // y
                                    temp = 2.0 * a1 * dip[iind + iym1][jind][0];
                                    if (m1) temp -= m1 * dip[iind - iym1][jind][0];
                                    buffer_[center_i + (1 * size) + ao12] -= temp * over_pf * moments_[cc + off];

                                    temp = 2.0 * a2 * dip[iind][jind + jym1][0];
                                    if (m2) temp -= m2 * dip[iind][jind - jym1][0];
                                    buffer_[center_j + (1 * size) + ao12] -= temp * over_pf * moments_[cc + off];

                                    // z
                                    temp = 2.0 * a1 * dip[iind + izm1][jind][0];
                                    if (n1) temp -= n1 * dip[iind - izm1][jind][0];
                                    buffer_[center_i + (2 * size) + ao12] -= temp * over_pf * moments_[cc + off];

                                    temp = 2.0 * a2 * dip[iind][jind + jzm1][0];
                                    if (n2) temp -= n2 * dip[iind][jind - jzm1][0];
                                    buffer_[center_j + (2 * size) + ao12] -= temp * over_pf * moments_[cc + off];
                                }
                            }
                            // could also be done in one go (charges, dipoles, quadrupoles etc.)
                            if (have_quadrupole) {
                                off = 4;
                                const std::vector<double ***> quads{xx, xy, xz, yy, yz, zz};
                                for (int cc = 0; cc < 6; ++cc) {
                                    double*** quad = quads[cc];
                                    // x
                                    temp = 2.0 * a1 * quad[iind + ixm1][jind][0];
                                    if (l1) temp -= l1 * quad[iind - ixm1][jind][0];
                                    buffer_[center_i + (0 * size) + ao12] -= temp * over_pf * moments_[cc + off];

                                    temp = 2.0 * a2 * quad[iind][jind + jxm1][0];
                                    if (l2) temp -= l2 * quad[iind][jind - jxm1][0];
                                    buffer_[center_j + (0 * size) + ao12] -= temp * over_pf * moments_[cc + off];

                                    // y
                                    temp = 2.0 * a1 * quad[iind + iym1][jind][0];
                                    if (m1) temp -= m1 * quad[iind - iym1][jind][0];
                                    buffer_[center_i + (1 * size) + ao12] -= temp * over_pf * moments_[cc + off];

                                    temp = 2.0 * a2 * quad[iind][jind + jym1][0];
                                    if (m2) temp -= m2 * quad[iind][jind - jym1][0];
                                    buffer_[center_j + (1 * size) + ao12] -= temp * over_pf * moments_[cc + off];

                                    // z
                                    temp = 2.0 * a1 * quad[iind + izm1][jind][0];
                                    if (n1) temp -= n1 * quad[iind - izm1][jind][0];
                                    buffer_[center_i + (2 * size) + ao12] -= temp * over_pf * moments_[cc + off];

                                    temp = 2.0 * a2 * quad[iind][jind + jzm1][0];
                                    if (n2) temp -= n2 * quad[iind][jind - jzm1][0];
                                    buffer_[center_j + (2 * size) + ao12] -= temp * over_pf * moments_[cc + off];
                                }
                            }
                            ao12++;
                        }
                    }
                }
            }
        }
    }
}


void MultipolePotentialInt::compute_shell_deriv1(int sh1, int sh2) {
    const GaussianShell &s1 = bs1_->shell(sh1);
    const GaussianShell &s2 = bs2_->shell(sh2);

    // Call the child's compute_pair method, results better be in buffer_.
    compute_pair_deriv1(s1, s2);

    // Normalize for angular momentum
    normalize_am(s1, s2, nchunk_);

    // Pure angular momentum (6d->5d, ...) transformation
    if (!force_cartesian_) pure_transform(s1, s2, nchunk_);
}