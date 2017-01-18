/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/libmints/efpmultipolepotential.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))

;
using namespace psi;
using namespace std;

EFPMultipolePotentialInt::EFPMultipolePotentialInt(vector<SphericalTransform>& spherical_transforms, std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2, int nderiv) :
    OneBodyAOInt(spherical_transforms, bs1, bs2, nderiv),
    mvi_recur_(bs1->max_am(), bs2->max_am())
{
    int maxam1 = bs1_->max_am();
    int maxam2 = bs2_->max_am();

    int maxnao1 = INT_NCART(maxam1);
    int maxnao2 = INT_NCART(maxam2);

    if (nderiv == 0) {
        buffer_ = new double[20*maxnao1*maxnao2];
        set_chunks(20);
    }
    else
        throw FeatureNotImplemented("LibMints", "MultipolePotentialInts called with deriv > 0",  __FILE__, __LINE__);
}

EFPMultipolePotentialInt::~EFPMultipolePotentialInt()
{
    delete[] buffer_;
}


void EFPMultipolePotentialInt::compute_pair(const GaussianShell& s1,
                                         const GaussianShell& s2)
{
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
    int size =  INT_NCART(am1) * INT_NCART(am2);

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    memset(buffer_, 0, 20 * size * sizeof(double));

    double ***q   = mvi_recur_.q();
    double ***x   = mvi_recur_.x();
    double ***y   = mvi_recur_.y();
    double ***z   = mvi_recur_.z();
    double ***xx  = mvi_recur_.xx();
    double ***yy  = mvi_recur_.yy();
    double ***zz  = mvi_recur_.zz();
    double ***xy  = mvi_recur_.xy();
    double ***xz  = mvi_recur_.xz();
    double ***yz  = mvi_recur_.yz();
    double ***xxx = mvi_recur_.xxx();
    double ***yyy = mvi_recur_.yyy();
    double ***zzz = mvi_recur_.zzz();
    double ***xxy = mvi_recur_.xxy();
    double ***xxz = mvi_recur_.xxz();
    double ***xyy = mvi_recur_.xyy();
    double ***yyz = mvi_recur_.yyz();
    double ***xzz = mvi_recur_.xzz();
    double ***yzz = mvi_recur_.yzz();
    double ***xyz = mvi_recur_.xyz();

    double Cx = origin_[0];
    double Cy = origin_[1];
    double Cz = origin_[2];

    for (int p1=0; p1<nprim1; ++p1) {
        double a1 = s1.exp(p1);
        double c1 = s1.coef(p1);
        for (int p2=0; p2<nprim2; ++p2) {
            double a2 = s2.exp(p2);
            double c2 = s2.coef(p2);
            double gamma = a1 + a2;
            double oog = 1.0 / gamma;

            double PA[3], PB[3];
            double P[3];

            P[0] = (a1*A[0] + a2*B[0])*oog;
            P[1] = (a1*A[1] + a2*B[1])*oog;
            P[2] = (a1*A[2] + a2*B[2])*oog;
            PA[0] = P[0] - A[0];
            PA[1] = P[1] - A[1];
            PA[2] = P[2] - A[2];
            PB[0] = P[0] - B[0];
            PB[1] = P[1] - B[1];
            PB[2] = P[2] - B[2];

            double over_pf = exp(-a1*a2*AB2*oog) * sqrt(M_PI*oog) * M_PI * oog * c1 * c2;
            double PC[3];

            PC[0] = P[0] - Cx;
            PC[1] = P[1] - Cy;
            PC[2] = P[2] - Cz;

            // Get recursive
            mvi_recur_.compute(PA, PB, PC, gamma, am1, am2);

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

                            buffer_[ao12 + 0 * size]  += q[iind][jind][0] * over_pf;
                            buffer_[ao12 + 1 * size]  += x[iind][jind][0] * over_pf;
                            buffer_[ao12 + 2 * size]  += y[iind][jind][0] * over_pf;
                            buffer_[ao12 + 3 * size]  += z[iind][jind][0] * over_pf;
                            buffer_[ao12 + 4 * size]  += xx[iind][jind][0] * over_pf;
                            buffer_[ao12 + 5 * size]  += yy[iind][jind][0] * over_pf;
                            buffer_[ao12 + 6 * size]  += zz[iind][jind][0] * over_pf;
                            buffer_[ao12 + 7 * size]  += xy[iind][jind][0] * over_pf;
                            buffer_[ao12 + 8 * size]  += xz[iind][jind][0] * over_pf;
                            buffer_[ao12 + 9 * size]  += yz[iind][jind][0] * over_pf;
                            buffer_[ao12 + 10 * size] += xxx[iind][jind][0] * over_pf;
                            buffer_[ao12 + 11 * size] += yyy[iind][jind][0] * over_pf;
                            buffer_[ao12 + 12 * size] += zzz[iind][jind][0] * over_pf;
                            buffer_[ao12 + 13 * size] += xxy[iind][jind][0] * over_pf;
                            buffer_[ao12 + 14 * size] += xxz[iind][jind][0] * over_pf;
                            buffer_[ao12 + 15 * size] += xyy[iind][jind][0] * over_pf;
                            buffer_[ao12 + 16 * size] += yyz[iind][jind][0] * over_pf;
                            buffer_[ao12 + 17 * size] += xzz[iind][jind][0] * over_pf;
                            buffer_[ao12 + 18 * size] += yzz[iind][jind][0] * over_pf;
                            buffer_[ao12 + 19 * size] += xyz[iind][jind][0] * over_pf;

                            ao12++;
                        }
                    }
                }
            }
        }
    }
}
