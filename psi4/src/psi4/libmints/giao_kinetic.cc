/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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

#include <stdexcept>
#include "psi4/libciomr/libciomr.h"
#include "psi4/physconst.h"
#include "psi4/libmints/giao_kinetic.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))

;
using namespace psi;

GiaoKineticInt::GiaoKineticInt(std::vector<SphericalTransform>& st, std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2, int deriv) :
    OneBodyAOInt(st, bs1, bs2, deriv), overlap_recur_(bs1->max_am()+1, bs2->max_am()+2)
{
    int maxam1 = bs1_->max_am();
    int maxam2 = bs2_->max_am();

    int maxnao1 = (maxam1+1)*(maxam1+2)/2;
    int maxnao2 = (maxam2+1)*(maxam2+2)/2;

    // Increase buffer size to handle x, y, and z components
    if (deriv_ == 0) {
        buffer_ = new double[3*maxnao1*maxnao2];
        set_chunks(3);
    }

    if (deriv_ > 0) {
        throw std::runtime_error("GiaoKineticInt: does not support any derivatives");
    }
}

GiaoKineticInt::~GiaoKineticInt()
{
    delete[] buffer_;
}

 /*  The engine only supports segmented basis sets */

void GiaoKineticInt::compute_pair(const GaussianShell& s1, const GaussianShell& s2)
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

    int ydisp = INT_NCART(am1) * INT_NCART(am2);
    int zdisp = ydisp + INT_NCART(am1) * INT_NCART(am2);

    // compute intermediates
    double AB[3];
    AB[0] = A[0]-B[0];
    AB[1] = A[1]-B[1];
    AB[2] = A[2]-B[2];

    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    memset(buffer_, 0, s1.ncartesian() * s2.ncartesian() * sizeof(double));

    double **x = overlap_recur_.x();
    double **y = overlap_recur_.y();
    double **z = overlap_recur_.z();

    for (int p1=0; p1<nprim1; ++p1) {
        double a1 = s1.exp(p1);
        double c1 = s1.coef(p1);
        for (int p2=0; p2<nprim2; ++p2) {
            double a2 = s2.exp(p2);
            double c2 = s2.coef(p2);
            double gamma = a1 + a2;
            double oog = 1.0/gamma;

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
            double dB_pfac = 0.5*over_pf;

            // Do recursion
            overlap_recur_.compute(PA, PB, gamma, am1+1, am2+2);

            ao12 = 0;
            for(int ii = 0; ii <= am1; ii++) {
                int l1 = am1 - ii;
                for(int jj = 0; jj <= ii; jj++) {
                    int m1 = ii - jj;
                    int n1 = jj;
                    /*--- create all am components of sj ---*/
                    for(int kk = 0; kk <= am2; kk++) {
                        int l2 = am2 - kk;
                        for(int ll = 0; ll <= kk; ll++) {
                            int m2 = kk - ll;
                            int n2 = ll;


                            double x00 = x[l1][l2];
                            double y00 = y[m1][m2];
                            double z00 = z[n1][n2];
                            double x01 = x[l1][l2+1];
                            double y01 = y[m1][m2+1];
                            double z01 = z[n1][n2+1];
                            double x1  = x01 + x00*B[0];
                            double y1  = y01 + y00*B[1];
                            double z1  = z01 + z00*B[2];

                            double tx00 = a2*(2*l2+1)*x[l1][l2] - 2*a2*a2*x[l1][l2+2];
                            if (l2 >= 2)
                              tx00 -= 0.5*l2*(l2-1)*x[l1][l2-2];
                            double ty00 = a2*(2*m2+1)*y[m1][m2] - 2*a2*a2*y[m1][m2+2];
                            if (m2 >= 2)
                              ty00 -= 0.5*m2*(m2-1)*y[m1][m2-2];
                            double tz00 = a2*(2*n2+1)*z[n1][n2] - 2*a2*a2*z[n1][n2+2];
                            if (n2 >= 2)
                              tz00 -= 0.5*n2*(n2-1)*z[n1][n2-2];

                            double tx10 = a2*(2*l2+1)*(x[l1+1][l2] + A[0]*x[l1][l2]) - 2*a2*a2*(x[l1+1][l2+2] + A[0]*x[l1][l2+2]);
                            if (l2 >= 2)
                              tx10 -= 0.5*l2*(l2-1)*(x[l1+1][l2-2] + A[0]*x[l1][l2-2]);
                            double ty10 = a2*(2*m2+1)*(y[m1+1][m2] + A[1]*y[m1][m2]) - 2*a2*a2*(y[m1+1][m2+2] + A[1]*y[m1][m2+2]);
                            if (m2 >= 2)
                              ty10 -= 0.5*m2*(m2-1)*(y[m1+1][m2-2] + A[1]*y[m1][m2-2]);
                            double tz10 = a2*(2*n2+1)*(z[n1+1][n2] + A[2]*z[n1][n2]) - 2*a2*a2*(z[n1+1][n2+2] + A[2]*z[n1][n2+2]);
                            if (n2 >= 2)
                              tz10 -= 0.5*n2*(n2-1)*(z[n1+1][n2-2] + A[2]*z[n1][n2-2]);


                            double xT = tx10*y00*z00 + x1*ty00*z00  + x1*y00*tz00;
                            double yT = tx00*y1*z00  + x00*ty10*z00 + x00*y1*tz00;
                            double zT = tx00*y00*z1  + x00*ty00*z1  + x00*y00*tz10;

                            buffer_[ao12]       += dB_pfac * ( AB[1] * zT - AB[2] * yT); 
                            buffer_[ao12+ydisp] += dB_pfac * ( AB[2] * xT - AB[0] * zT); 
                            buffer_[ao12+zdisp] += dB_pfac * ( AB[0] * yT - AB[1] * xT);
                            ao12++;
                        }
                    }
                }
            }
        }
    }
}
