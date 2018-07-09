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
#include "psi4/libmints/giao_potential.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))

;
using namespace psi;

GiaoPotentialInt::GiaoPotentialInt(std::vector<SphericalTransform>& st, std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2, int deriv) :
    OneBodyAOInt(st, bs1, bs2, deriv)
{
    int maxam1 = bs1_->max_am();
    int maxam2 = bs2_->max_am();

    int maxnao1 = (maxam1+1)*(maxam1+2)/2;
    int maxnao2 = (maxam2+1)*(maxam2+2)/2;

    // Increase buffer size to handle x, y, and z components
    if (deriv_ == 0) {
        potential_recur_ = new ObaraSaikaTwoCenterVIRecursion(bs1->max_am()+1, bs2->max_am()+1);
        buffer_ = new double[3*maxnao1*maxnao2];
        set_chunks(3);
    }

    if (deriv_ > 0) {
        throw std::runtime_error("GiaoPotentialInt: does not support any derivatives");
    }

    // Setup the initial field of partial charges
    Zxyz_ = std::make_shared<Matrix>("Partial Charge Field (Z,x,y,z)", bs1_->molecule()->natom(), 4);
    double** Zxyzp = Zxyz_->pointer();

    for (int A = 0; A < bs1_->molecule()->natom(); A++) {
        Zxyzp[A][0] = (double) bs1_->molecule()->Z(A);
        Zxyzp[A][1] = bs1_->molecule()->x(A);
        Zxyzp[A][2] = bs1_->molecule()->y(A);
        Zxyzp[A][3] = bs1_->molecule()->z(A);
    }
}

GiaoPotentialInt::~GiaoPotentialInt()
{
    delete[] buffer_;
    delete potential_recur_;
}

 /*  The engine only supports segmented basis sets */

void GiaoPotentialInt::compute_pair(const GaussianShell& s1, const GaussianShell& s2)
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
    B[3] = s2.center()[2];

    int ydisp = INT_NCART(am1) * INT_NCART(am2);
    int zdisp = ydisp + INT_NCART(am1) * INT_NCART(am2);

    int izm = 1;
    int iym = am1 + 2;
                           // Note +2 instead of +1 here -- we'll be computing nuclear attraction
                          // integrals with extra quantum in bra!
    int ixm = iym * iym;
    int jzm = 1;
    int jym = am2 + 1;
    int jxm = jym * jym;

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

    double ***vi = potential_recur_->vi();

    double** Zxyzp = Zxyz_->pointer();
    int ncharge = Zxyz_->rowspi()[0];

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

            // Loop over atoms of basis set 1 (only works if bs1_ and bs2_ are on the same
            // molecule)
            for (int atom=0; atom<ncharge; ++atom) {
                double PC[3];

                double Z = Zxyzp[atom][0];
                double dB_pfac = 0.5 * over_pf * Z;

                PC[0] = P[0] - Zxyzp[atom][1];
                PC[1] = P[1] - Zxyzp[atom][2];
                PC[2] = P[2] - Zxyzp[atom][3];

                // Do recursion
                potential_recur_->compute(PA, PB, PC, gamma, am1, am2);

                ao12 = 0;
                for(int ii = 0; ii <= am1; ii++) {
                    int l1 = am1 - ii;
                    for(int jj = 0; jj <= ii; jj++) {
                        int m1 = ii - jj;
                        int n1 = jj;
                        int iind_000 = n1*izm + m1*iym + l1*ixm;
                        int iind_100 = n1*izm + m1*iym + (l1+1)*ixm;
                        int iind_010 = n1*izm + (m1+1)*iym + l1*ixm;
                        int iind_001 = (n1+1)*izm + m1*iym + l1*ixm;
                        /*--- create all am components of sj ---*/
                        for(int kk = 0; kk <= am2; kk++) {
                            int l2 = am2 - kk;
                            for(int ll = 0; ll <= kk; ll++) {
                                int m2 = kk - ll;
                                int n2 = ll;

                                // Compute location in the recursion
                                int jind = l2 * jxm + m2 * jym + n2 * jzm;

                                double V  = vi[iind_000][jind][0];
                                double xV = vi[iind_100][jind][0] + A[0] * V;
                                double yV = vi[iind_010][jind][0] + A[1] * V;
                                double zV = vi[iind_001][jind][0] + A[2] * V;

                                buffer_[ao12] += dB_pfac * ( AB[1] * zV - AB[2] * yV);
                                buffer_[ao12+ydisp] += dB_pfac * ( AB[2] * xV - AB[0] * zV);
                                buffer_[ao12+zdisp] += dB_pfac * ( AB[0] * yV - AB[1] * xV);
                                ao12++;

                            }
                        }
                    }
                }
            }
        }
    }
}
