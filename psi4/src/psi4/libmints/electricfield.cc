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
#include "psi4/libmints/electricfield.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include <stdexcept>
#include <vector>
#include "psi4/libciomr/libciomr.h"
#include "psi4/physconst.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))

;
using namespace psi;

ElectricFieldInt::ElectricFieldInt(std::vector<SphericalTransform>& spherical_transforms, std::shared_ptr<BasisSet> bs1,
                                   std::shared_ptr<BasisSet> bs2, int nderiv)
    : OneBodyAOInt(spherical_transforms, bs1, bs2, nderiv),
      efield_recur_(bs1->max_am() + 2, bs2->max_am() + 2),
      natom_(bs1->molecule()->natom()) {
    int maxam1 = bs1_->max_am();
    int maxam2 = bs2_->max_am();

    int maxnao1 = INT_NCART(maxam1);
    int maxnao2 = INT_NCART(maxam2);

    if (nderiv == 0) {
        buffer_ = new double[3 * maxnao1 * maxnao2];
        set_chunks(3);
    } else if (nderiv == 1) {
        throw PSIEXCEPTION("Field derivative ints not coded yet!");
        buffer_ = new double[6 * maxnao1 * maxnao2];
        set_chunks(6);
    } else
        throw FeatureNotImplemented("LibMints", "ElectricFieldInts called with deriv > 1", __FILE__, __LINE__);
}

ElectricFieldInt::~ElectricFieldInt() { delete[] buffer_; }

Vector3 ElectricFieldInt::nuclear_contribution(const Vector3& origin, std::shared_ptr<Molecule> mol) {
    int natom = mol->natom();

    double Ex = 0.0;
    double Ey = 0.0;
    double Ez = 0.0;
    for (int i = 0; i < natom; ++i) {
        double x = origin[0] - mol->x(i);
        double y = origin[1] - mol->y(i);
        double z = origin[2] - mol->z(i);
        double r2 = x * x + y * y + z * z;
        double r = sqrt(r2);
        // This allows us to compute the field at the nuclei correctly
        if (r < 1.0E-8) continue;

        Ex += mol->Z(i) * x / (r * r2);
        Ey += mol->Z(i) * y / (r * r2);
        Ez += mol->Z(i) * z / (r * r2);
    }

    Vector3 result(Ex, Ey, Ez);

    return result;
}

SharedMatrix ElectricFieldInt::nuclear_contribution_to_gradient(const Vector3& /*origin*/,
                                                                std::shared_ptr<Molecule> /*mol*/) {
    //    int natom = mol->natom();
    //    auto result = std::make_shared<Matrix>("Nuclear contribution to electric field gradient", 1, 6);

    //    for (int i=0; i<natom; ++i) {
    //        double x = mol->x(i) - origin[0];
    //        double y = mol->y(i) - origin[1];
    //        double z = mol->z(i) - origin[2];
    //        double r2 = x*x + y*y + z*z;
    //        double r = sqrt(r2);

    //        result->add(0, 0, -mol->Z(i) * (3*x*x-r2) / (r*r2*r2));
    //        result->add(0, 1, -mol->Z(i) * 3*x*y / (r*r2*r2));
    //        result->add(0, 2, -mol->Z(i) * 3*x*z / (r*r2*r2));
    //        result->add(0, 3, -mol->Z(i) * (3*y*y-r2) / (r*r2*r2));
    //        result->add(0, 4, -mol->Z(i) * 3*y*z / (r*r2*r2));
    //        result->add(0, 5, -mol->Z(i) * (3*z*z-r2) / (r*r2*r2));
    //    }
    //    return result;
    throw NOT_IMPLEMENTED_EXCEPTION();
}

void ElectricFieldInt::compute_pair(const GaussianShell& s1, const GaussianShell& s2) {
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
    int ydisp = size;
    int zdisp = ydisp + size;

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    memset(buffer_, 0, 3 * size * sizeof(double));

    double*** ex = efield_recur_.x();
    double*** ey = efield_recur_.y();
    double*** ez = efield_recur_.z();

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

            int atomoffset = 0;

            double PC[3];
            Vector3 C = origin_;

            PC[0] = P[0] - C[0];
            PC[1] = P[1] - C[1];
            PC[2] = P[2] - C[2];

            // Get recursive
            efield_recur_.compute(PA, PB, PC, gamma, am1, am2);
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

                            buffer_[atomoffset + ao12] += ex[iind][jind][0] * over_pf;
                            buffer_[atomoffset + ao12 + ydisp] += ey[iind][jind][0] * over_pf;
                            buffer_[atomoffset + ao12 + zdisp] += ez[iind][jind][0] * over_pf;

                            ao12++;
                        }
                    }
                }
            }
        }
    }
}

void ElectricFieldInt::compute_pair_deriv1(const GaussianShell& /*s1*/, const GaussianShell& /*s2*/) {
    throw NOT_IMPLEMENTED_EXCEPTION();
    // NOT IMPLEMENTED!!!!!!
    //    int ao12;
    //    int am1 = s1.am();
    //    int am2 = s2.am();
    //    int nprim1 = s1.nprimitive();
    //    int nprim2 = s2.nprimitive();
    //    double A[3], B[3];
    //    A[0] = s1.center()[0];
    //    A[1] = s1.center()[1];
    //    A[2] = s1.center()[2];
    //    B[0] = s2.center()[0];
    //    B[1] = s2.center()[1];
    //    B[2] = s2.center()[2];

    //    int izm = 1;
    //    int iym = am1 + 1;
    //    int ixm = iym * iym;
    //    int jzm = 1;
    //    int jym = am2 + 1;
    //    int jxm = jym * jym;

    //    // Not sure if these are needed.
    //    int size =  INT_NCART(am1) * INT_NCART(am2);
    //    int xydisp = size;
    //    int xzdisp = xydisp + size;
    //    int yydisp = xzdisp + size;
    //    int yzdisp = yydisp + size;
    //    int zzdisp = yzdisp + size;

    //    // compute intermediates
    //    double AB2 = 0.0;
    //    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    //    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    //    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    //    memset(buffer_, 0, 6 * size * sizeof(double));

    //    double ***exx = efield_recur_.vxx();
    //    double ***exy = efield_recur_.vxy();
    //    double ***exz = efield_recur_.vxz();
    //    double ***eyy = efield_recur_.vyy();
    //    double ***eyz = efield_recur_.vyz();
    //    double ***ezz = efield_recur_.vzz();

    //    for (int p1=0; p1<nprim1; ++p1) {
    //        double a1 = s1.exp(p1);
    //        double c1 = s1.coef(p1);
    //        for (int p2=0; p2<nprim2; ++p2) {
    //            double a2 = s2.exp(p2);
    //            double c2 = s2.coef(p2);
    //            double gamma = a1 + a2;
    //            double oog = 1.0 / gamma;

    //            double PA[3], PB[3];
    //            double P[3];

    //            P[0] = (a1*A[0] + a2*B[0])*oog;
    //            P[1] = (a1*A[1] + a2*B[1])*oog;
    //            P[2] = (a1*A[2] + a2*B[2])*oog;
    //            PA[0] = P[0] - A[0];
    //            PA[1] = P[1] - A[1];
    //            PA[2] = P[2] - A[2];
    //            PB[0] = P[0] - B[0];
    //            PB[1] = P[1] - B[1];
    //            PB[2] = P[2] - B[2];

    //            double over_pf = exp(-a1*a2*AB2*oog) * sqrt(M_PI*oog) * M_PI * oog * c1 * c2;

    //            int atomoffset = 0;

    //            double PC[3];
    //            Vector3 C = origin_;

    //            PC[0] = P[0] - C[0];
    //            PC[1] = P[1] - C[1];
    //            PC[2] = P[2] - C[2];

    //            // Get recursive
    //            efield_recur_.compute(PA, PB, PC, gamma, am1+2, am2+2);

    //            // Gather contributions.
    //            ao12 = 0;
    //            for (int ii = 0; ii <= am1; ++ii) {
    //                int l1 = am1 - ii;
    //                for (int jj = 0; jj <= ii; ++jj) {
    //                    int m1 = ii - jj;
    //                    int n1 = jj;

    //                    for (int kk = 0; kk <= am2; ++kk) {
    //                        int l2 = am2 - kk;
    //                        for (int ll = 0; ll <= kk; ++ll) {
    //                            int m2 = kk - ll;
    //                            int n2 = ll;

    //                            // Compute location in the recursion
    //                            int iind = l1 * ixm + m1 * iym + n1 * izm;
    //                            int jind = l2 * jxm + m2 * jym + n2 * jzm;

    //                            buffer_[atomoffset+ao12]        += exx[iind][jind][0] * over_pf;
    //                            buffer_[atomoffset+ao12+xydisp] += exy[iind][jind][0] * over_pf;
    //                            buffer_[atomoffset+ao12+xzdisp] += exz[iind][jind][0] * over_pf;
    //                            buffer_[atomoffset+ao12+yydisp] += eyy[iind][jind][0] * over_pf;
    //                            buffer_[atomoffset+ao12+yzdisp] += eyz[iind][jind][0] * over_pf;
    //                            buffer_[atomoffset+ao12+zzdisp] += ezz[iind][jind][0] * over_pf;

    //                            ao12++;
    //                        }
    //                    }
    //                }
    //            }
    //        }
    //    }
}
