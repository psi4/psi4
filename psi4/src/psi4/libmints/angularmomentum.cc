/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
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

#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/angularmomentum.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/integral.h"
#include "psi4/physconst.h"

#include <memory>
#include <stdexcept>
#include <iostream>

using namespace psi;

// Initialize overlap_recur_ to +1 basis set angular momentum, +1 on each center is sufficient
// to compute the dipole derivatives
AngularMomentumInt::AngularMomentumInt(std::vector<SphericalTransform>& spherical_transforms,
                                       std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2, int nderiv)
    : OneBodyAOInt(spherical_transforms, bs1, bs2, nderiv), overlap_recur_(bs1->max_am() + 2, bs2->max_am() + 2) {
    int maxam1 = bs1_->max_am();
    int maxam2 = bs2_->max_am();

    int maxnao1 = (maxam1 + 1) * (maxam1 + 2) / 2;
    int maxnao2 = (maxam2 + 1) * (maxam2 + 2) / 2;

    // Increase buffer size to handle x, y, and z components
    if (deriv_ == 0) {
        buffer_ = new double[3 * maxnao1 * maxnao2];
        set_chunks(3);
    } else if (deriv_ == 1) {
        natom_ = bs1_->molecule()->natom();
        //buffer_ = new double[3 * 3 * natom_ * maxnao1 * maxnao2];  // 3 * 3 * N * maxnao1 * maxnao2
        //set_chunks(3 * 3 * natom_);
        buffer_ = new double[6 * 3 * maxnao1 * maxnao2];
        set_chunks(6 * 3);
    }
}

AngularMomentumInt::~AngularMomentumInt() { delete[] buffer_; }

// The engine only supports segmented basis sets
void AngularMomentumInt::compute_pair(const GaussianShell& s1, const GaussianShell& s2) {
    int ao12;
    int am1 = s1.am();
    int am2 = s2.am();
    int nprim1 = s1.nprimitive();
    int nprim2 = s2.nprimitive();
    double A[3], B[3];
    double Sxy1, Sxy2, Sxz1, Sxz2;
    double Syx1, Syx2, Syz1, Syz2;
    double Szx1, Szx2, Szy1, Szy2;
    double S0x1, S0x2, S0y1, S0y2, S0z1, S0z2;
    double muxy1, muxy2, muxz1, muxz2;
    double muyx1, muyx2, muyz1, muyz2;
    double muzx1, muzx2, muzy1, muzy2;

    A[0] = s1.center()[0];
    A[1] = s1.center()[1];
    A[2] = s1.center()[2];
    B[0] = s2.center()[0];
    B[1] = s2.center()[1];
    B[2] = s2.center()[2];

    int ydisp = INT_NCART(am1) * INT_NCART(am2);
    int zdisp = ydisp + INT_NCART(am1) * INT_NCART(am2);

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    memset(buffer_, 0, 3 * INT_NCART(am1) * INT_NCART(am2) * sizeof(double));

    double** x = overlap_recur_.x();
    double** y = overlap_recur_.y();
    double** z = overlap_recur_.z();

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

            // Do recursion
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

                            double x00 = x[l1][l2], y00 = y[m1][m2], z00 = z[n1][n2];
                            double x10 = x[l1 + 1][l2], y10 = y[m1 + 1][m2], z10 = z[n1 + 1][n2];
                            double x01 = x[l1][l2 + 1], y01 = y[m1][m2 + 1], z01 = z[n1][n2 + 1];

                            Sxy1 = Sxy2 = Sxz1 = Sxz2 = 0.0;

                            //
                            // Overlaps
                            //

                            // (a+1x|b+1y)
                            Sxy1 = x10 * y01 * z00 * over_pf;
                            // (a+1x|b-1y)
                            if (m2) Sxy2 = x10 * y[m1][m2 - 1] * z00 * over_pf;
                            // (a+1x|b+1z)
                            Sxz1 = x10 * y00 * z01 * over_pf;
                            // (a+1x|b-1z)
                            if (n2) Sxz2 = x10 * y00 * z[n1][n2 - 1] * over_pf;

                            //                            outfile->Printf( "Sxy1 %f Sxy2 %f Sxz1 %f Sxz2 %f\n", Sxy1,
                            //                            Sxy2, Sxz1, Sxz2);

                            Syx1 = Syx2 = Syz1 = Syz2 = 0.0;

                            // (a+1y|b+1x)
                            Syx1 = x01 * y10 * z00 * over_pf;
                            // (a+1y|b-1x)
                            if (l2) Syx2 = x[l1][l2 - 1] * y10 * z00 * over_pf;
                            // (a+1y|b+1z)
                            Syz1 = x00 * y10 * z01 * over_pf;
                            // (a+1y|b-1z)
                            if (n2) Syz2 = x00 * y10 * z[n1][n2 - 1] * over_pf;

                            //                            outfile->Printf( "Syx1 %f Syx2 %f Syz1 %f Syz2 %f\n", Syx1,
                            //                            Syx2, Syz1, Syz2);

                            Szx1 = Szx2 = Szy1 = Szy2 = 0.0;

                            // (a+1z|b+1x)
                            Szx1 = x01 * y00 * z10 * over_pf;
                            // (a+1z|b-1x)
                            if (l2) Szx2 = x[l1][l2 - 1] * y00 * z10 * over_pf;
                            // (a+1z|b+1y)
                            Szy1 = x00 * y01 * z10 * over_pf;
                            // (a+1z|b-1y)
                            if (m2) Szy2 = x00 * y[m1][m2 - 1] * z10 * over_pf;

                            //                            outfile->Printf( "Szx1 %f Szx2 %f Szy1 %f Szy2 %f\n", Szx1,
                            //                            Szx2, Szy1, Szy2);

                            S0x1 = S0x2 = S0y1 = S0y2 = S0z1 = S0z2 = 0.0;

                            // (a|b+1x)
                            S0x1 = x01 * y00 * z00 * over_pf;
                            // (a|b-1x)
                            if (l2) S0x2 = x[l1][l2 - 1] * y00 * z00 * over_pf;
                            // (a|b+1y)
                            S0y1 = x00 * y01 * z00 * over_pf;
                            // (a|b-1y)
                            if (m2) S0y2 = x00 * y[m1][m2 - 1] * z00 * over_pf;
                            // (a|b+1z)
                            S0z1 = x00 * y00 * z01 * over_pf;
                            // (a|b-1z)
                            if (n2) S0z2 = x00 * y00 * z[n1][n2 - 1] * over_pf;

                            //                            outfile->Printf( "S0x1 %f S0x2 %f S0y1 %f S0y2 %f S0z1
                            //                            S0z2\n", S0x1, S0x2, S0y1, S0y2, S0z1, S0z2);

                            //
                            // Moment integrals
                            //

                            muxy1 = muxy2 = muxz1 = muxz2 = 0.0;
                            // (a|x|b+1y)
                            muxy1 = Sxy1 + (A[0] - origin_[0]) * S0y1;
                            // (a|x|b+1z)
                            muxz1 = Sxz1 + (A[0] - origin_[0]) * S0z1;
                            // (a|x|b-1y)
                            muxy2 = Sxy2 + (A[0] - origin_[0]) * S0y2;
                            // (a|x|b-1z)
                            muxz2 = Sxz2 + (A[0] - origin_[0]) * S0z2;

                            muyx1 = muyx2 = muyz1 = muyz2 = 0.0;
                            // (a|y|b+1x)
                            muyx1 = Syx1 + (A[1] - origin_[1]) * S0x1;
                            // (a|y|b+1z)
                            muyz1 = Syz1 + (A[1] - origin_[1]) * S0z1;
                            // (a|y|b-1x)
                            muyx2 = Syx2 + (A[1] - origin_[1]) * S0x2;
                            // (a|y|b-1z)
                            muyz2 = Syz2 + (A[1] - origin_[1]) * S0z2;

                            muzx1 = muzx2 = muzy1 = muzy2 = 0.0;
                            // (a|z|b+1x)
                            muzx1 = Szx1 + (A[2] - origin_[2]) * S0x1;
                            // (a|z|b+1y)
                            muzy1 = Szy1 + (A[2] - origin_[2]) * S0y1;
                            // (a|z|b+1x)
                            muzx2 = Szx2 + (A[2] - origin_[2]) * S0x2;
                            // (a|z|b+1y)
                            muzy2 = Szy2 + (A[2] - origin_[2]) * S0y2;

                            /* (a|Lx|b) = 2 a2 * (a|(y-Cy)|b+1z) - B.z * (a|(y-Cy)|b-1z)
                               - 2 a2 * (a|(z-Cz)|b+1y) + B.y * (a|(z-Cz)|b-1y) */
//                            muyz1 = muyz2 = muzy2 = 0.0;
                            double Lx = (2.0 * a2 * muyz1 - n2 * muyz2 - 2.0 * a2 * muzy1 + m2 * muzy2) /* * over_pf*/;
//                            double Lx = 0.0;

                            /* (a|Ly|b) = 2 a2 * (a|(z-Cz)|b+1x) - B.x * (a|(z-Cz)|b-1x)
                               - 2 a2 * (a|(x-Cx)|b+1z) + B.z * (a|(x-Cx)|b-1z) */
//                            muzx1 = muzx2 = muxz1 = 0.0;
                            double Ly = (2.0 * a2 * muzx1 - l2 * muzx2 - 2.0 * a2 * muxz1 + n2 * muxz2) /* * over_pf*/;
//                            double Ly = 0.0;

                            /* (a|Lz|b) = 2 a2 * (a|(x-Cx)|b+1y) - B.y * (a|(x-Cx)|b-1y)
                               - 2 a2 * (a|(y-Cy)|b+1x) + B.x * (a|(y-Cy)|b-1x) */
//                            muxy1 = muxy2 = muyx1 = 0.0;
                            double Lz = (2.0 * a2 * muxy1 - m2 * muxy2 - 2.0 * a2 * muyx1 + l2 * muyx2) /* * over_pf*/;
//                            double Lz = 0.0;

                            // Electrons have a negative charge
                            buffer_[ao12] += (Lx);
                            buffer_[ao12 + ydisp] += (Ly);
                            buffer_[ao12 + zdisp] += (Lz);

                            ao12++;
                        }
                    }
                }
            }
        }
    }
}

// The engine only supports segmented basis sets
void AngularMomentumInt::compute_pair_deriv1(const GaussianShell& s1, const GaussianShell& s2)
{
    int ao12;
    int am1 = s1.am();
    int am2 = s2.am();
    int at1 = s1.ncenter();
    int at2 = s2.ncenter();
    int nprim1 = s1.nprimitive();
    int nprim2 = s2.nprimitive();
    size_t length = INT_NCART(am1) * INT_NCART(am2);
    double A[3], B[3], C[3];

    A[0] = s1.center()[0];
    A[1] = s1.center()[1];
    A[2] = s1.center()[2];
    B[0] = s2.center()[0];
    B[1] = s2.center()[1];
    B[2] = s2.center()[2];
    C[0] = origin_[0];
    C[1] = origin_[1];
    C[2] = origin_[2];

    //size_t xaxdisp = at1 * length * 9;
    //size_t xaydisp = xaxdisp + length;
    //size_t xazdisp = xaydisp + length;
    //size_t yaxdisp = xazdisp + length;
    //size_t yaydisp = yaxdisp + length;
    //size_t yazdisp = yaydisp + length;
    //size_t zaxdisp = yazdisp + length;
    //size_t zaydisp = zaxdisp + length;
    //size_t zazdisp = zaydisp + length;

    //size_t xbxdisp = at2 * length * 9;
    //size_t xbydisp = xbxdisp + length;
    //size_t xbzdisp = xbydisp + length;
    //size_t ybxdisp = xbzdisp + length;
    //size_t ybydisp = ybxdisp + length;
    //size_t ybzdisp = ybydisp + length;
    //size_t zbxdisp = ybzdisp + length;
    //size_t zbydisp = zbxdisp + length;
    //size_t zbzdisp = zbydisp + length;

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

    //std::cout << "at1 = " << at1 << std::endl;
    //std::cout << "length = " << length << std::endl;
    //std::cout << "xaxdisp = " << xaxdisp << std::endl;
    //std::cout << "xaydisp = " << xaydisp << std::endl;
    //std::cout << "xazdisp = " << xazdisp << std::endl;
    //std::cout << "yaxdisp = " << yaxdisp << std::endl;
    //std::cout << "yaydisp = " << yaydisp << std::endl;
    //std::cout << "yazdisp = " << yazdisp << std::endl;
    //std::cout << "zaxdisp = " << zaxdisp << std::endl;
    //std::cout << "zaydisp = " << zaydisp << std::endl;
    //std::cout << "zazdisp = " << zazdisp << std::endl;

    //std::cout << "at2 = " << at2 << std::endl;
    //std::cout << "length = " << length << std::endl;
    //std::cout << "xbxdisp = " << xbxdisp << std::endl;
    //std::cout << "xbydisp = " << xbydisp << std::endl;
    //std::cout << "xbzdisp = " << xbzdisp << std::endl;
    //std::cout << "ybxdisp = " << ybxdisp << std::endl;
    //std::cout << "ybydisp = " << ybydisp << std::endl;
    //std::cout << "ybzdisp = " << ybzdisp << std::endl;
    //std::cout << "zbxdisp = " << zbxdisp << std::endl;
    //std::cout << "zbydisp = " << zbydisp << std::endl;
    //std::cout << "zbzdisp = " << zbzdisp << std::endl;

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    // Zero out the buffer
    //memset(buffer_, 0, 3 * 3 * natom_ * length * sizeof(double));
    memset(buffer_, 0, 6 * 3 * length * sizeof(double));

    double **x = overlap_recur_.x();
    double **y = overlap_recur_.y();
    double **z = overlap_recur_.z();
    double v1, v2, v3, v4, v5;      // temporary value storage
    double v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16, v17, v18, v19, v20;      // temporary value storage

    //std::printf("%5s \t %5s \t %5s \t %5s \t %5s \t %5s \t %5s \t %5s \t %5s \t %10s \t %10s \n", "p1", "p2", "ao12", "am1", "am2", "ii", "jj", "kk", "ll", "xaxdisp", "xaydisp");
    //std::printf("--------------------------------------------------------------------------------------------------------- \n");
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
            //overlap_recur_.compute(PA, PB, gamma, am1+1, am2+1);
            overlap_recur_.compute(PA, PB, gamma, am1+2, am2+2);

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

                            //std::cout << "am1 = " << am1 << "\t am2 = " << am2 << "\t ao12 = " << ao12 << std::endl;
                            //std::printf("%5d \t %5d \t %5d \t %5d \t %5d \t %5d \t %5d \t %5d \t %5d \n", p1, p2, ao12, am1, am2, ii, jj, kk, ll);

//                            //
//                            // Ax derivatives with Lx
//                            //
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_y+1_x|Lx|b+1_z)
//                            v1 = x[l1+1][l2] * y[m1+1][m2] * z[n1][n2+1];
//                            // (a+1_x|Lx|b+1_z)
//                            v2 = x[l1+1][l2] * y[m1][m2] * z[n1][n2+1];
//                            if (l1) {
//                                // (a+1_y-1_x|Lx|b+1_z)
//                                v3 = x[l1-1][l2] * y[m1+1][m2] * z[n1][n2+1];
//                                // (a-1_x|Lx|b+1_z)
//                                v4 = x[l1-1][l2] * y[m1][m2] * z[n1][n2+1];
//                            }
//                            // (a|Lx|b+1_z)
//                            //v5 = x[l1][l2] * y[m1][m2] * z[n1][n2+1]; // because kronecker_delta(j,l) = (y,x) = 0
//                            buffer_[ao12+xaxdisp] += 2.0 * a2 * (2.0 * a1 * (v1 + (A[1] - C[1]) * v2) - l1 * (v3 + (A[1] - C[1]) * v4) + v5) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            if (n2) {
//                                // (a+1_y+1_x|Lx|b-1_z)
//                                v1 = x[l1+1][l2] * y[m1+1][m2] * z[n1][n2-1];
//                                // (a+1_x|Lx|b-1_z)
//                                v2 = x[l1+1][l2] * y[m1][m2] * z[n1][n2-1];
//                                if (l1) {
//                                    // (a+1_y-1_x|Lx|b-1_z)
//                                    v3 = x[l1-1][l2] * y[m1+1][m2] * z[n1][n2-1];
//                                    // (a-1_x|Lx|b-1_z)
//                                    v4 = x[l1-1][l2] * y[m1][m2] * z[n1][n2-1];
//                                }
//                                // (a|Lx|b-1_z)
//                                //v5 = x[l1][l2] * y[m1][m2] * z[n1][n2-1]; // because kronecker_delta(j,l) = (y,x) = 0
//                            }
//                            buffer_[ao12+xaxdisp] += -1.0 * n2 * (2.0 * a1 * (v1 + (A[1] - C[1]) * v2) - l1 * (v3 + (A[1] - C[1]) * v4) + v5) * over_pf;
//                            
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_z+1_x|Lx|b+1_y)
//                            v1 = x[l1+1][l2] * y[m1][m2+1] * z[n1+1][n2];
//                            // (a+1_x|Lx|b+1_y)
//                            v2 = x[l1+1][l2] * y[m1][m2+1] * z[n1][n2];
//                            if (l1) {
//                                // (a+1_z-1_x|Lx|b+1_y)
//                                v3 = x[l1-1][l2] * y[m1+1][m2+1] * z[n1+1][n2];
//                                // (a-1_x|Lx|b+1_y)
//                                v4 = x[l1-1][l2] * y[m1][m2+1] * z[n1][n2];
//                            }
//                            // (a|Lx|b+1_y)
//                            //v5 = x[l1][l2] * y[m1][m2+1] * z[n1][n2]; // because kronecker_delta(k,l) = (z,x) = 0
//                            buffer_[ao12+xaxdisp] += -2.0 * a2 * (2.0 * a1 * (v1 + (A[2] - C[2]) * v2) - l1 * (v3 + (A[2] - C[2]) * v4) + v5) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            if (m2) {
//                                // (a+1_z+1_x|Lx|b-1_y)
//                                v1 = x[l1+1][l2] * y[m1][m2-1] * z[n1+1][n2];
//                                // (a+1_x|Lx|b-1_y)
//                                v2 = x[l1+1][l2] * y[m1][m2-1] * z[n1][n2];
//                                if (l1) {
//                                    // (a+1_z-1_x|Lx|b-1_y)
//                                    v3 = x[l1-1][l2] * y[m1+1][m2-1] * z[n1+1][n2];
//                                    // (a-1_x|Lx|b-1_y)
//                                    v4 = x[l1-1][l2] * y[m1][m2-1] * z[n1][n2];
//                                }
//                                // (a|Lx|b-1_y)
//                                //v5 = x[l1][l2] * y[m1][m2-1] * z[n1][n2]; // because kronecker_delta(k,l) = (z,x) = 0
//                            }
//                            buffer_[ao12+xaxdisp] += 1.0 * m2 * (2.0 * a1 * (v1 + (A[2] - C[2]) * v2) - l1 * (v3 + (A[2] - C[2]) * v4) + v5) * over_pf;
//
//                            //
//                            // Ay derivatives with Lx
//                            //
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_y+1_y|Lx|b+1_z)
//                            v1 = x[l1][l2] * y[m1+2][m2] * z[n1][n2+1];
//                            // (a+1_y|Lx|b+1_z)
//                            v2 = x[l1][l2] * y[m1+1][m2] * z[n1][n2+1];
//                            // (a+1_y-1_y|Lx|b+1_z)
//                            v3 = x[l1][l2] * y[m1][m2] * z[n1][n2+1];
//                            if (m1) {
//                                // (a-1_y|Lx|b+1_z)
//                                v4 = x[l1][l2] * y[m1-1][m2] * z[n1][n2+1];
//                            }
//                            // (a|Lx|b+1_z)
//                            v5 = x[l1][l2] * y[m1][m2] * z[n1][n2+1]; // because kronecker_delta(j,l) = (y,y) = 1
//                            buffer_[ao12+xaydisp] += 2.0 * a2 * (2.0 * a1 * (v1 + (A[1] - C[1]) * v2) - m1 * (v3 + (A[1] - C[1]) * v4) + v5) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            if (n2) {
//                                // (a+1_y+1_y|Lx|b-1_z)
//                                v1 = x[l1][l2] * y[m1+2][m2] * z[n1][n2-1];
//                                // (a+1_y|Lx|b-1_z)
//                                v2 = x[l1][l2] * y[m1+1][m2] * z[n1][n2-1];
//                                // (a+1_y-1_y|Lx|b-1_z)
//                                v3 = x[l1][l2] * y[m1][m2] * z[n1][n2-1];
//                                if (m1) {
//                                    // (a-1_y|Lx|b-1_z)
//                                    v4 = x[l1][l2] * y[m1-1][m2] * z[n1][n2-1];
//                                }
//                                // (a|Lx|b-1_z)
//                                v5 = x[l1][l2] * y[m1][m2] * z[n1][n2-1]; // because kronecker_delta(j,l) = (y,y) = 1
//                            }
//                            buffer_[ao12+xaydisp] += -1.0 * n2 * (2.0 * a1 * (v1 + (A[1] - C[1]) * v2) - m1 * (v3 + (A[1] - C[1]) * v4) + v5) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_z+1_y|Lx|b+1_y)
//                            v1 = x[l1][l2] * y[m1+1][m2+1] * z[n1+1][n2];
//                            // (a+1_y|Lx|b+1_y)
//                            v2 = x[l1][l2] * y[m1+1][m2+1] * z[n1][n2];
//                            if (m1) {
//                                // (a+1_z-1_y|Lx|b+1_y)
//                                v3 = x[l1][l2] * y[m1-1][m2+1] * z[n1+1][n2];
//                                // (a-1_y|Lx|b+1_y)
//                                v4 = x[l1][l2] * y[m1-1][m2+1] * z[n1][n2];
//                            }
//                            // (a|Lx|b+1_y)
//                            //v5 = x[l1][l2] * y[m1][m2+1] * z[n1][n2]; // because kronecker_delta(k,l) = (z,y) = 0
//                            buffer_[ao12+xaydisp] += -2.0 * a2 * (2.0 * a1 * (v1 + (A[2] - C[2]) * v2) - m1 * (v3 + (A[2] - C[2]) * v4) + v5) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            if (m2) {
//                                // (a+1_z+1_y|Lx|b-1_y)
//                                v1 = x[l1][l2] * y[m1+1][m2-1] * z[n1+1][n2];
//                                // (a+1_y|Lx|b-1_y)
//                                v2 = x[l1][l2] * y[m1+1][m2-1] * z[n1][n2];
//                                if (m1) {
//                                    // (a+1_z-1_y|Lx|b-1_y)
//                                    v3 = x[l1][l2] * y[m1-1][m2-1] * z[n1+1][n2];
//                                    // (a-1_y|Lx|b-1_y)
//                                    v4 = x[l1][l2] * y[m1-1][m2-1] * z[n1][n2];
//                                }
//                                // (a|Lx|b-1_y)
//                                //v5 = x[l1][l2] * y[m1][m2-1] * z[n1][n2]; // because kronecker_delta(k,l) = (z,y) = 0
//                            }
//                            buffer_[ao12+xaydisp] += 1.0 * m2 * (2.0 * a1 * (v1 + (A[2] - C[2]) * v2) - m1 * (v3 + (A[2] - C[2]) * v4) + v5) * over_pf;
                            
                            //
                            // Az derivatives with Lx
                            //

                            v1 = v2 = v3 = v4 = v5 = 0.0;
                            // (a+1_y+1_z|Lx|b+1_z)
                            v1 = x[l1][l2] * y[m1+1][m2] * z[n1+1][n2+1];
                            // (a+1_z|Lx|b+1_z)
                            v2 = x[l1][l2] * y[m1][m2] * z[n1+1][n2+1];
                            if (n1) {
                                // (a+1_y-1_z|Lx|b+1_z)
                                v3 = x[l1][l2] * y[m1+1][m2] * z[n1-1][n2+1];
                                // (a-1_z|Lx|b+1_z)
                                v4 = x[l1][l2] * y[m1][m2] * z[n1-1][n2+1];
                            }
                            // (a|Lx|b+1_z)
                            //v5 = x[l1][l2] * y[m1][m2] * z[n1][n2+1]; // because kronecker_delta(j,l) = (y,z) = 0
                            buffer_[ao12+xazdisp] += 2.0 * a2 * (2.0 * a1 * (v1 + (A[1] - C[1]) * v2) - n1 * (v3 + (A[1] - C[1]) * v4) + v5) * over_pf;

                            v1 = v2 = v3 = v4 = v5 = 0.0;
                            if (n2) {
                                // (a+1_y+1_z|Lx|b-1_z)
                                v1 = x[l1][l2] * y[m1+1][m2] * z[n1+1][n2-1];
                                // (a+1_z|Lx|b-1_z)
                                v2 = x[l1][l2] * y[m1][m2] * z[n1+1][n2-1];
                                if (n1) {
                                    // (a+1_y-1_z|Lx|b-1_z)
                                    v3 = x[l1][l2] * y[m1+1][m2] * z[n1-1][n2-1];
                                    // (a-1_z|Lx|b-1_z)
                                    v4 = x[l1][l2] * y[m1][m2] * z[n1-1][n2-1];
                                }
                                // (a|Lx|b-1_z)
                                //v5 = x[l1][l2] * y[m1][m2] * z[n1][n2-1]; // because kronecker_delta(j,l) = (y,z) = 0
                            }
                            buffer_[ao12+xazdisp] += -1.0 * n2 * (2.0 * a1 * (v1 + (A[1] - C[1]) * v2) - n1 * (v3 + (A[1] - C[1]) * v4) + v5) * over_pf;

                            v1 = v2 = v3 = v4 = v5 = 0.0;
                            // (a+1_z+1_z|Lx|b+1_y)
                            v1 = x[l1][l2] * y[m1][m2+1] * z[n1+2][n2];
                            // (a+1_z|Lx|b+1_y)
                            v2 = x[l1][l2] * y[m1][m2+1] * z[n1+1][n2];
                            // (a+1_z-1_z|Lx|b+1_y)
                            v3 = x[l1][l2] * y[m1][m2+1] * z[n1][n2];
                            if (n1) {
                                // (a-1_z|Lx|b+1_y)
                                v4 = x[l1][l2] * y[m1][m2+1] * z[n1-1][n2];
                            }
                            // (a|Lx|b+1_y)
                            v5 = x[l1][l2] * y[m1][m2+1] * z[n1][n2]; // because kronecker_delta(k,l) = (z,z) = 1
                            buffer_[ao12+xazdisp] += -2.0 * a2 * (2.0 * a1 * (v1 + (A[2] - C[2]) * v2) - n1 * (v3 + (A[2] - C[2]) * v4) + v5) * over_pf;

                            v1 = v2 = v3 = v4 = v5 = 0.0;
                            if (m2) {
                                // (a+1_z+1_z|Lx|b-1_y)
                                v1 = x[l1][l2] * y[m1][m2-1] * z[n1+2][n2];
                                // (a+1_z|Lx|b-1_y)
                                v2 = x[l1][l2] * y[m1][m2-1] * z[n1+1][n2];
                                // (a+1_z-1_z|Lx|b-1_y)
                                v3 = x[l1][l2] * y[m1][m2-1] * z[n1][n2];
                                if (n1) {
                                    // (a-1_z|Lx|b-1_y)
                                    v4 = x[l1][l2] * y[m1][m2-1] * z[n1-1][n2];
                                }
                                // (a|Lx|b-1_y)
                                v5 = x[l1][l2] * y[m1][m2-1] * z[n1][n2]; // because kronecker_delta(k,l) = (z,z) = 1
                            }
                            buffer_[ao12+xazdisp] += 1.0 * m2 * (2.0 * a1 * (v1 + (A[2] - C[2]) * v2) - n1 * (v3 + (A[2] - C[2]) * v4) + v5) * over_pf;

//                            //
//                            // Bx derivatives with Lx
//                            //
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_y|Lx|b+1_z+1_x)
//                            v1 = x[l1][l2+1] * y[m1+1][m2] * z[n1][n2+1];
//                            // (a|Lx|b+1_z+1_x)
//                            v2 = x[l1][l2+1] * y[m1][m2] * z[n1][n2+1];
//                            if (l2) {
//                                // (a+1_y|Lx|b+1_z-1_x)
//                                v3 = x[l1][l2-1] * y[m1+1][m2] * z[n1][n2+1];
//                                // (a|Lx|b+1_z-1_x)
//                                v4 = x[l1][l2-1] * y[m1][m2] * z[n1][n2+1];
//                            }
//                            buffer_[ao12+xbxdisp] += 2.0 * a2 * (2.0 * a2 * (v1 + (A[1] - C[1]) * v2) - l2 * (v3 + (A[1] - C[1]) * v4)) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            if (n2) {
//                                // (a+1_y|Lx|b-1_z+1_x)
//                                v1 = x[l1][l2+1] * y[m1+1][m2] * z[n1][n2-1];
//                                // (a|Lx|b-1_z+1_x)
//                                v2 = x[l1][l2+1] * y[m1][m2] * z[n1][n2-1];
//                                if (l2) {
//                                    // (a+1_y|Lx|b-1_z-1_x)
//                                    v3 = x[l1][l2-1] * y[m1+1][m2] * z[n1][n2-1];
//                                    // (a|Lx|b-1_z-1_x)
//                                    v4 = x[l1][l2-1] * y[m1][m2] * z[n1][n2-1];
//                                }
//                            }
//                            buffer_[ao12+xbxdisp] += -1.0 * n2 * (2.0 * a2 * (v1 + (A[1] - C[1]) * v2) - l2 * (v3 + (A[1] - C[1]) * v4)) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_z|Lx|b+1_y+1_x)
//                            v1 = x[l1][l2+1] * y[m1][m2+1] * z[n1+1][n2];
//                            // (a|Lx|b+1_y+1_x)
//                            v2 = x[l1][l2+1] * y[m1][m2+1] * z[n1][n2];
//                            if (l2) {
//                                // (a+1_z|Lx|b+1_y-1_x)
//                                v3 = x[l1][l2-1] * y[m1][m2+1] * z[n1+1][n2];
//                                // (a|Lx|b+1_y-1_x)
//                                v4 = x[l1][l2-1] * y[m1][m2+1] * z[n1][n2];
//                            }
//                            buffer_[ao12+xbxdisp] += -2.0 * a2 * (2.0 * a2 * (v1 + (A[2] - C[2]) * v2) - l2 * (v3 + (A[2] - C[2]) * v4)) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            if (m2) {
//                                // (a+1_z|Lx|b-1_y+1_x)
//                                v1 = x[l1][l2+1] * y[m1][m2-1] * z[n1+1][n2];
//                                // (a|Lx|b-1_y+1_x)
//                                v2 = x[l1][l2+1] * y[m1][m2-1] * z[n1][n2];
//                                if (l2) {
//                                    // (a+1_z|Lx|b-1_y-1_x)
//                                    v3 = x[l1][l2-1] * y[m1][m2-1] * z[n1+1][n2];
//                                    // (a|Lx|b-1_y-1_x)
//                                    v4 = x[l1][l2-1] * y[m1][m2-1] * z[n1][n2];
//                                }
//                            }
//                            buffer_[ao12+xbxdisp] += 1.0 * m2 * (2.0 * a2 * (v1 + (A[2] - C[2]) * v2) - l2 * (v3 + (A[2] - C[2]) * v4)) * over_pf;
//
//                            //
//                            // By derivatives with Lx
//                            //
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_y|Lx|b+1_z+1_y)
//                            v1 = x[l1][l2] * y[m1+1][m2+1] * z[n1][n2+1];
//                            // (a|Lx|b+1_z+1_y)
//                            v2 = x[l1][l2] * y[m1][m2+1] * z[n1][n2+1];
//                            if (m2) {
//                                // (a+1_y|Lx|b+1_z-1_y)
//                                v3 = x[l1][l2] * y[m1+1][m2-1] * z[n1][n2+1];
//                                // (a|Lx|b+1_z-1_y)
//                                v4 = x[l1][l2] * y[m1][m2-1] * z[n1][n2+1];
//                            }
//                            buffer_[ao12+xbydisp] += 2.0 * a2 * (2.0 * a2 * (v1 + (A[1] - C[1]) * v2) - m2 * (v3 + (A[1] - C[1]) * v4)) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            if (n2) {
//                                // (a+1_y|Lx|b-1_z+1_y)
//                                v1 = x[l1][l2] * y[m1+1][m2+1] * z[n1][n2-1];
//                                // (a|Lx|b-1_z+1_y)
//                                v2 = x[l1][l2] * y[m1][m2+1] * z[n1][n2-1];
//                                if (m2) {
//                                    // (a+1_y|Lx|b-1_z-1_y)
//                                    v3 = x[l1][l2] * y[m1+1][m2-1] * z[n1][n2-1];
//                                    // (a|Lx|b-1_z-1_y)
//                                    v4 = x[l1][l2] * y[m1][m2-1] * z[n1][n2-1];
//                                }
//                            }
//                            buffer_[ao12+xbydisp] += -1.0 * n2 * (2.0 * a2 * (v1 + (A[1] - C[1]) * v2) - m2 * (v3 + (A[1] - C[1]) * v4)) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_z|Lx|b+1_y+1_y)
//                            v1 = x[l1][l2] * y[m1][m2+2] * z[n1+1][n2];
//                            // (a|Lx|b+1_y+1_y)
//                            v2 = x[l1][l2] * y[m1][m2+2] * z[n1][n2];
//                            // (a+1_z|Lx|b+1_y-1_y)
//                            v3 = x[l1][l2] * y[m1][m2] * z[n1+1][n2];
//                            // (a|Lx|b+1_y-1_y)
//                            v4 = x[l1][l2] * y[m1][m2] * z[n1][n2];
//                            buffer_[ao12+xbydisp] += -2.0 * a2 * (2.0 * a2 * (v1 + (A[2] - C[2]) * v2) - m2 * (v3 + (A[2] - C[2]) * v4)) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_z|Lx|b-1_y+1_y)
//                            v1 = x[l1][l2] * y[m1][m2] * z[n1+1][n2];
//                            // (a|Lx|b-1_y+1_y)
//                            v2 = x[l1][l2] * y[m1][m2] * z[n1][n2];
//                            if (m2 || m2-1) {
//                                // (a+1_z|Lx|b-1_y-1_y)
//                                v3 = x[l1][l2] * y[m1][m2-2] * z[n1+1][n2];
//                                // (a|Lx|b-1_y-1_y)
//                                v4 = x[l1][l2] * y[m1][m2-2] * z[n1][n2];
//                            }
//                            buffer_[ao12+xbydisp] += 1.0 * m2 * (2.0 * a2 * (v1 + (A[2] - C[2]) * v2) - m2 * (v3 + (A[2] - C[2]) * v4)) * over_pf;

                            //
                            // Bz derivatives with Lx
                            //

                            v1 = v2 = v3 = v4 = v5 = 0.0;
                            // (a+1_y|Lx|b+1_z+1_z)
                            v1 = x[l1][l2] * y[m1+1][m2] * z[n1][n2+2];
                            // (a|Lx|b+1_z+1_z)
                            v2 = x[l1][l2] * y[m1][m2] * z[n1][n2+2];
                            // (a+1_y|Lx|b+1_z-1_z)
                            v3 = x[l1][l2] * y[m1+1][m2] * z[n1][n2];
                            // (a|Lx|b+1_z-1_z)
                            v4 = x[l1][l2] * y[m1][m2] * z[n1][n2];
                            buffer_[ao12+xbzdisp] += 2.0 * a2 * (2.0 * a2 * (v1 + (A[1] - C[1]) * v2) - n2 * (v3 + (A[1] - C[1]) * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = v5 = 0.0;
                            // (a+1_y|Lx|b-1_z+1_z)
                            v1 = x[l1][l2] * y[m1+1][m2] * z[n1][n2];
                            // (a|Lx|b-1_z+1_z)
                            v2 = x[l1][l2] * y[m1][m2] * z[n1][n2];
                            if (n2 || n2-1) {
                                // (a+1_y|Lx|b-1_z-1_z)
                                v3 = x[l1][l2] * y[m1+1][m2] * z[n1][n2-2];
                                // (a|Lx|b-1_z-1_z)
                                v4 = x[l1][l2] * y[m1][m2] * z[n1][n2-2];
                            }
                            buffer_[ao12+xbzdisp] += -1.0 * n2 * (2.0 * a2 * (v1 + (A[1] - C[1]) * v2) - n2 * (v3 + (A[1] - C[1]) * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = v5 = 0.0;
                            // (a+1_z|Lx|b+1_y+1_z)
                            v1 = x[l1][l2] * y[m1][m2+1] * z[n1+1][n2+1];
                            // (a|Lx|b+1_y+1_z)
                            v2 = x[l1][l2] * y[m1][m2+1] * z[n1][n2+1];
                            if (n2) {
                                // (a+1_z|Lx|b+1_y-1_z)
                                v3 = x[l1][l2] * y[m1][m2+1] * z[n1+1][n2-1];
                                // (a|Lx|b+1_y-1_z)
                                v4 = x[l1][l2] * y[m1][m2+1] * z[n1][n2-1];
                            }
                            buffer_[ao12+xbzdisp] += -2.0 * a2 * (2.0 * a2 * (v1 + (A[2] - C[2]) * v2) - n2 * (v3 + (A[2] - C[2]) * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = v5 = 0.0;
                            if (m2) {
                                // (a+1_z|Lx|b-1_y+1_z)
                                v1 = x[l1][l2] * y[m1][m2-1] * z[n1+1][n2+1];
                                // (a|Lx|b-1_y+1_z)
                                v2 = x[l1][l2] * y[m1][m2-1] * z[n1][n2+1];
                                if (n2) {
                                    // (a+1_z|Lx|b-1_y-1_z)
                                    v3 = x[l1][l2] * y[m1][m2-1] * z[n1+1][n2-1];
                                    // (a|Lx|b-1_y-1_z)
                                    v4 = x[l1][l2] * y[m1][m2-1] * z[n1][n2-1];
                                }
                            }
                            buffer_[ao12+xbzdisp] += 1.0 * m2 * (2.0 * a2 * (v1 + (A[2] - C[2]) * v2) - n2 * (v3 + (A[2] - C[2]) * v4)) * over_pf;

//                            //
//                            // Ax derivatives with Ly
//                            //
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_z+1_x|Ly|b+1_x)
//                            v1 = x[l1+1][l2+1] * y[m1][m2] * z[n1+1][n2];
//                            // (a+1_x|Ly|b+1_x)
//                            v2 = x[l1+1][l2+1] * y[m1][m2] * z[n1][n2];
//                            if (l1) {
//                                // (a+1_z-1_x|Ly|b+1_x)
//                                v3 = x[l1-1][l2+1] * y[m1][m2] * z[n1+1][n2];
//                                // (a-1_x|Ly|b+1_x)
//                                v4 = x[l1-1][l2+1] * y[m1][m2] * z[n1][n2];
//                            }
//                            // (a|Ly|b+1_x)
//                            //v5 = x[l1][l2+1] * y[m1][m2] * z[n1][n2]; // because kronecker_delta(j,l) = (z,x) = 0
//                            buffer_[ao12+yaxdisp] += 2.0 * a2 * (2.0 * a1 * (v1 + (A[2] - C[2]) * v2) - l1 * (v3 + (A[2] - C[2]) * v4) + v5) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            if (l2) {
//                                // (a+1_z+1_x|Ly|b-1_x)
//                                v1 = x[l1+1][l2-1] * y[m1][m2] * z[n1+1][n2];
//                                // (a+1_x|Ly|b-1_x)
//                                v2 = x[l1+1][l2-1] * y[m1][m2] * z[n1][n2];
//                                if (l1) {
//                                    // (a+1_z-1_x|Ly|b-1_x)
//                                    v3 = x[l1-1][l2-1] * y[m1][m2] * z[n1+1][n2];
//                                    // (a-1_x|Ly|b-1_x)
//                                    v4 = x[l1-1][l2-1] * y[m1][m2] * z[n1][n2];
//                                }
//                                // (a|Ly|b-1_x)
//                                //v5 = x[l1][l2-1] * y[m1][m2] * z[n1][n2]; // because kronecker_delta(j,l) = (z,x) = 0
//                            }
//                            buffer_[ao12+yaxdisp] += -1.0 * l2 * (2.0 * a1 * (v1 + (A[2] - C[2]) * v2) - l1 * (v3 + (A[2] - C[2]) * v4) + v5) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_x+1_x|Ly|b+1_z)
//                            v1 = x[l1+2][l2] * y[m1][m2] * z[n1][n2+1];
//                            // (a+1_x|Ly|b+1_z)
//                            v2 = x[l1+1][l2] * y[m1][m2] * z[n1][n2+1];
//                            // (a+1_x-1_x|Ly|b+1_z)
//                            v3 = x[l1][l2] * y[m1][m2] * z[n1][n2+1];
//                            if (l1) {
//                                // (a-1_x|Ly|b+1_z)
//                                v4 = x[l1-1][l2] * y[m1][m2] * z[n1][n2+1];
//                            }
//                            // (a|Ly|b+1_z)
//                            v5 = x[l1][l2] * y[m1][m2] * z[n1][n2+1]; // because kronecker_delta(k,l) = (x,x) = 1
//                            buffer_[ao12+yaxdisp] += -2.0 * a2 * (2.0 * a1 * (v1 + (A[0] - C[0]) * v2) - l1 * (v3 + (A[0] - C[0]) * v4) + v5) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            if (n2) {
//                                // (a+1_x+1_x|Ly|b-1_z)
//                                v1 = x[l1+2][l2] * y[m1][m2] * z[n1][n2-1];
//                                // (a+1_x|Ly|b-1_z)
//                                v2 = x[l1+1][l2] * y[m1][m2] * z[n1][n2-1];
//                                // (a+1_x-1_x|Ly|b-1_z)
//                                v3 = x[l1][l2] * y[m1][m2] * z[n1][n2-1];
//                                if (l1) {
//                                    // (a-1_x|Ly|b-1_z)
//                                    v4 = x[l1-1][l2] * y[m1][m2] * z[n1][n2-1];
//                                }
//                                // (a|Ly|b-1_z)
//                                v5 = x[l1][l2] * y[m1][m2] * z[n1][n2-1]; // because kronecker_delta(k,l) = (x,x) = 1
//                            }
//                            buffer_[ao12+yaxdisp] += 1.0 * n2 * (2.0 * a1 * (v1 + (A[0] - C[0]) * v2) - l1 * (v3 + (A[0] - C[0]) * v4) + v5) * over_pf;
//
//                            //
//                            // Ay derivatives with Ly
//                            //
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_z+1_y|Ly|b+1_x)
//                            v1 = x[l1][l2+1] * y[m1+1][m2] * z[n1+1][n2];
//                            // (a+1_y|Ly|b+1_x)
//                            v2 = x[l1][l2+1] * y[m1+1][m2] * z[n1][n2];
//                            if (m1) {
//                                // (a+1_z-1_y|Ly|b+1_x)
//                                v3 = x[l1][l2+1] * y[m1-1][m2] * z[n1+1][n2];
//                                // (a-1_y|Ly|b+1_x)
//                                v4 = x[l1][l2+1] * y[m1-1][m2] * z[n1][n2];
//                            }
//                            // (a|Ly|b+1_x)
//                            //v5 = x[l1][l2+1] * y[m1][m2] * z[n1][n2]; // because kronecker_delta(j,l) = (z,y) = 0
//                            buffer_[ao12+yaydisp] += 2.0 * a2 * (2.0 * a1 * (v1 + (A[2] - C[2]) * v2) - m1 * (v3 + (A[2] - C[2]) * v4) + v5) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            if (l2) {
//                                // (a+1_z+1_y|Ly|b-1_x)
//                                v1 = x[l1][l2-1] * y[m1+1][m2] * z[n1+1][n2];
//                                // (a+1_y|Ly|b-1_x)
//                                v2 = x[l1][l2-1] * y[m1+1][m2] * z[n1][n2];
//                                if (m1) {
//                                    // (a+1_z-1_y|Ly|b-1_x)
//                                    v3 = x[l1][l2-1] * y[m1-1][m2] * z[n1+1][n2];
//                                    // (a-1_y|Ly|b-1_x)
//                                    v4 = x[l1][l2-1] * y[m1-1][m2] * z[n1][n2];
//                                }
//                                // (a|Ly|b-1_x)
//                                //v5 = x[l1][l2-1] * y[m1][m2] * z[n1][n2]; // because kronecker_delta(j,l) = (z,y) = 0
//                            }
//                            buffer_[ao12+yaydisp] += -1.0 * l2 * (2.0 * a1 * (v1 + (A[2] - C[2]) * v2) - m1 * (v3 + (A[2] - C[2]) * v4) + v5) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_x+1_y|Ly|b+1_z)
//                            v1 = x[l1+1][l2] * y[m1+1][m2] * z[n1][n2+1];
//                            // (a+1_y|Ly|b+1_z)
//                            v2 = x[l1][l2] * y[m1+1][m2] * z[n1][n2+1];
//                            if (m1) {
//                                // (a+1_x-1_y|Ly|b+1_z)
//                                v3 = x[l1+1][l2] * y[m1-1][m2] * z[n1][n2+1];
//                                // (a-1_y|Ly|b+1_z)
//                                v4 = x[l1][l2] * y[m1-1][m2] * z[n1][n2+1];
//                            }
//                            // (a|Ly|b+1_z)
//                            //v5 = x[l1][l2] * y[m1][m2] * z[n1][n2+1]; // because kronecker_delta(k,l) = (x,y) = 0
//                            buffer_[ao12+yaydisp] += -2.0 * a2 * (2.0 * a1 * (v1 + (A[0] - C[0]) * v2) - m1 * (v3 + (A[0] - C[0]) * v4) + v5) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            if (n2) {
//                                // (a+1_x+1_y|Ly|b-1_z)
//                                v1 = x[l1+1][l2] * y[m1+1][m2] * z[n1][n2-1];
//                                // (a+1_y|Ly|b-1_z)
//                                v2 = x[l1][l2] * y[m1+1][m2] * z[n1][n2-1];
//                                if (m1) {
//                                    // (a+1_x-1_y|Ly|b-1_z)
//                                    v3 = x[l1+1][l2] * y[m1-1][m2] * z[n1][n2-1];
//                                    // (a-1_y|Ly|b-1_z)
//                                    v4 = x[l1][l2] * y[m1-1][m2] * z[n1][n2-1];
//                                }
//                                // (a|Ly|b-1_z)
//                                //v5 = x[l1][l2] * y[m1][m2] * z[n1][n2-1]; // because kronecker_delta(k,l) = (x,y) = 0
//                            }
//                            buffer_[ao12+yaydisp] += 1.0 * n2 * (2.0 * a1 * (v1 + (A[0] - C[0]) * v2) - m1 * (v3 + (A[0] - C[0]) * v4) + v5) * over_pf;

                            //
                            // Az derivatives with Ly
                            //

                            v1 = v2 = v3 = v4 = v5 = 0.0;
                            // (a+1_z+1_z|Ly|b+1_x)
                            v1 = x[l1][l2+1] * y[m1][m2] * z[n1+2][n2];
                            // (a+1_z|Ly|b+1_x)
                            v2 = x[l1][l2+1] * y[m1][m2] * z[n1+1][n2];
                            // (a+1_z-1_z|Ly|b+1_x)
                            v3 = x[l1][l2+1] * y[m1][m2] * z[n1][n2];
                            if (n1) {
                                // (a-1_z|Ly|b+1_x)
                                v4 = x[l1][l2+1] * y[m1][m2] * z[n1-1][n2];
                            }
                            // (a|Ly|b+1_x)
                            v5 = x[l1][l2+1] * y[m1][m2] * z[n1][n2]; // because kronecker_delta(j,l) = (z,z) = 1
                            buffer_[ao12+yazdisp] += 2.0 * a2 * (2.0 * a1 * (v1 + (A[2] - C[2]) * v2) - n1 * (v3 + (A[2] - C[2]) * v4) + v5) * over_pf;

                            v1 = v2 = v3 = v4 = v5 = 0.0;
                            if (l2) {
                                // (a+1_z+1_z|Ly|b-1_x)
                                v1 = x[l1][l2-1] * y[m1][m2] * z[n1+2][n2];
                                // (a+1_z|Ly|b-1_x)
                                v2 = x[l1][l2-1] * y[m1][m2] * z[n1+1][n2];
                                // (a+1_z-1_z|Ly|b-1_x)
                                v3 = x[l1][l2-1] * y[m1][m2] * z[n1][n2];
                                if (n1) {
                                    // (a-1_z|Ly|b-1_x)
                                    v4 = x[l1][l2-1] * y[m1][m2] * z[n1-1][n2];
                                }
                                // (a|Ly|b-1_x)
                                v5 = x[l1][l2-1] * y[m1][m2] * z[n1][n2]; // because kronecker_delta(j,l) = (z,z) = 1
                            }
                            buffer_[ao12+yazdisp] += -1.0 * l2 * (2.0 * a1 * (v1 + (A[2] - C[2]) * v2) - n1 * (v3 + (A[2] - C[2]) * v4) + v5) * over_pf;

                            v1 = v2 = v3 = v4 = v5 = 0.0;
                            // (a+1_x+1_z|Ly|b+1_z)
                            v1 = x[l1+1][l2] * y[m1][m2] * z[n1+1][n2+1];
                            // (a+1_z|Ly|b+1_z)
                            v2 = x[l1][l2] * y[m1][m2] * z[n1+1][n2+1];
                            if (n1) {
                                // (a+1_x-1_z|Ly|b+1_z)
                                v3 = x[l1+1][l2] * y[m1][m2] * z[n1-1][n2+1];
                                // (a-1_z|Ly|b+1_z)
                                v4 = x[l1][l2] * y[m1][m2] * z[n1-1][n2+1];
                            }
                            // (a|Ly|b+1_z)
                            //v5 = x[l1][l2] * y[m1][m2] * z[n1][n2+1]; // because kronecker_delta(k,l) = (x,z) = 0
                            buffer_[ao12+yazdisp] += -2.0 * a2 * (2.0 * a1 * (v1 + (A[0] - C[0]) * v2) - n1 * (v3 + (A[0] - C[0]) * v4) + v5) * over_pf;

                            v1 = v2 = v3 = v4 = v5 = 0.0;
                            if (n2) {
                                // (a+1_x+1_z|Ly|b-1_z)
                                v1 = x[l1+1][l2] * y[m1][m2] * z[n1+1][n2-1];
                                // (a+1_z|Ly|b-1_z)
                                v2 = x[l1][l2] * y[m1][m2] * z[n1+1][n2-1];
                                if (n1) {
                                    // (a+1_x-1_z|Ly|b-1_z)
                                    v3 = x[l1+1][l2] * y[m1][m2] * z[n1-1][n2-1];
                                    // (a-1_z|Ly|b-1_z)
                                    v4 = x[l1][l2] * y[m1][m2] * z[n1-1][n2-1];
                                }
                                // (a|Ly|b-1_z)
                                //v5 = x[l1][l2] * y[m1][m2] * z[n1][n2-1]; // because kronecker_delta(k,l) = (x,z) = 0
                            }
                            buffer_[ao12+yazdisp] += 1.0 * n2 * (2.0 * a1 * (v1 + (A[0] - C[0]) * v2) - n1 * (v3 + (A[0] - C[0]) * v4) + v5) * over_pf;

//                            //
//                            // Bx derivatives with Ly
//                            //
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_z|Ly|b+1_x+1_x)
//                            v1 = x[l1][l2+2] * y[m1][m2] * z[n1+1][n2];
//                            // (a|Ly|b+1_x+1_x)
//                            v2 = x[l1][l2+2] * y[m1][m2] * z[n1][n2];
//                            // (a+1_z|Ly|b+1_x-1_x)
//                            v3 = x[l1][l2] * y[m1][m2] * z[n1+1][n2];
//                            // (a|Ly|b+1_x-1_x)
//                            v4 = x[l1][l2] * y[m1][m2] * z[n1][n2];
//                            buffer_[ao12+ybxdisp] += 2.0 * a2 * (2.0 * a2 * (v1 + (A[2] - C[2]) * v2) - l2 * (v3 + (A[2] - C[2]) * v4)) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_z|Ly|b-1_x+1_x)
//                            v1 = x[l1][l2] * y[m1][m2] * z[n1+1][n2];
//                            // (a|Ly|b-1_x+1_x)
//                            v2 = x[l1][l2] * y[m1][m2] * z[n1][n2];
//                            if (l2 || l2-1) {
//                                // (a+1_z|Ly|b-1_x-1_x)
//                                v3 = x[l1][l2-2] * y[m1][m2] * z[n1+1][n2];
//                                // (a|Ly|b-1_x-1_x)
//                                v4 = x[l1][l2-2] * y[m1][m2] * z[n1][n2];
//                            }
//                            buffer_[ao12+ybxdisp] += -1.0 * l2 * (2.0 * a2 * (v1 + (A[2] - C[2]) * v2) - l2 * (v3 + (A[2] - C[2]) * v4)) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_x|Ly|b+1_z+1_x)
//                            v1 = x[l1+1][l2+1] * y[m1][m2] * z[n1][n2+1];
//                            // (a|Ly|b+1_z+1_x)
//                            v2 = x[l1][l2+1] * y[m1][m2] * z[n1][n2+1];
//                            if (l2) {
//                                // (a+1_x|Ly|b+1_z-1_x)
//                                v3 = x[l1+1][l2-1] * y[m1][m2] * z[n1][n2+1];
//                                // (a|Ly|b+1_z-1_x)
//                                v4 = x[l1][l2-1] * y[m1][m2] * z[n1][n2+1];
//                            }
//                            buffer_[ao12+ybxdisp] += -2.0 * a2 * (2.0 * a2 * (v1 + (A[0] - C[0]) * v2) - l2 * (v3 + (A[0] - C[0]) * v4)) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            if (n2) {
//                                // (a+1_x|Ly|b-1_z+1_x)
//                                v1 = x[l1+1][l2+1] * y[m1][m2] * z[n1][n2-1];
//                                // (a|Ly|b-1_z+1_x)
//                                v2 = x[l1][l2+1] * y[m1][m2] * z[n1][n2-1];
//                                if (l2) {
//                                    // (a+1_x|Ly|b-1_z-1_x)
//                                    v3 = x[l1+1][l2-1] * y[m1][m2] * z[n1][n2-1];
//                                    // (a|Ly|b-1_z-1_x)
//                                    v4 = x[l1][l2-1] * y[m1][m2] * z[n1][n2-1];
//                                }
//                            }
//                            buffer_[ao12+ybxdisp] += 1.0 * n2 * (2.0 * a2 * (v1 + (A[0] - C[0]) * v2) - l2 * (v3 + (A[0] - C[0]) * v4)) * over_pf;
//
//                            //
//                            // By derivatives with Ly
//                            //
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_z|Ly|b+1_x+1_y)
//                            v1 = x[l1][l2+1] * y[m1][m2+1] * z[n1+1][n2];
//                            // (a|Ly|b+1_x+1_y)
//                            v2 = x[l1][l2+1] * y[m1][m2+1] * z[n1][n2];
//                            if (m2) {
//                                // (a+1_z|Ly|b+1_x-1_y)
//                                v3 = x[l1][l2+1] * y[m1][m2-1] * z[n1+1][n2];
//                                // (a|Ly|b+1_x-1_y)
//                                v4 = x[l1][l2+1] * y[m1][m2-1] * z[n1][n2];
//                            }
//                            buffer_[ao12+ybydisp] += 2.0 * a2 * (2.0 * a2 * (v1 + (A[2] - C[2]) * v2) - m2 * (v3 + (A[2] - C[2]) * v4)) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            if (l2) {
//                                // (a+1_z|Ly|b-1_x+1_y)
//                                v1 = x[l1][l2-1] * y[m1][m2+1] * z[n1+1][n2];
//                                // (a|Ly|b-1_x+1_y)
//                                v2 = x[l1][l2-1] * y[m1][m2+1] * z[n1][n2];
//                                if (m2) {
//                                    // (a+1_z|Ly|b-1_x-1_y)
//                                    v3 = x[l1][l2-1] * y[m1][m2-1] * z[n1+1][n2];
//                                    // (a|Ly|b-1_x-1_y)
//                                    v4 = x[l1][l2-1] * y[m1][m2-1] * z[n1][n2];
//                                }
//                            }
//                            buffer_[ao12+ybydisp] += -1.0 * l2 * (2.0 * a2 * (v1 + (A[2] - C[2]) * v2) - m2 * (v3 + (A[2] - C[2]) * v4)) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_x|Ly|b+1_z+1_y)
//                            v1 = x[l1+1][l2] * y[m1][m2+1] * z[n1][n2+1];
//                            // (a|Ly|b+1_z+1_y)
//                            v2 = x[l1][l2] * y[m1][m2+1] * z[n1][n2+1];
//                            if (m2) {
//                                // (a+1_x|Ly|b+1_z-1_y)
//                                v3 = x[l1+1][l2] * y[m1][m2-1] * z[n1][n2+1];
//                                // (a|Ly|b+1_z-1_y)
//                                v4 = x[l1][l2] * y[m1][m2-1] * z[n1][n2+1];
//                            }
//                            buffer_[ao12+ybydisp] += -2.0 * a2 * (2.0 * a2 * (v1 + (A[0] - C[0]) * v2) - m2 * (v3 + (A[0] - C[0]) * v4)) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            if (n2) {
//                                // (a+1_x|Ly|b-1_z+1_y)
//                                v1 = x[l1+1][l2] * y[m1][m2+1] * z[n1][n2-1];
//                                // (a|Ly|b-1_z+1_y)
//                                v2 = x[l1][l2] * y[m1][m2+1] * z[n1][n2-1];
//                                if (m2) {
//                                    // (a+1_x|Ly|b-1_z-1_y)
//                                    v3 = x[l1+1][l2] * y[m1][m2-1] * z[n1][n2-1];
//                                    // (a|Ly|b-1_z-1_y)
//                                    v4 = x[l1][l2] * y[m1][m2-1] * z[n1][n2-1];
//                                }
//                            }
//                            buffer_[ao12+ybydisp] += 1.0 * n2 * (2.0 * a2 * (v1 + (A[0] - C[0]) * v2) - m2 * (v3 + (A[0] - C[0]) * v4)) * over_pf;

                            //
                            // Bz derivatives with Ly
                            //

                            v1 = v2 = v3 = v4 = v5 = 0.0;
                            // (a+1_z|Ly|b+1_x+1_z)
                            v1 = x[l1][l2+1] * y[m1][m2] * z[n1+1][n2+1];
                            // (a|Ly|b+1_x+1_z)
                            v2 = x[l1][l2+1] * y[m1][m2] * z[n1][n2+1];
                            if (n2) {
                                // (a+1_z|Ly|b+1_x-1_z)
                                v3 = x[l1][l2+1] * y[m1][m2] * z[n1+1][n2-1];
                                // (a|Ly|b+1_x-1_z)
                                v4 = x[l1][l2+1] * y[m1][m2] * z[n1][n2-1];
                            }
                            buffer_[ao12+ybzdisp] += 2.0 * a2 * (2.0 * a2 * (v1 + (A[2] - C[2]) * v2) - n2 * (v3 + (A[2] - C[2]) * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = v5 = 0.0;
                            if (l2) {
                                // (a+1_z|Ly|b-1_x+1_z)
                                v1 = x[l1][l2-1] * y[m1][m2] * z[n1+1][n2+1];
                                // (a|Ly|b-1_x+1_z)
                                v2 = x[l1][l2-1] * y[m1][m2] * z[n1][n2+1];
                                if (n2) {
                                    // (a+1_z|Ly|b-1_x-1_z)
                                    v3 = x[l1][l2-1] * y[m1][m2] * z[n1+1][n2-1];
                                    // (a|Ly|b-1_x-1_z)
                                    v4 = x[l1][l2-1] * y[m1][m2] * z[n1][n2-1];
                                }
                            }
                            buffer_[ao12+ybzdisp] += -1.0 * l2 * (2.0 * a2 * (v1 + (A[2] - C[2]) * v2) - n2 * (v3 + (A[2] - C[2]) * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = v5 = 0.0;
                            // (a+1_x|Ly|b+1_z+1_z)
                            v1 = x[l1+1][l2] * y[m1][m2] * z[n1][n2+2];
                            // (a|Ly|b+1_z+1_z)
                            v2 = x[l1][l2] * y[m1][m2] * z[n1][n2+2];
                            // (a+1_x|Ly|b+1_z-1_z)
                            v3 = x[l1+1][l2] * y[m1][m2] * z[n1][n2];
                            // (a|Ly|b+1_z-1_z)
                            v4 = x[l1][l2] * y[m1][m2] * z[n1][n2];
                            buffer_[ao12+ybzdisp] += -2.0 * a2 * (2.0 * a2 * (v1 + (A[0] - C[0]) * v2) - n2 * (v3 + (A[0] - C[0]) * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = v5 = 0.0;
                            // (a+1_x|Ly|b-1_z+1_z)
                            v1 = x[l1+1][l2] * y[m1][m2] * z[n1][n2];
                            // (a|Ly|b-1_z+1_z)
                            v2 = x[l1][l2] * y[m1][m2] * z[n1][n2];
                            if (n2 || n2-1) {
                                // (a+1_x|Ly|b-1_z-1_z)
                                v3 = x[l1+1][l2] * y[m1][m2] * z[n1][n2-2];
                                // (a|Ly|b-1_z-1_z)
                                v4 = x[l1][l2] * y[m1][m2] * z[n1][n2-2];
                            }
                            buffer_[ao12+ybzdisp] += 1.0 * n2 * (2.0 * a2 * (v1 + (A[0] - C[0]) * v2) - n2 * (v3 + (A[0] - C[0]) * v4)) * over_pf;

//                            //
//                            // Ax derivatives with Lz
//                            //
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_x+1_x|Lz|b+1_y)
//                            v1 = x[l1+2][l2] * y[m1][m2+1] * z[n1][n2];
//                            // (a+1_x|Lz|b+1_y)
//                            v2 = x[l1+1][l2] * y[m1][m2+1] * z[n1][n2];
//                            // (a+1_x-1_x|Lz|b+1_y)
//                            v3 = x[l1][l2] * y[m1][m2+1] * z[n1][n2];
//                            if (l1) {
//                                // (a-1_x|Lz|b+1_y)
//                                v4 = x[l1-1][l2] * y[m1][m2+1] * z[n1][n2];
//                            }
//                            // (a|Lz|b+1_y)
//                            v5 = x[l1][l2] * y[m1][m2+1] * z[n1][n2]; // because kronecker_delta(j,l) = (x,x) = 1
//                            buffer_[ao12+zaxdisp] += 2.0 * a2 * (2.0 * a1 * (v1 + (A[0] - C[0]) * v2) - l1 * (v3 + (A[0] - C[0]) * v4) + v5) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            if (m2) {
//                                // (a+1_x+1_x|Lz|b-1_y)
//                                v1 = x[l1+2][l2] * y[m1][m2-1] * z[n1][n2];
//                                // (a+1_x|Lz|b-1_y)
//                                v2 = x[l1+1][l2] * y[m1][m2-1] * z[n1][n2];
//                                // (a+1_x-1_x|Lz|b-1_y)
//                                v3 = x[l1][l2] * y[m1][m2-1] * z[n1][n2];
//                                if (l1) {
//                                    // (a-1_x|Lz|b-1_y)
//                                    v4 = x[l1-1][l2] * y[m1][m2-1] * z[n1][n2];
//                                }
//                                // (a|Lz|b-1_y)
//                                v5 = x[l1][l2] * y[m1][m2-1] * z[n1][n2]; // because kronecker_delta(j,l) = (x,x) = 1
//                            }
//                            buffer_[ao12+zaxdisp] += -1.0 * m2 * (2.0 * a1 * (v1 + (A[0] - C[0]) * v2) - l1 * (v3 + (A[0] - C[0]) * v4) + v5) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_y+1_x|Lz|b+1_x)
//                            v1 = x[l1+1][l2+1] * y[m1+1][m2] * z[n1][n2];
//                            // (a+1_x|Lz|b+1_x)
//                            v2 = x[l1+1][l2+1] * y[m1][m2] * z[n1][n2];
//                            if (l1) {
//                                // (a+1_y-1_x|Lz|b+1_x)
//                                v3 = x[l1-1][l2+1] * y[m1+1][m2] * z[n1][n2];
//                                // (a-1_x|Lz|b+1_x)
//                                v4 = x[l1-1][l2+1] * y[m1][m2] * z[n1][n2];
//                            }
//                            // (a|Lz|b+1_x)
//                            //v5 = x[l1][l2+1] * y[m1][m2] * z[n1][n2]; // because kronecker_delta(k,l) = (y,x) = 0
//                            buffer_[ao12+zaxdisp] += -2.0 * a2 * (2.0 * a1 * (v1 + (A[1] - C[1]) * v2) - l1 * (v3 + (A[1] - C[1]) * v4) + v5) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            if (l2) {
//                                // (a+1_y+1_x|Lz|b-1_x)
//                                v1 = x[l1+1][l2-1] * y[m1+1][m2] * z[n1][n2];
//                                // (a+1_x|Lz|b-1_x)
//                                v2 = x[l1+1][l2-1] * y[m1][m2] * z[n1][n2];
//                                if (l1) {
//                                    // (a+1_y-1_x|Lz|b-1_x)
//                                    v3 = x[l1-1][l2-1] * y[m1+1][m2] * z[n1][n2];
//                                    // (a-1_x|Lz|b-1_x)
//                                    v4 = x[l1-1][l2-1] * y[m1][m2] * z[n1][n2];
//                                }
//                                // (a|Lz|b-1_x)
//                                //v5 = x[l1][l2-1] * y[m1][m2] * z[n1][n2]; // because kronecker_delta(k,l) = (y,x) = 0
//                            }
//                            buffer_[ao12+zaxdisp] += 1.0 * l2 * (2.0 * a1 * (v1 + (A[1] - C[1]) * v2) - l1 * (v3 + (A[1] - C[1]) * v4) + v5) * over_pf;
//
//                            //
//                            // Ay derivatives with Lz
//                            //
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_x+1_y|Lz|b+1_y)
//                            v1 = x[l1+1][l2] * y[m1+1][m2+1] * z[n1][n2];
//                            // (a+1_y|Lz|b+1_y)
//                            v2 = x[l1][l2] * y[m1+1][m2+1] * z[n1][n2];
//                            if (m1) {
//                                // (a+1_x-1_y|Lz|b+1_y)
//                                v3 = x[l1+1][l2] * y[m1-1][m2+1] * z[n1][n2];
//                                // (a-1_y|Lz|b+1_y)
//                                v4 = x[l1][l2] * y[m1-1][m2+1] * z[n1][n2];
//                            }
//                            // (a|Lz|b+1_y)
//                            //v5 = x[l1][l2] * y[m1][m2+1] * z[n1][n2]; // because kronecker_delta(j,l) = (x,y) = 0
//                            buffer_[ao12+zaydisp] += 2.0 * a2 * (2.0 * a1 * (v1 + (A[0] - C[0]) * v2) - m1 * (v3 + (A[0] - C[0]) * v4) + v5) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            if (m2) {
//                                // (a+1_x+1_y|Lz|b-1_y)
//                                v1 = x[l1+1][l2] * y[m1+1][m2-1] * z[n1][n2];
//                                // (a+1_y|Lz|b-1_y)
//                                v2 = x[l1][l2] * y[m1+1][m2-1] * z[n1][n2];
//                                if (m1) {
//                                    // (a+1_x-1_y|Lz|b-1_y)
//                                    v3 = x[l1+1][l2] * y[m1-1][m2-1] * z[n1][n2];
//                                    // (a-1_y|Lz|b-1_y)
//                                    v4 = x[l1][l2] * y[m1-1][m2-1] * z[n1][n2];
//                                }
//                                // (a|Lz|b-1_y)
//                                //v5 = x[l1][l2] * y[m1][m2-1] * z[n1][n2]; // because kronecker_delta(j,l) = (x,y) = 0
//                            }
//                            buffer_[ao12+zaydisp] += -1.0 * m2 * (2.0 * a1 * (v1 + (A[0] - C[0]) * v2) - m1 * (v3 + (A[0] - C[0]) * v4) + v5) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_y+1_y|Lz|b+1_x)
//                            v1 = x[l1][l2+1] * y[m1+2][m2] * z[n1][n2];
//                            // (a+1_y|Lz|b+1_x)
//                            v2 = x[l1][l2+1] * y[m1+1][m2] * z[n1][n2];
//                            // (a+1_y-1_y|Lz|b+1_x)
//                            v3 = x[l1][l2+1] * y[m1][m2] * z[n1][n2];
//                            if (m1) {
//                                // (a-1_y|Lz|b+1_x)
//                                v4 = x[l1][l2+1] * y[m1-1][m2] * z[n1][n2];
//                            }
//                            // (a|Lz|b+1_x)
//                            v5 = x[l1][l2+1] * y[m1][m2] * z[n1][n2]; // because kronecker_delta(k,l) = (y,y) = 1
//                            buffer_[ao12+zaydisp] += -2.0 * a2 * (2.0 * a1 * (v1 + (A[1] - C[1]) * v2) - m1 * (v3 + (A[1] - C[1]) * v4) + v5) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            if (l2) {
//                                // (a+1_y+1_y|Lz|b-1_x)
//                                v1 = x[l1][l2-1] * y[m1+2][m2] * z[n1][n2];
//                                // (a+1_y|Lz|b-1_x)
//                                v2 = x[l1][l2-1] * y[m1+1][m2] * z[n1][n2];
//                                // (a+1_y-1_y|Lz|b-1_x)
//                                v3 = x[l1][l2-1] * y[m1][m2] * z[n1][n2];
//                                if (m1) {
//                                    // (a-1_y|Lz|b-1_x)
//                                    v4 = x[l1][l2-1] * y[m1-1][m2] * z[n1][n2];
//                                }
//                                // (a|Lz|b-1_x)
//                                v5 = x[l1][l2-1] * y[m1][m2] * z[n1][n2]; // because kronecker_delta(k,l) = (y,y) = 1
//                            }
//                            buffer_[ao12+zaydisp] += 1.0 * l2 * (2.0 * a1 * (v1 + (A[1] - C[1]) * v2) - m1 * (v3 + (A[1] - C[1]) * v4) + v5) * over_pf;

                            ////
                            //// Az derivatives with Lz
                            ////

                            //v1 = v2 = v3 = v4 = v5 = 0.0;
                            //// (a+1_x+1_z|Lz|b+1_y)
                            //v1 = x[l1+1][l2] * y[m1][m2+1] * z[n1+1][n2];
                            //// (a+1_z|Lz|b+1_y)
                            //v2 = x[l1][l2] * y[m1][m2+1] * z[n1+1][n2];
                            //if (n1) {
                            //    // (a+1_x-1_z|Lz|b+1_y)
                            //    v3 = x[l1+1][l2] * y[m1][m2+1] * z[n1-1][n2];
                            //    // (a-1_z|Lz|b+1_y)
                            //    v4 = x[l1][l2] * y[m1][m2+1] * z[n1-1][n2];
                            //}
                            //// (a|Lz|b+1_y)
                            ////v5 = x[l1][l2] * y[m1][m2+1] * z[n1][n2]; // because kronecker_delta(j,l) = (x,z) = 0
                            //buffer_[ao12+zazdisp] += 2.0 * a2 * (2.0 * a1 * (v1 + (A[0] - C[0]) * v2) - n1 * (v3 + (A[0] - C[0]) * v4) + v5) * over_pf;

                            //v1 = v2 = v3 = v4 = v5 = 0.0;
                            //if (m2) {
                            //    // (a+1_x+1_z|Lz|b-1_y)
                            //    v1 = x[l1+1][l2] * y[m1][m2-1] * z[n1+1][n2];
                            //    // (a+1_z|Lz|b-1_y)
                            //    v2 = x[l1][l2] * y[m1][m2-1] * z[n1+1][n2];
                            //    if (n1) {
                            //        // (a+1_x-1_z|Lz|b-1_y)
                            //        v3 = x[l1+1][l2] * y[m1][m2-1] * z[n1-1][n2];
                            //        // (a-1_z|Lz|b-1_y)
                            //        v4 = x[l1][l2] * y[m1][m2-1] * z[n1-1][n2];
                            //    }
                            //    // (a|Lz|b-1_y)
                            //    //v5 = x[l1][l2] * y[m1][m2-1] * z[n1][n2]; // because kronecker_delta(j,l) = (x,z) = 0
                            //}
                            //buffer_[ao12+zazdisp] += -1.0 * m2 * (2.0 * a1 * (v1 + (A[0] - C[0]) * v2) - n1 * (v3 + (A[0] - C[0]) * v4) + v5) * over_pf;

                            //v1 = v2 = v3 = v4 = v5 = 0.0;
                            //// (a+1_y+1_z|Lz|b+1_x)
                            //v1 = x[l1][l2+1] * y[m1+1][m2] * z[n1+1][n2];
                            //// (a+1_z|Lz|b+1_x)
                            //v2 = x[l1][l2+1] * y[m1][m2] * z[n1+1][n2];
                            //if (n1) {
                            //    // (a+1_y-1_z|Lz|b+1_x)
                            //    v3 = x[l1][l2+1] * y[m1+1][m2] * z[n1-1][n2];
                            //    // (a-1_z|Lz|b+1_x)
                            //    v4 = x[l1][l2+1] * y[m1][m2] * z[n1-1][n2];
                            //}
                            //// (a|Lz|b+1_x)
                            ////v5 = x[l1][l2+1] * y[m1][m2] * z[n1][n2]; // because kronecker_delta(k,l) = (y,z) = 0
                            //buffer_[ao12+zazdisp] += -2.0 * a2 * (2.0 * a1 * (v1 + (A[1] - C[1]) * v2) - n1 * (v3 + (A[1] - C[1]) * v4) + v5) * over_pf;

                            //v1 = v2 = v3 = v4 = v5 = 0.0;
                            //if (l2) {
                            //    // (a+1_y+1_z|Lz|b-1_x)
                            //    v1 = x[l1][l2-1] * y[m1+1][m2] * z[n1+1][n2];
                            //    // (a+1_z|Lz|b-1_x)
                            //    v2 = x[l1][l2-1] * y[m1][m2] * z[n1+1][n2];
                            //    if (n1) {
                            //        // (a+1_y-1_z|Lz|b-1_x)
                            //        v3 = x[l1][l2-1] * y[m1+1][m2] * z[n1-1][n2];
                            //        // (a-1_z|Lz|b-1_x)
                            //        v4 = x[l1][l2-1] * y[m1][m2] * z[n1-1][n2];
                            //    }
                            //    // (a|Lz|b-1_x)
                            //    //v5 = x[l1][l2-1] * y[m1][m2] * z[n1][n2]; // because kronecker_delta(k,l) = (y,z) = 0
                            //}
                            //buffer_[ao12+zazdisp] += 1.0 * l2 * (2.0 * a1 * (v1 + (A[1] - C[1]) * v2) - n1 * (v3 + (A[1] - C[1]) * v4) + v5) * over_pf;
                            
//                            //
//                            // Bx derivatives with Lz
//                            //
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_x|Lz|b+1_y+1_x)
//                            v1 = x[l1+1][l2+1] * y[m1][m2+1] * z[n1][n2];
//                            // (a|Lz|b+1_y+1_x)
//                            v2 = x[l1][l2+1] * y[m1][m2+1] * z[n1][n2];
//                            if (l2) {
//                                // (a+1_x|Lz|b+1_y-1_x)
//                                v3 = x[l1+1][l2-1] * y[m1][m2+1] * z[n1][n2];
//                                // (a|Lz|b+1_y-1_x)
//                                v4 = x[l1][l2-1] * y[m1][m2+1] * z[n1][n2];
//                            }
//                            buffer_[ao12+zbxdisp] += 2.0 * a2 * (2.0 * a2 * (v1 + (A[0] - C[0]) * v2) - l2 * (v3 + (A[0] - C[0]) * v4)) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            if (m2) {
//                                // (a+1_x|Lz|b-1_y+1_x)
//                                v1 = x[l1+1][l2+1] * y[m1][m2-1] * z[n1][n2];
//                                // (a|Lz|b-1_y+1_x)
//                                v2 = x[l1][l2+1] * y[m1][m2-1] * z[n1][n2];
//                                if (l2) {
//                                    // (a+1_x|Lz|b-1_y-1_x)
//                                    v3 = x[l1+1][l2-1] * y[m1][m2-1] * z[n1][n2];
//                                    // (a|Lz|b-1_y-1_x)
//                                    v4 = x[l1][l2-1] * y[m1][m2-1] * z[n1][n2];
//                                }
//                            }
//                            buffer_[ao12+zbxdisp] += -1.0 * m2 * (2.0 * a2 * (v1 + (A[0] - C[0]) * v2) - l2 * (v3 + (A[0] - C[0]) * v4)) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_y|Lz|b+1_x+1_x)
//                            v1 = x[l1][l2+2] * y[m1+1][m2] * z[n1][n2];
//                            // (a|Lz|b+1_x+1_x)
//                            v2 = x[l1][l2+2] * y[m1][m2] * z[n1][n2];
//                            // (a+1_y|Lz|b+1_x-1_x)
//                            v3 = x[l1][l2] * y[m1+1][m2] * z[n1][n2];
//                            // (a|Lz|b+1_x-1_x)
//                            v4 = x[l1][l2] * y[m1][m2] * z[n1][n2];
//                            buffer_[ao12+zbxdisp] += -2.0 * a2 * (2.0 * a2 * (v1 + (A[1] - C[1]) * v2) - l2 * (v3 + (A[1] - C[1]) * v4)) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_y|Lz|b-1_x+1_x)
//                            v1 = x[l1][l2] * y[m1+1][m2] * z[n1][n2];
//                            // (a|Lz|b-1_x+1_x)
//                            v2 = x[l1][l2] * y[m1][m2] * z[n1][n2];
//                            if (l2 || l2-1) {
//                                // (a+1_y|Lz|b-1_x-1_x)
//                                v3 = x[l1][l2-2] * y[m1+1][m2] * z[n1][n2];
//                                // (a|Lz|b-1_x-1_x)
//                                v4 = x[l1][l2-2] * y[m1][m2] * z[n1][n2];
//                            }
//                            buffer_[ao12+zbxdisp] += 1.0 * l2 * (2.0 * a2 * (v1 + (A[1] - C[1]) * v2) - l2 * (v3 + (A[1] - C[1]) * v4)) * over_pf;
//
//                            //
//                            // By derivatives with Lz
//                            //
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_x|Lz|b+1_y+1_y)
//                            v1 = x[l1+1][l2] * y[m1][m2+2] * z[n1][n2];
//                            // (a|Lz|b+1_y+1_y)
//                            v2 = x[l1][l2] * y[m1][m2+2] * z[n1][n2];
//                            // (a+1_x|Lz|b+1_y-1_y)
//                            v3 = x[l1+1][l2] * y[m1][m2] * z[n1][n2];
//                            // (a|Lz|b+1_y-1_y)
//                            v4 = x[l1][l2] * y[m1][m2] * z[n1][n2];
//                            buffer_[ao12+zbydisp] += 2.0 * a2 * (2.0 * a2 * (v1 + (A[0] - C[0]) * v2) - m2 * (v3 + (A[0] - C[0]) * v4)) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_x|Lz|b-1_y+1_y)
//                            v1 = x[l1+1][l2] * y[m1][m2] * z[n1][n2];
//                            // (a|Lz|b-1_y+1_y)
//                            v2 = x[l1][l2] * y[m1][m2] * z[n1][n2];
//                            if (m2 || m2-1) {
//                                // (a+1_x|Lz|b-1_y-1_y)
//                                v3 = x[l1+1][l2] * y[m1][m2-2] * z[n1][n2];
//                                // (a|Lz|b-1_y-1_y)
//                                v4 = x[l1][l2] * y[m1][m2-2] * z[n1][n2];
//                            }
//                            buffer_[ao12+zbydisp] += -1.0 * m2 * (2.0 * a2 * (v1 + (A[0] - C[0]) * v2) - m2 * (v3 + (A[0] - C[0]) * v4)) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            // (a+1_y|Lz|b+1_x+1_y)
//                            v1 = x[l1][l2+1] * y[m1+1][m2+1] * z[n1][n2];
//                            // (a|Lz|b+1_x+1_y)
//                            v2 = x[l1][l2+1] * y[m1][m2+1] * z[n1][n2];
//                            if (m2) {
//                                // (a+1_y|Lz|b+1_x-1_y)
//                                v3 = x[l1][l2+1] * y[m1+1][m2-1] * z[n1][n2];
//                                // (a|Lz|b+1_x-1_y)
//                                v4 = x[l1][l2+1] * y[m1][m2-1] * z[n1][n2];
//                            }
//                            buffer_[ao12+zbydisp] += -2.0 * a2 * (2.0 * a2 * (v1 + (A[1] - C[1]) * v2) - m2 * (v3 + (A[1] - C[1]) * v4)) * over_pf;
//
//                            v1 = v2 = v3 = v4 = v5 = 0.0;
//                            if (l2) {
//                                // (a+1_y|Lz|b-1_x+1_y)
//                                v1 = x[l1][l2-1] * y[m1+1][m2+1] * z[n1][n2];
//                                // (a|Lz|b-1_x+1_y)
//                                v2 = x[l1][l2-1] * y[m1][m2+1] * z[n1][n2];
//                                if (m2) {
//                                    // (a+1_y|Lz|b-1_x-1_y)
//                                    v3 = x[l1][l2-1] * y[m1+1][m2-1] * z[n1][n2];
//                                    // (a|Lz|b-1_x-1_y)
//                                    v4 = x[l1][l2-1] * y[m1][m2-1] * z[n1][n2];
//                                }
//                            }
//                            buffer_[ao12+zbydisp] += 1.0 * l2 * (2.0 * a2 * (v1 + (A[1] - C[1]) * v2) - m2 * (v3 + (A[1] - C[1]) * v4)) * over_pf;

                            ////
                            //// Bz derivatives with Lz
                            ////

                            //v1 = v2 = v3 = v4 = v5 = 0.0;
                            //// (a+1_x|Lz|b+1_y+1_z)
                            //v1 = x[l1+1][l2] * y[m1][m2+1] * z[n1][n2+1];
                            //// (a|Lz|b+1_y+1_z)
                            //v2 = x[l1][l2] * y[m1][m2+1] * z[n1][n2+1];
                            //if (n2) {
                            //    // (a+1_x|Lz|b+1_y-1_z)
                            //    v3 = x[l1+1][l2] * y[m1][m2+1] * z[n1][n2-1];
                            //    // (a|Lz|b+1_y-1_z)
                            //    v4 = x[l1][l2] * y[m1][m2+1] * z[n1][n2-1];
                            //}
                            //buffer_[ao12+zbzdisp] += 2.0 * a2 * (2.0 * a2 * (v1 + (A[0] - C[0]) * v2) - n2 * (v3 + (A[0] - C[0]) * v4)) * over_pf;

                            //v1 = v2 = v3 = v4 = v5 = 0.0;
                            //if (m2) {
                            //    // (a+1_x|Lz|b-1_y+1_z)
                            //    v1 = x[l1+1][l2] * y[m1][m2-1] * z[n1][n2+1];
                            //    // (a|Lz|b-1_y+1_z)
                            //    v2 = x[l1][l2] * y[m1][m2-1] * z[n1][n2+1];
                            //    if (n2) {
                            //        // (a+1_x|Lz|b-1_y-1_z)
                            //        v3 = x[l1+1][l2] * y[m1][m2-1] * z[n1][n2-1];
                            //        // (a|Lz|b-1_y-1_z)
                            //        v4 = x[l1][l2] * y[m1][m2-1] * z[n1][n2-1];
                            //    }
                            //}
                            //buffer_[ao12+zbzdisp] += -1.0 * m2 * (2.0 * a2 * (v1 + (A[0] - C[0]) * v2) - n2 * (v3 + (A[0] - C[0]) * v4)) * over_pf;

                            //v1 = v2 = v3 = v4 = v5 = 0.0;
                            //// (a+1_y|Lz|b+1_x+1_z)
                            //v1 = x[l1][l2+1] * y[m1+1][m2] * z[n1][n2+1];
                            //// (a|Lz|b+1_x+1_z)
                            //v2 = x[l1][l2+1] * y[m1][m2] * z[n1][n2+1];
                            //if (n2) {
                            //    // (a+1_y|Lz|b+1_x-1_z)
                            //    v3 = x[l1][l2+1] * y[m1+1][m2] * z[n1][n2-1];
                            //    // (a|Lz|b+1_x-1_z)
                            //    v4 = x[l1][l2+1] * y[m1][m2] * z[n1][n2-1];
                            //}
                            //buffer_[ao12+zbzdisp] += -2.0 * a2 * (2.0 * a2 * (v1 + (A[1] - C[1]) * v2) - n2 * (v3 + (A[1] - C[1]) * v4)) * over_pf;

                            //v1 = v2 = v3 = v4 = v5 = 0.0;
                            //if (l2) {
                            //    // (a+1_y|Lz|b-1_x+1_z)
                            //    v1 = x[l1][l2-1] * y[m1+1][m2] * z[n1][n2+1];
                            //    // (a|Lz|b-1_x+1_z)
                            //    v2 = x[l1][l2-1] * y[m1][m2] * z[n1][n2+1];
                            //    if (n2) {
                            //        // (a+1_y|Lz|b-1_x-1_z)
                            //        v3 = x[l1][l2-1] * y[m1+1][m2] * z[n1][n2-1];
                            //        // (a|Lz|b-1_x-1_z)
                            //        v4 = x[l1][l2-1] * y[m1][m2] * z[n1][n2-1];
                            //    }
                            //}
                            //buffer_[ao12+zbzdisp] += 1.0 * l2 * (2.0 * a2 * (v1 + (A[1] - C[1]) * v2) - n2 * (v3 + (A[1] - C[1]) * v4)) * over_pf;




                            //TEST

                            //
                            // Az derivatives with Lz
                            //

                            v1 = v2 = v3 = v4 = v5 = v6 = v7 = v8 = v9 = v10 = v11 = v12 = v13 = v14 = v15 = v16 = v17 = v18 = v19 = v20 = 0.0;

                            // (a+1_x+1_z|Lz|b+1_y)
                            v1 = x[l1+1][l2] * y[m1][m2+1] * z[n1+1][n2];
                            // (a+1_z|Lz|b+1_y)
                            v2 = x[l1][l2] * y[m1][m2+1] * z[n1+1][n2];
                            if (n1) {
                                // (a+1_x-1_z|Lz|b+1_y)
                                v3 = x[l1+1][l2] * y[m1][m2+1] * z[n1-1][n2];
                                // (a-1_z|Lz|b+1_y)
                                v4 = x[l1][l2] * y[m1][m2+1] * z[n1-1][n2];
                            }
                            // (a|Lz|b+1_y)
                            //v5 = x[l1][l2] * y[m1][m2+1] * z[n1][n2]; // because kronecker_delta(j,l) = (x,z) = 0

                            if (m2) {
                                // (a+1_x+1_z|Lz|b-1_y)
                                v6 = x[l1+1][l2] * y[m1][m2-1] * z[n1+1][n2];
                                // (a+1_z|Lz|b-1_y)
                                v7 = x[l1][l2] * y[m1][m2-1] * z[n1+1][n2];
                                if (n1) {
                                    // (a+1_x-1_z|Lz|b-1_y)
                                    v8 = x[l1+1][l2] * y[m1][m2-1] * z[n1-1][n2];
                                    // (a-1_z|Lz|b-1_y)
                                    v9 = x[l1][l2] * y[m1][m2-1] * z[n1-1][n2];
                                }
                                // (a|Lz|b-1_y)
                                //v10 = x[l1][l2] * y[m1][m2-1] * z[n1][n2]; // because kronecker_delta(j,l) = (x,z) = 0
                            }

                            // (a+1_y+1_z|Lz|b+1_x)
                            v11 = x[l1][l2+1] * y[m1+1][m2] * z[n1+1][n2];
                            // (a+1_z|Lz|b+1_x)
                            v12 = x[l1][l2+1] * y[m1][m2] * z[n1+1][n2];
                            if (n1) {
                                // (a+1_y-1_z|Lz|b+1_x)
                                v13 = x[l1][l2+1] * y[m1+1][m2] * z[n1-1][n2];
                                // (a-1_z|Lz|b+1_x)
                                v14 = x[l1][l2+1] * y[m1][m2] * z[n1-1][n2];
                            }
                            // (a|Lz|b+1_x)
                            //v15 = x[l1][l2+1] * y[m1][m2] * z[n1][n2]; // because kronecker_delta(k,l) = (y,z) = 0

                            if (l2) {
                                // (a+1_y+1_z|Lz|b-1_x)
                                v16 = x[l1][l2-1] * y[m1+1][m2] * z[n1+1][n2];
                                // (a+1_z|Lz|b-1_x)
                                v17 = x[l1][l2-1] * y[m1][m2] * z[n1+1][n2];
                                if (n1) {
                                    // (a+1_y-1_z|Lz|b-1_x)
                                    v18 = x[l1][l2-1] * y[m1+1][m2] * z[n1-1][n2];
                                    // (a-1_z|Lz|b-1_x)
                                    v19 = x[l1][l2-1] * y[m1][m2] * z[n1-1][n2];
                                }
                                // (a|Lz|b-1_x)
                                //v20 = x[l1][l2-1] * y[m1][m2] * z[n1][n2]; // because kronecker_delta(k,l) = (y,z) = 0
                            }
                            //std::cout << "Az derivatives with Lz:" << std::endl;
                            buffer_[ao12+zazdisp] += 2.0 * a2 * (2.0 * a1 * (v1 + (A[0] - C[0]) * v2) - n1 * (v3 + (A[0] - C[0]) * v4) + v5) * over_pf;
                            //std::cout << "buffer_[ao12+zazdisp] += " << 2.0 * a2 * (2.0 * a1 * (v1 + (A[0] - C[0]) * v2) - n1 * (v3 + (A[0] - C[0]) * v4) + v5) * over_pf << std::endl;
                            buffer_[ao12+zazdisp] += -1.0 * m2 * (2.0 * a1 * (v6 + (A[0] - C[0]) * v7) - n1 * (v3 + (A[0] - C[0]) * v9) + v10) * over_pf;
                            //std::cout << "buffer_[ao12+zazdisp] += " << -1.0 * m2 * (2.0 * a1 * (v6 + (A[0] - C[0]) * v7) - n1 * (v3 + (A[0] - C[0]) * v9) + v10) * over_pf << std::endl;
                            buffer_[ao12+zazdisp] += -2.0 * a2 * (2.0 * a1 * (v11 + (A[1] - C[1]) * v12) - n1 * (v13 + (A[1] - C[1]) * v14) + v15) * over_pf;
                            //std::cout << "buffer_[ao12+zazdisp] += " << -2.0 * a2 * (2.0 * a1 * (v11 + (A[1] - C[1]) * v12) - n1 * (v13 + (A[1] - C[1]) * v14) + v15) * over_pf << std::endl;
                            buffer_[ao12+zazdisp] += 1.0 * l2 * (2.0 * a1 * (v16 + (A[1] - C[1]) * v17) - n1 * (v18 + (A[1] - C[1]) * v19) + v20) * over_pf;
                            //std::cout << "buffer_[ao12+zazdisp] += " << 1.0 * l2 * (2.0 * a1 * (v16 + (A[1] - C[1]) * v17) - n1 * (v18 + (A[1] - C[1]) * v19) + v20) * over_pf << std::endl;

                            //
                            // Bz derivatives with Lz
                            //

                            v1 = v2 = v3 = v4 = v5 = v6 = v7 = v8 = v9 = v10 = v11 = v12 = v13 = v14 = v15 = v16 = v17 = v18 = v19 = v20 = 0.0;

                            // (a+1_x|Lz|b+1_y+1_z)
                            v1 = x[l1+1][l2] * y[m1][m2+1] * z[n1][n2+1];
                            // (a|Lz|b+1_y+1_z)
                            v2 = x[l1][l2] * y[m1][m2+1] * z[n1][n2+1];
                            if (n2) {
                                // (a+1_x|Lz|b+1_y-1_z)
                                v3 = x[l1+1][l2] * y[m1][m2+1] * z[n1][n2-1];
                                // (a|Lz|b+1_y-1_z)
                                v4 = x[l1][l2] * y[m1][m2+1] * z[n1][n2-1];
                            }

                            if (m2) {
                                // (a+1_x|Lz|b-1_y+1_z)
                                v6 = x[l1+1][l2] * y[m1][m2-1] * z[n1][n2+1];
                                // (a|Lz|b-1_y+1_z)
                                v7 = x[l1][l2] * y[m1][m2-1] * z[n1][n2+1];
                                if (n2) {
                                    // (a+1_x|Lz|b-1_y-1_z)
                                    v8 = x[l1+1][l2] * y[m1][m2-1] * z[n1][n2-1];
                                    // (a|Lz|b-1_y-1_z)
                                    v9 = x[l1][l2] * y[m1][m2-1] * z[n1][n2-1];
                                }
                            }

                            // (a+1_y|Lz|b+1_x+1_z)
                            v11 = x[l1][l2+1] * y[m1+1][m2] * z[n1][n2+1];
                            // (a|Lz|b+1_x+1_z)
                            v12 = x[l1][l2+1] * y[m1][m2] * z[n1][n2+1];
                            if (n2) {
                                // (a+1_y|Lz|b+1_x-1_z)
                                v13 = x[l1][l2+1] * y[m1+1][m2] * z[n1][n2-1];
                                // (a|Lz|b+1_x-1_z)
                                v14 = x[l1][l2+1] * y[m1][m2] * z[n1][n2-1];
                            }

                            if (l2) {
                                // (a+1_y|Lz|b-1_x+1_z)
                                v16 = x[l1][l2-1] * y[m1+1][m2] * z[n1][n2+1];
                                // (a|Lz|b-1_x+1_z)
                                v17 = x[l1][l2-1] * y[m1][m2] * z[n1][n2+1];
                                if (n2) {
                                    // (a+1_y|Lz|b-1_x-1_z)
                                    v18 = x[l1][l2-1] * y[m1+1][m2] * z[n1][n2-1];
                                    // (a|Lz|b-1_x-1_z)
                                    v19 = x[l1][l2-1] * y[m1][m2] * z[n1][n2-1];
                                }
                            }
                            //std::cout << "Bz derivatives with Lz:" << std::endl;
                            buffer_[ao12+zbzdisp] += 2.0 * a2 * (2.0 * a2 * (v1 + (A[0] - C[0]) * v2) - n2 * (v3 + (A[0] - C[0]) * v4)) * over_pf;
                            //std::cout << "buffer_[ao12+zbzdisp] += " << 2.0 * a2 * (2.0 * a2 * (v1 + (A[0] - C[0]) * v2) - n2 * (v3 + (A[0] - C[0]) * v4)) * over_pf << std::endl;
                            buffer_[ao12+zbzdisp] += -1.0 * m2 * (2.0 * a2 * (v6 + (A[0] - C[0]) * v7) - n2 * (v8 + (A[0] - C[0]) * v9)) * over_pf;
                            //std::cout << "buffer_[ao12+zbzdisp] += " << -1.0 * m2 * (2.0 * a2 * (v6 + (A[0] - C[0]) * v7) - n2 * (v8 + (A[0] - C[0]) * v9)) * over_pf << std::endl;
                            buffer_[ao12+zbzdisp] += -2.0 * a2 * (2.0 * a2 * (v11 + (A[1] - C[1]) * v12) - n2 * (v13 + (A[1] - C[1]) * v14)) * over_pf;
                            //std::cout << "buffer_[ao12+zbzdisp] += " << -2.0 * a2 * (2.0 * a2 * (v11 + (A[1] - C[1]) * v12) - n2 * (v13 + (A[1] - C[1]) * v14)) * over_pf << std::endl;
                            buffer_[ao12+zbzdisp] += 1.0 * l2 * (2.0 * a2 * (v16 + (A[1] - C[1]) * v17) - n2 * (v18 + (A[1] - C[1]) * v19)) * over_pf;
                            //std::cout << "buffer_[ao12+zbzdisp] += " << 1.0 * l2 * (2.0 * a2 * (v16 + (A[1] - C[1]) * v17) - n2 * (v18 + (A[1] - C[1]) * v19)) * over_pf << std::endl;
                            //std::printf("\n");

                            ao12++;
                        }
                    }
                }
            }
            //std::printf("\n\n\n");
        }
    }
    //std::printf("\n\n\n\n\n");
}
