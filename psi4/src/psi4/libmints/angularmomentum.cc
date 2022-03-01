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

#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/angularmomentum.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/integral.h"

#include <memory>
#include <stdexcept>

#include <libint2/shell.h>

using namespace psi;

// Initialize overlap_recur_ to +1 basis set angular momentum, +1 on each center is sufficient
// to compute the dipole derivatives
AngularMomentumInt::AngularMomentumInt(std::vector<SphericalTransform>& spherical_transforms,
                                       std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2, int nderiv)
    : OneBodyAOInt(spherical_transforms, bs1, bs2, nderiv), overlap_recur_(bs1->max_am() + 1, bs2->max_am() + 1) {
    int maxam1 = bs1_->max_am();
    int maxam2 = bs2_->max_am();

    int maxnao1 = (maxam1 + 1) * (maxam1 + 2) / 2;
    int maxnao2 = (maxam2 + 1) * (maxam2 + 2) / 2;

    // Increase buffer size to handle x, y, and z components
    if (deriv_ == 0) {
        buffer_ = new double[3 * maxnao1 * maxnao2];
        set_chunks(3);
    } else {
        throw PSIEXCEPTION("AngularMomentumInt does not provide derivatives");
    }
    buffers_.resize(nchunk_);
}

AngularMomentumInt::~AngularMomentumInt() { delete[] buffer_; }


void AngularMomentumInt::compute_pair(const libint2::Shell &s1, const libint2::Shell &s2) {
    int ao12;
    int am1 = s1.contr[0].l;
    int am2 = s2.contr[0].l;
    int nprim1 = s1.nprim();
    int nprim2 = s2.nprim();
    double A[3], B[3];
    double Sxy1, Sxy2, Sxz1, Sxz2;
    double Syx1, Syx2, Syz1, Syz2;
    double Szx1, Szx2, Szy1, Szy2;
    double S0x1, S0x2, S0y1, S0y2, S0z1, S0z2;
    double muxy1, muxy2, muxz1, muxz2;
    double muyx1, muyx2, muyz1, muyz2;
    double muzx1, muzx2, muzy1, muzy2;

    A[0] = s1.O[0];
    A[1] = s1.O[1];
    A[2] = s1.O[2];
    B[0] = s2.O[0];
    B[1] = s2.O[1];
    B[2] = s2.O[2];

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
                            double Lx = (2.0 * a2 * muyz1 - n2 * muyz2 - 2.0 * a2 * muzy1 + m2 * muzy2) /* * over_pf*/;

                            /* (a|Ly|b) = 2 a2 * (a|(z-Cz)|b+1x) - B.x * (a|(z-Cz)|b-1x)
                               - 2 a2 * (a|(x-Cx)|b+1z) + B.z * (a|(x-Cx)|b-1z) */
                            double Ly = (2.0 * a2 * muzx1 - l2 * muzx2 - 2.0 * a2 * muxz1 + n2 * muxz2) /* * over_pf*/;

                            /* (a|Lz|b) = 2 a2 * (a|(x-Cx)|b+1y) - B.y * (a|(x-Cx)|b-1y)
                               - 2 a2 * (a|(y-Cy)|b+1x) + B.x * (a|(y-Cy)|b-1x) */
                            double Lz = (2.0 * a2 * muxy1 - m2 * muxy2 - 2.0 * a2 * muyx1 + l2 * muyx2) /* * over_pf*/;

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
    pure_transform(s1, s2, 3);
    buffers_[0] = buffer_ + 0 * s1.size() * s2.size();
    buffers_[1] = buffer_ + 1 * s1.size() * s2.size();
    buffers_[2] = buffer_ + 2 * s1.size() * s2.size();
}
