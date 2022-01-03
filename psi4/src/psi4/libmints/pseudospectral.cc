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
#include "psi4/libmints/pseudospectral.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/physconst.h"
#include "libint2/shell.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))

;
using namespace psi;

// Initialize potential_recur_ to +1 basis set angular momentum
PseudospectralInt::PseudospectralInt(std::vector<SphericalTransform>& st, std::shared_ptr<BasisSet> bs1,
                                     std::shared_ptr<BasisSet> bs2, int deriv)
    : OneBodyAOInt(st, bs1, bs2, deriv),
      potential_recur_(bs1->max_am() + 1, bs2->max_am() + 1) {
    int maxam1 = bs1_->max_am();
    int maxam2 = bs2_->max_am();

    use_omega_ = false;
    omega_ = 0.0;
    C_[0] = 0.0;
    C_[1] = 0.0;
    C_[2] = 0.0;
    int maxnao1 = INT_NCART(maxam1);
    int maxnao2 = INT_NCART(maxam2);

    buffer_ = new double[maxnao1 * maxnao2];
}

PseudospectralInt::~PseudospectralInt() { delete[] buffer_; }


// The engine only supports segmented basis sets
void PseudospectralInt::compute_pair(const libint2::Shell& s1, const libint2::Shell& s2) {
    int ao12;
    int am1 = s1.contr[0].l;
    int am2 = s2.contr[0].l;
    int nprim1 = s1.nprim();
    int nprim2 = s2.nprim();
    double A[3], B[3];
    A[0] = s1.O[0];
    A[1] = s1.O[1];
    A[2] = s1.O[2];
    B[0] = s2.O[0];
    B[1] = s2.O[1];
    B[2] = s2.O[2];

    int izm = 1;
    int iym = am1 + 1;
    int ixm = iym * iym;
    int jzm = 1;
    int jym = am2 + 1;
    int jxm = jym * jym;

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    memset(buffer_, 0, s1.cartesian_size() * s2.cartesian_size() * sizeof(double));

    double*** vi = potential_recur_.vi();

    for (int p1 = 0; p1 < nprim1; ++p1) {
        double a1 = s1.alpha[p1];
        double c1 = s1.contr[0].coeff[p1];
        for (int p2 = 0; p2 < nprim2; ++p2) {
            double a2 = s2.alpha[p2];
            double c2 = s2.contr[0].coeff[p2];

            // The GPT gamma of the two cartesian functions.
            // This is used for all GPT operations here, such as over_pf and P,
            // And for OS relations AFTER the modified (0|A|0)^(m) integrals are built
            double gamma0 = a1 + a2;
            double oog = 1.0 / gamma0;

            // An effective gamma if range-separation is to be used.
            // This gamma is only for use in the generation of auxiliary integrals,
            // particularly the (0|A|0)^(m) auxiliary integrals as built in potential_recur_.compute().
            double gamma = gamma0;
            if (use_omega_) {
                gamma = gamma0 * omega_ * omega_ / (gamma0 + omega_ * omega_);
            }

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

            // Loop over atoms of basis set 1 (only works if bs1_ and bs2_ are on the same
            // molecule)
            double PC[3];

            // C_ is the pseudospectral grid point
            PC[0] = P[0] - C_[0];
            PC[1] = P[1] - C_[1];
            PC[2] = P[2] - C_[2];

            // Do recursion
            potential_recur_.compute_erf(PA, PB, PC, gamma0, am1, am2, gamma);

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

                            // Compute location in the recursion
                            int iind = l1 * ixm + m1 * iym + n1 * izm;
                            int jind = l2 * jxm + m2 * jym + n2 * jzm;

                            buffer_[ao12++] += vi[iind][jind][0] * over_pf;

                            //                            outfile->Printf( "ao12=%d, vi[%d][%d][0] = %20.14f, over_pf =
                            //                            %20.14f, Z = %f\n", ao12-1, iind, jind, vi[iind][jind][0],
                            //                            over_pf, Z);
                        }
                    }
                }
            }
        }
    }
    pure_transform(s1, s2, 1);
    buffers_[0] = buffer_;
}
