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

#include <cmath>
#include <stdexcept>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/wavefunction.h"   // for df
#include "psi4/libmints/osrecur.h"
#include "psi4/libpsi4util/exception.h"

using namespace psi;

double ***init_box(int a, int b, int c)
{
    int i,j;
    double ***box;

    box = (double ***) malloc(sizeof(double **)*a);
    for(i=0;i<a;i++)
        box[i] = (double **) malloc(sizeof(double *)*b);
    for(i=0;i<a;i++) {
        for(j=0;j<b;j++) {
            box[i][j] = (double *) malloc(sizeof(double)*c);
            memset((void *)box[i][j], '\0', sizeof(double)*c);
        }
    }
    return box;
}

void zero_box(double ***box, int a, int b, int c)
{
    int i, j;
    for (i=0; i<a; ++i) {
        for (j=0; j<b; ++j) {
            memset((void*)box[i][j], 0, sizeof(double)*c);
        }
    }
}

void free_box(double ***box, int a, int b)
{
    int i,j;

    for(i=0;i<a;i++)
        for(j=0;j<b;j++)
            free(box[i][j]);

    for(i=0;i<a;i++)
        free(box[i]);

    free(box);
}

ObaraSaikaTwoCenterMIRecursion::ObaraSaikaTwoCenterMIRecursion(int max_am1, int max_am2, int max_m):
    max_am1_(max_am1), max_am2_(max_am2), max_m_(max_m)
{
    if (max_am1 < 0)
        throw SanityCheckError("ObaraSaikaTwoCenterMIRecursion -- max_am1 must be nonnegative", __FILE__, __LINE__);
    if (max_am2 < 0)
        throw SanityCheckError("ERROR: ObaraSaikaTwoCenterMIRecursion -- max_am2 must be nonnegative", __FILE__, __LINE__);
//    if (max_m > 3)
//        throw SanityCheckError("ERROR: ObaraSaikaTwoCenterMIRecursion -- max_m must be nonnegative and less than 4", __FILE__, __LINE__);

    x_ = init_box(max_am1+1, max_am2+1, max_m+1);
    y_ = init_box(max_am1+1, max_am2+1, max_m+1);
    z_ = init_box(max_am1+1, max_am2+1, max_m+1);
}

ObaraSaikaTwoCenterMIRecursion::~ObaraSaikaTwoCenterMIRecursion()
{
    free_box(x_, max_am1_+1, max_am2_+1);
    free_box(y_, max_am1_+1, max_am2_+1);
    free_box(z_, max_am1_+1, max_am2_+1);
}

void ObaraSaikaTwoCenterMIRecursion::compute(double PA[3], double PB[3], double gamma, int am1, int am2)
{
    if (am1 < 0 || am1 > max_am1_)
        throw SanityCheckError("ERROR: ObaraSaikaTwoCenterMIRecursion::compute -- am1 out of bounds", __FILE__, __LINE__);
    if (am2 < 0 || am2 > max_am2_)
        throw SanityCheckError("ERROR: ObaraSaikaTwoCenterMIRecursion::compute -- am2 out of bounds", __FILE__, __LINE__);

    int i, j, k;
    double oog = 1.0 / (2.0 * gamma);


    // Generate the fundamental integrals.  N.B. Only even parity terms survive!
    x_[0][0][0] = y_[0][0][0] = z_[0][0][0] = 1.0;
    for(int n = 1; n < max_m_; n += 2) {
        x_[0][0][n+1] = n * oog * x_[0][0][n-1];
        y_[0][0][n+1] = n * oog * y_[0][0][n-1];
        z_[0][0][n+1] = n * oog * z_[0][0][n-1];
    }

    // Upward recursion in j for i=0
    for (j=0; j<am2; ++j) {
        for (k=0; k<=max_m_; ++k) {
            x_[0][j+1][k] = PB[0] * x_[0][j][k];
            y_[0][j+1][k] = PB[1] * y_[0][j][k];
            z_[0][j+1][k] = PB[2] * z_[0][j][k];

            if (j > 0) {
                x_[0][j+1][k] += j * oog * x_[0][j-1][k];
                y_[0][j+1][k] += j * oog * y_[0][j-1][k];
                z_[0][j+1][k] += j * oog * z_[0][j-1][k];
            }

            if (k > 0) {
                x_[0][j+1][k] += k * oog * x_[0][j][k-1];
                y_[0][j+1][k] += k * oog * y_[0][j][k-1];
                z_[0][j+1][k] += k * oog * z_[0][j][k-1];
            }
        }
    }

    // Upward recursion in i for all j's
    for (i=0; i<am1; ++i) {
        for (j=0; j<=am2; ++j) {
            for (k=0; k<=max_m_; ++k) {
                x_[i+1][j][k] = PA[0] * x_[i][j][k];
                y_[i+1][j][k] = PA[1] * y_[i][j][k];
                z_[i+1][j][k] = PA[2] * z_[i][j][k];

                if (i > 0) {
                    x_[i+1][j][k] += i * oog * x_[i-1][j][k];
                    y_[i+1][j][k] += i * oog * y_[i-1][j][k];
                    z_[i+1][j][k] += i * oog * z_[i-1][j][k];
                }

                if (j > 0) {
                    x_[i+1][j][k] += j * oog * x_[i][j-1][k];
                    y_[i+1][j][k] += j * oog * y_[i][j-1][k];
                    z_[i+1][j][k] += j * oog * z_[i][j-1][k];
                }

                if (k > 0) {
                    x_[i+1][j][k] += k * oog * x_[i][j][k-1];
                    y_[i+1][j][k] += k * oog * y_[i][j][k-1];
                    z_[i+1][j][k] += k * oog * z_[i][j][k-1];
                }
            }
        }
    }
}

ObaraSaikaTwoCenterEFPRecursion::ObaraSaikaTwoCenterEFPRecursion(int max_am1, int max_am2):
    max_am1_(max_am1), max_am2_(max_am2)
{
    if (max_am1 < 0)
        throw SanityCheckError("ERROR: ObaraSaikaTwoCenterMVIRecursion -- max_am1 must be nonnegative", __FILE__, __LINE__);
    if (max_am2 < 0)
        throw SanityCheckError("ERROR: ObaraSaikaTwoCenterMVIRecursion -- max_am2 must be nonnegative", __FILE__, __LINE__);

    size_ = max_am1 > max_am2 ? max_am1 : max_am2;
    size_ += 1;
    size_ = (size_-1)*size_*(size_+1)+1;
    q_   = init_box(size_, size_, max_am1_ + max_am2_ + 4);
    x_   = init_box(size_, size_, max_am1_ + max_am2_ + 3);
    y_   = init_box(size_, size_, max_am1_ + max_am2_ + 3);
    z_   = init_box(size_, size_, max_am1_ + max_am2_ + 3);
    xx_  = init_box(size_, size_, max_am1_ + max_am2_ + 2);
    yy_  = init_box(size_, size_, max_am1_ + max_am2_ + 2);
    zz_  = init_box(size_, size_, max_am1_ + max_am2_ + 2);
    xy_  = init_box(size_, size_, max_am1_ + max_am2_ + 2);
    xz_  = init_box(size_, size_, max_am1_ + max_am2_ + 2);
    yz_  = init_box(size_, size_, max_am1_ + max_am2_ + 2);
    xxx_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    yyy_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    zzz_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    xxy_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    xxz_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    xyy_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    yyz_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    xzz_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    yzz_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    xyz_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
}

ObaraSaikaTwoCenterEFPRecursion::~ObaraSaikaTwoCenterEFPRecursion()
{
    free_box(q_  , size_, size_);
    free_box(x_  , size_, size_);
    free_box(y_  , size_, size_);
    free_box(z_  , size_, size_);
    free_box(xx_ , size_, size_);
    free_box(yy_ , size_, size_);
    free_box(zz_ , size_, size_);
    free_box(xy_ , size_, size_);
    free_box(xz_ , size_, size_);
    free_box(yz_ , size_, size_);
    free_box(xxx_, size_, size_);
    free_box(yyy_, size_, size_);
    free_box(zzz_, size_, size_);
    free_box(xxy_, size_, size_);
    free_box(xxz_, size_, size_);
    free_box(xyy_, size_, size_);
    free_box(yyz_, size_, size_);
    free_box(xzz_, size_, size_);
    free_box(yzz_, size_, size_);
    free_box(xyz_, size_, size_);
}

#define EPS 1.0e-17

void ObaraSaikaTwoCenterEFPRecursion::calculate_f(double *F, int n, double t)
{
    int i, m;
    int m2;
    double t2;
    double num;
    double sum;
    double term1;
    static double K = 1.0/M_2_SQRTPI;
    double et;


    if (t>20.0){
        t2 = 2*t;
        et = exp(-t);
        t = sqrt(t);
        F[0] = K*erf(t)/t;
        for(m=0; m<=n-1; m++){
            F[m+1] = ((2*m + 1)*F[m] - et)/(t2);
        }
    }
    else {
        et = exp(-t);
        t2 = 2*t;
        m2 = 2*n;
        num = df[m2];
        i=0;
        sum = 1.0/(m2+1);
        do{
            i++;
            num = num*t2;
            term1 = num/df[m2+2*i+2];
            sum += term1;
        } while (fabs(term1) > EPS && i < MAX_FAC);
        F[n] = sum*et;
        for(m=n-1;m>=0;m--){
            F[m] = (t2*F[m+1] + et)/(2*m+1);
        }
    }
}

void ObaraSaikaTwoCenterEFPRecursion::compute(double PA[3], double PB[3], double PC[3], double zeta, int am1, int am2)
{
    int a, b, m;
    int azm = 1;
    int aym = am1 + 1;
    int axm = aym * aym;
    int bzm = 1;
    int bym = am2 + 1;
    int bxm = bym * bym;
    int ax, ay, az, bx, by, bz;
    int aind, bind;
    double ooz = 1.0/(2.0 * zeta);
    int mmax = am1 + am2 + 3;

    // Prefactor from A20
    double tmp = sqrt(zeta) * M_2_SQRTPI;
    // U from A21
    double u = zeta * (PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2]);
    double *F = new double[mmax+1]; // TODO: Move this allocation into constructor

    // Zero out F
    memset(F, 0, sizeof(double) * (mmax+1));

    // Form Fm(U) from A20
    calculate_f(F, mmax, u);

    // Perform recursion in m for (a|A(0)|s) using A20
    for (m=0; m<=mmax; ++m) {
        q_[0][0][m] = tmp * F[m];
    }
    for (m=0; m<=mmax-1; ++m) {
        x_[0][0][m] = 2.0*zeta*PC[0]*q_[0][0][m+1];
        y_[0][0][m] = 2.0*zeta*PC[1]*q_[0][0][m+1];
        z_[0][0][m] = 2.0*zeta*PC[2]*q_[0][0][m+1];
    }
    for (m=0; m<=mmax-2; ++m) {
        xx_[0][0][m] = 4.0*zeta*zeta*PC[0]*PC[0]*q_[0][0][m+2] - 2.0*zeta*q_[0][0][m+1];
        yy_[0][0][m] = 4.0*zeta*zeta*PC[1]*PC[1]*q_[0][0][m+2] - 2.0*zeta*q_[0][0][m+1];
        zz_[0][0][m] = 4.0*zeta*zeta*PC[2]*PC[2]*q_[0][0][m+2] - 2.0*zeta*q_[0][0][m+1];
        xy_[0][0][m] = 4.0*zeta*zeta*PC[0]*PC[1]*q_[0][0][m+2];
        xz_[0][0][m] = 4.0*zeta*zeta*PC[0]*PC[2]*q_[0][0][m+2];
        yz_[0][0][m] = 4.0*zeta*zeta*PC[1]*PC[2]*q_[0][0][m+2];
    }
    for (m=0; m<=mmax-3; ++m) {
        xxx_[0][0][m] = 8.0*zeta*zeta*zeta*PC[0]*PC[0]*PC[0]*q_[0][0][m+3] - 12.0*zeta*zeta*PC[0]*q_[0][0][m+2];
        yyy_[0][0][m] = 8.0*zeta*zeta*zeta*PC[1]*PC[1]*PC[1]*q_[0][0][m+3] - 12.0*zeta*zeta*PC[1]*q_[0][0][m+2];
        zzz_[0][0][m] = 8.0*zeta*zeta*zeta*PC[2]*PC[2]*PC[2]*q_[0][0][m+3] - 12.0*zeta*zeta*PC[2]*q_[0][0][m+2];
        xxy_[0][0][m] = 8.0*zeta*zeta*zeta*PC[0]*PC[0]*PC[1]*q_[0][0][m+3] -  4.0*zeta*zeta*PC[1]*q_[0][0][m+2];
        xxz_[0][0][m] = 8.0*zeta*zeta*zeta*PC[0]*PC[0]*PC[2]*q_[0][0][m+3] -  4.0*zeta*zeta*PC[2]*q_[0][0][m+2];
        xyy_[0][0][m] = 8.0*zeta*zeta*zeta*PC[0]*PC[1]*PC[1]*q_[0][0][m+3] -  4.0*zeta*zeta*PC[0]*q_[0][0][m+2];
        yyz_[0][0][m] = 8.0*zeta*zeta*zeta*PC[1]*PC[1]*PC[2]*q_[0][0][m+3] -  4.0*zeta*zeta*PC[2]*q_[0][0][m+2];
        xzz_[0][0][m] = 8.0*zeta*zeta*zeta*PC[0]*PC[2]*PC[2]*q_[0][0][m+3] -  4.0*zeta*zeta*PC[0]*q_[0][0][m+2];
        yzz_[0][0][m] = 8.0*zeta*zeta*zeta*PC[1]*PC[2]*PC[2]*q_[0][0][m+3] -  4.0*zeta*zeta*PC[1]*q_[0][0][m+2];
        xyz_[0][0][m] = 8.0*zeta*zeta*zeta*PC[0]*PC[1]*PC[2]*q_[0][0][m+3];
    }

    // Perform recursion in b with a=0
    //  subset of A19
    for (b=1; b<=am2; ++b) {
        for (bx=0; bx<=b; ++bx) {
            for (by=0; by<=b-bx; ++by) {
                bz = b-bx-by;

                // Compute the index into VI for bx,by,bz
                bind = bx*bxm + by*bym + bz*bzm;

                // Compute each x, y, z contribution
                if (bz > 0) {
                    for (m=0; m<=mmax-b; ++m) { /* Electrostatic potential integrals */
                        q_[0][bind][m] = PB[2] * q_[0][bind-bzm][m] - PC[2] * q_[0][bind-bzm][m+1];
                    }
                    for (m=0; m<=mmax-b-1; ++m) { /* Electric field integrals */
                        x_[0][bind][m] = PB[2] * x_[0][bind-bzm][m] - PC[2] * x_[0][bind-bzm][m+1];
                        y_[0][bind][m] = PB[2] * y_[0][bind-bzm][m] - PC[2] * y_[0][bind-bzm][m+1];
                        z_[0][bind][m] = PB[2] * z_[0][bind-bzm][m] - PC[2] * z_[0][bind-bzm][m+1] + q_[0][bind-bzm][m+1];
                    }
                    for (m=0; m<=mmax-b-2; ++m) { /* Gradients of the electric field */
                        xx_[0][bind][m] = PB[2]*xx_[0][bind-bzm][m] - PC[2]*xx_[0][bind-bzm][m+1];
                        yy_[0][bind][m] = PB[2]*yy_[0][bind-bzm][m] - PC[2]*yy_[0][bind-bzm][m+1];
                        zz_[0][bind][m] = PB[2]*zz_[0][bind-bzm][m] - PC[2]*zz_[0][bind-bzm][m+1] + 2*z_[0][bind-bzm][m+1];
                        xy_[0][bind][m] = PB[2]*xy_[0][bind-bzm][m] - PC[2]*xy_[0][bind-bzm][m+1];
                        xz_[0][bind][m] = PB[2]*xz_[0][bind-bzm][m] - PC[2]*xz_[0][bind-bzm][m+1] + x_[0][bind-bzm][m+1];
                        yz_[0][bind][m] = PB[2]*yz_[0][bind-bzm][m] - PC[2]*yz_[0][bind-bzm][m+1] + y_[0][bind-bzm][m+1];
                    }
                    for (m=0; m <=mmax-b-3; ++m) { /* Hessians of the electric field */
                        xxx_[0][bind][m] = PB[2]*xxx_[0][bind-bzm][m] - PC[2] * xxx_[0][bind-bzm][m+1];
                        yyy_[0][bind][m] = PB[2]*yyy_[0][bind-bzm][m] - PC[2] * yyy_[0][bind-bzm][m+1];
                        zzz_[0][bind][m] = PB[2]*zzz_[0][bind-bzm][m] - PC[2] * zzz_[0][bind-bzm][m+1] + 3 * zz_[0][bind-bzm][m+1];
                        xxy_[0][bind][m] = PB[2]*xxy_[0][bind-bzm][m] - PC[2] * xxy_[0][bind-bzm][m+1];
                        xxz_[0][bind][m] = PB[2]*xxz_[0][bind-bzm][m] - PC[2] * xxz_[0][bind-bzm][m+1] + 1 * xx_[0][bind-bzm][m+1];
                        xyy_[0][bind][m] = PB[2]*xyy_[0][bind-bzm][m] - PC[2] * xyy_[0][bind-bzm][m+1];
                        yyz_[0][bind][m] = PB[2]*yyz_[0][bind-bzm][m] - PC[2] * yyz_[0][bind-bzm][m+1] + 1 * yy_[0][bind-bzm][m+1];
                        xzz_[0][bind][m] = PB[2]*xzz_[0][bind-bzm][m] - PC[2] * xzz_[0][bind-bzm][m+1] + 2 * xz_[0][bind-bzm][m+1];
                        yzz_[0][bind][m] = PB[2]*yzz_[0][bind-bzm][m] - PC[2] * yzz_[0][bind-bzm][m+1] + 2 * yz_[0][bind-bzm][m+1];
                        xyz_[0][bind][m] = PB[2]*xyz_[0][bind-bzm][m] - PC[2] * xyz_[0][bind-bzm][m+1] + 1 * xy_[0][bind-bzm][m+1];
                    }
                    if (bz > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            q_[0][bind][m] += ooz * (bz-1) * (q_[0][bind-2*bzm][m] - q_[0][bind-2*bzm][m+1]);
                        }
                        for (m=0; m<=mmax-b-1; ++m) {
                            x_[0][bind][m] += ooz * (bz-1) * (x_[0][bind-2*bzm][m] - x_[0][bind-2*bzm][m+1]);
                            y_[0][bind][m] += ooz * (bz-1) * (y_[0][bind-2*bzm][m] - y_[0][bind-2*bzm][m+1]);
                            z_[0][bind][m] += ooz * (bz-1) * (z_[0][bind-2*bzm][m] - z_[0][bind-2*bzm][m+1]);
                        }
                        for (m=0; m<=mmax-b-2; m++) {
                            xx_[0][bind][m] += ooz*(bz-1)*(xx_[0][bind-2*bzm][m] - xx_[0][bind-2*bzm][m+1]);
                            yy_[0][bind][m] += ooz*(bz-1)*(yy_[0][bind-2*bzm][m] - yy_[0][bind-2*bzm][m+1]);
                            zz_[0][bind][m] += ooz*(bz-1)*(zz_[0][bind-2*bzm][m] - zz_[0][bind-2*bzm][m+1]);
                            xy_[0][bind][m] += ooz*(bz-1)*(xy_[0][bind-2*bzm][m] - xy_[0][bind-2*bzm][m+1]);
                            xz_[0][bind][m] += ooz*(bz-1)*(xz_[0][bind-2*bzm][m] - xz_[0][bind-2*bzm][m+1]);
                            yz_[0][bind][m] += ooz*(bz-1)*(yz_[0][bind-2*bzm][m] - yz_[0][bind-2*bzm][m+1]);
                        }
                        for (m=0; m<=mmax-b-3; m++) {
                            xxx_[0][bind][m] += ooz*(bz-1)*(xxx_[0][bind-2*bzm][m] - xxx_[0][bind-2*bzm][m+1]);
                            yyy_[0][bind][m] += ooz*(bz-1)*(yyy_[0][bind-2*bzm][m] - yyy_[0][bind-2*bzm][m+1]);
                            zzz_[0][bind][m] += ooz*(bz-1)*(zzz_[0][bind-2*bzm][m] - zzz_[0][bind-2*bzm][m+1]);
                            xxy_[0][bind][m] += ooz*(bz-1)*(xxy_[0][bind-2*bzm][m] - xxy_[0][bind-2*bzm][m+1]);
                            xxz_[0][bind][m] += ooz*(bz-1)*(xxz_[0][bind-2*bzm][m] - xxz_[0][bind-2*bzm][m+1]);
                            xyy_[0][bind][m] += ooz*(bz-1)*(xyy_[0][bind-2*bzm][m] - xyy_[0][bind-2*bzm][m+1]);
                            yyz_[0][bind][m] += ooz*(bz-1)*(yyz_[0][bind-2*bzm][m] - yyz_[0][bind-2*bzm][m+1]);
                            xzz_[0][bind][m] += ooz*(bz-1)*(xzz_[0][bind-2*bzm][m] - xzz_[0][bind-2*bzm][m+1]);
                            yzz_[0][bind][m] += ooz*(bz-1)*(yzz_[0][bind-2*bzm][m] - yzz_[0][bind-2*bzm][m+1]);
                            xyz_[0][bind][m] += ooz*(bz-1)*(xyz_[0][bind-2*bzm][m] - xyz_[0][bind-2*bzm][m+1]);
                        }
                    }
                }
                else if (by > 0) {
                    for (m=0; m<=mmax-b; ++m) {
                        q_[0][bind][m] = PB[1] * q_[0][bind-bym][m] - PC[1] * q_[0][bind-bym][m+1];
                    }
                    for (m=0; m<=mmax-b-1; ++m) {
                        x_[0][bind][m] = PB[1] * x_[0][bind-bym][m] - PC[1] * x_[0][bind-bym][m+1];
                        y_[0][bind][m] = PB[1] * y_[0][bind-bym][m] - PC[1] * y_[0][bind-bym][m+1] + q_[0][bind-bym][m+1];
                        z_[0][bind][m] = PB[1] * z_[0][bind-bym][m] - PC[1] * z_[0][bind-bym][m+1];
                    }
                    for(m=0;m<=mmax-b-2;m++) {
                        xx_[0][bind][m] = PB[1]*xx_[0][bind-bym][m] - PC[1]*xx_[0][bind-bym][m+1];
                        yy_[0][bind][m] = PB[1]*yy_[0][bind-bym][m] - PC[1]*yy_[0][bind-bym][m+1] + 2*y_[0][bind-bym][m+1];
                        zz_[0][bind][m] = PB[1]*zz_[0][bind-bym][m] - PC[1]*zz_[0][bind-bym][m+1];
                        xy_[0][bind][m] = PB[1]*xy_[0][bind-bym][m] - PC[1]*xy_[0][bind-bym][m+1] + x_[0][bind-bym][m+1];
                        xz_[0][bind][m] = PB[1]*xz_[0][bind-bym][m] - PC[1]*xz_[0][bind-bym][m+1];
                        yz_[0][bind][m] = PB[1]*yz_[0][bind-bym][m] - PC[1]*yz_[0][bind-bym][m+1] + z_[0][bind-bym][m+1];
                    }
                    for (m=0; m <=mmax-b-3; ++m) {
                        xxx_[0][bind][m] = PB[1]*xxx_[0][bind-bym][m] - PC[1] * xxx_[0][bind-bym][m+1];
                        yyy_[0][bind][m] = PB[1]*yyy_[0][bind-bym][m] - PC[1] * yyy_[0][bind-bym][m+1] + 3 * yy_[0][bind-bym][m+1];
                        zzz_[0][bind][m] = PB[1]*zzz_[0][bind-bym][m] - PC[1] * zzz_[0][bind-bym][m+1];
                        xxy_[0][bind][m] = PB[1]*xxy_[0][bind-bym][m] - PC[1] * xxy_[0][bind-bym][m+1] + 1 * xx_[0][bind-bym][m+1];
                        xxz_[0][bind][m] = PB[1]*xxz_[0][bind-bym][m] - PC[1] * xxz_[0][bind-bym][m+1];
                        xyy_[0][bind][m] = PB[1]*xyy_[0][bind-bym][m] - PC[1] * xyy_[0][bind-bym][m+1] + 2 * xy_[0][bind-bym][m+1];
                        yyz_[0][bind][m] = PB[1]*yyz_[0][bind-bym][m] - PC[1] * yyz_[0][bind-bym][m+1] + 2 * yz_[0][bind-bym][m+1];
                        xzz_[0][bind][m] = PB[1]*xzz_[0][bind-bym][m] - PC[1] * xzz_[0][bind-bym][m+1];
                        yzz_[0][bind][m] = PB[1]*yzz_[0][bind-bym][m] - PC[1] * yzz_[0][bind-bym][m+1] + 1 * zz_[0][bind-bym][m+1];
                        xyz_[0][bind][m] = PB[1]*xyz_[0][bind-bym][m] - PC[1] * xyz_[0][bind-bym][m+1] + 1 * xz_[0][bind-bym][m+1];
                    }

                    if (by > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            q_[0][bind][m] += ooz * (by-1) * (q_[0][bind-2*bym][m] - q_[0][bind-2*bym][m+1]);
                        }
                        for (m=0; m<=mmax-b-1; ++m) {
                            x_[0][bind][m] += ooz * (by-1) * (x_[0][bind-2*bym][m] - x_[0][bind-2*bym][m+1]);
                            y_[0][bind][m] += ooz * (by-1) * (y_[0][bind-2*bym][m] - y_[0][bind-2*bym][m+1]);
                            z_[0][bind][m] += ooz * (by-1) * (z_[0][bind-2*bym][m] - z_[0][bind-2*bym][m+1]);
                        }
                        for(m=0;m<=mmax-b-2;m++) {
                            xx_[0][bind][m] += ooz*(by-1)*(xx_[0][bind-2*bym][m] - xx_[0][bind-2*bym][m+1]);
                            yy_[0][bind][m] += ooz*(by-1)*(yy_[0][bind-2*bym][m] - yy_[0][bind-2*bym][m+1]);
                            zz_[0][bind][m] += ooz*(by-1)*(zz_[0][bind-2*bym][m] - zz_[0][bind-2*bym][m+1]);
                            xy_[0][bind][m] += ooz*(by-1)*(xy_[0][bind-2*bym][m] - xy_[0][bind-2*bym][m+1]);
                            xz_[0][bind][m] += ooz*(by-1)*(xz_[0][bind-2*bym][m] - xz_[0][bind-2*bym][m+1]);
                            yz_[0][bind][m] += ooz*(by-1)*(yz_[0][bind-2*bym][m] - yz_[0][bind-2*bym][m+1]);
                        }
                        for (m=0; m<=mmax-b-3; m++) {
                            xxx_[0][bind][m] += ooz*(by-1)*(xxx_[0][bind-2*bym][m] - xxx_[0][bind-2*bym][m+1]);
                            yyy_[0][bind][m] += ooz*(by-1)*(yyy_[0][bind-2*bym][m] - yyy_[0][bind-2*bym][m+1]);
                            zzz_[0][bind][m] += ooz*(by-1)*(zzz_[0][bind-2*bym][m] - zzz_[0][bind-2*bym][m+1]);
                            xxy_[0][bind][m] += ooz*(by-1)*(xxy_[0][bind-2*bym][m] - xxy_[0][bind-2*bym][m+1]);
                            xxz_[0][bind][m] += ooz*(by-1)*(xxz_[0][bind-2*bym][m] - xxz_[0][bind-2*bym][m+1]);
                            xyy_[0][bind][m] += ooz*(by-1)*(xyy_[0][bind-2*bym][m] - xyy_[0][bind-2*bym][m+1]);
                            yyz_[0][bind][m] += ooz*(by-1)*(yyz_[0][bind-2*bym][m] - yyz_[0][bind-2*bym][m+1]);
                            xzz_[0][bind][m] += ooz*(by-1)*(xzz_[0][bind-2*bym][m] - xzz_[0][bind-2*bym][m+1]);
                            yzz_[0][bind][m] += ooz*(by-1)*(yzz_[0][bind-2*bym][m] - yzz_[0][bind-2*bym][m+1]);
                            xyz_[0][bind][m] += ooz*(by-1)*(xyz_[0][bind-2*bym][m] - xyz_[0][bind-2*bym][m+1]);
                        }
                    }
                }
                else if (bx > 0) {
                    for (m=0; m<=mmax-b; ++m) {
                        q_[0][bind][m] = PB[0] * q_[0][bind-bxm][m] - PC[0] * q_[0][bind-bxm][m+1];
                    }
                    for (m=0; m<=mmax-b-1; ++m) {
                        x_[0][bind][m] = PB[0] * x_[0][bind-bxm][m] - PC[0] * x_[0][bind-bxm][m+1] + q_[0][bind-bxm][m+1];
                        y_[0][bind][m] = PB[0] * y_[0][bind-bxm][m] - PC[0] * y_[0][bind-bxm][m+1];
                        z_[0][bind][m] = PB[0] * z_[0][bind-bxm][m] - PC[0] * z_[0][bind-bxm][m+1];
                    }
                    for(m=0;m<=mmax-b-2;m++) {
                        xx_[0][bind][m] = PB[0]*xx_[0][bind-bxm][m] - PC[0]*xx_[0][bind-bxm][m+1] + 2*x_[0][bind-bxm][m+1];
                        yy_[0][bind][m] = PB[0]*yy_[0][bind-bxm][m] - PC[0]*yy_[0][bind-bxm][m+1];
                        zz_[0][bind][m] = PB[0]*zz_[0][bind-bxm][m] - PC[0]*zz_[0][bind-bxm][m+1];
                        xy_[0][bind][m] = PB[0]*xy_[0][bind-bxm][m] - PC[0]*xy_[0][bind-bxm][m+1] + y_[0][bind-bxm][m+1];
                        xz_[0][bind][m] = PB[0]*xz_[0][bind-bxm][m] - PC[0]*xz_[0][bind-bxm][m+1] + z_[0][bind-bxm][m+1];
                        yz_[0][bind][m] = PB[0]*yz_[0][bind-bxm][m] - PC[0]*yz_[0][bind-bxm][m+1];
                    }
                    for (m=0; m <=mmax-b-3; ++m) {
                        xxx_[0][bind][m] = PB[0]*xxx_[0][bind-bxm][m] - PC[0] * xxx_[0][bind-bxm][m+1] + 3 * xx_[0][bind-bxm][m+1];
                        yyy_[0][bind][m] = PB[0]*yyy_[0][bind-bxm][m] - PC[0] * yyy_[0][bind-bxm][m+1];
                        zzz_[0][bind][m] = PB[0]*zzz_[0][bind-bxm][m] - PC[0] * zzz_[0][bind-bxm][m+1];
                        xxy_[0][bind][m] = PB[0]*xxy_[0][bind-bxm][m] - PC[0] * xxy_[0][bind-bxm][m+1] + 2 * xy_[0][bind-bxm][m+1];
                        xxz_[0][bind][m] = PB[0]*xxz_[0][bind-bxm][m] - PC[0] * xxz_[0][bind-bxm][m+1] + 2 * xz_[0][bind-bxm][m+1];
                        xyy_[0][bind][m] = PB[0]*xyy_[0][bind-bxm][m] - PC[0] * xyy_[0][bind-bxm][m+1] + 1 * yy_[0][bind-bxm][m+1];
                        yyz_[0][bind][m] = PB[0]*yyz_[0][bind-bxm][m] - PC[0] * yyz_[0][bind-bxm][m+1];
                        xzz_[0][bind][m] = PB[0]*xzz_[0][bind-bxm][m] - PC[0] * xzz_[0][bind-bxm][m+1] + 1 * zz_[0][bind-bxm][m+1];
                        yzz_[0][bind][m] = PB[0]*yzz_[0][bind-bxm][m] - PC[0] * yzz_[0][bind-bxm][m+1];
                        xyz_[0][bind][m] = PB[0]*xyz_[0][bind-bxm][m] - PC[0] * xyz_[0][bind-bxm][m+1] + 1 * yz_[0][bind-bxm][m+1];
                    }

                    if (bx > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            q_[0][bind][m] += ooz * (bx-1) * (q_[0][bind-2*bxm][m] - q_[0][bind-2*bxm][m+1]);
                        }
                        for (m=0; m<=mmax-b-1; ++m) {
                            x_[0][bind][m] += ooz * (bx-1) * (x_[0][bind-2*bxm][m] - x_[0][bind-2*bxm][m+1]);
                            y_[0][bind][m] += ooz * (bx-1) * (y_[0][bind-2*bxm][m] - y_[0][bind-2*bxm][m+1]);
                            z_[0][bind][m] += ooz * (bx-1) * (z_[0][bind-2*bxm][m] - z_[0][bind-2*bxm][m+1]);
                        }
                        for(m=0;m<=mmax-b-2;m++) {
                            xx_[0][bind][m] += ooz*(bx-1)*(xx_[0][bind-2*bxm][m] - xx_[0][bind-2*bxm][m+1]);
                            yy_[0][bind][m] += ooz*(bx-1)*(yy_[0][bind-2*bxm][m] - yy_[0][bind-2*bxm][m+1]);
                            zz_[0][bind][m] += ooz*(bx-1)*(zz_[0][bind-2*bxm][m] - zz_[0][bind-2*bxm][m+1]);
                            xy_[0][bind][m] += ooz*(bx-1)*(xy_[0][bind-2*bxm][m] - xy_[0][bind-2*bxm][m+1]);
                            xz_[0][bind][m] += ooz*(bx-1)*(xz_[0][bind-2*bxm][m] - xz_[0][bind-2*bxm][m+1]);
                            yz_[0][bind][m] += ooz*(bx-1)*(yz_[0][bind-2*bxm][m] - yz_[0][bind-2*bxm][m+1]);
                        }
                        for (m=0; m<=mmax-b-3; m++) {
                            xxx_[0][bind][m] += ooz*(bx-1)*(xxx_[0][bind-2*bxm][m] - xxx_[0][bind-2*bxm][m+1]);
                            yyy_[0][bind][m] += ooz*(bx-1)*(yyy_[0][bind-2*bxm][m] - yyy_[0][bind-2*bxm][m+1]);
                            zzz_[0][bind][m] += ooz*(bx-1)*(zzz_[0][bind-2*bxm][m] - zzz_[0][bind-2*bxm][m+1]);
                            xxy_[0][bind][m] += ooz*(bx-1)*(xxy_[0][bind-2*bxm][m] - xxy_[0][bind-2*bxm][m+1]);
                            xxz_[0][bind][m] += ooz*(bx-1)*(xxz_[0][bind-2*bxm][m] - xxz_[0][bind-2*bxm][m+1]);
                            xyy_[0][bind][m] += ooz*(bx-1)*(xyy_[0][bind-2*bxm][m] - xyy_[0][bind-2*bxm][m+1]);
                            yyz_[0][bind][m] += ooz*(bx-1)*(yyz_[0][bind-2*bxm][m] - yyz_[0][bind-2*bxm][m+1]);
                            xzz_[0][bind][m] += ooz*(bx-1)*(xzz_[0][bind-2*bxm][m] - xzz_[0][bind-2*bxm][m+1]);
                            yzz_[0][bind][m] += ooz*(bx-1)*(yzz_[0][bind-2*bxm][m] - yzz_[0][bind-2*bxm][m+1]);
                            xyz_[0][bind][m] += ooz*(bx-1)*(xyz_[0][bind-2*bxm][m] - xyz_[0][bind-2*bxm][m+1]);
                        }
                    }
                }
            }
        }
    }

    // Perform upward recursion in a with all b's
    for (b=0; b<=am2; b++) {
        for (bx=0; bx<=b; bx++) {
            for (by=0; by<=b-bx;by++) {
                bz = b-bx-by;
                bind = bx*bxm + by*bym + bz*bzm;

                for (a=1; a<=am1; a++) {
                    // This next for loop was for (ax=0; ax<=b; ax++)
                    // this could explain why dx2 was not being computed.
                    // change for for(ax=0; ax<a; ax++) on 2005-09-15 4:11pm
                    for (ax=0; ax<=a; ax++) {
                        for (ay=0; ay<=a-ax; ay++) {
                            az = a-ax-ay;
                            aind = ax*axm + ay*aym + az*azm;

                            if (az > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    q_[aind][bind][m] = PA[2] * q_[aind-azm][bind][m] - PC[2] * q_[aind-azm][bind][m+1];
                                }
                                for (m=0; m<=mmax-a-b-1; ++m) {
                                    x_[aind][bind][m] = PA[2] * x_[aind-azm][bind][m] - PC[2] * x_[aind-azm][bind][m+1];
                                    y_[aind][bind][m] = PA[2] * y_[aind-azm][bind][m] - PC[2] * y_[aind-azm][bind][m+1];
                                    z_[aind][bind][m] = PA[2] * z_[aind-azm][bind][m] - PC[2] * z_[aind-azm][bind][m+1] + q_[aind-azm][bind][m+1];
                                }
                                for(m=0;m<=mmax-a-b-2;m++) {
                                    xx_[aind][bind][m] = PA[2]*xx_[aind-azm][bind][m] - PC[2]*xx_[aind-azm][bind][m+1];
                                    yy_[aind][bind][m] = PA[2]*yy_[aind-azm][bind][m] - PC[2]*yy_[aind-azm][bind][m+1];
                                    zz_[aind][bind][m] = PA[2]*zz_[aind-azm][bind][m] - PC[2]*zz_[aind-azm][bind][m+1] + 2*z_[aind-azm][bind][m+1];
                                    xy_[aind][bind][m] = PA[2]*xy_[aind-azm][bind][m] - PC[2]*xy_[aind-azm][bind][m+1];
                                    xz_[aind][bind][m] = PA[2]*xz_[aind-azm][bind][m] - PC[2]*xz_[aind-azm][bind][m+1] + x_[aind-azm][bind][m+1];
                                    yz_[aind][bind][m] = PA[2]*yz_[aind-azm][bind][m] - PC[2]*yz_[aind-azm][bind][m+1] + y_[aind-azm][bind][m+1];
                                }
                                for(m=0;m<=mmax-a-b-3;m++) {
                                    xxx_[aind][bind][m] = PA[2]*xxx_[aind-azm][bind][m] - PC[2]*xxx_[aind-azm][bind][m+1];
                                    yyy_[aind][bind][m] = PA[2]*yyy_[aind-azm][bind][m] - PC[2]*yyy_[aind-azm][bind][m+1];
                                    zzz_[aind][bind][m] = PA[2]*zzz_[aind-azm][bind][m] - PC[2]*zzz_[aind-azm][bind][m+1] + 3*zz_[aind-azm][bind][m+1];
                                    xxy_[aind][bind][m] = PA[2]*xxy_[aind-azm][bind][m] - PC[2]*xxy_[aind-azm][bind][m+1];
                                    xxz_[aind][bind][m] = PA[2]*xxz_[aind-azm][bind][m] - PC[2]*xxz_[aind-azm][bind][m+1] + 1*xx_[aind-azm][bind][m+1];
                                    xyy_[aind][bind][m] = PA[2]*xyy_[aind-azm][bind][m] - PC[2]*xyy_[aind-azm][bind][m+1];
                                    yyz_[aind][bind][m] = PA[2]*yyz_[aind-azm][bind][m] - PC[2]*yyz_[aind-azm][bind][m+1] + 1*yy_[aind-azm][bind][m+1];
                                    xzz_[aind][bind][m] = PA[2]*xzz_[aind-azm][bind][m] - PC[2]*xzz_[aind-azm][bind][m+1] + 2*xz_[aind-azm][bind][m+1];
                                    yzz_[aind][bind][m] = PA[2]*yzz_[aind-azm][bind][m] - PC[2]*yzz_[aind-azm][bind][m+1] + 2*yz_[aind-azm][bind][m+1];
                                    xyz_[aind][bind][m] = PA[2]*xyz_[aind-azm][bind][m] - PC[2]*xyz_[aind-azm][bind][m+1] + 1*xy_[aind-azm][bind][m+1];
                                }

                                if (az > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        q_[aind][bind][m] += ooz * (az-1) * (q_[aind-2*azm][bind][m] - q_[aind-2*azm][bind][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        x_[aind][bind][m] += ooz * (az-1) * (x_[aind-2*azm][bind][m] - x_[aind-2*azm][bind][m+1]);
                                        y_[aind][bind][m] += ooz * (az-1) * (y_[aind-2*azm][bind][m] - y_[aind-2*azm][bind][m+1]);
                                        z_[aind][bind][m] += ooz * (az-1) * (z_[aind-2*azm][bind][m] - z_[aind-2*azm][bind][m+1]);
                                    }
                                    for(m=0;m<=mmax-a-b-2;m++) {
                                        xx_[aind][bind][m] += ooz*(az-1)*(xx_[aind-2*azm][bind][m] - xx_[aind-2*azm][bind][m+1]);
                                        yy_[aind][bind][m] += ooz*(az-1)*(yy_[aind-2*azm][bind][m] - yy_[aind-2*azm][bind][m+1]);
                                        zz_[aind][bind][m] += ooz*(az-1)*(zz_[aind-2*azm][bind][m] - zz_[aind-2*azm][bind][m+1]);
                                        xy_[aind][bind][m] += ooz*(az-1)*(xy_[aind-2*azm][bind][m] - xy_[aind-2*azm][bind][m+1]);
                                        xz_[aind][bind][m] += ooz*(az-1)*(xz_[aind-2*azm][bind][m] - xz_[aind-2*azm][bind][m+1]);
                                        yz_[aind][bind][m] += ooz*(az-1)*(yz_[aind-2*azm][bind][m] - yz_[aind-2*azm][bind][m+1]);
                                    }
                                    for(m=0;m<=mmax-a-b-3;m++) {
                                        xxx_[aind][bind][m] += ooz*(az-1)*(xxx_[aind-2*azm][bind][m] - xxx_[aind-2*azm][bind][m+1]);
                                        yyy_[aind][bind][m] += ooz*(az-1)*(yyy_[aind-2*azm][bind][m] - yyy_[aind-2*azm][bind][m+1]);
                                        zzz_[aind][bind][m] += ooz*(az-1)*(zzz_[aind-2*azm][bind][m] - zzz_[aind-2*azm][bind][m+1]);
                                        xxy_[aind][bind][m] += ooz*(az-1)*(xxy_[aind-2*azm][bind][m] - xxy_[aind-2*azm][bind][m+1]);
                                        xxz_[aind][bind][m] += ooz*(az-1)*(xxz_[aind-2*azm][bind][m] - xxz_[aind-2*azm][bind][m+1]);
                                        xyy_[aind][bind][m] += ooz*(az-1)*(xyy_[aind-2*azm][bind][m] - xyy_[aind-2*azm][bind][m+1]);
                                        yyz_[aind][bind][m] += ooz*(az-1)*(yyz_[aind-2*azm][bind][m] - yyz_[aind-2*azm][bind][m+1]);
                                        xzz_[aind][bind][m] += ooz*(az-1)*(xzz_[aind-2*azm][bind][m] - xzz_[aind-2*azm][bind][m+1]);
                                        yzz_[aind][bind][m] += ooz*(az-1)*(yzz_[aind-2*azm][bind][m] - yzz_[aind-2*azm][bind][m+1]);
                                        xyz_[aind][bind][m] += ooz*(az-1)*(xyz_[aind-2*azm][bind][m] - xyz_[aind-2*azm][bind][m+1]);
                                    }
                                }
                                if (bz > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        q_[aind][bind][m] += ooz * bz * (q_[aind-azm][bind-bzm][m] - q_[aind-azm][bind-bzm][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        x_[aind][bind][m] += ooz * bz * (x_[aind-azm][bind-bzm][m] - x_[aind-azm][bind-bzm][m+1]);
                                        y_[aind][bind][m] += ooz * bz * (y_[aind-azm][bind-bzm][m] - y_[aind-azm][bind-bzm][m+1]);
                                        z_[aind][bind][m] += ooz * bz * (z_[aind-azm][bind-bzm][m] - z_[aind-azm][bind-bzm][m+1]);
                                    }
                                    for(m=0;m<=mmax-a-b-2;m++) {
                                        xx_[aind][bind][m] += ooz*bz*(xx_[aind-azm][bind-bzm][m] - xx_[aind-azm][bind-bzm][m+1]);
                                        yy_[aind][bind][m] += ooz*bz*(yy_[aind-azm][bind-bzm][m] - yy_[aind-azm][bind-bzm][m+1]);
                                        zz_[aind][bind][m] += ooz*bz*(zz_[aind-azm][bind-bzm][m] - zz_[aind-azm][bind-bzm][m+1]);
                                        xy_[aind][bind][m] += ooz*bz*(xy_[aind-azm][bind-bzm][m] - xy_[aind-azm][bind-bzm][m+1]);
                                        xz_[aind][bind][m] += ooz*bz*(xz_[aind-azm][bind-bzm][m] - xz_[aind-azm][bind-bzm][m+1]);
                                        yz_[aind][bind][m] += ooz*bz*(yz_[aind-azm][bind-bzm][m] - yz_[aind-azm][bind-bzm][m+1]);
                                    }
                                    for(m=0;m<=mmax-a-b-3;m++) {
                                        xxx_[aind][bind][m] += ooz*bz*(xxx_[aind-azm][bind-bzm][m] - xxx_[aind-azm][bind-bzm][m+1]);
                                        yyy_[aind][bind][m] += ooz*bz*(yyy_[aind-azm][bind-bzm][m] - yyy_[aind-azm][bind-bzm][m+1]);
                                        zzz_[aind][bind][m] += ooz*bz*(zzz_[aind-azm][bind-bzm][m] - zzz_[aind-azm][bind-bzm][m+1]);
                                        xxy_[aind][bind][m] += ooz*bz*(xxy_[aind-azm][bind-bzm][m] - xxy_[aind-azm][bind-bzm][m+1]);
                                        xxz_[aind][bind][m] += ooz*bz*(xxz_[aind-azm][bind-bzm][m] - xxz_[aind-azm][bind-bzm][m+1]);
                                        xyy_[aind][bind][m] += ooz*bz*(xyy_[aind-azm][bind-bzm][m] - xyy_[aind-azm][bind-bzm][m+1]);
                                        yyz_[aind][bind][m] += ooz*bz*(yyz_[aind-azm][bind-bzm][m] - yyz_[aind-azm][bind-bzm][m+1]);
                                        xzz_[aind][bind][m] += ooz*bz*(xzz_[aind-azm][bind-bzm][m] - xzz_[aind-azm][bind-bzm][m+1]);
                                        yzz_[aind][bind][m] += ooz*bz*(yzz_[aind-azm][bind-bzm][m] - yzz_[aind-azm][bind-bzm][m+1]);
                                        xyz_[aind][bind][m] += ooz*bz*(xyz_[aind-azm][bind-bzm][m] - xyz_[aind-azm][bind-bzm][m+1]);
                                    }
                                }
                            }
                            else if (ay > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    q_[aind][bind][m] = PA[1] * q_[aind-aym][bind][m] - PC[1] * q_[aind-aym][bind][m+1];
                                }
                                for (m=0; m<=mmax-a-b-1; ++m) {
                                    x_[aind][bind][m] = PA[1] * x_[aind-aym][bind][m] - PC[1] * x_[aind-aym][bind][m+1];
                                    y_[aind][bind][m] = PA[1] * y_[aind-aym][bind][m] - PC[1] * y_[aind-aym][bind][m+1] + q_[aind-aym][bind][m+1];
                                    z_[aind][bind][m] = PA[1] * z_[aind-aym][bind][m] - PC[1] * z_[aind-aym][bind][m+1];
                                }
                                for(m=0;m<=mmax-a-b-2;m++) {
                                    xx_[aind][bind][m] = PA[1]*xx_[aind-aym][bind][m] - PC[1]*xx_[aind-aym][bind][m+1];
                                    yy_[aind][bind][m] = PA[1]*yy_[aind-aym][bind][m] - PC[1]*yy_[aind-aym][bind][m+1] + 2*y_[aind-aym][bind][m+1];
                                    zz_[aind][bind][m] = PA[1]*zz_[aind-aym][bind][m] - PC[1]*zz_[aind-aym][bind][m+1];
                                    xy_[aind][bind][m] = PA[1]*xy_[aind-aym][bind][m] - PC[1]*xy_[aind-aym][bind][m+1] + x_[aind-aym][bind][m+1];
                                    xz_[aind][bind][m] = PA[1]*xz_[aind-aym][bind][m] - PC[1]*xz_[aind-aym][bind][m+1];
                                    yz_[aind][bind][m] = PA[1]*yz_[aind-aym][bind][m] - PC[1]*yz_[aind-aym][bind][m+1] + z_[aind-aym][bind][m+1];
                                }
                                for(m=0;m<=mmax-a-b-3;m++) {
                                    xxx_[aind][bind][m] = PA[1]*xxx_[aind-aym][bind][m] - PC[1]*xxx_[aind-aym][bind][m+1];
                                    yyy_[aind][bind][m] = PA[1]*yyy_[aind-aym][bind][m] - PC[1]*yyy_[aind-aym][bind][m+1] + 3*yy_[aind-aym][bind][m+1];
                                    zzz_[aind][bind][m] = PA[1]*zzz_[aind-aym][bind][m] - PC[1]*zzz_[aind-aym][bind][m+1];
                                    xxy_[aind][bind][m] = PA[1]*xxy_[aind-aym][bind][m] - PC[1]*xxy_[aind-aym][bind][m+1] + 1*xx_[aind-aym][bind][m+1];
                                    xxz_[aind][bind][m] = PA[1]*xxz_[aind-aym][bind][m] - PC[1]*xxz_[aind-aym][bind][m+1];
                                    xyy_[aind][bind][m] = PA[1]*xyy_[aind-aym][bind][m] - PC[1]*xyy_[aind-aym][bind][m+1] + 2*xy_[aind-aym][bind][m+1];
                                    yyz_[aind][bind][m] = PA[1]*yyz_[aind-aym][bind][m] - PC[1]*yyz_[aind-aym][bind][m+1] + 2*yz_[aind-aym][bind][m+1];
                                    xzz_[aind][bind][m] = PA[1]*xzz_[aind-aym][bind][m] - PC[1]*xzz_[aind-aym][bind][m+1];
                                    yzz_[aind][bind][m] = PA[1]*yzz_[aind-aym][bind][m] - PC[1]*yzz_[aind-aym][bind][m+1] + 1*zz_[aind-aym][bind][m+1];
                                    xyz_[aind][bind][m] = PA[1]*xyz_[aind-aym][bind][m] - PC[1]*xyz_[aind-aym][bind][m+1] + 1*xz_[aind-aym][bind][m+1];
                                }
                                if (ay > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        q_[aind][bind][m] += ooz * (ay-1) * (q_[aind-2*aym][bind][m] - q_[aind-2*aym][bind][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        x_[aind][bind][m] += ooz * (ay-1) * (x_[aind-2*aym][bind][m] - x_[aind-2*aym][bind][m+1]);
                                        y_[aind][bind][m] += ooz * (ay-1) * (y_[aind-2*aym][bind][m] - y_[aind-2*aym][bind][m+1]);
                                        z_[aind][bind][m] += ooz * (ay-1) * (z_[aind-2*aym][bind][m] - z_[aind-2*aym][bind][m+1]);
                                    }
                                    for(m=0;m<=mmax-a-b-2;m++) {
                                        xx_[aind][bind][m] += ooz*(ay-1)*(xx_[aind-2*aym][bind][m] - xx_[aind-2*aym][bind][m+1]);
                                        yy_[aind][bind][m] += ooz*(ay-1)*(yy_[aind-2*aym][bind][m] - yy_[aind-2*aym][bind][m+1]);
                                        zz_[aind][bind][m] += ooz*(ay-1)*(zz_[aind-2*aym][bind][m] - zz_[aind-2*aym][bind][m+1]);
                                        xy_[aind][bind][m] += ooz*(ay-1)*(xy_[aind-2*aym][bind][m] - xy_[aind-2*aym][bind][m+1]);
                                        xz_[aind][bind][m] += ooz*(ay-1)*(xz_[aind-2*aym][bind][m] - xz_[aind-2*aym][bind][m+1]);
                                        yz_[aind][bind][m] += ooz*(ay-1)*(yz_[aind-2*aym][bind][m] - yz_[aind-2*aym][bind][m+1]);
                                    }
                                    for(m=0;m<=mmax-a-b-3;m++) {
                                        xxx_[aind][bind][m] += ooz*(ay-1)*(xxx_[aind-2*aym][bind][m] - xxx_[aind-2*aym][bind][m+1]);
                                        yyy_[aind][bind][m] += ooz*(ay-1)*(yyy_[aind-2*aym][bind][m] - yyy_[aind-2*aym][bind][m+1]);
                                        zzz_[aind][bind][m] += ooz*(ay-1)*(zzz_[aind-2*aym][bind][m] - zzz_[aind-2*aym][bind][m+1]);
                                        xxy_[aind][bind][m] += ooz*(ay-1)*(xxy_[aind-2*aym][bind][m] - xxy_[aind-2*aym][bind][m+1]);
                                        xxz_[aind][bind][m] += ooz*(ay-1)*(xxz_[aind-2*aym][bind][m] - xxz_[aind-2*aym][bind][m+1]);
                                        xyy_[aind][bind][m] += ooz*(ay-1)*(xyy_[aind-2*aym][bind][m] - xyy_[aind-2*aym][bind][m+1]);
                                        yyz_[aind][bind][m] += ooz*(ay-1)*(yyz_[aind-2*aym][bind][m] - yyz_[aind-2*aym][bind][m+1]);
                                        xzz_[aind][bind][m] += ooz*(ay-1)*(xzz_[aind-2*aym][bind][m] - xzz_[aind-2*aym][bind][m+1]);
                                        yzz_[aind][bind][m] += ooz*(ay-1)*(yzz_[aind-2*aym][bind][m] - yzz_[aind-2*aym][bind][m+1]);
                                        xyz_[aind][bind][m] += ooz*(ay-1)*(xyz_[aind-2*aym][bind][m] - xyz_[aind-2*aym][bind][m+1]);
                                    }
                                }
                                if (by > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        q_[aind][bind][m] += ooz * by * (q_[aind-aym][bind-bym][m] - q_[aind-aym][bind-bym][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        x_[aind][bind][m] += ooz * by * (x_[aind-aym][bind-bym][m] - x_[aind-aym][bind-bym][m+1]);
                                        y_[aind][bind][m] += ooz * by * (y_[aind-aym][bind-bym][m] - y_[aind-aym][bind-bym][m+1]);
                                        z_[aind][bind][m] += ooz * by * (z_[aind-aym][bind-bym][m] - z_[aind-aym][bind-bym][m+1]);
                                    }
                                    for(m=0;m<=mmax-a-b-2;m++) {
                                        xx_[aind][bind][m] += ooz*by*(xx_[aind-aym][bind-bym][m] - xx_[aind-aym][bind-bym][m+1]);
                                        yy_[aind][bind][m] += ooz*by*(yy_[aind-aym][bind-bym][m] - yy_[aind-aym][bind-bym][m+1]);
                                        zz_[aind][bind][m] += ooz*by*(zz_[aind-aym][bind-bym][m] - zz_[aind-aym][bind-bym][m+1]);
                                        xy_[aind][bind][m] += ooz*by*(xy_[aind-aym][bind-bym][m] - xy_[aind-aym][bind-bym][m+1]);
                                        xz_[aind][bind][m] += ooz*by*(xz_[aind-aym][bind-bym][m] - xz_[aind-aym][bind-bym][m+1]);
                                        yz_[aind][bind][m] += ooz*by*(yz_[aind-aym][bind-bym][m] - yz_[aind-aym][bind-bym][m+1]);
                                    }
                                    for(m=0;m<=mmax-a-b-3;m++) {
                                        xxx_[aind][bind][m] += ooz*by*(xxx_[aind-aym][bind-bym][m] - xxx_[aind-aym][bind-bym][m+1]);
                                        yyy_[aind][bind][m] += ooz*by*(yyy_[aind-aym][bind-bym][m] - yyy_[aind-aym][bind-bym][m+1]);
                                        zzz_[aind][bind][m] += ooz*by*(zzz_[aind-aym][bind-bym][m] - zzz_[aind-aym][bind-bym][m+1]);
                                        xxy_[aind][bind][m] += ooz*by*(xxy_[aind-aym][bind-bym][m] - xxy_[aind-aym][bind-bym][m+1]);
                                        xxz_[aind][bind][m] += ooz*by*(xxz_[aind-aym][bind-bym][m] - xxz_[aind-aym][bind-bym][m+1]);
                                        xyy_[aind][bind][m] += ooz*by*(xyy_[aind-aym][bind-bym][m] - xyy_[aind-aym][bind-bym][m+1]);
                                        yyz_[aind][bind][m] += ooz*by*(yyz_[aind-aym][bind-bym][m] - yyz_[aind-aym][bind-bym][m+1]);
                                        xzz_[aind][bind][m] += ooz*by*(xzz_[aind-aym][bind-bym][m] - xzz_[aind-aym][bind-bym][m+1]);
                                        yzz_[aind][bind][m] += ooz*by*(yzz_[aind-aym][bind-bym][m] - yzz_[aind-aym][bind-bym][m+1]);
                                        xyz_[aind][bind][m] += ooz*by*(xyz_[aind-aym][bind-bym][m] - xyz_[aind-aym][bind-bym][m+1]);
                                    }
                                }
                            }
                            else if (ax > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    q_[aind][bind][m] = PA[0] * q_[aind-axm][bind][m] - PC[0] * q_[aind-axm][bind][m+1];
                                }
                                for (m=0; m<=mmax-a-b-1; ++m) {
                                    x_[aind][bind][m] = PA[0] * x_[aind-axm][bind][m] - PC[0] * x_[aind-axm][bind][m+1] + q_[aind-axm][bind][m+1];
                                    y_[aind][bind][m] = PA[0] * y_[aind-axm][bind][m] - PC[0] * y_[aind-axm][bind][m+1];
                                    z_[aind][bind][m] = PA[0] * z_[aind-axm][bind][m] - PC[0] * z_[aind-axm][bind][m+1];
                                }
                                for(m=0;m<=mmax-a-b-2;m++) {  /* Gradients of the electric field */
                                    xx_[aind][bind][m] = PA[0]*xx_[aind-axm][bind][m] - PC[0]*xx_[aind-axm][bind][m+1] + 2*x_[aind-axm][bind][m+1];
                                    yy_[aind][bind][m] = PA[0]*yy_[aind-axm][bind][m] - PC[0]*yy_[aind-axm][bind][m+1];
                                    zz_[aind][bind][m] = PA[0]*zz_[aind-axm][bind][m] - PC[0]*zz_[aind-axm][bind][m+1];
                                    xy_[aind][bind][m] = PA[0]*xy_[aind-axm][bind][m] - PC[0]*xy_[aind-axm][bind][m+1] + y_[aind-axm][bind][m+1];
                                    xz_[aind][bind][m] = PA[0]*xz_[aind-axm][bind][m] - PC[0]*xz_[aind-axm][bind][m+1] + z_[aind-axm][bind][m+1];
                                    yz_[aind][bind][m] = PA[0]*yz_[aind-axm][bind][m] - PC[0]*yz_[aind-axm][bind][m+1];
                                }
                                for(m=0;m<=mmax-a-b-3;m++) {
                                    xxx_[aind][bind][m] = PA[0]*xxx_[aind-axm][bind][m] - PC[0]*xxx_[aind-axm][bind][m+1] + 3*xx_[aind-axm][bind][m+1];
                                    yyy_[aind][bind][m] = PA[0]*yyy_[aind-axm][bind][m] - PC[0]*yyy_[aind-axm][bind][m+1];
                                    zzz_[aind][bind][m] = PA[0]*zzz_[aind-axm][bind][m] - PC[0]*zzz_[aind-axm][bind][m+1];
                                    xxy_[aind][bind][m] = PA[0]*xxy_[aind-axm][bind][m] - PC[0]*xxy_[aind-axm][bind][m+1] + 2*xy_[aind-axm][bind][m+1];
                                    xxz_[aind][bind][m] = PA[0]*xxz_[aind-axm][bind][m] - PC[0]*xxz_[aind-axm][bind][m+1] + 2*xz_[aind-axm][bind][m+1];
                                    xyy_[aind][bind][m] = PA[0]*xyy_[aind-axm][bind][m] - PC[0]*xyy_[aind-axm][bind][m+1] + 1*yy_[aind-axm][bind][m+1];
                                    yyz_[aind][bind][m] = PA[0]*yyz_[aind-axm][bind][m] - PC[0]*yyz_[aind-axm][bind][m+1];
                                    xzz_[aind][bind][m] = PA[0]*xzz_[aind-axm][bind][m] - PC[0]*xzz_[aind-axm][bind][m+1] + 1*zz_[aind-axm][bind][m+1];
                                    yzz_[aind][bind][m] = PA[0]*yzz_[aind-axm][bind][m] - PC[0]*yzz_[aind-axm][bind][m+1];
                                    xyz_[aind][bind][m] = PA[0]*xyz_[aind-axm][bind][m] - PC[0]*xyz_[aind-axm][bind][m+1] + 1*yz_[aind-axm][bind][m+1];
                                }

                                if (ax > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        q_[aind][bind][m] += ooz * (ax-1) * (q_[aind-2*axm][bind][m] - q_[aind-2*axm][bind][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        x_[aind][bind][m] += ooz * (ax-1) * (x_[aind-2*axm][bind][m] - x_[aind-2*axm][bind][m+1]);
                                        y_[aind][bind][m] += ooz * (ax-1) * (y_[aind-2*axm][bind][m] - y_[aind-2*axm][bind][m+1]);
                                        z_[aind][bind][m] += ooz * (ax-1) * (z_[aind-2*axm][bind][m] - z_[aind-2*axm][bind][m+1]);
                                    }
                                    for(m=0;m<=mmax-a-b-2;m++) {
                                        xx_[aind][bind][m] += ooz*(ax-1)*(xx_[aind-2*axm][bind][m] - xx_[aind-2*axm][bind][m+1]);
                                        yy_[aind][bind][m] += ooz*(ax-1)*(yy_[aind-2*axm][bind][m] - yy_[aind-2*axm][bind][m+1]);
                                        zz_[aind][bind][m] += ooz*(ax-1)*(zz_[aind-2*axm][bind][m] - zz_[aind-2*axm][bind][m+1]);
                                        xy_[aind][bind][m] += ooz*(ax-1)*(xy_[aind-2*axm][bind][m] - xy_[aind-2*axm][bind][m+1]);
                                        xz_[aind][bind][m] += ooz*(ax-1)*(xz_[aind-2*axm][bind][m] - xz_[aind-2*axm][bind][m+1]);
                                        yz_[aind][bind][m] += ooz*(ax-1)*(yz_[aind-2*axm][bind][m] - yz_[aind-2*axm][bind][m+1]);
                                    }
                                    for(m=0;m<=mmax-a-b-3;m++) {
                                        xxx_[aind][bind][m] += ooz*(ax-1)*(xxx_[aind-2*axm][bind][m] - xxx_[aind-2*axm][bind][m+1]);
                                        yyy_[aind][bind][m] += ooz*(ax-1)*(yyy_[aind-2*axm][bind][m] - yyy_[aind-2*axm][bind][m+1]);
                                        zzz_[aind][bind][m] += ooz*(ax-1)*(zzz_[aind-2*axm][bind][m] - zzz_[aind-2*axm][bind][m+1]);
                                        xxy_[aind][bind][m] += ooz*(ax-1)*(xxy_[aind-2*axm][bind][m] - xxy_[aind-2*axm][bind][m+1]);
                                        xxz_[aind][bind][m] += ooz*(ax-1)*(xxz_[aind-2*axm][bind][m] - xxz_[aind-2*axm][bind][m+1]);
                                        xyy_[aind][bind][m] += ooz*(ax-1)*(xyy_[aind-2*axm][bind][m] - xyy_[aind-2*axm][bind][m+1]);
                                        yyz_[aind][bind][m] += ooz*(ax-1)*(yyz_[aind-2*axm][bind][m] - yyz_[aind-2*axm][bind][m+1]);
                                        xzz_[aind][bind][m] += ooz*(ax-1)*(xzz_[aind-2*axm][bind][m] - xzz_[aind-2*axm][bind][m+1]);
                                        yzz_[aind][bind][m] += ooz*(ax-1)*(yzz_[aind-2*axm][bind][m] - yzz_[aind-2*axm][bind][m+1]);
                                        xyz_[aind][bind][m] += ooz*(ax-1)*(xyz_[aind-2*axm][bind][m] - xyz_[aind-2*axm][bind][m+1]);
                                    }
                                }
                                if (bx > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        q_[aind][bind][m] += ooz * bx * (q_[aind-axm][bind-bxm][m] - q_[aind-axm][bind-bxm][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        x_[aind][bind][m] += ooz * bx * (x_[aind-axm][bind-bxm][m] - x_[aind-axm][bind-bxm][m+1]);
                                        y_[aind][bind][m] += ooz * bx * (y_[aind-axm][bind-bxm][m] - y_[aind-axm][bind-bxm][m+1]);
                                        z_[aind][bind][m] += ooz * bx * (z_[aind-axm][bind-bxm][m] - z_[aind-axm][bind-bxm][m+1]);
                                    }
                                    for(m=0;m<=mmax-a-b-2;m++) {
                                        xx_[aind][bind][m] += ooz*bx*(xx_[aind-axm][bind-bxm][m] - xx_[aind-axm][bind-bxm][m+1]);
                                        yy_[aind][bind][m] += ooz*bx*(yy_[aind-axm][bind-bxm][m] - yy_[aind-axm][bind-bxm][m+1]);
                                        zz_[aind][bind][m] += ooz*bx*(zz_[aind-axm][bind-bxm][m] - zz_[aind-axm][bind-bxm][m+1]);
                                        xy_[aind][bind][m] += ooz*bx*(xy_[aind-axm][bind-bxm][m] - xy_[aind-axm][bind-bxm][m+1]);
                                        xz_[aind][bind][m] += ooz*bx*(xz_[aind-axm][bind-bxm][m] - xz_[aind-axm][bind-bxm][m+1]);
                                        yz_[aind][bind][m] += ooz*bx*(yz_[aind-axm][bind-bxm][m] - yz_[aind-axm][bind-bxm][m+1]);
                                    }
                                    for(m=0;m<=mmax-a-b-3;m++) {
                                        xxx_[aind][bind][m] += ooz*bx*(xxx_[aind-axm][bind-bxm][m] - xxx_[aind-axm][bind-bxm][m+1]);
                                        yyy_[aind][bind][m] += ooz*bx*(yyy_[aind-axm][bind-bxm][m] - yyy_[aind-axm][bind-bxm][m+1]);
                                        zzz_[aind][bind][m] += ooz*bx*(zzz_[aind-axm][bind-bxm][m] - zzz_[aind-axm][bind-bxm][m+1]);
                                        xxy_[aind][bind][m] += ooz*bx*(xxy_[aind-axm][bind-bxm][m] - xxy_[aind-axm][bind-bxm][m+1]);
                                        xxz_[aind][bind][m] += ooz*bx*(xxz_[aind-axm][bind-bxm][m] - xxz_[aind-axm][bind-bxm][m+1]);
                                        xyy_[aind][bind][m] += ooz*bx*(xyy_[aind-axm][bind-bxm][m] - xyy_[aind-axm][bind-bxm][m+1]);
                                        yyz_[aind][bind][m] += ooz*bx*(yyz_[aind-axm][bind-bxm][m] - yyz_[aind-axm][bind-bxm][m+1]);
                                        xzz_[aind][bind][m] += ooz*bx*(xzz_[aind-axm][bind-bxm][m] - xzz_[aind-axm][bind-bxm][m+1]);
                                        yzz_[aind][bind][m] += ooz*bx*(yzz_[aind-axm][bind-bxm][m] - yzz_[aind-axm][bind-bxm][m+1]);
                                        xyz_[aind][bind][m] += ooz*bx*(xyz_[aind-axm][bind-bxm][m] - xyz_[aind-axm][bind-bxm][m+1]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    delete[] F;
}


ObaraSaikaTwoCenterRecursion::ObaraSaikaTwoCenterRecursion(int max_am1, int max_am2):
    max_am1_(max_am1), max_am2_(max_am2)
{
    if (max_am1 < 0)
        throw SanityCheckError("ERROR: ObaraSaikaTwoCenterRecursion -- max_am1 must be nonnegative", __FILE__, __LINE__);
    if (max_am2 < 0)
        throw SanityCheckError("ERROR: ObaraSaikaTwoCenterRecursion -- max_am2 must be nonnegative", __FILE__, __LINE__);

    int am1 = max_am1;
    int am2 = max_am2;
    if (max_am1_ == 0)
        am1 = 1;
    if (max_am2_ == 0)
        am2 = 1;

    x_ = block_matrix(am1+1, am2+1);
    y_ = block_matrix(am1+1, am2+1);
    z_ = block_matrix(am1+1, am2+1);
}

ObaraSaikaTwoCenterRecursion::~ObaraSaikaTwoCenterRecursion()
{
    free_block(x_);
    free_block(y_);
    free_block(z_);
}

void ObaraSaikaTwoCenterRecursion::compute(double PA[3], double PB[3], double gamma, int am1, int am2)
{
    if (am1 < 0 || am1 > max_am1_)
        throw SanityCheckError("ERROR: ObaraSaikaTwoCenterRecursion::compute -- am1 out of bounds", __FILE__, __LINE__);
    if (am2 < 0 || am2 > max_am2_)
        throw SanityCheckError("ERROR: ObaraSaikaTwoCenterRecursion::compute -- am2 out of bounds", __FILE__, __LINE__);

    int i,j;
    double pp = 1/(2*gamma);
    int lmaxi = am1;
    int lmaxj = am2;

    // Try zeroing out the matrices to fix dipoles
    memset(x_[0], 0, sizeof(double) * (max_am1_+1) * (max_am2_+1));
    memset(y_[0], 0, sizeof(double) * (max_am1_+1) * (max_am2_+1));
    memset(z_[0], 0, sizeof(double) * (max_am1_+1) * (max_am2_+1));

    x_[0][0] = y_[0][0] = z_[0][0] = 1.0;

    /* Upward recursion in j for i=0 */

    x_[0][1] = PB[0];
    y_[0][1] = PB[1];
    z_[0][1] = PB[2];

    for(j=1;j<lmaxj;j++) {
        x_[0][j+1] = PB[0]*x_[0][j];
        y_[0][j+1] = PB[1]*y_[0][j];
        z_[0][j+1] = PB[2]*z_[0][j];
        x_[0][j+1] += j*pp*x_[0][j-1];
        y_[0][j+1] += j*pp*y_[0][j-1];
        z_[0][j+1] += j*pp*z_[0][j-1];
    }

    /* Upward recursion in i for all j's */
    if (lmaxi > 0) {
        x_[1][0] = PA[0];
        y_[1][0] = PA[1];
        z_[1][0] = PA[2];
        for(j=1;j<=lmaxj;j++) {
            x_[1][j] = PA[0]*x_[0][j];
            y_[1][j] = PA[1]*y_[0][j];
            z_[1][j] = PA[2]*z_[0][j];
            x_[1][j] += j*pp*x_[0][j-1];
            y_[1][j] += j*pp*y_[0][j-1];
            z_[1][j] += j*pp*z_[0][j-1];
        }
        for(i=1;i<lmaxi;i++) {
            x_[i+1][0] = PA[0]*x_[i][0];
            y_[i+1][0] = PA[1]*y_[i][0];
            z_[i+1][0] = PA[2]*z_[i][0];
            x_[i+1][0] += i*pp*x_[i-1][0];
            y_[i+1][0] += i*pp*y_[i-1][0];
            z_[i+1][0] += i*pp*z_[i-1][0];
            for(j=1;j<=lmaxj;j++) {
                x_[i+1][j] = PA[0]*x_[i][j];
                y_[i+1][j] = PA[1]*y_[i][j];
                z_[i+1][j] = PA[2]*z_[i][j];
                x_[i+1][j] += i*pp*x_[i-1][j];
                y_[i+1][j] += i*pp*y_[i-1][j];
                z_[i+1][j] += i*pp*z_[i-1][j];
                x_[i+1][j] += j*pp*x_[i][j-1];
                y_[i+1][j] += j*pp*y_[i][j-1];
                z_[i+1][j] += j*pp*z_[i][j-1];
            }
        }
    }
}

ObaraSaikaTwoCenterVIRecursion::ObaraSaikaTwoCenterVIRecursion(int max_am1, int max_am2):
    max_am1_(max_am1), max_am2_(max_am2)
{
    if (max_am1 < 0)
        throw SanityCheckError("ERROR: ObaraSaikaTwoCenterVIRecursion -- max_am1 must be nonnegative", __FILE__, __LINE__);
    if (max_am2 < 0)
        throw SanityCheckError("ERROR: ObaraSaikaTwoCenterVIRecursion -- max_am2 must be nonnegative", __FILE__, __LINE__);

    size_ = max_am1 > max_am2 ? max_am1 : max_am2;
    size_ += 1;
    size_ = (size_-1)*size_*(size_+1)+1;
    vi_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
}

ObaraSaikaTwoCenterVIRecursion::~ObaraSaikaTwoCenterVIRecursion()
{
    free_box(vi_, size_, size_);
}

#define EPS 1.0e-17

void ObaraSaikaTwoCenterVIRecursion::calculate_f(double *F, int n, double t)
{
    int i, m;
    int m2;
    double t2;
    double num;
    double sum;
    double term1;
    static double K = 1.0/M_2_SQRTPI;
    double et;


    if (t>20.0){
        t2 = 2*t;
        et = exp(-t);
        t = sqrt(t);
        F[0] = K*erf(t)/t;
        for(m=0; m<=n-1; m++){
            F[m+1] = ((2*m + 1)*F[m] - et)/(t2);
        }
    }
    else {
        et = exp(-t);
        t2 = 2*t;
        m2 = 2*n;
        num = df[m2];
        i=0;
        sum = 1.0/(m2+1);
        do{
            i++;
            num = num*t2;
            term1 = num/df[m2+2*i+2];
            sum += term1;
        } while (fabs(term1) > EPS && i < MAX_FAC);
        F[n] = sum*et;
        for(m=n-1;m>=0;m--){
            F[m] = (t2*F[m+1] + et)/(2*m+1);
        }
    }
}

void ObaraSaikaTwoCenterVIRecursion::compute(double PA[3], double PB[3], double PC[3], double zeta, int am1, int am2)
{
    int a, b, m;
    int azm = 1;
    int aym = am1 + 1;
    int axm = aym * aym;
    int bzm = 1;
    int bym = am2 + 1;
    int bxm = bym * bym;
    int ax, ay, az, bx, by, bz;
    int aind, bind;
    double ooz = 1.0/(2.0 * zeta);
    int mmax = max_am1_ + max_am2_;

    // Prefactor from A20
    double tmp = sqrt(zeta) * M_2_SQRTPI;
    // U from A21
    double u = zeta * (PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2]);
    double *F = new double[mmax+1];

    // Form Fm(U) from A20
    calculate_f(F, mmax, u);

    // Think we're having problems with values being left over.
    //zero_box(vi_, size_, size_, mmax + 1);

    // Perform recursion in m for (a|A(0)|s) using A20
    for (m=0; m<=mmax; ++m) {
        vi_[0][0][m] = tmp * F[m];
    }

    // Perform recursion in b with a=0
    //  subset of A19
    for (b=1; b<=am2; ++b) {
        for (bx=0; bx<=b; ++bx) {
            for (by=0; by<=b-bx; ++by) {
                bz = b-bx-by;

                // Compute the index into VI for bx,by,bz
                bind = bx*bxm + by*bym + bz*bzm;

                // Compute each x, y, z contribution
                if (bz > 0) {
                    for (m=0; m<=mmax-b; ++m) {
                        vi_[0][bind][m] = PB[2] * vi_[0][bind-bzm][m] - PC[2] * vi_[0][bind-bzm][m+1];
                    }
                    if (bz > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            vi_[0][bind][m] += ooz * (bz-1) * (vi_[0][bind-2*bzm][m] - vi_[0][bind-2*bzm][m+1]);
                        }
                    }
                }
                else if (by > 0) {
                    for (m=0; m<=mmax-b; ++m) {
                        vi_[0][bind][m] = PB[1] * vi_[0][bind-bym][m] - PC[1] * vi_[0][bind-bym][m+1];
                    }
                    if (by > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            vi_[0][bind][m] += ooz * (by-1) * (vi_[0][bind-2*bym][m] - vi_[0][bind-2*bym][m+1]);
                        }
                    }
                }
                else if (bx > 0) {
                    for (m=0; m<=mmax-b; ++m) {
                        vi_[0][bind][m] = PB[0] * vi_[0][bind-bxm][m] - PC[0] * vi_[0][bind-bxm][m+1];
                    }
                    if (bx > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            vi_[0][bind][m] += ooz * (bx-1) * (vi_[0][bind-2*bxm][m] - vi_[0][bind-2*bxm][m+1]);
                        }
                    }
                }
            }
        }
    }

    // Perform upward recursion in a with all b's
    for (b=0; b<=am2; b++) {
        for (bx=0; bx<=b; bx++) {
            for (by=0; by<=b-bx;by++) {
                bz = b-bx-by;
                bind = bx*bxm + by*bym + bz*bzm;

                for (a=1; a<=am1; a++) {
                    // This next for loop was for (ax=0; ax<=b; ax++)
                    // this could explain why dx2 was not being computed.
                    // change for for(ax=0; ax<a; ax++) on 2005-09-15 4:11pm
                    for (ax=0; ax<=a; ax++) {
                        for (ay=0; ay<=a-ax; ay++) {
                            az = a-ax-ay;
                            aind = ax*axm + ay*aym + az*azm;

                            if (az > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    vi_[aind][bind][m] = PA[2] * vi_[aind-azm][bind][m] - PC[2] * vi_[aind-azm][bind][m+1];
                                }

                                if (az > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * (az-1) * (vi_[aind-2*azm][bind][m] - vi_[aind-2*azm][bind][m+1]);
                                    }
                                }
                                if (bz > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * bz * (vi_[aind-azm][bind-bzm][m] - vi_[aind-azm][bind-bzm][m+1]);
                                    }
                                }
                            }
                            else if (ay > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    vi_[aind][bind][m] = PA[1] * vi_[aind-aym][bind][m] - PC[1] * vi_[aind-aym][bind][m+1];
                                }
                                if (ay > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * (ay-1) * (vi_[aind-2*aym][bind][m] - vi_[aind-2*aym][bind][m+1]);
                                    }
                                }
                                if (by > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * by * (vi_[aind-aym][bind-bym][m] - vi_[aind-aym][bind-bym][m+1]);
                                    }
                                }
                            }
                            else if (ax > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    vi_[aind][bind][m] = PA[0] * vi_[aind-axm][bind][m] - PC[0] * vi_[aind-axm][bind][m+1];
                                }

                                if (ax > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * (ax-1) * (vi_[aind-2*axm][bind][m] - vi_[aind-2*axm][bind][m+1]);
                                    }
                                }
                                if (bx > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * bx * (vi_[aind-axm][bind-bxm][m] - vi_[aind-axm][bind-bxm][m+1]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    delete[] F;

}

void ObaraSaikaTwoCenterVIRecursion::compute_erf(double PA[3], double PB[3], double PC[3], double zeta, int am1, int am2, double zetam)
{
    int a, b, m;
    int azm = 1;
    int aym = am1 + 1;
    int axm = aym * aym;
    int bzm = 1;
    int bym = am2 + 1;
    int bxm = bym * bym;
    int ax, ay, az, bx, by, bz;
    int aind, bind;
    double ooz = 1.0/(2.0 * zeta);
    int mmax = max_am1_ + max_am2_;

    // Prefactor from A20
    double tmp = sqrt(zetam) * M_2_SQRTPI;
    // U from A21
    double u = zetam * (PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2]);
    double *F = new double[mmax+1];

    // Form Fm(U) from A20
    calculate_f(F, mmax, u);

    // Think we're having problems with values being left over.
    //zero_box(vi_, size_, size_, mmax + 1);

    // Backsub, to account for the fact that Jet/OS roll the \zeta^m in later
    double F_prefac = 1.0;
    double T_prefac = zetam / zeta;

    // Perform recursion in m for (a|A(0)|s) using A20
    for (m=0; m<=mmax; ++m) {
        vi_[0][0][m] = tmp * F_prefac * F[m];
        F_prefac *= T_prefac;
    }

    // Perform recursion in b with a=0
    //  subset of A19
    for (b=1; b<=am2; ++b) {
        for (bx=0; bx<=b; ++bx) {
            for (by=0; by<=b-bx; ++by) {
                bz = b-bx-by;

                // Compute the index into VI for bx,by,bz
                bind = bx*bxm + by*bym + bz*bzm;

                // Compute each x, y, z contribution
                if (bz > 0) {
                    for (m=0; m<=mmax-b; ++m) {
                        vi_[0][bind][m] = PB[2] * vi_[0][bind-bzm][m] - PC[2] * vi_[0][bind-bzm][m+1];
                    }
                    if (bz > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            vi_[0][bind][m] += ooz * (bz-1) * (vi_[0][bind-2*bzm][m] - vi_[0][bind-2*bzm][m+1]);
                        }
                    }
                }
                else if (by > 0) {
                    for (m=0; m<=mmax-b; ++m) {
                        vi_[0][bind][m] = PB[1] * vi_[0][bind-bym][m] - PC[1] * vi_[0][bind-bym][m+1];
                    }
                    if (by > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            vi_[0][bind][m] += ooz * (by-1) * (vi_[0][bind-2*bym][m] - vi_[0][bind-2*bym][m+1]);
                        }
                    }
                }
                else if (bx > 0) {
                    for (m=0; m<=mmax-b; ++m) {
                        vi_[0][bind][m] = PB[0] * vi_[0][bind-bxm][m] - PC[0] * vi_[0][bind-bxm][m+1];
                    }
                    if (bx > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            vi_[0][bind][m] += ooz * (bx-1) * (vi_[0][bind-2*bxm][m] - vi_[0][bind-2*bxm][m+1]);
                        }
                    }
                }
            }
        }
    }

    // Perform upward recursion in a with all b's
    for (b=0; b<=am2; b++) {
        for (bx=0; bx<=b; bx++) {
            for (by=0; by<=b-bx;by++) {
                bz = b-bx-by;
                bind = bx*bxm + by*bym + bz*bzm;

                for (a=1; a<=am1; a++) {
                    // This next for loop was for (ax=0; ax<=b; ax++)
                    // this could explain why dx2 was not being computed.
                    // change for for(ax=0; ax<a; ax++) on 2005-09-15 4:11pm
                    for (ax=0; ax<=a; ax++) {
                        for (ay=0; ay<=a-ax; ay++) {
                            az = a-ax-ay;
                            aind = ax*axm + ay*aym + az*azm;

                            if (az > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    vi_[aind][bind][m] = PA[2] * vi_[aind-azm][bind][m] - PC[2] * vi_[aind-azm][bind][m+1];
                                }

                                if (az > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * (az-1) * (vi_[aind-2*azm][bind][m] - vi_[aind-2*azm][bind][m+1]);
                                    }
                                }
                                if (bz > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * bz * (vi_[aind-azm][bind-bzm][m] - vi_[aind-azm][bind-bzm][m+1]);
                                    }
                                }
                            }
                            else if (ay > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    vi_[aind][bind][m] = PA[1] * vi_[aind-aym][bind][m] - PC[1] * vi_[aind-aym][bind][m+1];
                                }
                                if (ay > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * (ay-1) * (vi_[aind-2*aym][bind][m] - vi_[aind-2*aym][bind][m+1]);
                                    }
                                }
                                if (by > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * by * (vi_[aind-aym][bind-bym][m] - vi_[aind-aym][bind-bym][m+1]);
                                    }
                                }
                            }
                            else if (ax > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    vi_[aind][bind][m] = PA[0] * vi_[aind-axm][bind][m] - PC[0] * vi_[aind-axm][bind][m+1];
                                }

                                if (ax > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * (ax-1) * (vi_[aind-2*axm][bind][m] - vi_[aind-2*axm][bind][m+1]);
                                    }
                                }
                                if (bx > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * bx * (vi_[aind-axm][bind-bxm][m] - vi_[aind-axm][bind-bxm][m+1]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    delete[] F;

}

ObaraSaikaTwoCenterVIDerivRecursion::ObaraSaikaTwoCenterVIDerivRecursion(int max_am1, int max_am2)
    : ObaraSaikaTwoCenterVIRecursion(max_am1+1, max_am2+1)
{
    vx_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    vy_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    vz_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
}

ObaraSaikaTwoCenterVIDerivRecursion::~ObaraSaikaTwoCenterVIDerivRecursion()
{
    free_box(vx_, size_, size_);
    free_box(vy_, size_, size_);
    free_box(vz_, size_, size_);
}

void ObaraSaikaTwoCenterVIDerivRecursion::compute(double PA[3], double PB[3], double PC[3], double zeta, int am1, int am2)
{
    int a, b, m;
    int azm = 1;
    int aym = am1 + 1;
    int axm = aym * aym;
    int bzm = 1;
    int bym = am2 + 1;
    int bxm = bym * bym;
    int ax, ay, az, bx, by, bz;
    int aind, bind;
    double ooz = 1.0/(2.0 * zeta);
    int mmax = am1 + am2;

    // Prefactor from A20
    double tmp = sqrt(zeta) * M_2_SQRTPI;
    // U from A21
    double u = zeta * (PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2]);
    double *F = new double[mmax+1];

    // Zero out F
    memset(F, 0, sizeof(double) * (mmax+1));

    // Form Fm(U) from A20
    calculate_f(F, mmax, u);

    // Perform recursion in m for (a|A(0)|s) using A20
    for (m=0; m<=mmax; ++m) {
        vi_[0][0][m] = tmp * F[m];
    }
    for (m=0; m<=mmax-1; ++m) {
        vx_[0][0][m] = 2.0*zeta*PC[0]*vi_[0][0][m+1];
        vy_[0][0][m] = 2.0*zeta*PC[1]*vi_[0][0][m+1];
        vz_[0][0][m] = 2.0*zeta*PC[2]*vi_[0][0][m+1];
    }

    // Perform recursion in b with a=0
    //  subset of A19
    for (b=1; b<=am2; ++b) {
        for (bx=0; bx<=b; ++bx) {
            for (by=0; by<=b-bx; ++by) {
                bz = b-bx-by;

                // Compute the index into VI for bx,by,bz
                bind = bx*bxm + by*bym + bz*bzm;

                // Compute each x, y, z contribution
                if (bz > 0) {
                    for (m=0; m<=mmax-b; ++m) {
                        vi_[0][bind][m] = PB[2] * vi_[0][bind-bzm][m] - PC[2] * vi_[0][bind-bzm][m+1];
                    }
                    for (m=0; m<=mmax-b-1; ++m) {
                        vx_[0][bind][m] = PB[2] * vx_[0][bind-bzm][m] - PC[2] * vx_[0][bind-bzm][m+1];
                        vy_[0][bind][m] = PB[2] * vy_[0][bind-bzm][m] - PC[2] * vy_[0][bind-bzm][m+1];
                        vz_[0][bind][m] = PB[2] * vz_[0][bind-bzm][m] - PC[2] * vz_[0][bind-bzm][m+1] + vi_[0][bind-bzm][m+1];
                    }
                    if (bz > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            vi_[0][bind][m] += ooz * (bz-1) * (vi_[0][bind-2*bzm][m] - vi_[0][bind-2*bzm][m+1]);
                        }
                        for (m=0; m<=mmax-b-1; ++m) {
                            vx_[0][bind][m] += ooz * (bz-1) * (vx_[0][bind-2*bzm][m] - vx_[0][bind-2*bzm][m+1]);
                            vy_[0][bind][m] += ooz * (bz-1) * (vy_[0][bind-2*bzm][m] - vy_[0][bind-2*bzm][m+1]);
                            vz_[0][bind][m] += ooz * (bz-1) * (vz_[0][bind-2*bzm][m] - vz_[0][bind-2*bzm][m+1]);
                        }
                    }
                }
                else if (by > 0) {
                    for (m=0; m<=mmax-b; ++m) {
                        vi_[0][bind][m] = PB[1] * vi_[0][bind-bym][m] - PC[1] * vi_[0][bind-bym][m+1];
                    }
                    for (m=0; m<=mmax-b-1; ++m) {
                        vx_[0][bind][m] = PB[1] * vx_[0][bind-bym][m] - PC[1] * vx_[0][bind-bym][m+1];
                        vy_[0][bind][m] = PB[1] * vy_[0][bind-bym][m] - PC[1] * vy_[0][bind-bym][m+1] + vi_[0][bind-bym][m+1];
                        vz_[0][bind][m] = PB[1] * vz_[0][bind-bym][m] - PC[1] * vz_[0][bind-bym][m+1];
                    }
                    if (by > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            vi_[0][bind][m] += ooz * (by-1) * (vi_[0][bind-2*bym][m] - vi_[0][bind-2*bym][m+1]);
                        }
                        for (m=0; m<=mmax-b-1; ++m) {
                            vx_[0][bind][m] += ooz * (by-1) * (vx_[0][bind-2*bym][m] - vx_[0][bind-2*bym][m+1]);
                            vy_[0][bind][m] += ooz * (by-1) * (vy_[0][bind-2*bym][m] - vy_[0][bind-2*bym][m+1]);
                            vz_[0][bind][m] += ooz * (by-1) * (vz_[0][bind-2*bym][m] - vz_[0][bind-2*bym][m+1]);
                        }
                    }
                }
                else if (bx > 0) {
                    for (m=0; m<=mmax-b; ++m) {
                        vi_[0][bind][m] = PB[0] * vi_[0][bind-bxm][m] - PC[0] * vi_[0][bind-bxm][m+1];
                    }
                    for (m=0; m<=mmax-b-1; ++m) {
                        vx_[0][bind][m] = PB[0] * vx_[0][bind-bxm][m] - PC[0] * vx_[0][bind-bxm][m+1] + vi_[0][bind-bxm][m+1];
                        vy_[0][bind][m] = PB[0] * vy_[0][bind-bxm][m] - PC[0] * vy_[0][bind-bxm][m+1];
                        vz_[0][bind][m] = PB[0] * vz_[0][bind-bxm][m] - PC[0] * vz_[0][bind-bxm][m+1];
                    }
                    if (bx > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            vi_[0][bind][m] += ooz * (bx-1) * (vi_[0][bind-2*bxm][m] - vi_[0][bind-2*bxm][m+1]);
                        }
                        for (m=0; m<=mmax-b-1; ++m) {
                            vx_[0][bind][m] += ooz * (bx-1) * (vx_[0][bind-2*bxm][m] - vx_[0][bind-2*bxm][m+1]);
                            vy_[0][bind][m] += ooz * (bx-1) * (vy_[0][bind-2*bxm][m] - vy_[0][bind-2*bxm][m+1]);
                            vz_[0][bind][m] += ooz * (bx-1) * (vz_[0][bind-2*bxm][m] - vz_[0][bind-2*bxm][m+1]);
                        }
                    }
                }
            }
        }
    }

    // Perform upward recursion in a with all b's
    for (b=0; b<=am2; b++) {
        for (bx=0; bx<=b; bx++) {
            for (by=0; by<=b-bx;by++) {
                bz = b-bx-by;
                bind = bx*bxm + by*bym + bz*bzm;

                for (a=1; a<=am1; a++) {
                    // This next for loop was for (ax=0; ax<=b; ax++)
                    // this could explain why dx2 was not being computed.
                    // change for for(ax=0; ax<a; ax++) on 2005-09-15 4:11pm
                    for (ax=0; ax<=a; ax++) {
                        for (ay=0; ay<=a-ax; ay++) {
                            az = a-ax-ay;
                            aind = ax*axm + ay*aym + az*azm;

                            if (az > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    vi_[aind][bind][m] = PA[2] * vi_[aind-azm][bind][m] - PC[2] * vi_[aind-azm][bind][m+1];
                                }
                                for (m=0; m<=mmax-a-b-1; ++m) {
                                    vx_[aind][bind][m] = PA[2] * vx_[aind-azm][bind][m] - PC[2] * vx_[aind-azm][bind][m+1];
                                    vy_[aind][bind][m] = PA[2] * vy_[aind-azm][bind][m] - PC[2] * vy_[aind-azm][bind][m+1];
                                    vz_[aind][bind][m] = PA[2] * vz_[aind-azm][bind][m] - PC[2] * vz_[aind-azm][bind][m+1] + vi_[aind-azm][bind][m+1];
                                }

                                if (az > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * (az-1) * (vi_[aind-2*azm][bind][m] - vi_[aind-2*azm][bind][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        vx_[aind][bind][m] += ooz * (az-1) * (vx_[aind-2*azm][bind][m] - vx_[aind-2*azm][bind][m+1]);
                                        vy_[aind][bind][m] += ooz * (az-1) * (vy_[aind-2*azm][bind][m] - vy_[aind-2*azm][bind][m+1]);
                                        vz_[aind][bind][m] += ooz * (az-1) * (vz_[aind-2*azm][bind][m] - vz_[aind-2*azm][bind][m+1]);
                                    }
                                }
                                if (bz > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * bz * (vi_[aind-azm][bind-bzm][m] - vi_[aind-azm][bind-bzm][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        vx_[aind][bind][m] += ooz * bz * (vx_[aind-azm][bind-bzm][m] - vx_[aind-azm][bind-bzm][m+1]);
                                        vy_[aind][bind][m] += ooz * bz * (vy_[aind-azm][bind-bzm][m] - vy_[aind-azm][bind-bzm][m+1]);
                                        vz_[aind][bind][m] += ooz * bz * (vz_[aind-azm][bind-bzm][m] - vz_[aind-azm][bind-bzm][m+1]);
                                    }
                                }
                            }
                            else if (ay > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    vi_[aind][bind][m] = PA[1] * vi_[aind-aym][bind][m] - PC[1] * vi_[aind-aym][bind][m+1];
                                }
                                for (m=0; m<=mmax-a-b-1; ++m) {
                                    vx_[aind][bind][m] = PA[1] * vx_[aind-aym][bind][m] - PC[1] * vx_[aind-aym][bind][m+1];
                                    vy_[aind][bind][m] = PA[1] * vy_[aind-aym][bind][m] - PC[1] * vy_[aind-aym][bind][m+1] + vi_[aind-aym][bind][m+1];
                                    vz_[aind][bind][m] = PA[1] * vz_[aind-aym][bind][m] - PC[1] * vz_[aind-aym][bind][m+1];
                                }
                                if (ay > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * (ay-1) * (vi_[aind-2*aym][bind][m] - vi_[aind-2*aym][bind][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        vx_[aind][bind][m] += ooz * (ay-1) * (vx_[aind-2*aym][bind][m] - vx_[aind-2*aym][bind][m+1]);
                                        vy_[aind][bind][m] += ooz * (ay-1) * (vy_[aind-2*aym][bind][m] - vy_[aind-2*aym][bind][m+1]);
                                        vz_[aind][bind][m] += ooz * (ay-1) * (vz_[aind-2*aym][bind][m] - vz_[aind-2*aym][bind][m+1]);
                                    }
                                }
                                if (by > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * by * (vi_[aind-aym][bind-bym][m] - vi_[aind-aym][bind-bym][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        vx_[aind][bind][m] += ooz * by * (vx_[aind-aym][bind-bym][m] - vx_[aind-aym][bind-bym][m+1]);
                                        vy_[aind][bind][m] += ooz * by * (vy_[aind-aym][bind-bym][m] - vy_[aind-aym][bind-bym][m+1]);
                                        vz_[aind][bind][m] += ooz * by * (vz_[aind-aym][bind-bym][m] - vz_[aind-aym][bind-bym][m+1]);
                                    }
                                }
                            }
                            else if (ax > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    vi_[aind][bind][m] = PA[0] * vi_[aind-axm][bind][m] - PC[0] * vi_[aind-axm][bind][m+1];
                                }
                                for (m=0; m<=mmax-a-b-1; ++m) {
                                    vx_[aind][bind][m] = PA[0] * vx_[aind-axm][bind][m] - PC[0] * vx_[aind-axm][bind][m+1] + vi_[aind-axm][bind][m+1];
                                    vy_[aind][bind][m] = PA[0] * vy_[aind-axm][bind][m] - PC[0] * vy_[aind-axm][bind][m+1];
                                    vz_[aind][bind][m] = PA[0] * vz_[aind-axm][bind][m] - PC[0] * vz_[aind-axm][bind][m+1];
                                }

                                if (ax > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * (ax-1) * (vi_[aind-2*axm][bind][m] - vi_[aind-2*axm][bind][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        vx_[aind][bind][m] += ooz * (ax-1) * (vx_[aind-2*axm][bind][m] - vx_[aind-2*axm][bind][m+1]);
                                        vy_[aind][bind][m] += ooz * (ax-1) * (vy_[aind-2*axm][bind][m] - vy_[aind-2*axm][bind][m+1]);
                                        vz_[aind][bind][m] += ooz * (ax-1) * (vz_[aind-2*axm][bind][m] - vz_[aind-2*axm][bind][m+1]);
                                    }
                                }
                                if (bx > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * bx * (vi_[aind-axm][bind-bxm][m] - vi_[aind-axm][bind-bxm][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        vx_[aind][bind][m] += ooz * bx * (vx_[aind-axm][bind-bxm][m] - vx_[aind-axm][bind-bxm][m+1]);
                                        vy_[aind][bind][m] += ooz * bx * (vy_[aind-axm][bind-bxm][m] - vy_[aind-axm][bind-bxm][m+1]);
                                        vz_[aind][bind][m] += ooz * bx * (vz_[aind-axm][bind-bxm][m] - vz_[aind-axm][bind-bxm][m+1]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    delete[] F;
}

ObaraSaikaTwoCenterVIDeriv2Recursion::ObaraSaikaTwoCenterVIDeriv2Recursion(int max_am1, int max_am2)
    : ObaraSaikaTwoCenterVIDerivRecursion(max_am1+1, max_am2+1)
{
    int max_am = 2*(std::max(max_am1, max_am2) + 2) + 1;

    vxx_ = init_box(size_, size_, max_am);
    vxy_ = init_box(size_, size_, max_am);
    vxz_ = init_box(size_, size_, max_am);
    vyy_ = init_box(size_, size_, max_am);
    vyz_ = init_box(size_, size_, max_am);
    vzz_ = init_box(size_, size_, max_am);
}

ObaraSaikaTwoCenterVIDeriv2Recursion::~ObaraSaikaTwoCenterVIDeriv2Recursion()
{
    free_box(vxx_, size_, size_);
    free_box(vxy_, size_, size_);
    free_box(vxz_, size_, size_);
    free_box(vyy_, size_, size_);
    free_box(vyz_, size_, size_);
    free_box(vzz_, size_, size_);
}

void ObaraSaikaTwoCenterVIDeriv2Recursion::compute(double PA[3], double PB[3], double PC[3], double zeta, int am1, int am2)
{
    int a, b, m;
    int azm = 1;
    int aym = am1 + 1;
    int axm = aym * aym;
    int bzm = 1;
    int bym = am2 + 1;
    int bxm = bym * bym;
    int ax, ay, az, bx, by, bz;
    int aind, bind;
    double ooz = 1.0/(2.0 * zeta);
    int mmax = am1 + am2;

    // Prefactor from A20
    double tmp = sqrt(zeta) * M_2_SQRTPI;
    // U from A21
    double u = zeta * (PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2]);
    double *F = new double[mmax+1]; // TODO: Move this allocation into constructor

    // Zero out F
    memset(F, 0, sizeof(double) * (mmax+1));

    // Form Fm(U) from A20
    calculate_f(F, mmax, u);

    // Perform recursion in m for (a|A(0)|s) using A20
    for (m=0; m<=mmax; ++m) {
        vi_[0][0][m] = tmp * F[m];
    }
    for (m=0; m<=mmax-1; ++m) {
        vx_[0][0][m] = 2.0*zeta*PC[0]*vi_[0][0][m+1];
        vy_[0][0][m] = 2.0*zeta*PC[1]*vi_[0][0][m+1];
        vz_[0][0][m] = 2.0*zeta*PC[2]*vi_[0][0][m+1];
    }
    for (m=0; m<=mmax-2; ++m) {
        vxx_[0][0][m] = 4*zeta*zeta*PC[0]*PC[0]*vi_[0][0][m+2] - 2.0*zeta*vi_[0][0][m+1];
        vyy_[0][0][m] = 4*zeta*zeta*PC[1]*PC[1]*vi_[0][0][m+2] - 2.0*zeta*vi_[0][0][m+1];
        vzz_[0][0][m] = 4*zeta*zeta*PC[2]*PC[2]*vi_[0][0][m+2] - 2.0*zeta*vi_[0][0][m+1];
        vxy_[0][0][m] = 4*zeta*zeta*PC[0]*PC[1]*vi_[0][0][m+2];
        vxz_[0][0][m] = 4*zeta*zeta*PC[0]*PC[2]*vi_[0][0][m+2];
        vyz_[0][0][m] = 4*zeta*zeta*PC[1]*PC[2]*vi_[0][0][m+2];
    }

    // Perform recursion in b with a=0
    //  subset of A19
    for (b=1; b<=am2; ++b) {
        for (bx=0; bx<=b; ++bx) {
            for (by=0; by<=b-bx; ++by) {
                bz = b-bx-by;

                // Compute the index into VI for bx,by,bz
                bind = bx*bxm + by*bym + bz*bzm;

                // Compute each x, y, z contribution
                if (bz > 0) {
                    for (m=0; m<=mmax-b; ++m) { /* Electrostatic potential integrals */
                        vi_[0][bind][m] = PB[2] * vi_[0][bind-bzm][m] - PC[2] * vi_[0][bind-bzm][m+1];
                    }
                    for (m=0; m<=mmax-b-1; ++m) { /* Electric field integrals */
                        vx_[0][bind][m] = PB[2] * vx_[0][bind-bzm][m] - PC[2] * vx_[0][bind-bzm][m+1];
                        vy_[0][bind][m] = PB[2] * vy_[0][bind-bzm][m] - PC[2] * vy_[0][bind-bzm][m+1];
                        vz_[0][bind][m] = PB[2] * vz_[0][bind-bzm][m] - PC[2] * vz_[0][bind-bzm][m+1] + vi_[0][bind-bzm][m+1];
                    }
                    for (m=0; m<=mmax-b-2; ++m) { /* Gradients of the electric field */
                        vxx_[0][bind][m] = PB[2]*vxx_[0][bind-bzm][m] - PC[2]*vxx_[0][bind-bzm][m+1];
                        vyy_[0][bind][m] = PB[2]*vyy_[0][bind-bzm][m] - PC[2]*vyy_[0][bind-bzm][m+1];
                        vzz_[0][bind][m] = PB[2]*vzz_[0][bind-bzm][m] - PC[2]*vzz_[0][bind-bzm][m+1] + 2*vz_[0][bind-bzm][m+1];
                        vxy_[0][bind][m] = PB[2]*vxy_[0][bind-bzm][m] - PC[2]*vxy_[0][bind-bzm][m+1];
                        vxz_[0][bind][m] = PB[2]*vxz_[0][bind-bzm][m] - PC[2]*vxz_[0][bind-bzm][m+1] + vx_[0][bind-bzm][m+1];
                        vyz_[0][bind][m] = PB[2]*vyz_[0][bind-bzm][m] - PC[2]*vyz_[0][bind-bzm][m+1] + vy_[0][bind-bzm][m+1];
                    }
                    if (bz > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            vi_[0][bind][m] += ooz * (bz-1) * (vi_[0][bind-2*bzm][m] - vi_[0][bind-2*bzm][m+1]);
                        }
                        for (m=0; m<=mmax-b-1; ++m) {
                            vx_[0][bind][m] += ooz * (bz-1) * (vx_[0][bind-2*bzm][m] - vx_[0][bind-2*bzm][m+1]);
                            vy_[0][bind][m] += ooz * (bz-1) * (vy_[0][bind-2*bzm][m] - vy_[0][bind-2*bzm][m+1]);
                            vz_[0][bind][m] += ooz * (bz-1) * (vz_[0][bind-2*bzm][m] - vz_[0][bind-2*bzm][m+1]);
                        }
                        for (m=0; m<=mmax-b-2; m++) {
                            vxx_[0][bind][m] += ooz*(bz-1)*(vxx_[0][bind-2*bzm][m] - vxx_[0][bind-2*bzm][m+1]);
                            vyy_[0][bind][m] += ooz*(bz-1)*(vyy_[0][bind-2*bzm][m] - vyy_[0][bind-2*bzm][m+1]);
                            vzz_[0][bind][m] += ooz*(bz-1)*(vzz_[0][bind-2*bzm][m] - vzz_[0][bind-2*bzm][m+1]);
                            vxy_[0][bind][m] += ooz*(bz-1)*(vxy_[0][bind-2*bzm][m] - vxy_[0][bind-2*bzm][m+1]);
                            vxz_[0][bind][m] += ooz*(bz-1)*(vxz_[0][bind-2*bzm][m] - vxz_[0][bind-2*bzm][m+1]);
                            vyz_[0][bind][m] += ooz*(bz-1)*(vyz_[0][bind-2*bzm][m] - vyz_[0][bind-2*bzm][m+1]);
                        }
                    }
                }
                else if (by > 0) {
                    for (m=0; m<=mmax-b; ++m) {
                        vi_[0][bind][m] = PB[1] * vi_[0][bind-bym][m] - PC[1] * vi_[0][bind-bym][m+1];
                    }
                    for (m=0; m<=mmax-b-1; ++m) {
                        vx_[0][bind][m] = PB[1] * vx_[0][bind-bym][m] - PC[1] * vx_[0][bind-bym][m+1];
                        vy_[0][bind][m] = PB[1] * vy_[0][bind-bym][m] - PC[1] * vy_[0][bind-bym][m+1] + vi_[0][bind-bym][m+1];
                        vz_[0][bind][m] = PB[1] * vz_[0][bind-bym][m] - PC[1] * vz_[0][bind-bym][m+1];
                    }
                    for(m=0;m<=mmax-b-2;m++) {
                        vxx_[0][bind][m] = PB[1]*vxx_[0][bind-bym][m] - PC[1]*vxx_[0][bind-bym][m+1];
                        vyy_[0][bind][m] = PB[1]*vyy_[0][bind-bym][m] - PC[1]*vyy_[0][bind-bym][m+1] + 2*vy_[0][bind-bym][m+1];
                        vzz_[0][bind][m] = PB[1]*vzz_[0][bind-bym][m] - PC[1]*vzz_[0][bind-bym][m+1];
                        vxy_[0][bind][m] = PB[1]*vxy_[0][bind-bym][m] - PC[1]*vxy_[0][bind-bym][m+1] + vx_[0][bind-bym][m+1];
                        vxz_[0][bind][m] = PB[1]*vxz_[0][bind-bym][m] - PC[1]*vxz_[0][bind-bym][m+1];
                        vyz_[0][bind][m] = PB[1]*vyz_[0][bind-bym][m] - PC[1]*vyz_[0][bind-bym][m+1] + vz_[0][bind-bym][m+1];
                    }

                    if (by > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            vi_[0][bind][m] += ooz * (by-1) * (vi_[0][bind-2*bym][m] - vi_[0][bind-2*bym][m+1]);
                        }
                        for (m=0; m<=mmax-b-1; ++m) {
                            vx_[0][bind][m] += ooz * (by-1) * (vx_[0][bind-2*bym][m] - vx_[0][bind-2*bym][m+1]);
                            vy_[0][bind][m] += ooz * (by-1) * (vy_[0][bind-2*bym][m] - vy_[0][bind-2*bym][m+1]);
                            vz_[0][bind][m] += ooz * (by-1) * (vz_[0][bind-2*bym][m] - vz_[0][bind-2*bym][m+1]);
                        }
                        for(m=0;m<=mmax-b-2;m++) {
                            vxx_[0][bind][m] += ooz*(by-1)*(vxx_[0][bind-2*bym][m] - vxx_[0][bind-2*bym][m+1]);
                            vyy_[0][bind][m] += ooz*(by-1)*(vyy_[0][bind-2*bym][m] - vyy_[0][bind-2*bym][m+1]);
                            vzz_[0][bind][m] += ooz*(by-1)*(vzz_[0][bind-2*bym][m] - vzz_[0][bind-2*bym][m+1]);
                            vxy_[0][bind][m] += ooz*(by-1)*(vxy_[0][bind-2*bym][m] - vxy_[0][bind-2*bym][m+1]);
                            vxz_[0][bind][m] += ooz*(by-1)*(vxz_[0][bind-2*bym][m] - vxz_[0][bind-2*bym][m+1]);
                            vyz_[0][bind][m] += ooz*(by-1)*(vyz_[0][bind-2*bym][m] - vyz_[0][bind-2*bym][m+1]);
                        }
                    }
                }
                else if (bx > 0) {
                    for (m=0; m<=mmax-b; ++m) {
                        vi_[0][bind][m] = PB[0] * vi_[0][bind-bxm][m] - PC[0] * vi_[0][bind-bxm][m+1];
                    }
                    for (m=0; m<=mmax-b-1; ++m) {
                        vx_[0][bind][m] = PB[0] * vx_[0][bind-bxm][m] - PC[0] * vx_[0][bind-bxm][m+1] + vi_[0][bind-bxm][m+1];
                        vy_[0][bind][m] = PB[0] * vy_[0][bind-bxm][m] - PC[0] * vy_[0][bind-bxm][m+1];
                        vz_[0][bind][m] = PB[0] * vz_[0][bind-bxm][m] - PC[0] * vz_[0][bind-bxm][m+1];
                    }
                    for(m=0;m<=mmax-b-2;m++) {
                        vxx_[0][bind][m] = PB[0]*vxx_[0][bind-bxm][m] - PC[0]*vxx_[0][bind-bxm][m+1] + 2*vx_[0][bind-bxm][m+1];
                        vyy_[0][bind][m] = PB[0]*vyy_[0][bind-bxm][m] - PC[0]*vyy_[0][bind-bxm][m+1];
                        vzz_[0][bind][m] = PB[0]*vzz_[0][bind-bxm][m] - PC[0]*vzz_[0][bind-bxm][m+1];
                        vxy_[0][bind][m] = PB[0]*vxy_[0][bind-bxm][m] - PC[0]*vxy_[0][bind-bxm][m+1] + vy_[0][bind-bxm][m+1];
                        vxz_[0][bind][m] = PB[0]*vxz_[0][bind-bxm][m] - PC[0]*vxz_[0][bind-bxm][m+1] + vz_[0][bind-bxm][m+1];
                        vyz_[0][bind][m] = PB[0]*vyz_[0][bind-bxm][m] - PC[0]*vyz_[0][bind-bxm][m+1];
                    }

                    if (bx > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            vi_[0][bind][m] += ooz * (bx-1) * (vi_[0][bind-2*bxm][m] - vi_[0][bind-2*bxm][m+1]);
                        }
                        for (m=0; m<=mmax-b-1; ++m) {
                            vx_[0][bind][m] += ooz * (bx-1) * (vx_[0][bind-2*bxm][m] - vx_[0][bind-2*bxm][m+1]);
                            vy_[0][bind][m] += ooz * (bx-1) * (vy_[0][bind-2*bxm][m] - vy_[0][bind-2*bxm][m+1]);
                            vz_[0][bind][m] += ooz * (bx-1) * (vz_[0][bind-2*bxm][m] - vz_[0][bind-2*bxm][m+1]);
                        }
                        for(m=0;m<=mmax-b-2;m++) {
                            vxx_[0][bind][m] += ooz*(bx-1)*(vxx_[0][bind-2*bxm][m] - vxx_[0][bind-2*bxm][m+1]);
                            vyy_[0][bind][m] += ooz*(bx-1)*(vyy_[0][bind-2*bxm][m] - vyy_[0][bind-2*bxm][m+1]);
                            vzz_[0][bind][m] += ooz*(bx-1)*(vzz_[0][bind-2*bxm][m] - vzz_[0][bind-2*bxm][m+1]);
                            vxy_[0][bind][m] += ooz*(bx-1)*(vxy_[0][bind-2*bxm][m] - vxy_[0][bind-2*bxm][m+1]);
                            vxz_[0][bind][m] += ooz*(bx-1)*(vxz_[0][bind-2*bxm][m] - vxz_[0][bind-2*bxm][m+1]);
                            vyz_[0][bind][m] += ooz*(bx-1)*(vyz_[0][bind-2*bxm][m] - vyz_[0][bind-2*bxm][m+1]);
                        }
                    }
                }
            }
        }
    }

    // Perform upward recursion in a with all b's
    for (b=0; b<=am2; b++) {
        for (bx=0; bx<=b; bx++) {
            for (by=0; by<=b-bx;by++) {
                bz = b-bx-by;
                bind = bx*bxm + by*bym + bz*bzm;

                for (a=1; a<=am1; a++) {
                    // This next for loop was for (ax=0; ax<=b; ax++)
                    // this could explain why dx2 was not being computed.
                    // change for for(ax=0; ax<a; ax++) on 2005-09-15 4:11pm
                    for (ax=0; ax<=a; ax++) {
                        for (ay=0; ay<=a-ax; ay++) {
                            az = a-ax-ay;
                            aind = ax*axm + ay*aym + az*azm;

                            if (az > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    vi_[aind][bind][m] = PA[2] * vi_[aind-azm][bind][m] - PC[2] * vi_[aind-azm][bind][m+1];
                                }
                                for (m=0; m<=mmax-a-b-1; ++m) {
                                    vx_[aind][bind][m] = PA[2] * vx_[aind-azm][bind][m] - PC[2] * vx_[aind-azm][bind][m+1];
                                    vy_[aind][bind][m] = PA[2] * vy_[aind-azm][bind][m] - PC[2] * vy_[aind-azm][bind][m+1];
                                    vz_[aind][bind][m] = PA[2] * vz_[aind-azm][bind][m] - PC[2] * vz_[aind-azm][bind][m+1] + vi_[aind-azm][bind][m+1];
                                }
                                for(m=0;m<=mmax-a-b-2;m++) {  /* Gradients of the electric field */
                                    vxx_[aind][bind][m] = PA[2]*vxx_[aind-azm][bind][m] - PC[2]*vxx_[aind-azm][bind][m+1];
                                    vyy_[aind][bind][m] = PA[2]*vyy_[aind-azm][bind][m] - PC[2]*vyy_[aind-azm][bind][m+1];
                                    vzz_[aind][bind][m] = PA[2]*vzz_[aind-azm][bind][m] - PC[2]*vzz_[aind-azm][bind][m+1] + 2*vz_[aind-azm][bind][m+1];
                                    vxy_[aind][bind][m] = PA[2]*vxy_[aind-azm][bind][m] - PC[2]*vxy_[aind-azm][bind][m+1];
                                    vxz_[aind][bind][m] = PA[2]*vxz_[aind-azm][bind][m] - PC[2]*vxz_[aind-azm][bind][m+1] + vx_[aind-azm][bind][m+1];
                                    vyz_[aind][bind][m] = PA[2]*vyz_[aind-azm][bind][m] - PC[2]*vyz_[aind-azm][bind][m+1] + vy_[aind-azm][bind][m+1];
                                }

                                if (az > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * (az-1) * (vi_[aind-2*azm][bind][m] - vi_[aind-2*azm][bind][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        vx_[aind][bind][m] += ooz * (az-1) * (vx_[aind-2*azm][bind][m] - vx_[aind-2*azm][bind][m+1]);
                                        vy_[aind][bind][m] += ooz * (az-1) * (vy_[aind-2*azm][bind][m] - vy_[aind-2*azm][bind][m+1]);
                                        vz_[aind][bind][m] += ooz * (az-1) * (vz_[aind-2*azm][bind][m] - vz_[aind-2*azm][bind][m+1]);
                                    }
                                    for(m=0;m<=mmax-a-b-2;m++) {
                                        vxx_[aind][bind][m] += ooz*(az-1)*(vxx_[aind-2*azm][bind][m] - vxx_[aind-2*azm][bind][m+1]);
                                        vyy_[aind][bind][m] += ooz*(az-1)*(vyy_[aind-2*azm][bind][m] - vyy_[aind-2*azm][bind][m+1]);
                                        vzz_[aind][bind][m] += ooz*(az-1)*(vzz_[aind-2*azm][bind][m] - vzz_[aind-2*azm][bind][m+1]);
                                        vxy_[aind][bind][m] += ooz*(az-1)*(vxy_[aind-2*azm][bind][m] - vxy_[aind-2*azm][bind][m+1]);
                                        vxz_[aind][bind][m] += ooz*(az-1)*(vxz_[aind-2*azm][bind][m] - vxz_[aind-2*azm][bind][m+1]);
                                        vyz_[aind][bind][m] += ooz*(az-1)*(vyz_[aind-2*azm][bind][m] - vyz_[aind-2*azm][bind][m+1]);
                                    }
                                }
                                if (bz > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * bz * (vi_[aind-azm][bind-bzm][m] - vi_[aind-azm][bind-bzm][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        vx_[aind][bind][m] += ooz * bz * (vx_[aind-azm][bind-bzm][m] - vx_[aind-azm][bind-bzm][m+1]);
                                        vy_[aind][bind][m] += ooz * bz * (vy_[aind-azm][bind-bzm][m] - vy_[aind-azm][bind-bzm][m+1]);
                                        vz_[aind][bind][m] += ooz * bz * (vz_[aind-azm][bind-bzm][m] - vz_[aind-azm][bind-bzm][m+1]);
                                    }
                                    for(m=0;m<=mmax-a-b-2;m++) {
                                        vxx_[aind][bind][m] += ooz*bz*(vxx_[aind-azm][bind-bzm][m] - vxx_[aind-azm][bind-bzm][m+1]);
                                        vyy_[aind][bind][m] += ooz*bz*(vyy_[aind-azm][bind-bzm][m] - vyy_[aind-azm][bind-bzm][m+1]);
                                        vzz_[aind][bind][m] += ooz*bz*(vzz_[aind-azm][bind-bzm][m] - vzz_[aind-azm][bind-bzm][m+1]);
                                        vxy_[aind][bind][m] += ooz*bz*(vxy_[aind-azm][bind-bzm][m] - vxy_[aind-azm][bind-bzm][m+1]);
                                        vxz_[aind][bind][m] += ooz*bz*(vxz_[aind-azm][bind-bzm][m] - vxz_[aind-azm][bind-bzm][m+1]);
                                        vyz_[aind][bind][m] += ooz*bz*(vyz_[aind-azm][bind-bzm][m] - vyz_[aind-azm][bind-bzm][m+1]);
                                    }
                                }
                            }
                            else if (ay > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    vi_[aind][bind][m] = PA[1] * vi_[aind-aym][bind][m] - PC[1] * vi_[aind-aym][bind][m+1];
                                }
                                for (m=0; m<=mmax-a-b-1; ++m) {
                                    vx_[aind][bind][m] = PA[1] * vx_[aind-aym][bind][m] - PC[1] * vx_[aind-aym][bind][m+1];
                                    vy_[aind][bind][m] = PA[1] * vy_[aind-aym][bind][m] - PC[1] * vy_[aind-aym][bind][m+1] + vi_[aind-aym][bind][m+1];
                                    vz_[aind][bind][m] = PA[1] * vz_[aind-aym][bind][m] - PC[1] * vz_[aind-aym][bind][m+1];
                                }
                                for(m=0;m<=mmax-a-b-2;m++) {
                                    vxx_[aind][bind][m] = PA[1]*vxx_[aind-aym][bind][m] - PC[1]*vxx_[aind-aym][bind][m+1];
                                    vyy_[aind][bind][m] = PA[1]*vyy_[aind-aym][bind][m] - PC[1]*vyy_[aind-aym][bind][m+1] + 2*vy_[aind-aym][bind][m+1];
                                    vzz_[aind][bind][m] = PA[1]*vzz_[aind-aym][bind][m] - PC[1]*vzz_[aind-aym][bind][m+1];
                                    vxy_[aind][bind][m] = PA[1]*vxy_[aind-aym][bind][m] - PC[1]*vxy_[aind-aym][bind][m+1] + vx_[aind-aym][bind][m+1];
                                    vxz_[aind][bind][m] = PA[1]*vxz_[aind-aym][bind][m] - PC[1]*vxz_[aind-aym][bind][m+1];
                                    vyz_[aind][bind][m] = PA[1]*vyz_[aind-aym][bind][m] - PC[1]*vyz_[aind-aym][bind][m+1] + vz_[aind-aym][bind][m+1];
                                }
                                if (ay > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * (ay-1) * (vi_[aind-2*aym][bind][m] - vi_[aind-2*aym][bind][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        vx_[aind][bind][m] += ooz * (ay-1) * (vx_[aind-2*aym][bind][m] - vx_[aind-2*aym][bind][m+1]);
                                        vy_[aind][bind][m] += ooz * (ay-1) * (vy_[aind-2*aym][bind][m] - vy_[aind-2*aym][bind][m+1]);
                                        vz_[aind][bind][m] += ooz * (ay-1) * (vz_[aind-2*aym][bind][m] - vz_[aind-2*aym][bind][m+1]);
                                    }
                                    for(m=0;m<=mmax-a-b-2;m++) {
                                        vxx_[aind][bind][m] += ooz*(ay-1)*(vxx_[aind-2*aym][bind][m] - vxx_[aind-2*aym][bind][m+1]);
                                        vyy_[aind][bind][m] += ooz*(ay-1)*(vyy_[aind-2*aym][bind][m] - vyy_[aind-2*aym][bind][m+1]);
                                        vzz_[aind][bind][m] += ooz*(ay-1)*(vzz_[aind-2*aym][bind][m] - vzz_[aind-2*aym][bind][m+1]);
                                        vxy_[aind][bind][m] += ooz*(ay-1)*(vxy_[aind-2*aym][bind][m] - vxy_[aind-2*aym][bind][m+1]);
                                        vxz_[aind][bind][m] += ooz*(ay-1)*(vxz_[aind-2*aym][bind][m] - vxz_[aind-2*aym][bind][m+1]);
                                        vyz_[aind][bind][m] += ooz*(ay-1)*(vyz_[aind-2*aym][bind][m] - vyz_[aind-2*aym][bind][m+1]);
                                    }
                                }
                                if (by > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * by * (vi_[aind-aym][bind-bym][m] - vi_[aind-aym][bind-bym][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        vx_[aind][bind][m] += ooz * by * (vx_[aind-aym][bind-bym][m] - vx_[aind-aym][bind-bym][m+1]);
                                        vy_[aind][bind][m] += ooz * by * (vy_[aind-aym][bind-bym][m] - vy_[aind-aym][bind-bym][m+1]);
                                        vz_[aind][bind][m] += ooz * by * (vz_[aind-aym][bind-bym][m] - vz_[aind-aym][bind-bym][m+1]);
                                    }
                                    for(m=0;m<=mmax-a-b-2;m++) {
                                        vxx_[aind][bind][m] += ooz*by*(vxx_[aind-aym][bind-bym][m] - vxx_[aind-aym][bind-bym][m+1]);
                                        vyy_[aind][bind][m] += ooz*by*(vyy_[aind-aym][bind-bym][m] - vyy_[aind-aym][bind-bym][m+1]);
                                        vzz_[aind][bind][m] += ooz*by*(vzz_[aind-aym][bind-bym][m] - vzz_[aind-aym][bind-bym][m+1]);
                                        vxy_[aind][bind][m] += ooz*by*(vxy_[aind-aym][bind-bym][m] - vxy_[aind-aym][bind-bym][m+1]);
                                        vxz_[aind][bind][m] += ooz*by*(vxz_[aind-aym][bind-bym][m] - vxz_[aind-aym][bind-bym][m+1]);
                                        vyz_[aind][bind][m] += ooz*by*(vyz_[aind-aym][bind-bym][m] - vyz_[aind-aym][bind-bym][m+1]);
                                    }
                                }
                            }
                            else if (ax > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    vi_[aind][bind][m] = PA[0] * vi_[aind-axm][bind][m] - PC[0] * vi_[aind-axm][bind][m+1];
                                }
                                for (m=0; m<=mmax-a-b-1; ++m) {
                                    vx_[aind][bind][m] = PA[0] * vx_[aind-axm][bind][m] - PC[0] * vx_[aind-axm][bind][m+1] + vi_[aind-axm][bind][m+1];
                                    vy_[aind][bind][m] = PA[0] * vy_[aind-axm][bind][m] - PC[0] * vy_[aind-axm][bind][m+1];
                                    vz_[aind][bind][m] = PA[0] * vz_[aind-axm][bind][m] - PC[0] * vz_[aind-axm][bind][m+1];
                                }
                                for(m=0;m<=mmax-a-b-2;m++) {  /* Gradients of the electric field */
                                    vxx_[aind][bind][m] = PA[0]*vxx_[aind-axm][bind][m] - PC[0]*vxx_[aind-axm][bind][m+1] + 2*vx_[aind-axm][bind][m+1];
                                    vyy_[aind][bind][m] = PA[0]*vyy_[aind-axm][bind][m] - PC[0]*vyy_[aind-axm][bind][m+1];
                                    vzz_[aind][bind][m] = PA[0]*vzz_[aind-axm][bind][m] - PC[0]*vzz_[aind-axm][bind][m+1];
                                    vxy_[aind][bind][m] = PA[0]*vxy_[aind-axm][bind][m] - PC[0]*vxy_[aind-axm][bind][m+1] + vy_[aind-axm][bind][m+1];
                                    vxz_[aind][bind][m] = PA[0]*vxz_[aind-axm][bind][m] - PC[0]*vxz_[aind-axm][bind][m+1] + vz_[aind-axm][bind][m+1];
                                    vyz_[aind][bind][m] = PA[0]*vyz_[aind-axm][bind][m] - PC[0]*vyz_[aind-axm][bind][m+1];
                                }


                                if (ax > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * (ax-1) * (vi_[aind-2*axm][bind][m] - vi_[aind-2*axm][bind][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        vx_[aind][bind][m] += ooz * (ax-1) * (vx_[aind-2*axm][bind][m] - vx_[aind-2*axm][bind][m+1]);
                                        vy_[aind][bind][m] += ooz * (ax-1) * (vy_[aind-2*axm][bind][m] - vy_[aind-2*axm][bind][m+1]);
                                        vz_[aind][bind][m] += ooz * (ax-1) * (vz_[aind-2*axm][bind][m] - vz_[aind-2*axm][bind][m+1]);
                                    }
                                    for(m=0;m<=mmax-a-b-2;m++) {
                                        vxx_[aind][bind][m] += ooz*(ax-1)*(vxx_[aind-2*axm][bind][m] - vxx_[aind-2*axm][bind][m+1]);
                                        vyy_[aind][bind][m] += ooz*(ax-1)*(vyy_[aind-2*axm][bind][m] - vyy_[aind-2*axm][bind][m+1]);
                                        vzz_[aind][bind][m] += ooz*(ax-1)*(vzz_[aind-2*axm][bind][m] - vzz_[aind-2*axm][bind][m+1]);
                                        vxy_[aind][bind][m] += ooz*(ax-1)*(vxy_[aind-2*axm][bind][m] - vxy_[aind-2*axm][bind][m+1]);
                                        vxz_[aind][bind][m] += ooz*(ax-1)*(vxz_[aind-2*axm][bind][m] - vxz_[aind-2*axm][bind][m+1]);
                                        vyz_[aind][bind][m] += ooz*(ax-1)*(vyz_[aind-2*axm][bind][m] - vyz_[aind-2*axm][bind][m+1]);
                                    }
                                }
                                if (bx > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * bx * (vi_[aind-axm][bind-bxm][m] - vi_[aind-axm][bind-bxm][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        vx_[aind][bind][m] += ooz * bx * (vx_[aind-axm][bind-bxm][m] - vx_[aind-axm][bind-bxm][m+1]);
                                        vy_[aind][bind][m] += ooz * bx * (vy_[aind-axm][bind-bxm][m] - vy_[aind-axm][bind-bxm][m+1]);
                                        vz_[aind][bind][m] += ooz * bx * (vz_[aind-axm][bind-bxm][m] - vz_[aind-axm][bind-bxm][m+1]);
                                    }
                                    for(m=0;m<=mmax-a-b-2;m++) {
                                        vxx_[aind][bind][m] += ooz*bx*(vxx_[aind-axm][bind-bxm][m] - vxx_[aind-axm][bind-bxm][m+1]);
                                        vyy_[aind][bind][m] += ooz*bx*(vyy_[aind-axm][bind-bxm][m] - vyy_[aind-axm][bind-bxm][m+1]);
                                        vzz_[aind][bind][m] += ooz*bx*(vzz_[aind-axm][bind-bxm][m] - vzz_[aind-axm][bind-bxm][m+1]);
                                        vxy_[aind][bind][m] += ooz*bx*(vxy_[aind-axm][bind-bxm][m] - vxy_[aind-axm][bind-bxm][m+1]);
                                        vxz_[aind][bind][m] += ooz*bx*(vxz_[aind-axm][bind-bxm][m] - vxz_[aind-axm][bind-bxm][m+1]);
                                        vyz_[aind][bind][m] += ooz*bx*(vyz_[aind-axm][bind-bxm][m] - vyz_[aind-axm][bind-bxm][m+1]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    delete[] F;
}

ObaraSaikaTwoCenterElectricField::ObaraSaikaTwoCenterElectricField(int max_am1, int max_am2)
    : ObaraSaikaTwoCenterVIRecursion(max_am1, max_am2) // This is not used, because the elec deriv (3rd arg) needed is different
{
    q_   = init_box(size_, size_, max_am1_ + max_am2_ + 2);
    x_   = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    y_   = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    z_   = init_box(size_, size_, max_am1_ + max_am2_ + 1);
}

ObaraSaikaTwoCenterElectricField::~ObaraSaikaTwoCenterElectricField()
{
    free_box(x_, size_, size_);
    free_box(y_, size_, size_);
    free_box(z_, size_, size_);
}

void ObaraSaikaTwoCenterElectricField::compute(double PA[3], double PB[3], double PC[3], double zeta, int am1, int am2)
{
    int a, b, m;
    int azm = 1;
    int aym = am1 + 1;
    int axm = aym * aym;
    int bzm = 1;
    int bym = am2 + 1;
    int bxm = bym * bym;
    int ax, ay, az, bx, by, bz;
    int aind, bind;
    double ooz = 1.0/(2.0 * zeta);
    int mmax = am1 + am2 + 1;

    // Prefactor from A20
    double tmp = sqrt(zeta) * M_2_SQRTPI;
    // U from A21
    double u = zeta * (PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2]);
    double *F = new double[mmax+1]; // TODO: Move this allocation into constructor

    // Zero out F
    memset(F, 0, sizeof(double) * (mmax+1));

    // Form Fm(U) from A20
    calculate_f(F, mmax, u);

    // Perform recursion in m for (a|A(0)|s) using A20
    for (m=0; m<=mmax; ++m) {
        q_[0][0][m] = tmp * F[m];
    }
    for (m=0; m<=mmax-1; ++m) {
        x_[0][0][m] = 2.0*zeta*PC[0]*q_[0][0][m+1];
        y_[0][0][m] = 2.0*zeta*PC[1]*q_[0][0][m+1];
        z_[0][0][m] = 2.0*zeta*PC[2]*q_[0][0][m+1];
    }

    // Perform recursion in b with a=0
    //  subset of A19
    for (b=1; b<=am2; ++b) {
        for (bx=0; bx<=b; ++bx) {
            for (by=0; by<=b-bx; ++by) {
                bz = b-bx-by;

                // Compute the index into VI for bx,by,bz
                bind = bx*bxm + by*bym + bz*bzm;

                // Compute each x, y, z contribution
                if (bz > 0) {
                    for (m=0; m<=mmax-b; ++m) { /* Electrostatic potential integrals */
                        q_[0][bind][m] = PB[2] * q_[0][bind-bzm][m] - PC[2] * q_[0][bind-bzm][m+1];
                    }
                    for (m=0; m<=mmax-b-1; ++m) { /* Electric field integrals */
                        x_[0][bind][m] = PB[2] * x_[0][bind-bzm][m] - PC[2] * x_[0][bind-bzm][m+1];
                        y_[0][bind][m] = PB[2] * y_[0][bind-bzm][m] - PC[2] * y_[0][bind-bzm][m+1];
                        z_[0][bind][m] = PB[2] * z_[0][bind-bzm][m] - PC[2] * z_[0][bind-bzm][m+1] + q_[0][bind-bzm][m+1];
                    }
                    if (bz > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            q_[0][bind][m] += ooz * (bz-1) * (q_[0][bind-2*bzm][m] - q_[0][bind-2*bzm][m+1]);
                        }
                        for (m=0; m<=mmax-b-1; ++m) {
                            x_[0][bind][m] += ooz * (bz-1) * (x_[0][bind-2*bzm][m] - x_[0][bind-2*bzm][m+1]);
                            y_[0][bind][m] += ooz * (bz-1) * (y_[0][bind-2*bzm][m] - y_[0][bind-2*bzm][m+1]);
                            z_[0][bind][m] += ooz * (bz-1) * (z_[0][bind-2*bzm][m] - z_[0][bind-2*bzm][m+1]);
                        }
                    }
                }
                else if (by > 0) {
                    for (m=0; m<=mmax-b; ++m) {
                        q_[0][bind][m] = PB[1] * q_[0][bind-bym][m] - PC[1] * q_[0][bind-bym][m+1];
                    }
                    for (m=0; m<=mmax-b-1; ++m) {
                        x_[0][bind][m] = PB[1] * x_[0][bind-bym][m] - PC[1] * x_[0][bind-bym][m+1];
                        y_[0][bind][m] = PB[1] * y_[0][bind-bym][m] - PC[1] * y_[0][bind-bym][m+1] + q_[0][bind-bym][m+1];
                        z_[0][bind][m] = PB[1] * z_[0][bind-bym][m] - PC[1] * z_[0][bind-bym][m+1];
                    }

                    if (by > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            q_[0][bind][m] += ooz * (by-1) * (q_[0][bind-2*bym][m] - q_[0][bind-2*bym][m+1]);
                        }
                        for (m=0; m<=mmax-b-1; ++m) {
                            x_[0][bind][m] += ooz * (by-1) * (x_[0][bind-2*bym][m] - x_[0][bind-2*bym][m+1]);
                            y_[0][bind][m] += ooz * (by-1) * (y_[0][bind-2*bym][m] - y_[0][bind-2*bym][m+1]);
                            z_[0][bind][m] += ooz * (by-1) * (z_[0][bind-2*bym][m] - z_[0][bind-2*bym][m+1]);
                        }
                    }
                }
                else if (bx > 0) {
                    for (m=0; m<=mmax-b; ++m) {
                        q_[0][bind][m] = PB[0] * q_[0][bind-bxm][m] - PC[0] * q_[0][bind-bxm][m+1];
                    }
                    for (m=0; m<=mmax-b-1; ++m) {
                        x_[0][bind][m] = PB[0] * x_[0][bind-bxm][m] - PC[0] * x_[0][bind-bxm][m+1] + q_[0][bind-bxm][m+1];
                        y_[0][bind][m] = PB[0] * y_[0][bind-bxm][m] - PC[0] * y_[0][bind-bxm][m+1];
                        z_[0][bind][m] = PB[0] * z_[0][bind-bxm][m] - PC[0] * z_[0][bind-bxm][m+1];
                    }

                    if (bx > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            q_[0][bind][m] += ooz * (bx-1) * (q_[0][bind-2*bxm][m] - q_[0][bind-2*bxm][m+1]);
                        }
                        for (m=0; m<=mmax-b-1; ++m) {
                            x_[0][bind][m] += ooz * (bx-1) * (x_[0][bind-2*bxm][m] - x_[0][bind-2*bxm][m+1]);
                            y_[0][bind][m] += ooz * (bx-1) * (y_[0][bind-2*bxm][m] - y_[0][bind-2*bxm][m+1]);
                            z_[0][bind][m] += ooz * (bx-1) * (z_[0][bind-2*bxm][m] - z_[0][bind-2*bxm][m+1]);
                        }
                    }
                }
            }
        }
    }

    // Perform upward recursion in a with all b's
    for (b=0; b<=am2; b++) {
        for (bx=0; bx<=b; bx++) {
            for (by=0; by<=b-bx;by++) {
                bz = b-bx-by;
                bind = bx*bxm + by*bym + bz*bzm;

                for (a=1; a<=am1; a++) {
                    // This next for loop was for (ax=0; ax<=b; ax++)
                    // this could explain why dx2 was not being computed.
                    // change for for(ax=0; ax<a; ax++) on 2005-09-15 4:11pm
                    for (ax=0; ax<=a; ax++) {
                        for (ay=0; ay<=a-ax; ay++) {
                            az = a-ax-ay;
                            aind = ax*axm + ay*aym + az*azm;

                            if (az > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    q_[aind][bind][m] = PA[2] * q_[aind-azm][bind][m] - PC[2] * q_[aind-azm][bind][m+1];
                                }
                                for (m=0; m<=mmax-a-b-1; ++m) {
                                    x_[aind][bind][m] = PA[2] * x_[aind-azm][bind][m] - PC[2] * x_[aind-azm][bind][m+1];
                                    y_[aind][bind][m] = PA[2] * y_[aind-azm][bind][m] - PC[2] * y_[aind-azm][bind][m+1];
                                    z_[aind][bind][m] = PA[2] * z_[aind-azm][bind][m] - PC[2] * z_[aind-azm][bind][m+1] + q_[aind-azm][bind][m+1];
                                }

                                if (az > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        q_[aind][bind][m] += ooz * (az-1) * (q_[aind-2*azm][bind][m] - q_[aind-2*azm][bind][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        x_[aind][bind][m] += ooz * (az-1) * (x_[aind-2*azm][bind][m] - x_[aind-2*azm][bind][m+1]);
                                        y_[aind][bind][m] += ooz * (az-1) * (y_[aind-2*azm][bind][m] - y_[aind-2*azm][bind][m+1]);
                                        z_[aind][bind][m] += ooz * (az-1) * (z_[aind-2*azm][bind][m] - z_[aind-2*azm][bind][m+1]);
                                    }
                                }
                                if (bz > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        q_[aind][bind][m] += ooz * bz * (q_[aind-azm][bind-bzm][m] - q_[aind-azm][bind-bzm][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        x_[aind][bind][m] += ooz * bz * (x_[aind-azm][bind-bzm][m] - x_[aind-azm][bind-bzm][m+1]);
                                        y_[aind][bind][m] += ooz * bz * (y_[aind-azm][bind-bzm][m] - y_[aind-azm][bind-bzm][m+1]);
                                        z_[aind][bind][m] += ooz * bz * (z_[aind-azm][bind-bzm][m] - z_[aind-azm][bind-bzm][m+1]);
                                    }
                                }
                            }
                            else if (ay > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    q_[aind][bind][m] = PA[1] * q_[aind-aym][bind][m] - PC[1] * q_[aind-aym][bind][m+1];
                                }
                                for (m=0; m<=mmax-a-b-1; ++m) {
                                    x_[aind][bind][m] = PA[1] * x_[aind-aym][bind][m] - PC[1] * x_[aind-aym][bind][m+1];
                                    y_[aind][bind][m] = PA[1] * y_[aind-aym][bind][m] - PC[1] * y_[aind-aym][bind][m+1] + q_[aind-aym][bind][m+1];
                                    z_[aind][bind][m] = PA[1] * z_[aind-aym][bind][m] - PC[1] * z_[aind-aym][bind][m+1];
                                }
                                if (ay > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        q_[aind][bind][m] += ooz * (ay-1) * (q_[aind-2*aym][bind][m] - q_[aind-2*aym][bind][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        x_[aind][bind][m] += ooz * (ay-1) * (x_[aind-2*aym][bind][m] - x_[aind-2*aym][bind][m+1]);
                                        y_[aind][bind][m] += ooz * (ay-1) * (y_[aind-2*aym][bind][m] - y_[aind-2*aym][bind][m+1]);
                                        z_[aind][bind][m] += ooz * (ay-1) * (z_[aind-2*aym][bind][m] - z_[aind-2*aym][bind][m+1]);
                                    }
                                }
                                if (by > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        q_[aind][bind][m] += ooz * by * (q_[aind-aym][bind-bym][m] - q_[aind-aym][bind-bym][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        x_[aind][bind][m] += ooz * by * (x_[aind-aym][bind-bym][m] - x_[aind-aym][bind-bym][m+1]);
                                        y_[aind][bind][m] += ooz * by * (y_[aind-aym][bind-bym][m] - y_[aind-aym][bind-bym][m+1]);
                                        z_[aind][bind][m] += ooz * by * (z_[aind-aym][bind-bym][m] - z_[aind-aym][bind-bym][m+1]);
                                    }
                                }
                            }
                            else if (ax > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    q_[aind][bind][m] = PA[0] * q_[aind-axm][bind][m] - PC[0] * q_[aind-axm][bind][m+1];
                                }
                                for (m=0; m<=mmax-a-b-1; ++m) {
                                    x_[aind][bind][m] = PA[0] * x_[aind-axm][bind][m] - PC[0] * x_[aind-axm][bind][m+1] + q_[aind-axm][bind][m+1];
                                    y_[aind][bind][m] = PA[0] * y_[aind-axm][bind][m] - PC[0] * y_[aind-axm][bind][m+1];
                                    z_[aind][bind][m] = PA[0] * z_[aind-axm][bind][m] - PC[0] * z_[aind-axm][bind][m+1];
                                }

                                if (ax > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        q_[aind][bind][m] += ooz * (ax-1) * (q_[aind-2*axm][bind][m] - q_[aind-2*axm][bind][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        x_[aind][bind][m] += ooz * (ax-1) * (x_[aind-2*axm][bind][m] - x_[aind-2*axm][bind][m+1]);
                                        y_[aind][bind][m] += ooz * (ax-1) * (y_[aind-2*axm][bind][m] - y_[aind-2*axm][bind][m+1]);
                                        z_[aind][bind][m] += ooz * (ax-1) * (z_[aind-2*axm][bind][m] - z_[aind-2*axm][bind][m+1]);
                                    }
                                }
                                if (bx > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        q_[aind][bind][m] += ooz * bx * (q_[aind-axm][bind-bxm][m] - q_[aind-axm][bind-bxm][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        x_[aind][bind][m] += ooz * bx * (x_[aind-axm][bind-bxm][m] - x_[aind-axm][bind-bxm][m+1]);
                                        y_[aind][bind][m] += ooz * bx * (y_[aind-axm][bind-bxm][m] - y_[aind-axm][bind-bxm][m+1]);
                                        z_[aind][bind][m] += ooz * bx * (z_[aind-axm][bind-bxm][m] - z_[aind-axm][bind-bxm][m+1]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    delete[] F;
}

ObaraSaikaTwoCenterElectricFieldGradient::ObaraSaikaTwoCenterElectricFieldGradient(int max_am1, int max_am2)
    : ObaraSaikaTwoCenterElectricField(max_am1+1, max_am2+1)
{
    // Check memory allocation for 3rd parameter. I think it should be max_am1_ + max_am2_ - 2
    exx_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    eyy_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    ezz_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    exy_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    exz_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    eyz_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
}

ObaraSaikaTwoCenterElectricFieldGradient::~ObaraSaikaTwoCenterElectricFieldGradient()
{
    free_box(exx_, size_, size_);
    free_box(eyy_, size_, size_);
    free_box(ezz_, size_, size_);
    free_box(exy_, size_, size_);
    free_box(exz_, size_, size_);
    free_box(eyz_, size_, size_);
}

void ObaraSaikaTwoCenterElectricFieldGradient::compute(double PA[3], double PB[3], double PC[3], double zeta, int am1, int am2)
{
    int a, b, m;
    int azm = 1;
    int aym = am1 + 1;
    int axm = aym * aym;
    int bzm = 1;
    int bym = am2 + 1;
    int bxm = bym * bym;
    int ax, ay, az, bx, by, bz;
    int aind, bind;
    double ooz = 1.0/(2.0 * zeta);
    int mmax = am1 + am2;

    // Call our super class to compute electric field and potential integrals
    ObaraSaikaTwoCenterElectricField::compute(PA, PB, PC, zeta, am1+1, am2+1);

    // Compute starting integrals using A26
    for (m=0; m<=mmax-2; ++m) {
        exx_[0][0][m] = 4.0 * zeta * zeta * PC[0] * PC[0] * vi_[0][0][m+2] - 2.0 * zeta * vi_[0][0][m+1];
        eyy_[0][0][m] = 4.0 * zeta * zeta * PC[1] * PC[1] * vi_[0][0][m+2] - 2.0 * zeta * vi_[0][0][m+1];
        ezz_[0][0][m] = 4.0 * zeta * zeta * PC[2] * PC[2] * vi_[0][0][m+2] - 2.0 * zeta * vi_[0][0][m+1];
        exy_[0][0][m] = 4.0 * zeta * zeta * PC[0] * PC[1] * vi_[0][0][m+2];
        exz_[0][0][m] = 4.0 * zeta * zeta * PC[0] * PC[2] * vi_[0][0][m+2];
        eyz_[0][0][m] = 4.0 * zeta * zeta * PC[1] * PC[2] * vi_[0][0][m+2];
    }

    // Perform recursion in b with a=0
    //    subset of A24
    for (b=1; b<=am2; ++b) {
        for (bx=0; bx<=b; ++bx) {
            for (by=0; by<=b-bx; ++by) {
                bz = b-bx-by;

                // Compute the index into exx_, eyy_, ezz_, exy_, exz_, eyz
                bind = bx*bxm + by*bym + bz*bzm;

                // Compute each x, y, z contribution
                if (bz > 0) {
                    for (m=0; m<=mmax-b-2; ++m) {
                        exx_[0][bind][m] = PB[2] * exx_[0][bind-bzm][m] - PC[2] * exx_[0][bind-bzm][m+1];
                        eyy_[0][bind][m] = PB[2] * eyy_[0][bind-bzm][m] - PC[2] * eyy_[0][bind-bzm][m+1];
                        ezz_[0][bind][m] = PB[2] * ezz_[0][bind-bzm][m] - PC[2] * ezz_[0][bind-bzm][m+1] + 2.0 * z_[0][bind-bzm][m+1];
                        exy_[0][bind][m] = PB[2] * exy_[0][bind-bzm][m] - PC[2] * exy_[0][bind-bzm][m+1];
                        exz_[0][bind][m] = PB[2] * exz_[0][bind-bzm][m] - PC[2] * exz_[0][bind-bzm][m+1] + x_[0][bind-bzm][m+1];
                        eyz_[0][bind][m] = PB[2] * eyz_[0][bind-bzm][m] - PC[2] * eyz_[0][bind-bzm][m+1] + y_[0][bind-bzm][m+1];
                    }
                    if (bz > 1) {
                        for (m=0; m<=mmax-b-2; ++m) {
                            exx_[0][bind][m] += ooz * (bz-1) * (exx_[0][bind-2*bzm][m] - exx_[0][bind-2*bzm][m+1]);
                            eyy_[0][bind][m] += ooz * (bz-1) * (eyy_[0][bind-2*bzm][m] - eyy_[0][bind-2*bzm][m+1]);
                            ezz_[0][bind][m] += ooz * (bz-1) * (ezz_[0][bind-2*bzm][m] - ezz_[0][bind-2*bzm][m+1]);
                            exy_[0][bind][m] += ooz * (bz-1) * (exy_[0][bind-2*bzm][m] - exy_[0][bind-2*bzm][m+1]);
                            exz_[0][bind][m] += ooz * (bz-1) * (exz_[0][bind-2*bzm][m] - exz_[0][bind-2*bzm][m+1]);
                            eyz_[0][bind][m] += ooz * (bz-1) * (eyz_[0][bind-2*bzm][m] - eyz_[0][bind-2*bzm][m+1]);
                        }
                    }
                }
                else if (by > 0) {
                    for (m=0; m<=mmax-b-2; ++m) {
                        exx_[0][bind][m] = PB[1] * exx_[0][bind-bym][m] - PC[1] * exx_[0][bind-bym][m+1];
                        eyy_[0][bind][m] = PB[1] * eyy_[0][bind-bym][m] - PC[1] * eyy_[0][bind-bym][m+1] + 2.0 * y_[0][bind-bym][m+1];
                        ezz_[0][bind][m] = PB[1] * ezz_[0][bind-bym][m] - PC[1] * ezz_[0][bind-bym][m+1];
                        exy_[0][bind][m] = PB[1] * exy_[0][bind-bym][m] - PC[1] * exy_[0][bind-bym][m+1] - x_[0][bind-bym][m+1];
                        exz_[0][bind][m] = PB[1] * exz_[0][bind-bym][m] - PC[1] * exz_[0][bind-bym][m+1];
                        eyz_[0][bind][m] = PB[1] * eyz_[0][bind-bym][m] - PC[1] * eyz_[0][bind-bym][m+1] - z_[0][bind-bym][m+1];
                    }
                    if (by > 1) {
                        for (m=0; m<=mmax-b-2; ++m) {
                            exx_[0][bind][m] += ooz * (by-1) * (exx_[0][bind-2*bym][m] - exx_[0][bind-2*bym][m+1]);
                            eyy_[0][bind][m] += ooz * (by-1) * (eyy_[0][bind-2*bym][m] - eyy_[0][bind-2*bym][m+1]);
                            ezz_[0][bind][m] += ooz * (by-1) * (ezz_[0][bind-2*bym][m] - ezz_[0][bind-2*bym][m+1]);
                            exy_[0][bind][m] += ooz * (by-1) * (exy_[0][bind-2*bym][m] - exy_[0][bind-2*bym][m+1]);
                            exz_[0][bind][m] += ooz * (by-1) * (exz_[0][bind-2*bym][m] - exz_[0][bind-2*bym][m+1]);
                            eyz_[0][bind][m] += ooz * (by-1) * (eyz_[0][bind-2*bym][m] - eyz_[0][bind-2*bym][m+1]);
                        }
                    }
                }
                else if (bx > 0) {
                    for (m=0; m<=mmax-b-2; ++m) {
                        exx_[0][bind][m] = PB[0] * exx_[0][bind-bxm][m] - PC[0] * exx_[0][bind-bxm][m+1] + 2.0 * x_[0][bind-bxm][m+1];
                        eyy_[0][bind][m] = PB[0] * eyy_[0][bind-bxm][m] - PC[0] * eyy_[0][bind-bxm][m+1];
                        ezz_[0][bind][m] = PB[0] * ezz_[0][bind-bxm][m] - PC[0] * ezz_[0][bind-bxm][m+1];
                        exy_[0][bind][m] = PB[0] * exy_[0][bind-bxm][m] - PC[0] * exy_[0][bind-bxm][m+1] - y_[0][bind-bxm][m+1];
                        exz_[0][bind][m] = PB[0] * exz_[0][bind-bxm][m] - PC[0] * exz_[0][bind-bxm][m+1] - z_[0][bind-bxm][m+1];
                        eyz_[0][bind][m] = PB[0] * eyz_[0][bind-bxm][m] - PC[0] * eyz_[0][bind-bxm][m+1];
                    }
                    if (bx > 1) {
                        for (m=0; m<=mmax-b-2; ++m) {
                            exx_[0][bind][m] += ooz * (bx-1) * (exx_[0][bind-2*bxm][m] - exx_[0][bind-2*bxm][m+1]);
                            eyy_[0][bind][m] += ooz * (bx-1) * (eyy_[0][bind-2*bxm][m] - eyy_[0][bind-2*bxm][m+1]);
                            ezz_[0][bind][m] += ooz * (bx-1) * (ezz_[0][bind-2*bxm][m] - ezz_[0][bind-2*bxm][m+1]);
                            exy_[0][bind][m] += ooz * (bx-1) * (exy_[0][bind-2*bxm][m] - exy_[0][bind-2*bxm][m+1]);
                            exz_[0][bind][m] += ooz * (bx-1) * (exz_[0][bind-2*bxm][m] - exz_[0][bind-2*bxm][m+1]);
                            eyz_[0][bind][m] += ooz * (bx-1) * (eyz_[0][bind-2*bxm][m] - eyz_[0][bind-2*bxm][m+1]);
                        }
                    }
                }
            }
        }
    }

    // Perform upward recursion in a with all b's
    for (b=0; b<=am2; ++b) {
        for (bx=0; bx<=b; bx++) {
            for (by=0; by<=b-bx;by++) {
                bz = b-bx-by;
                bind = bx*bxm + by*bym + bz*bzm;

                for (a=1; a<=am1; a++) {
                    for (ax=0; ax<=a; ax++) {
                        for (ay=0; ay<=a-ax; ay++) {
                            az = a-ax-ay;
                            aind = ax*axm + ay*aym + az*azm;

                            if (az > 0) {
                                for (m=0; m<=mmax-a-b-2; ++m) {
                                    exx_[aind][bind][m] = PA[2] * exx_[aind-azm][bind][m] - PC[2] * exx_[aind-azm][bind][m+1];
                                    eyy_[aind][bind][m] = PA[2] * eyy_[aind-azm][bind][m] - PC[2] * eyy_[aind-azm][bind][m+1];
                                    ezz_[aind][bind][m] = PA[2] * ezz_[aind-azm][bind][m] - PC[2] * ezz_[aind-azm][bind][m+1] - 2.0 * z_[aind-azm][bind][m+1];
                                    exy_[aind][bind][m] = PA[2] * exy_[aind-azm][bind][m] - PC[2] * exy_[aind-azm][bind][m+1];
                                    exz_[aind][bind][m] = PA[2] * exz_[aind-azm][bind][m] - PC[2] * exz_[aind-azm][bind][m+1] + x_[aind-azm][bind][m+1];
                                    eyz_[aind][bind][m] = PA[2] * eyz_[aind-azm][bind][m] - PC[2] * eyz_[aind-azm][bind][m+1] + y_[aind-azm][bind][m+1];
                                }
                                if (az > 1) {
                                    for (m=0; m<=mmax-a-b-2; ++m) {
                                        exx_[aind][bind][m] += ooz * (az-1) * (exx_[aind-2*azm][bind][m] - exx_[aind-2*azm][bind][m+1]);
                                        eyy_[aind][bind][m] += ooz * (az-1) * (eyy_[aind-2*azm][bind][m] - eyy_[aind-2*azm][bind][m+1]);
                                        ezz_[aind][bind][m] += ooz * (az-1) * (ezz_[aind-2*azm][bind][m] - ezz_[aind-2*azm][bind][m+1]);
                                        exy_[aind][bind][m] += ooz * (az-1) * (exy_[aind-2*azm][bind][m] - exy_[aind-2*azm][bind][m+1]);
                                        exz_[aind][bind][m] += ooz * (az-1) * (exz_[aind-2*azm][bind][m] - exz_[aind-2*azm][bind][m+1]);
                                        eyz_[aind][bind][m] += ooz * (az-1) * (eyz_[aind-2*azm][bind][m] - eyz_[aind-2*azm][bind][m+1]);
                                    }
                                }
                                if (bz > 0) {
                                    for (m=0; m<=mmax-a-b-2; ++m) {
                                        exx_[aind][bind][m] += ooz * bz * (exx_[aind-azm][bind-bzm][m] - exx_[aind-azm][bind-bzm][m+1]);
                                        eyy_[aind][bind][m] += ooz * bz * (eyy_[aind-azm][bind-bzm][m] - eyy_[aind-azm][bind-bzm][m+1]);
                                        ezz_[aind][bind][m] += ooz * bz * (ezz_[aind-azm][bind-bzm][m] - ezz_[aind-azm][bind-bzm][m+1]);
                                        exy_[aind][bind][m] += ooz * bz * (exy_[aind-azm][bind-bzm][m] - exy_[aind-azm][bind-bzm][m+1]);
                                        exz_[aind][bind][m] += ooz * bz * (exz_[aind-azm][bind-bzm][m] - exz_[aind-azm][bind-bzm][m+1]);
                                        eyz_[aind][bind][m] += ooz * bz * (eyz_[aind-azm][bind-bzm][m] - eyz_[aind-azm][bind-bzm][m+1]);
                                    }
                                }
                            }
                            else if (ay > 0) {
                                for (m=0; m<=mmax-a-b-2; ++m) {
                                    exx_[aind][bind][m] = PA[1] * exx_[aind-aym][bind][m] - PC[1] * exx_[aind-aym][bind][m+1];
                                    eyy_[aind][bind][m] = PA[1] * eyy_[aind-aym][bind][m] - PC[1] * eyy_[aind-aym][bind][m+1] - 2.0 * y_[aind-aym][bind][m+1];
                                    ezz_[aind][bind][m] = PA[1] * ezz_[aind-aym][bind][m] - PC[1] * ezz_[aind-aym][bind][m+1];
                                    exy_[aind][bind][m] = PA[1] * exy_[aind-aym][bind][m] - PC[1] * exy_[aind-aym][bind][m+1] + x_[aind-aym][bind][m+1];
                                    exz_[aind][bind][m] = PA[1] * exz_[aind-aym][bind][m] - PC[1] * exz_[aind-aym][bind][m+1];
                                    eyz_[aind][bind][m] = PA[1] * eyz_[aind-aym][bind][m] - PC[1] * eyz_[aind-aym][bind][m+1] + z_[aind-aym][bind][m+1];
                                }
                                if (ay > 1) {
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        exx_[aind][bind][m] += ooz * (ay-1) * (exx_[aind-2*aym][bind][m] - exx_[aind-2*aym][bind][m+1]);
                                        eyy_[aind][bind][m] += ooz * (ay-1) * (eyy_[aind-2*aym][bind][m] - eyy_[aind-2*aym][bind][m+1]);
                                        ezz_[aind][bind][m] += ooz * (ay-1) * (ezz_[aind-2*aym][bind][m] - ezz_[aind-2*aym][bind][m+1]);
                                        exy_[aind][bind][m] += ooz * (ay-1) * (exy_[aind-2*aym][bind][m] - exy_[aind-2*aym][bind][m+1]);
                                        exz_[aind][bind][m] += ooz * (ay-1) * (exz_[aind-2*aym][bind][m] - exz_[aind-2*aym][bind][m+1]);
                                        eyz_[aind][bind][m] += ooz * (ay-1) * (eyz_[aind-2*aym][bind][m] - eyz_[aind-2*aym][bind][m+1]);
                                    }
                                }
                                if (by > 0) {
                                    for (m=0; m<=mmax-a-b-2; ++m) {
                                        exx_[aind][bind][m] += ooz * by * (exx_[aind-aym][bind-bym][m] - exx_[aind-aym][bind-bym][m+1]);
                                        eyy_[aind][bind][m] += ooz * by * (eyy_[aind-aym][bind-bym][m] - eyy_[aind-aym][bind-bym][m+1]);
                                        ezz_[aind][bind][m] += ooz * by * (ezz_[aind-aym][bind-bym][m] - ezz_[aind-aym][bind-bym][m+1]);
                                        exy_[aind][bind][m] += ooz * by * (exy_[aind-aym][bind-bym][m] - exy_[aind-aym][bind-bym][m+1]);
                                        exz_[aind][bind][m] += ooz * by * (exz_[aind-aym][bind-bym][m] - exz_[aind-aym][bind-bym][m+1]);
                                        eyz_[aind][bind][m] += ooz * by * (eyz_[aind-aym][bind-bym][m] - eyz_[aind-aym][bind-bym][m+1]);
                                    }
                                }
                            }
                            else if (ax > 0) {
                                for (m=0; m<=mmax-a-b-2; ++m) {
                                    exx_[aind][bind][m] = PA[0] * exx_[aind-axm][bind][m] - PC[0] * exx_[aind-axm][bind][m+1] - 2.0 * x_[aind-axm][bind][m+1];
                                    eyy_[aind][bind][m] = PA[0] * eyy_[aind-axm][bind][m] - PC[0] * eyy_[aind-axm][bind][m+1];
                                    ezz_[aind][bind][m] = PA[0] * ezz_[aind-axm][bind][m] - PC[0] * ezz_[aind-axm][bind][m+1];
                                    exy_[aind][bind][m] = PA[0] * exy_[aind-axm][bind][m] - PC[0] * exy_[aind-axm][bind][m+1] + y_[aind-axm][bind][m+1];
                                    exz_[aind][bind][m] = PA[0] * exz_[aind-axm][bind][m] - PC[0] * exz_[aind-axm][bind][m+1] + z_[aind-axm][bind][m+1];
                                    eyz_[aind][bind][m] = PA[0] * eyz_[aind-axm][bind][m] - PC[0] * eyz_[aind-axm][bind][m+1];
                                }
                                if (ax > 1) {
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        exx_[aind][bind][m] += ooz * (ax-1) * (exx_[aind-2*axm][bind][m] - exx_[aind-2*axm][bind][m+1]);
                                        eyy_[aind][bind][m] += ooz * (ax-1) * (eyy_[aind-2*axm][bind][m] - eyy_[aind-2*axm][bind][m+1]);
                                        ezz_[aind][bind][m] += ooz * (ax-1) * (ezz_[aind-2*axm][bind][m] - ezz_[aind-2*axm][bind][m+1]);
                                        exy_[aind][bind][m] += ooz * (ax-1) * (exy_[aind-2*axm][bind][m] - exy_[aind-2*axm][bind][m+1]);
                                        exz_[aind][bind][m] += ooz * (ax-1) * (exz_[aind-2*axm][bind][m] - exz_[aind-2*axm][bind][m+1]);
                                        eyz_[aind][bind][m] += ooz * (ax-1) * (eyz_[aind-2*axm][bind][m] - eyz_[aind-2*axm][bind][m+1]);
                                    }
                                }
                                if (bx > 0) {
                                    for (m=0; m<=mmax-a-b-2; ++m) {
                                        exx_[aind][bind][m] += ooz * bx * (exx_[aind-axm][bind-bxm][m] - exx_[aind-axm][bind-bxm][m+1]);
                                        eyy_[aind][bind][m] += ooz * bx * (eyy_[aind-axm][bind-bxm][m] - eyy_[aind-axm][bind-bxm][m+1]);
                                        ezz_[aind][bind][m] += ooz * bx * (ezz_[aind-axm][bind-bxm][m] - ezz_[aind-axm][bind-bxm][m+1]);
                                        exy_[aind][bind][m] += ooz * bx * (exy_[aind-axm][bind-bxm][m] - exy_[aind-axm][bind-bxm][m+1]);
                                        exz_[aind][bind][m] += ooz * bx * (exz_[aind-axm][bind-bxm][m] - exz_[aind-axm][bind-bxm][m+1]);
                                        eyz_[aind][bind][m] += ooz * bx * (eyz_[aind-axm][bind-bxm][m] - eyz_[aind-axm][bind-bxm][m+1]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/***********************************************************
 *                                                         *
 * ObaraSaikaThreeCenterRecurion                           *
 *                                                         *
 **********************************************************/
ObaraSaikaThreeCenterRecursion::ObaraSaikaThreeCenterRecursion(int max_am1, int max_am2, int max_am3)
    : max_am1_(max_am1), max_am2_(max_am2), max_am3_(max_am3)
{
    if (max_am1 < 0)
        throw SanityCheckError("ERROR: ObaraSaikaThreeCenterRecursion -- max_am1 must be nonnegative", __FILE__, __LINE__);
    if (max_am2 < 0)
        throw SanityCheckError("ERROR: ObaraSaikaThreeCenterRecursion -- max_am2 must be nonnegative", __FILE__, __LINE__);
    if (max_am3 < 0)
        throw SanityCheckError("ERROR: ObaraSaikaThreeCenterRecursion -- max_am3 must be nonnegative", __FILE__, __LINE__);

    x_ = init_box(max_am1+1, max_am3+1, max_am2+1);
    y_ = init_box(max_am1+1, max_am3+1, max_am2+1);
    z_ = init_box(max_am1+1, max_am3+1, max_am2+1);
}

ObaraSaikaThreeCenterRecursion::~ObaraSaikaThreeCenterRecursion()
{
    free_box(x_, max_am1_+1, max_am3_+1);
    free_box(y_, max_am1_+1, max_am3_+1);
    free_box(z_, max_am1_+1, max_am3_+1);
}

void ObaraSaikaThreeCenterRecursion::compute(double GA[3], double GB[3], double GC[3], double gamma, int amA, int amB, int amC)
{
    if (amA < 0 || amA > max_am1_)
        throw SanityCheckError("ERROR: ObaraSaikaThreeCenterRecursion::compute -- am1 out of bounds", __FILE__, __LINE__);
    if (amB < 0 || amB > max_am2_)
        throw SanityCheckError("ERROR: ObaraSaikaThreeCenterRecursion::compute -- am2 out of bounds", __FILE__, __LINE__);
    if (amC < 0 || amC > max_am3_)
        throw SanityCheckError("ERROR: ObaraSaikaThreeCenterRecursion::compute -- am3 out of bounds", __FILE__, __LINE__);

    int a,b,c;
    double pp = 1.0/(2.0*gamma);

    // Starting point
    x_[0][0][0] = y_[0][0][0] = z_[0][0][0] = 1.0;

    // 1 up from the bottom.
    if (amA >= 1) {
        x_[1][0][0] = GA[0];
        y_[1][0][0] = GA[1];
        z_[1][0][0] = GA[2];
    }
    if (amB >= 1) {
        x_[0][0][1] = GB[0];
        y_[0][0][1] = GB[1];
        z_[0][0][1] = GB[2];
    }
    if (amC >= 1) {
        x_[0][1][0] = GC[0];
        y_[0][1][0] = GC[1];
        z_[0][1][0] = GC[2];
    }

    // Begin - Upward recursion in b for a=c=0
    for (b=1; b<amB; ++b) {
        x_[0][0][b+1] = GB[0] * x_[0][0][b];
        y_[0][0][b+1] = GB[1] * y_[0][0][b];
        z_[0][0][b+1] = GB[2] * z_[0][0][b];
        x_[0][0][b+1] += b * pp * x_[0][0][b-1];
        y_[0][0][b+1] += b * pp * y_[0][0][b-1];
        z_[0][0][b+1] += b * pp * z_[0][0][b-1];
    }
    // End - Upward recursion in b for a=c=0

    // Begin - Upward recursion in c for all b and a=0
    if (amC) {
        for (b=1; b<=amB; ++b) {
            x_[0][1][b] = GC[0] * x_[0][0][b];
            y_[0][1][b] = GC[1] * y_[0][0][b];
            z_[0][1][b] = GC[2] * z_[0][0][b];
            x_[0][1][b] += b * pp * x_[0][0][b-1];
            y_[0][1][b] += b * pp * y_[0][0][b-1];
            z_[0][1][b] += b * pp * z_[0][0][b-1];
        }
        for (c=1; c<amC; ++c) {
            x_[0][c+1][0] = GC[0] * x_[0][c][0];
            y_[0][c+1][0] = GC[1] * y_[0][c][0];
            z_[0][c+1][0] = GC[2] * z_[0][c][0];
            x_[0][c+1][0] += c * pp * x_[0][c-1][0];
            y_[0][c+1][0] += c * pp * y_[0][c-1][0];
            z_[0][c+1][0] += c * pp * z_[0][c-1][0];

            for (b=1; b<=amB; ++b) {
                x_[0][c+1][b] = GC[0] * x_[0][c][b];
                y_[0][c+1][b] = GC[1] * y_[0][c][b];
                z_[0][c+1][b] = GC[2] * z_[0][c][b];
                x_[0][c+1][b] += c * pp * x_[0][c-1][b];
                y_[0][c+1][b] += c * pp * y_[0][c-1][b];
                z_[0][c+1][b] += c * pp * z_[0][c-1][b];
                x_[0][c+1][b] += b * pp * x_[0][c][b-1];
                y_[0][c+1][b] += b * pp * y_[0][c][b-1];
                z_[0][c+1][b] += b * pp * z_[0][c][b-1];
            }
        }
    }
    // End - Upward recursion in c for all b and a = 0

    // Begin - Upward recursion in a for all b and c
    if (amA) {
        for (b=1; b<=amB; ++b) {
            x_[1][0][b] = GA[0] * x_[0][0][b];
            y_[1][0][b] = GA[1] * y_[0][0][b];
            z_[1][0][b] = GA[2] * z_[0][0][b];
            x_[1][0][b] += b * pp * x_[0][0][b-1];
            y_[1][0][b] += b * pp * y_[0][0][b-1];
            z_[1][0][b] += b * pp * z_[0][0][b-1];
        }
        for (c=1; c<=amC; ++c) {
            x_[1][c][0] = GA[0] * x_[0][c][0];
            y_[1][c][0] = GA[1] * y_[0][c][0];
            z_[1][c][0] = GA[2] * z_[0][c][0];
            x_[1][c][0] += c * pp * x_[0][c-1][0];
            y_[1][c][0] += c * pp * y_[0][c-1][0];
            z_[1][c][0] += c * pp * z_[0][c-1][0];
        }
        for (a=1; a<amA; ++a) {
            x_[a+1][0][0] = GA[0] * x_[a][0][0];
            y_[a+1][0][0] = GA[1] * y_[a][0][0];
            z_[a+1][0][0] = GA[2] * z_[a][0][0];
            x_[a+1][0][0] += a * pp * x_[a-1][0][0];
            y_[a+1][0][0] += a * pp * y_[a-1][0][0];
            z_[a+1][0][0] += a * pp * z_[a-1][0][0];

            for (b=1; b<=amB; ++b) {
                x_[a+1][0][b] = GA[0] * x_[a][0][b];
                y_[a+1][0][b] = GA[1] * y_[a][0][b];
                z_[a+1][0][b] = GA[2] * z_[a][0][b];
                x_[a+1][0][b] += a * pp * x_[a-1][0][b];
                y_[a+1][0][b] += a * pp * y_[a-1][0][b];
                z_[a+1][0][b] += a * pp * z_[a-1][0][b];
                x_[a+1][0][b] += b * pp * x_[a][0][b-1];
                y_[a+1][0][b] += b * pp * y_[a][0][b-1];
                z_[a+1][0][b] += b * pp * z_[a][0][b-1];
            }

            for (c=1; c<=amC; ++c) {
                x_[a+1][c][0] = GA[0] * x_[a][c][0];
                y_[a+1][c][0] = GA[1] * y_[a][c][0];
                z_[a+1][c][0] = GA[2] * z_[a][c][0];
                x_[a+1][c][0] += a * pp * x_[a-1][c][0];
                y_[a+1][c][0] += a * pp * y_[a-1][c][0];
                z_[a+1][c][0] += a * pp * z_[a-1][c][0];
                x_[a+1][c][0] += c * pp * x_[a][c-1][0];
                y_[a+1][c][0] += c * pp * y_[a][c-1][0];
                z_[a+1][c][0] += c * pp * z_[a][c-1][0];
            }
        }
    }
    // End - Upward recursion in a for all b and c = 0

    if (amA && amB && amC) {
        for (b=1; b<=amB; ++b) {
            for (c=1; c<=amC; ++c) {
                x_[1][c][b] = GA[0] * x_[0][c][b];
                y_[1][c][b] = GA[1] * y_[0][c][b];
                z_[1][c][b] = GA[2] * z_[0][c][b];
                x_[1][c][b] += b * pp * x_[0][c][b-1];
                y_[1][c][b] += b * pp * y_[0][c][b-1];
                z_[1][c][b] += b * pp * z_[0][c][b-1];
                x_[1][c][b] += c * pp * x_[0][c-1][b];
                y_[1][c][b] += c * pp * y_[0][c-1][b];
                z_[1][c][b] += c * pp * z_[0][c-1][b];
            }
        }

        // Begin - Bring everything together
        for (a=1; a<amA; ++a) {
            for (b=1; b<=amB; ++b) {
                for (c=1; c<=amC; ++c) {
                    x_[a+1][c][b] = GA[0] * x_[a][c][b];
                    y_[a+1][c][b] = GA[1] * y_[a][c][b];
                    z_[a+1][c][b] = GA[2] * z_[a][c][b];
                    x_[a+1][c][b] += a * pp * x_[a-1][c][b];
                    y_[a+1][c][b] += a * pp * y_[a-1][c][b];
                    z_[a+1][c][b] += a * pp * z_[a-1][c][b];
                    x_[a+1][c][b] += b * pp * x_[a][c][b-1];
                    y_[a+1][c][b] += b * pp * y_[a][c][b-1];
                    z_[a+1][c][b] += b * pp * z_[a][c][b-1];
                    x_[a+1][c][b] += c * pp * x_[a][c-1][b];
                    y_[a+1][c][b] += c * pp * y_[a][c-1][b];
                    z_[a+1][c][b] += c * pp * z_[a][c-1][b];
                }
            }
        }
    }
}
