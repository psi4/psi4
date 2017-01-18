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

#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/cdsalclist.h"
#include "psi4/libmints/rel_potential.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/osrecur.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/physconst.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define RELVDEBUG 0

;
using namespace psi;

// Initialize potential_recur_ to +1 basis set angular momentum
RelPotentialInt::RelPotentialInt(std::vector<SphericalTransform>& st, std::shared_ptr<BasisSet> bs1,
                                 std::shared_ptr<BasisSet> bs2, int deriv) :
        OneBodyAOInt(st, bs1, bs2, deriv)
{
    if (deriv == 0)
        potential_recur_ = new ObaraSaikaTwoCenterVIRecursion(bs1->max_am() + 2, bs2->max_am() + 2);
    else
        throw PSIEXCEPTION("RelPotentialInt: deriv > 0 is not supported.");

    const int maxam1 = bs1_->max_am();
    const int maxam2 = bs2_->max_am();

    int maxnao1 = INT_NCART(maxam1);
    int maxnao2 = INT_NCART(maxam2);

    buffer_ = new double[maxnao1 * maxnao2];

    // Setup the initial field of partial charges
    Zxyz_ = SharedMatrix(new Matrix("Partial Charge Field (Z,x,y,z)", bs1_->molecule()->natom(), 4));
    double **Zxyzp = Zxyz_->pointer();

    for (int A = 0; A < bs1_->molecule()->natom(); A++) {
        Zxyzp[A][0] = (double) bs1_->molecule()->Z(A);
        Zxyzp[A][1] = bs1_->molecule()->x(A);
        Zxyzp[A][2] = bs1_->molecule()->y(A);
        Zxyzp[A][3] = bs1_->molecule()->z(A);
    }
}

RelPotentialInt::~RelPotentialInt()
{
    delete[] buffer_;
    delete potential_recur_;
}

/* Prakash
   changing compute_pair to do <p mu |1/r-C | p nu >
*/

// The engine only supports segmented basis sets
void RelPotentialInt::compute_pair(const GaussianShell& s1,
                                   const GaussianShell& s2)
{
    int ao12;
    int am1 = s1.am();
    int am2 = s2.am();
    int nprim1 = s1.nprimitive();
    int nprim2 = s2.nprimitive();
    double A[3], B[3];
    double Ixx, Iyy, Izz, Itotal;

    A[0] = s1.center()[0];
    A[1] = s1.center()[1];
    A[2] = s1.center()[2];
    B[0] = s2.center()[0];
    B[1] = s2.center()[1];
    B[2] = s2.center()[2];

    int izm = 1;
    int iym = am1 + 2;
    int ixm = iym * iym;
    int jzm = 1;
    int jym = am2 + 2;
    int jxm = jym * jym;

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    memset(buffer_, 0, s1.ncartesian() * s2.ncartesian() * sizeof(double));

    double ***vi = potential_recur_->vi();

    double **Zxyzp = Zxyz_->pointer();
    int ncharge = Zxyz_->rowspi()[0];

#if RELVDEBUG
    outfile->Printf("\n s1= %d, s2= %d,  \n\n",nprim1,nprim2);
    outfile->Printf("\n am1= %d, am2= %d,  \n\n",am1,am2);
    outfile->Printf("izm= %d, jzm= %d\n", izm, jzm);
    outfile->Printf("iym= %d, jym= %d\n", am1+1, am2+1);
    outfile->Printf("ixm= %d, jxm= %d\n", (am1+1)*(am1+1), (am2+1)*(am2+1));
#endif

    for (int p1 = 0; p1 < nprim1; ++p1) {
        double a1 = s1.exp(p1);
#if RELVDEBUG
        outfile->Printf("alpha_a=%20.14f\n",a1);
#endif
        double c1 = s1.coef(p1);
        for (int p2 = 0; p2 < nprim2; ++p2) {
            double a2 = s2.exp(p2);
#if RELVDEBUG
            outfile->Printf("alpha_b=%20.14f\n",a2);
#endif
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

            // Loop over atoms of basis set 1 (only works if bs1_ and bs2_ are on the same
            // molecule)
            for (int atom = 0; atom < ncharge; ++atom) {
#if RELVDEBUG
                outfile->Printf("\n atoms=%d \n",atom);
#endif
                double PC[3];

                double Z = Zxyzp[atom][0];

                PC[0] = P[0] - Zxyzp[atom][1];
                PC[1] = P[1] - Zxyzp[atom][2];
                PC[2] = P[2] - Zxyzp[atom][3];

                // Do recursion and ask to compute one higher angular momentum
                potential_recur_->compute(PA, PB, PC, gamma, am1 + 1, am2 + 1);

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

                                int iind = l1 * ixm + m1 * iym + n1 * izm;
                                int jind = l2 * jxm + m2 * jym + n2 * jzm;

                                // <a+1|V|b+1> terms
                                Ixx = 4.0 * a1 * a2 * vi[iind + ixm][jind + jxm][0];
                                Iyy = 4.0 * a1 * a2 * vi[iind + iym][jind + jym][0];
                                Izz = 4.0 * a1 * a2 * vi[iind + izm][jind + jzm][0];

#if RELVDEBUG
                                outfile->Printf("l1= %d, m1= %d, n1= %d  \n",l1,m1,n1);
                                outfile->Printf("l2= %d, m2= %d, n2= %d \n",l2,m2,n2);
                                outfile->Printf("iind =%d, jind= %d \n", iind,jind);
                                outfile->Printf("Siind+ixm= %d, jind+jxm= %d \n", iind+ixm, jind+jxm);
                                outfile->Printf("Siind+iym= %d, jind+jym= %d \n", iind+iym, jind+jym);
                                outfile->Printf("Siind+izm= %d, jind+jzm= %d \n", iind+izm, jind+jzm);
#endif

                                // <a+1|V|b-1>, <a-1|V|b+1>, <a-1|V|b-1> terms
                                if (l1 && l2) {
                                    Ixx += l1 * l2 * vi[iind - ixm][jind - jxm][0];
#if RELVDEBUG
                                    outfile->Printf("iind-ixm= %d, jind-jxm= %d \n", iind-ixm, jind-jxm);
#endif
                                }
                                if (l1) {
                                    Ixx -= 2.0 * l1 * a2 * vi[iind - ixm][jind + jxm][0];
#if RELVDEBUG
                                    outfile->Printf("iind-ixm= %d, jind+jxm= %d \n", iind-ixm, jind+jxm);
#endif
                                }
                                if (l2) {
                                    Ixx -= 2.0 * l2 * a1 * vi[iind + ixm][jind - jxm][0];
#if RELVDEBUG
                                    outfile->Printf("iind+ixm= %d, jind-jxm= %d \n", iind+ixm, jind-jxm);
#endif
                                }
                                if (m1 && m2) {
                                    Iyy += m1 * m2 * vi[iind - iym][jind - jym][0];
#if RELVDEBUG
                                    outfile->Printf("iind-iym= %d, jind-jym= %d \n", iind-iym, jind-jym);
#endif
                                }
                                if (m1) {
                                    Iyy -= 2.0 * m1 * a2 * vi[iind - iym][jind + jym][0];
#if RELVDEBUG
                                    outfile->Printf("iind-iym= %d, jind+jym= %d \n", iind-iym, jind+jym);
#endif
                                }
                                if (m2) {
                                    Iyy -= 2.0 * m2 * a1 * vi[iind + iym][jind - jym][0];
#if RELVDEBUG
                                    outfile->Printf("iind+iym= %d, jind-jym= %d \n", iind+iym, jind-jym);
#endif
                                }

                                if (n1 && n2) {
                                    Izz += n1 * n2 * vi[iind - izm][jind - jzm][0];
#if RELVDEBUG
                                    outfile->Printf("iind-izm= %d, jind-jzm= %d \n", iind-izm, jind-jzm);
#endif
                                }

                                if (n1) {
                                    Izz -= 2.0 * n1 * a2 * vi[iind - izm][jind + jzm][0];
#if RELVDEBUG
                                    outfile->Printf("iind-izm= %d, jind+jzm= %d \n", iind-izm, jind+jzm);
#endif
                                }
                                if (n2) {
                                    Izz -= 2.0 * n2 * a1 * vi[iind + izm][jind - jzm][0];
#if RELVDEBUG
                                    outfile->Printf("iind+izm= %d, jind-jzm= %d \n", iind+izm, jind-jzm);
#endif
                                }
                                Itotal = Ixx + Iyy + Izz;
                                buffer_[ao12++] += -Itotal * over_pf * Z;
#if RELVDEBUG
                                outfile->Printf("Itota=%20.14f \n",Itotal*over_pf*Z);
                                outfile->Printf("ao12=%d \n", ao12-1);
                                outfile->Printf("ao12=%d, vi[%d][%d][0] = %20.14f, over_pf = %20.14f, Z = %f\n", ao12-1, iind+ixm, jind+jxm, vi[iind+ixm][jind+jxm][0], over_pf, Z);
#endif
                            }
                        }
                    }
                }
            }
        }
    }
}

void RelPotentialInt::compute_pair_deriv1(const GaussianShell&, const GaussianShell&)
{
    throw SanityCheckError("RelPotentialInt::compute_pair_deriv1(): not implemented.", __FILE__, __LINE__);
}

void RelPotentialInt::compute_pair_deriv2(const GaussianShell&, const GaussianShell&)
{
    throw SanityCheckError("RelPotentialInt::compute_pair_deriv2(): not implemented.", __FILE__, __LINE__);
}

void RelPotentialInt::compute_deriv1(std::vector<SharedMatrix>&)
{
    throw SanityCheckError("RelPotentialInt::compute_deriv1(): not implemented.", __FILE__, __LINE__);
}

void RelPotentialInt::compute_deriv2(std::vector<SharedMatrix>&)
{
    throw SanityCheckError("RelPotentialInt::compute_deriv2(): not implemented.", __FILE__, __LINE__);
}

RelPotentialSOInt::RelPotentialSOInt(const std::shared_ptr<OneBodyAOInt>& aoint,
                                     const std::shared_ptr<IntegralFactory>& fact)
        : OneBodySOInt(aoint, fact)
{
    natom_ = ob_->basis1()->molecule()->natom();
}

RelPotentialSOInt::RelPotentialSOInt(const std::shared_ptr<OneBodyAOInt>& aoint, const IntegralFactory *fact)
        : OneBodySOInt(aoint, fact)
{
    natom_ = ob_->basis1()->molecule()->natom();
}

void RelPotentialSOInt::compute_deriv1(std::vector<SharedMatrix>,
                                       const CdSalcList&)
{
    throw SanityCheckError("RelPotentialSOInt::compute_deriv1(): not implemented.", __FILE__, __LINE__);
}
