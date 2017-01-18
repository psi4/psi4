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

#ifndef _psi_src_lib_libmints_potentialint_h_
#define _psi_src_lib_libmints_potentialint_h_

#include "psi4/libmints/potential.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/osrecur.h"

namespace psi{

class GaussianShell;
class OneBodyAOInt;
class IntegralFactory;
class SphericalTransform;
class Vector3;

/**
 * This is a cheesy modification to PotentialInt, to allow the in-place handling of integrals to avoid storage
 * N.B. The integrals are computed directly in the Cartesian basis and are not transformed, for efficiency.  To
 * use this code, you should transform any matrices to be contracted with these integrals to the Cartesian basis first.
 *
 * By defining the compute function of integral to be a template class, we can write classes (functors)
 * that will be inlined into the innermost loops, allowing us to do different tasks without re-writing
 * the code or having to make function calls. (AS)
 *
 * NB: This code must be specified in the .h file in order for the compiler to properly in-line the functors. (TDC)
 */
class PCMPotentialInt : public PotentialInt
{
public:
    PCMPotentialInt(std::vector<SphericalTransform>&, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>, int deriv=0);
    /// Drives the loops over all shell pairs, to compute integrals
    template<typename PCMPotentialIntFunctor>
    void compute(PCMPotentialIntFunctor &functor);
};

template<typename PCMPotentialIntFunctor>
void PCMPotentialInt::compute(PCMPotentialIntFunctor &functor)
{
    // Do not worry about zeroing out result
    int ns1 = bs1_->nshell();
    int ns2 = bs2_->nshell();
    int bf1_offset = 0;
    for (int i=0; i<ns1; ++i) {
        const GaussianShell& s1 = bs1_->shell(i);
        int ni = s1.ncartesian();
        int bf2_offset = 0;
        for (int j=0; j<ns2; ++j) {
            const GaussianShell& s2 = bs2_->shell(j);
            int nj = s2.ncartesian();
            // Compute the shell

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

            // compute intermediates
            double AB2 = 0.0;
            AB2 += (A[0] - B[0]) * (A[0] - B[0]);
            AB2 += (A[1] - B[1]) * (A[1] - B[1]);
            AB2 += (A[2] - B[2]) * (A[2] - B[2]);


            double ***vi = potential_recur_->vi();

            double** Zxyzp = Zxyz_->pointer();
            int ncharge = Zxyz_->rowspi()[0];

            for (int atom=0; atom<ncharge; ++atom) {
                memset(buffer_, 0, s1.ncartesian() * s2.ncartesian() * sizeof(double));
                double PC[3];

                double Z = Zxyzp[atom][0];

                double C[3];
                C[0] = Zxyzp[atom][1];
                C[1] = Zxyzp[atom][2];
                C[2] = Zxyzp[atom][3];
                for (int p1=0; p1<nprim1; ++p1) {
                    double a1 = s1.exp(p1);
                    double c1 = s1.coef(p1);
                    for (int p2=0; p2<nprim2; ++p2) {
                        double a2 = s2.exp(p2);
                        double c2 = s2.coef(p2);
                        double gamma = a1 + a2;
                        double oog = 1.0/gamma;

                        double PA[3], PB[3], P[3];
                        P[0] = (a1*A[0] + a2*B[0])*oog;
                        P[1] = (a1*A[1] + a2*B[1])*oog;
                        P[2] = (a1*A[2] + a2*B[2])*oog;
                        PA[0] = P[0] - A[0];
                        PA[1] = P[1] - A[1];
                        PA[2] = P[2] - A[2];
                        PB[0] = P[0] - B[0];
                        PB[1] = P[1] - B[1];
                        PB[2] = P[2] - B[2];
                        PC[0] = P[0] - C[0];
                        PC[1] = P[1] - C[1];
                        PC[2] = P[2] - C[2];

                        double over_pf = exp(-a1*a2*AB2*oog) * sqrt(M_PI*oog) * M_PI * oog * c1 * c2;


                        // Do recursion
                        potential_recur_->compute(PA, PB, PC, gamma, am1, am2);

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

                                        // Compute location in the recursion and store the value
                                        int iind = l1 * ixm + m1 * iym + n1 * izm;
                                        int jind = l2 * jxm + m2 * jym + n2 * jzm;
                                        buffer_[ao12++] += -vi[iind][jind][0] * over_pf * Z;
                                    }
                                }
                            }
                        }
                    } // End loop over primitives of shell 2
                } // End loop over primitives of shell 1
                ao12 = 0;
                int ao1 = 0;
                for(int ii = 0; ii <= am1; ii++) {
                    for(int jj = 0; jj <= ii; jj++) {
                        /*--- create all am components of sj ---*/
                        int ao2 = 0;
                        for(int kk = 0; kk <= am2; kk++) {
                            for(int ll = 0; ll <= kk; ll++) {
                                // Compute location in the recursion
                                double val = buffer_[ao12++];
                                // Hand the work off to the functor
                                functor(ao1+bf1_offset, ao2+bf2_offset, atom, val);
                                ao2++;
                            }
                        }
                        ao1++;
                    }
                }
            } // End loop over points
            bf2_offset += nj;
        } // End loop over shell 2
        bf1_offset += ni;
    } // End loop over shell 1
}

class PrintIntegralsFunctor
{
    public:
 /**
  * A functor, to be used with PCMPotentialInt, that just prints the integrals out for debugging
  */
  void operator()(int bf1, int bf2, int center, double integral)
  {
    outfile->Printf( "bf1: %3d bf2 %3d center (%5d) integral %16.10f\n", bf1, bf2, center, integral);
  }
};


class ContractOverDensityFunctor
{
 /**
  * A functor, to be used with PCMPotentialInt, that just contracts potential integrals and the
  * density matrix, over the basis function indices, giving the charge expectation value.
  */
    protected:
        /// Pointer to the density matrix.
        double **pD_;
        /// The array of charges
        double *charges_;
    public:
        ContractOverDensityFunctor(size_t /*ncenters*/, double *charges, SharedMatrix D):
            pD_(D->pointer()),
            charges_(charges)
        {
        }
        void operator()(int bf1, int bf2, int center, double integral)
        {
            charges_[center] += pD_[bf1][bf2] * integral;
        }
};


class ContractOverChargesFunctor
{
 /**
  * A functor, to be used with PCMPotentialInt, that just contracts potential integrals over charges,
  * leaving a contribution to the Fock matrix
  */
    protected:
        /// Pointer to the matrix that will contribute to the 2e part of the Fock matrix
        double **pF_;
        /// The array of charges
        const double *charges_;
    public:
        ContractOverChargesFunctor(const double* charges, SharedMatrix F):
            pF_(F->pointer()),
            charges_(charges)
        {
            if(F->rowdim() != F->coldim())
                throw PSIEXCEPTION("Invalid Fock matrix in ContractOverCharges");
            int nbf = F->rowdim();
            ::memset(pF_[0], 0, nbf*nbf*sizeof(double));
        }

        void operator()(int bf1, int bf2, int center, double integral)
        {
            pF_[bf1][bf2] += integral * charges_[center];
        }
};

} //Namespace
#endif
