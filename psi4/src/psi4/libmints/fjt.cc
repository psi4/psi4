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

//
// fjt.cc
//
// Copyright (C) 2001 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include "integral.h"
#include "fjt.h"
#include "wavefunction.h"
#include "integralparameters.h"

#include "psi4/libciomr/libciomr.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <cmath>

using namespace psi;
;

const double oon[] = {0.0, 1.0, 1.0/2.0, 1.0/3.0, 1.0/4.0, 1.0/5.0, 1.0/6.0, 1.0/7.0, 1.0/8.0, 1.0/9.0, 1.0/10.0, 1.0/11.0};


// calculates a value of the boys function. Slow, but accurate
// n = J in other code
// t = T in other code
void calculate_f(double * F, int n, double t)
{
    const double eps = 1e-17;

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
        } while (fabs(term1) > eps && i < MAX_FAC);
        F[n] = sum*et;
        for(m=n-1;m>=0;m--){
            F[m] = (t2*F[m+1] + et)/(2*m+1);
        }
    }
}


Split_Fjt::Split_Fjt(unsigned int maxJ)
{
    // If we have a grid, it is big enough?
    // If not, free it
	/*
    if(initialized_ && max_m_ < mmax)
    {
        free_block(grid_);
        grid_ = nullptr;
        initialized_ = false;
        max_m_ = 0;
    }
	*/

	initialized_ = false;

    // initialize the grid if we have to
    if(!initialized_)
    {
        max_Tval_ = 43.0; // rough
        max_T_ = 430;
        max_J_ = maxJ;

        // initialize the table
        const int nrow = max_T_+1;
        const int ncol = max_J_+8+1;  // +8 for the higher orders required by the taylor series
        grid_ = block_matrix(nrow, ncol);


        for(int i = 0; i <= max_T_; i++) 
        {
            // we are using a 0.1-spaced grid
            const double Tval = 0.1*static_cast<double>(i);

            // Fill in all the values for this row
            calculate_f(grid_[i], max_J_+8, Tval);
        }
    }
}

Split_Fjt::~Split_Fjt()
{
    free_block(grid_);
}

void Split_Fjt::calculate(double * F, int J, double T)
{
    // Calculate the boys function for the highest n value
    if(T > max_Tval_)
    {
        // long range asymptotic formula
        const double p = -(2*J+1);
        const double T2 = pow(T, p);
        F[J] = boys_longfac[J] * sqrt(T2);
    }
    else
    {
        // lookup + taylor series
        // this calculates the index of the row of the lookup table
        // (ie, finds the closest value of T for which we have precomputed values)
        const int lookup_idx = (int)(10*(T+0.05));

        // this is the value of that closest value on the grid
        const double xi = ((double)lookup_idx) * 0.1;

        // delta x parameter for the taylor series
        // (but including a negative sign)
        const double dx = xi - T;

        const double * gridpts = &(grid_[lookup_idx][J]);

        F[J] = gridpts[0]
                + dx * (                  gridpts[1]
                + dx * ( (1.0/2.0   )   * gridpts[2]
                + dx * ( (1.0/6.0   )   * gridpts[3]
                + dx * ( (1.0/24.0  )   * gridpts[4]
                + dx * ( (1.0/120.0 )   * gridpts[5]
                + dx * ( (1.0/720.0 )   * gridpts[6]
                + dx * ( (1.0/5040.0)   * gridpts[7]
                )))))));

    }

    // calculate the others by downward recursion
    if(J > 0)
    {
        const double eT = exp(-T);
        for(int m = J-1; m >= 0; m--)
            F[m] = (2*T*F[m+1] + eT)/(2*m+1);
    }
}

/////////////////////////////////////////////////////////////////////////////

////////
// GaussianFundamental
////////

GaussianFundamental::GaussianFundamental(std::shared_ptr<CorrelationFactor> cf, int maxJ)
{
    cf_ = cf;

    // For now, set the rho vars to zero. They will be set by the compute_shell routine.
    rho_ = 0.0;
}

GaussianFundamental::~GaussianFundamental()
{
}

void GaussianFundamental::set_rho(double rho)
{
    rho_ = rho;
}

////////
// F12Fundamental
////////

F12Fundamental::F12Fundamental(std::shared_ptr<CorrelationFactor> cf, int max)
    : GaussianFundamental(cf, max)
{

}

F12Fundamental::~F12Fundamental()
{

}

void F12Fundamental::calculate(double * F, int J, double T)
{
    // because the current implementation is just a hack of the eri
    // routines, we have to undo the eri prefactor of 2pi/rho that
    // will be added later
    double* exps = cf_->exponent();
    double* coeffs = cf_->coeff();
    int nparam = cf_->nparam();

    // zero the values array
    for (int n=0; n<=J; ++n)
        F[n] = 0.0;

    double pfac, expterm, rhotilde, omega;
    double eri_correct = rho_ / 2 / M_PI;
    for (int i=0; i<nparam; ++i) {
        omega = exps[i];
        rhotilde = omega / (rho_ + omega);
        pfac = coeffs[i] * pow(M_PI/(rho_ + omega), 1.5) * eri_correct;
        expterm = exp(-rhotilde*T)*pfac;
        for (int n=0; n<=J; ++n) {
            F[n] += expterm;
            expterm *= rhotilde;
        }
    }
}

////////
// F12ScaledFundamental
////////

F12ScaledFundamental::F12ScaledFundamental(std::shared_ptr<CorrelationFactor> cf, int max)
: GaussianFundamental(cf, max)
{

}

F12ScaledFundamental::~F12ScaledFundamental()
{

}

void F12ScaledFundamental::calculate(double * F, int J, double T)
{
    // because the current implementation is just a hack of the eri
    // routines, we have to undo the eri prefactor of 2pi/rho that
    // will be added later
    double* exps = cf_->exponent();
    double* coeffs = cf_->coeff();
    int nparam = cf_->nparam();

    // zero the values array
    for (int n=0; n<=J; ++n)
        F[n] = 0.0;

    double pfac, expterm, rhotilde, omega;
    double eri_correct = rho_ / 2 / M_PI;
    eri_correct /= cf_->slater_exponent();
    for (int i=0; i<nparam; ++i) {
        omega = exps[i];
        rhotilde = omega / (rho_ + omega);
        pfac = coeffs[i] * pow(M_PI/(rho_ + omega), 1.5) * eri_correct;
        expterm = exp(-rhotilde*T)*pfac;
        for (int n=0; n<=J; ++n) {
            F[n] += expterm;
            expterm *= rhotilde;
        }
    }
}

////////
// F12SquaredFundamental
////////

F12SquaredFundamental::F12SquaredFundamental(std::shared_ptr<CorrelationFactor> cf, int max)
    : GaussianFundamental(cf, max)
{

}

F12SquaredFundamental::~F12SquaredFundamental()
{

}

void F12SquaredFundamental::calculate(double * F, int J, double T)
{
    double* exps = cf_->exponent();
    double* coeffs = cf_->coeff();
    int nparam = cf_->nparam();

    double pfac, expterm, rhotilde, omega;
    double eri_correct = rho_ / 2 / M_PI;

    // zero the values
    for (int n=0; n<=J; ++n)
        F[n] = 0.0;

    for (int i=0; i<nparam; ++i) {
        for (int j=0; j<nparam; ++j) {
            omega = exps[i] + exps[j];
            rhotilde = omega / (rho_ + omega);
            pfac = coeffs[i] * coeffs[j] * pow(M_PI/(rho_+omega), 1.5) * eri_correct;
            expterm = exp(-rhotilde * T) * pfac;
            for (int n=0; n<=J; ++n) {
                F[n] += expterm;
                expterm *= rhotilde;
            }
        }
    }
}

////////
// F12G12Fundamental
////////

F12G12Fundamental::F12G12Fundamental(std::shared_ptr<CorrelationFactor> cf, int max)
    : GaussianFundamental(cf, max), Fm_(max)
{
}

F12G12Fundamental::~F12G12Fundamental()
{
}

void F12G12Fundamental::calculate(double * F, int J, double T)
{
    double Fvals[J+1];

    double* exps = cf_->exponent();
    double* coeffs = cf_->coeff();
    int nparam = cf_->nparam();

    double pfac, expterm, rhotilde, omega, rhohat;
    double boysterm, rhotilde_term, rhohat_term;
    double eri_correct = rho_ / 2 / M_PI;
    double binom_term;

    // Zero the values
    for (int n=0; n<=J; ++n)
        F[n] = 0.0;

    for (int i=0; i<nparam; ++i) {
        omega = exps[i];
        rhotilde = omega / (rho_ + omega);
        rhohat = rho_ / (rho_ + omega);
        expterm = exp(-rhotilde * T);
        pfac = 2*M_PI / (rho_ + omega) * coeffs[i] * expterm * eri_correct;
        Fm_.calculate(Fvals, J, rhohat * T);
        for (int n=0; n<=J; ++n) {
            boysterm = 0.0;
            rhotilde_term = pow(rhotilde, n);
            rhohat_term = 1.0;
            for (int m=0; m<=n; ++m) {
                binom_term = bc[n][m];   // bc is formed by Wavefunction::initialize_singletons
                boysterm += binom_term * rhotilde_term * rhohat_term * Fvals[m];

                rhotilde_term /= rhotilde;
                rhohat_term *= rhohat;
            }
            F[n] += pfac * boysterm;
        }
    }
}

////////
// F12DoubleCommutatorFundamental
////////

F12DoubleCommutatorFundamental::F12DoubleCommutatorFundamental(std::shared_ptr<CorrelationFactor> cf, int max)
    : GaussianFundamental(cf, max)
{
}

F12DoubleCommutatorFundamental::~F12DoubleCommutatorFundamental()
{
}

void F12DoubleCommutatorFundamental::calculate(double * F, int J, double T)
{
    double *exps = cf_->exponent();
    double *coeffs = cf_->coeff();
    int nparam = cf_->nparam();

    double pfac, expterm, rhotilde, omega, sqrt_term, rhohat;
    double eri_correct = rho_ / 2 / M_PI;
    double term1, term2;

    // Zero the values
    for (int n=0; n<=J; ++n)
        F[n] = 0.0;

    for (int i=0; i<nparam; ++i) {
        for (int j=0; j<nparam; ++j) {
            omega = exps[i] + exps[j];
            rhotilde = omega / (rho_ + omega);
            rhohat = rho_ / (rho_ + omega);
            expterm = exp(-rhotilde * T);
            sqrt_term = sqrt(M_PI*M_PI*M_PI / pow(rho_ + omega, 5.0));
            pfac = 4.0*coeffs[i] * coeffs[j] * exps[i] * exps[j] * sqrt_term * eri_correct * expterm;

            term1 = 1.5*rhotilde + rhotilde*rhohat*T;
            term2 = 1.0/rhotilde*pfac;
            for (int n=0; n<=J; ++n) {
                F[n] += term1 * term2;
                term1 -= rhohat;
                term2 *= rhotilde;
            }
        }
    }
}

////////
// ErfFundamental
////////

ErfFundamental::ErfFundamental(double omega, int max)
    : GaussianFundamental(std::shared_ptr<CorrelationFactor>(), max),
      boys_(max)
{
    omega_ = omega;
    rho_ = 0;
}

ErfFundamental::~ErfFundamental()
{
}

void ErfFundamental::calculate(double * F, int J, double T)
{
    double Fvals[J+1];
    boys_.calculate(Fvals, J, T);

    for (int n=0; n<=J; ++n)
        F[n] = 0.0;

    // build the erf constants
    double omegasq = omega_ * omega_;
    double T_prefac = omegasq / (omegasq + rho_);
    double F_prefac = sqrt(T_prefac);
    double erf_T = T_prefac * T;

    boys_.calculate(Fvals, J, erf_T);
    for (int n=0; n<=J; ++n) {
        F[n] += Fvals[n] * F_prefac;
        F_prefac *= T_prefac;
    }
}

////////
// ErfComplementFundamental
////////

ErfComplementFundamental::ErfComplementFundamental(double omega, int max)
    : GaussianFundamental(std::shared_ptr<CorrelationFactor>(), max),
      boys_(max)
{
    omega_ = omega;
    rho_ = 0;
}

ErfComplementFundamental::~ErfComplementFundamental()
{
}

void ErfComplementFundamental::calculate(double * F, int J, double T)
{
    double Fvals[J+1];
    boys_.calculate(Fvals, J, T);

    for (int n=0; n<=J; ++n)
        F[n] = Fvals[n];

    // build the erf constants
    double omegasq = omega_ * omega_;
    double T_prefac = omegasq / (omegasq + rho_);
    double F_prefac = sqrt(T_prefac);
    double erf_T = T_prefac * T;

    boys_.calculate(Fvals, J, erf_T);
    for (int n=0; n<=J; ++n) {
        F[n] -= Fvals[n] * F_prefac;
        F_prefac *= T_prefac;
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
