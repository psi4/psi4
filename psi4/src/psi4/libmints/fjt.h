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
// fjt.h
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
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifndef _chemistry_qc_basis_fjt_h
#define _chemistry_qc_basis_fjt_h

// This table stores some common factors for the boys function
// asymptotic form
#define BOYS_LONGFAC_MAXN 50 // maximum value of nu stored
static const double boys_longfac[BOYS_LONGFAC_MAXN + 1] =
{
 /* n =  0 */  0.886226925452758014    , /* n =  1 */  0.443113462726379007    ,
 /* n =  2 */  0.66467019408956851     , /* n =  3 */  1.66167548522392128     ,
 /* n =  4 */  5.81586419828372446     , /* n =  5 */  26.1713888922767601     ,
 /* n =  6 */  143.94263890752218      , /* n =  7 */  935.627152898894173     ,
 /* n =  8 */  7017.2036467417063      , /* n =  9 */  59646.2309973045035     ,
 /* n = 10 */  566639.194474392784     , /* n = 11 */  5949711.54198112423     ,
 /* n = 12 */  68421682.7327829286     , /* n = 13 */  855271034.159786608     ,
 /* n = 14 */  11546158961.1571192     , /* n = 15 */  167419304936.778228     ,
 /* n = 16 */  2594999226520.06254     , /* n = 17 */  42817487237581.0319     ,
 /* n = 18 */  749306026657668.059     , /* n = 19 */  13862161493166859.1     ,
 /* n = 20 */  270312149116753752.0    , /* n = 21 */  5.54139905689345192e+18 ,
 /* n = 22 */  1.19140079723209216e+20 , /* n = 23 */  2.68065179377220737e+21 ,
 /* n = 24 */  6.29953171536468731e+22 , /* n = 25 */  1.54338527026434839e+24 ,
 /* n = 26 */  3.9356324391740884e+25  , /* n = 27 */  1.04294259638113343e+27 ,
 /* n = 28 */  2.86809214004811692e+28 , /* n = 29 */  8.17406259913713322e+29 ,
 /* n = 30 */  2.4113484667454543e+31  , /* n = 31 */  7.35461282357363562e+32 ,
 /* n = 32 */  2.31670303942569522e+34 , /* n = 33 */  7.52928487813350946e+35 ,
 /* n = 34 */  2.52231043417472567e+37 , /* n = 35 */  8.70197099790280356e+38 ,
 /* n = 36 */  3.08919970425549526e+40 , /* n = 37 */  1.12755789205325577e+42 ,
 /* n = 38 */  4.22834209519970914e+43 , /* n = 39 */  1.62791170665188802e+45 ,
 /* n = 40 */  6.43025124127495768e+46 , /* n = 41 */  2.60425175271635786e+48 ,
 /* n = 42 */  1.08076447737728851e+50 , /* n = 43 */  4.59324902885347618e+51 ,
 /* n = 44 */  1.99806332755126214e+53 , /* n = 45 */  8.89138180760311651e+54 ,
 /* n = 46 */  4.04557872245941801e+56 , /* n = 47 */  1.88119410594362937e+58 ,
 /* n = 48 */  8.93567200323223953e+59 , /* n = 49 */  4.33380092156763617e+61 ,
 /* n = 50 */  2.14523145617597991e+63 ,
};




namespace psi {

class CorrelationFactor;

/// Evaluates the Boys function F_j(T)
class Fjt {
public:
    Fjt();
    virtual ~Fjt();
    /** Computed F_j(T) for every 0 <= j <= J (total of J+1 doubles).
        The user may read/write these values.
        The values will be overwritten with the next call to this functions.
        The pointer will be invalidated after the call to ~Fjt. */
    virtual double *values(int J, double T) =0;
    virtual void set_rho(double /*rho*/) { }
};

/// Uses Taylor interpolation of up to 8-th order to compute the Boys function
class Split_Fjt : public Fjt {
public:
    Split_Fjt(unsigned int jmax);
    virtual ~Split_Fjt();
    /// Implements Fjt::values()
    double *values(int J, double T);
private:

    bool initialized_;  /* Has the table been initialized */
    double **grid_;     /* Table of "exact" Fm(T) values. Row index corresponds to
                           values of T (max_T_+1 rows), column index to values
                           of m (max_m+1 columns) */
	int max_T_;         /* Maximum T index stored (10*max T value) */
	double max_Tval_;   /* Maximum T value stored */
    int max_m_;         /* Maximum value of m in the table, depends on cutoff
                                  and the number of terms in Taylor interpolation */
    double *F_;         /* Here computed values of Fj(T) are stored */
};

class GaussianFundamental : public Fjt {
protected:
    std::shared_ptr<CorrelationFactor> cf_;
    double rho_;
    double* value_;

public:
    GaussianFundamental(std::shared_ptr<CorrelationFactor> cf, int max);
    virtual ~GaussianFundamental();

    virtual double* values(int J, double T) = 0;
    void set_rho(double rho);
};

    /**
     *  Solves \scp -\gamma r_{12}
     */
class F12Fundamental : public GaussianFundamental {
public:
    F12Fundamental(std::shared_ptr<CorrelationFactor> cf, int max);
    virtual ~F12Fundamental();
    double* values(int J, double T);
};

/**
 *  Solves \frac{\exp -\gamma r_{12}}{\gamma}.
 */
class F12ScaledFundamental : public GaussianFundamental {
public:
    F12ScaledFundamental(std::shared_ptr<CorrelationFactor> cf, int max);
    virtual ~F12ScaledFundamental();
    double* values(int J, double T);
};

class F12SquaredFundamental : public GaussianFundamental {
public:
    F12SquaredFundamental(std::shared_ptr<CorrelationFactor> cf, int max);
    virtual ~F12SquaredFundamental();
    double* values(int J, double T);
};

class F12G12Fundamental : public GaussianFundamental {
private:
    std::shared_ptr<Split_Fjt> Fm_;
public:
    F12G12Fundamental(std::shared_ptr<CorrelationFactor> cf, int max);
    virtual ~F12G12Fundamental();
    double* values(int J, double T);
};

class F12DoubleCommutatorFundamental : public GaussianFundamental {
public:
    F12DoubleCommutatorFundamental(std::shared_ptr<CorrelationFactor> cf, int max);
    virtual ~F12DoubleCommutatorFundamental();
    double* values(int J, double T);
};

class ErfFundamental : public GaussianFundamental {
private:
    double omega_;
    std::shared_ptr<Split_Fjt> boys_;
public:
    ErfFundamental(double omega, int max);
    virtual ~ErfFundamental();
    double* values(int J, double T);
    void setOmega(double omega) { omega_ = omega; }
};

class ErfComplementFundamental : public GaussianFundamental {
private:
    double omega_;
    std::shared_ptr<Split_Fjt> boys_;
public:
    ErfComplementFundamental(double omega, int max);
    virtual ~ErfComplementFundamental();
    double* values(int J, double T);
    void setOmega(double omega) { omega_ = omega; }
};

} // end of namespace sc

#endif // header guard

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
