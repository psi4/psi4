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

namespace boost {
template<class T> class shared_ptr;
}

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

#define TAYLOR_INTERPOLATION_ORDER 6
#define TAYLOR_INTERPOLATION_AND_RECURSION 0  // compute F_lmax(T) and then iterate down to F_0(T)? Else use interpolation only
/// Uses Taylor interpolation of up to 8-th order to compute the Boys function
class Taylor_Fjt : public Fjt {
    static double relative_zero_;
public:
    static const int max_interp_order = 8;

    Taylor_Fjt(unsigned int jmax, double accuracy);
    virtual ~Taylor_Fjt();
    /// Implements Fjt::values()
    double *values(int J, double T);
private:
    double **grid_;            /* Table of "exact" Fm(T) values. Row index corresponds to
                                  values of T (max_T+1 rows), column index to values
                                  of m (max_m+1 columns) */
    double delT_;              /* The step size for T, depends on cutoff */
    double oodelT_;            /* 1.0 / delT_, see above */
    double cutoff_;            /* Tolerance cutoff used in all computations of Fm(T) */
    int interp_order_;         /* Order of (Taylor) interpolation */
    int max_m_;                /* Maximum value of m in the table, depends on cutoff
                                  and the number of terms in Taylor interpolation */
    int max_T_;                /* Maximum index of T in the table, depends on cutoff
                                  and m */
    double *T_crit_;           /* Maximum T for each row, depends on cutoff;
                                  for a given m and T_idx <= max_T_idx[m] use Taylor interpolation,
                                  for a given m and T_idx > max_T_idx[m] use the asymptotic formula */
    double *F_;                /* Here computed values of Fj(T) are stored */

    class ExpensiveMath {
    public:
        ExpensiveMath(int ifac, int idf);
        ~ExpensiveMath();
        double *fac;
        double *df;
    };
    ExpensiveMath ExpMath_;
};

/// "Old" intv3 code from Curt
/// Computes F_j(T) using 6-th order Taylor interpolation
class FJT: public Fjt {
private:
    double **gtable;

    int maxj;
    double *denomarray;
    double wval_infinity;
    int itable_infinity;

    double *int_fjttable;

    int ngtable() const { return maxj + 7; }
public:
    FJT(int n);
    virtual ~FJT();
    /// implementation of Fjt::values()
    double *values(int J, double T);
};

class GaussianFundamental : public Fjt {
protected:
    boost::shared_ptr<CorrelationFactor> cf_;
    double rho_;
    double* value_;

public:
    GaussianFundamental(boost::shared_ptr<CorrelationFactor> cf, int max);
    virtual ~GaussianFundamental();

    virtual double* values(int J, double T) = 0;
    void set_rho(double rho);
};

class F12Fundamental : public GaussianFundamental {
public:
    F12Fundamental(boost::shared_ptr<CorrelationFactor> cf, int max);
    virtual ~F12Fundamental();
    double* values(int J, double T);
};

class F12SquaredFundamental : public GaussianFundamental {
public:
    F12SquaredFundamental(boost::shared_ptr<CorrelationFactor> cf, int max);
    virtual ~F12SquaredFundamental();
    double* values(int J, double T);
};

class F12G12Fundamental : public GaussianFundamental {
private:
    boost::shared_ptr<FJT> Fm_;
public:
    F12G12Fundamental(boost::shared_ptr<CorrelationFactor> cf, int max);
    virtual ~F12G12Fundamental();
    double* values(int J, double T);
};

class F12DoubleCommutatorFundamental : public GaussianFundamental {
public:
    F12DoubleCommutatorFundamental(boost::shared_ptr<CorrelationFactor> cf, int max);
    virtual ~F12DoubleCommutatorFundamental();
    double* values(int J, double T);
};

class ErfFundamental : public GaussianFundamental {
private:
    double omega_;
    boost::shared_ptr<FJT> boys_;
public:
    ErfFundamental(double omega, int max);
    virtual ~ErfFundamental();
    double* values(int J, double T);
    void setOmega(double omega) { omega_ = omega; }
};

class ErfComplementFundamental : public GaussianFundamental {
private:
    double omega_;
    boost::shared_ptr<FJT> boys_;
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
