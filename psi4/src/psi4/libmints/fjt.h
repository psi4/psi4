/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#include "psi4/pragma.h"

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
    virtual double* values(int J, double T) = 0;
    virtual void set_rho(double /*rho*/) {}
};

#define TAYLOR_INTERPOLATION_ORDER 6
#define TAYLOR_INTERPOLATION_AND_RECURSION \
    0  // compute F_lmax(T) and then iterate down to F_0(T)? Else use interpolation only
/// Uses Taylor interpolation of up to 8-th order to compute the Boys function
class Taylor_Fjt : public Fjt {
    static double relative_zero_;

   public:
    static const int max_interp_order = 8;

    Taylor_Fjt(size_t jmax, double accuracy);
    ~Taylor_Fjt() override;
    /// Implements Fjt::values()
    double* values(int J, double T) override;

   private:
    double** grid_;    /* Table of "exact" Fm(T) values. Row index corresponds to
                          values of T (max_T+1 rows), column index to values
                          of m (max_m+1 columns) */
    double delT_;      /* The step size for T, depends on cutoff */
    double oodelT_;    /* 1.0 / delT_, see above */
    double cutoff_;    /* Tolerance cutoff used in all computations of Fm(T) */
    int interp_order_; /* Order of (Taylor) interpolation */
    int max_m_;        /* Maximum value of m in the table, depends on cutoff
                          and the number of terms in Taylor interpolation */
    int max_T_;        /* Maximum index of T in the table, depends on cutoff
                          and m */
    double* T_crit_;   /* Maximum T for each row, depends on cutoff;
                          for a given m and T_idx <= max_T_idx[m] use Taylor interpolation,
                          for a given m and T_idx > max_T_idx[m] use the asymptotic formula */
    double* F_;        /* Here computed values of Fj(T) are stored */
};

/// "Old" intv3 code from Curt
/// Computes F_j(T) using 6-th order Taylor interpolation
class FJT : public Fjt {
   private:
    double** gtable;

    int maxj;
    double* denomarray;
    double wval_infinity;
    int itable_infinity;

    double* int_fjttable;

    int ngtable() const { return maxj + 7; }

   public:
    FJT(int n);
    ~FJT() override;
    /// implementation of Fjt::values()
    double* values(int J, double T) override;
};

class GaussianFundamental : public Fjt {
   protected:
    std::shared_ptr<CorrelationFactor> cf_;
    double rho_;
    double* value_;

   public:
    GaussianFundamental(std::shared_ptr<CorrelationFactor> cf, int max);
    ~GaussianFundamental() override;

    double* values(int J, double T) override = 0;
    void set_rho(double rho) override;
};

/**
 *  Solves \scp -\gamma r_{12}
 */
class F12Fundamental : public GaussianFundamental {
   public:
    F12Fundamental(std::shared_ptr<CorrelationFactor> cf, int max);
    ~F12Fundamental() override;
    double* values(int J, double T) override;
};

/**
 *  Solves \frac{\exp -\gamma r_{12}}{\gamma}.
 */
class F12ScaledFundamental : public GaussianFundamental {
   public:
    F12ScaledFundamental(std::shared_ptr<CorrelationFactor> cf, int max);
    ~F12ScaledFundamental() override;
    double* values(int J, double T) override;
};

class F12SquaredFundamental : public GaussianFundamental {
   public:
    F12SquaredFundamental(std::shared_ptr<CorrelationFactor> cf, int max);
    ~F12SquaredFundamental() override;
    double* values(int J, double T) override;
};

class F12G12Fundamental : public GaussianFundamental {
   private:
    std::shared_ptr<FJT> Fm_;

   public:
    F12G12Fundamental(std::shared_ptr<CorrelationFactor> cf, int max);
    ~F12G12Fundamental() override;
    double* values(int J, double T) override;
};

class F12DoubleCommutatorFundamental : public GaussianFundamental {
   public:
    F12DoubleCommutatorFundamental(std::shared_ptr<CorrelationFactor> cf, int max);
    ~F12DoubleCommutatorFundamental() override;
    double* values(int J, double T) override;
};

class ErfFundamental : public GaussianFundamental {
   private:
    double omega_;
    std::shared_ptr<FJT> boys_;

   public:
    ErfFundamental(double omega, int max);
    ~ErfFundamental() override;
    double* values(int J, double T) override;
    void setOmega(double omega) { omega_ = omega; }
};

class ErfComplementFundamental : public GaussianFundamental {
   private:
    double omega_;
    std::shared_ptr<FJT> boys_;

   public:
    ErfComplementFundamental(double omega, int max);
    ~ErfComplementFundamental() override;
    double* values(int J, double T) override;
    void setOmega(double omega) { omega_ = omega; }
};

}  // namespace psi

#endif  // header guard

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
