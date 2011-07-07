/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680

  $Id: solvers.h 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/
#ifndef MADNESS_LINALG_SOLVERS_H__INCLUDED
#define MADNESS_LINALG_SOLVERS_H__INCLUDED

#include <tensor/tensor.h>
#include <world/print.h>
#include <iostream>
#include <linalg/tensor_lapack.h>

/*!
  \file solvers.h
  \brief Defines interfaces for optimization and non-linear equation solvers
  \ingroup solvers
*/

namespace madness {

    /*!
      \brief Solves non-linear equation using KAIN (returns coefficients to compute next vector)

      \ingroup solvers

      The Krylov-accelerated inexact-Newton method employs directional
      derivatives to estimate the Jacobian in the subspace and
      separately computes updates in the subspace and its complement.

      We wish to solve the non-linear equations \f$ f(x)=0 \f$ where
      \f$ f \f$ and \f$ x \f$ are vectors of the same dimension (e.g.,
      consider both being MADNESS functions).

      Define the following matrices and vector (with \f$ i \f$ and \f$
      j \f$ denoting previous iterations in the Krylov subspace and
      \f$ m \f$ the current iteration):
      \f{eqnarray*}{
      Q_{i j} & = & \langle x_i \mid f_j \rangle \\
      A_{i j} & = & \langle x_i - x_m \mid f_j - f_m \rangle = Q_{i j} - Q_{m j} - Q_{i m} + Q{m m} \\
      b_i & =& -\langle x_i - x_m \mid f_m \rangle = -Q_{i m} + Q_{m m}
      \f}
      The subspace equation is of dimension \f$ m \f$ (assuming iterations
      are indexed from zero) and is given by
      \f[
      A c = b
      \f]
      The interior and exterior updates may be combined into one simple expression
      as follows. First, define an additional element of the solution vector
      \f[
      c_m = 1 - \sum_{i<m} c_i
      \f]
      and then the new vector (guess for next iteration) is given by
      \f[
      x_{m+1} = \sum_{i \le m}{c_i ( x_i - f_i)}
      \f]

      To employ the solver, each iteration
      -# Compute the additional row and column of the matrix \f$ Q \f$
      that is the inner product between solution vectors (\f$ x_i \f$) and residuals
      (\f$ f_j \f$).
      -# Call this routine to compute the coefficients \f$ c \f$ and from these
      compute the next solution vector
      -# Employ step restriction or line search as necessary to ensure stable/robust solution.

      @param[in] Q The matrix of inner products between subspace vectors and residuals.
      @param[in] rcond Threshold for discarding small singular values in the subspace equations.
      @return Vector for computing next solution vector
    */
    template <typename T>
    Tensor<T> KAIN(const Tensor<T>& Q, double rcond=1e-12) {
        const int nvec = Q.dim(0);
        const int m = nvec-1;

        if (nvec == 1) {
            Tensor<T> c(1);
            c(0L) = 1.0;
            return c;
        }

        Tensor<T> A(m,m);
        Tensor<T> b(m);
        for (long i=0; i<m; ++i) {
            b(i) = Q(m,m) - Q(i,m);
            for (long j=0; j<m; ++j) {
                A(i,j) = Q(i,j) - Q(m,j) - Q(i,m) + Q(m,m);
            }
        }

    //     print("Q");
    //     print(Q);
    //     print("A");
    //     print(A);
    //     print("b");
    //     print(b);

        Tensor<T> x;
        Tensor<double> s, sumsq;
        long rank;
        gelss(A, b, rcond, x, s, rank, sumsq);
//         print("singular values", s);
//         print("rank", rank);
//         print("solution", x);

        Tensor<T> c(nvec);
        T sumC = 0.0;
        for (long i=0; i<m; ++i) sumC += x(i);
        c(Slice(0,m-1)) = x;
//         print("SUMC", nvec, m, sumC);
        c(m) = 1.0 - sumC;

//         print("returned C", c);

        return c;
    }


    /// The interface to be provided by targets for non-linear equation solver

    /// \ingroup solvers
    struct SolverTargetInterface {
        /// Should return the resdiual (vector F(x))
        virtual Tensor<double> residual(const Tensor<double>& x) = 0;

        /// Override this to return \c true if the Jacobian is implemented
        virtual bool provides_jacobian() const {return false;}

        /// Some solvers require the jacobian or are faster if an analytic form is available

        /// J(i,j) = partial F[i] over partial x[j] where F(x) is the vector valued residual
        virtual Tensor<double> jacobian(const Tensor<double>& x) {
            throw "not implemented";
        }

        /// Implement this if advantageous to compute residual and jacobian simultaneously
        virtual void residual_and_jacobian(const Tensor<double>& x,
                                           Tensor<double>& residual, Tensor<double>& jacobian) {
            residual = this->residual(x);
            jacobian = this->jacobian(x);
        }

        virtual ~SolverTargetInterface() {}
    };


    /// The interface to be provided by functions to be optimized

    /// \ingroup solvers
    struct OptimizationTargetInterface {
        /// Should return the value of the objective function
        virtual double value(const Tensor<double>& x) = 0;

        /// Override this to return true if the derivative is implemented
        virtual bool provides_gradient() const {return false;}

        /// Should return the derivative of the function
        virtual Tensor<double> gradient(const Tensor<double>& x) {
            throw "not implemented";
        }

        /// Reimplement if more efficient to evaluate both value and gradient in one call
        virtual void value_and_gradient(const Tensor<double>& x,
                                        double& value,
                                        Tensor<double>& gradient) {
            value = this->value(x);
            gradient = this->gradient(x);
        }

        /// Numerical test of the derivative ... optionally prints to stdout, returns max abs error
        double test_gradient(Tensor<double>& x, double value_precision, bool doprint=true);

	virtual ~OptimizationTargetInterface(){}
    };


    /// The interface to be provided by solvers ... NOT USED ANYWHERE?

    /// \ingroup solvers
    struct SolverInterface {
        virtual bool solve(Tensor<double>& x) = 0;
        virtual bool converged() const = 0;
        virtual double residual_norm() const = 0;
	virtual ~SolverInterface() {}
    };

    /// The interface to be provided by optimizers

    /// \ingroup solvers
    struct OptimizerInterface {
        virtual bool optimize(Tensor<double>& x) = 0;
        virtual bool converged() const = 0;
        virtual double value() const = 0;
        virtual double gradient_norm() const = 0;
	virtual ~OptimizerInterface(){}
    };


    /// Unconstrained minimization via steepest descent

    /// \ingroup solvers
    class SteepestDescent : public OptimizerInterface {
        std::shared_ptr<OptimizationTargetInterface> target;
        const double tol;
        const double value_precision;  // Numerical precision of value
        const double gradient_precision; // Numerical precision of each element of residual
        double f;
        double gnorm;

    public:
        SteepestDescent(const std::shared_ptr<OptimizationTargetInterface>& tar,
                        double tol = 1e-6,
                        double value_precision = 1e-12,
                        double gradient_precision = 1e-12);

        bool optimize(Tensor<double>& x);

        bool converged() const;

        double gradient_norm() const;

        double value() const;

        virtual ~SteepestDescent() { }
    };


    /// Optimization via quasi-Newton (BFGS or SR1 update)

    /// \ingroup solvers
    /// This is presently not a low memory algorithm ... we really need one!
    class QuasiNewton : public OptimizerInterface {
    private:
        std::string update;              // One of BFGS or SR1
        std::shared_ptr<OptimizationTargetInterface> target;
        const double tol;
        const double value_precision;  // Numerical precision of value
        const double gradient_precision; // Numerical precision of each element of residual
        double f;
        double gnorm;
        Tensor<double> h;
        int n;

        double line_search(double a1, double f0, double dxgrad, const Tensor<double>& x, const Tensor<double>& dx);

        void hessian_update_sr1(const Tensor<double>& s, const Tensor<double>& y);

        void hessian_update_bfgs(const Tensor<double>& dx,
                                 const Tensor<double>& dg);

        Tensor<double> new_search_direction(const Tensor<double>& g);

    public:
        QuasiNewton(const std::shared_ptr<OptimizationTargetInterface>& tar,
                    double tol = 1e-6,
                    double value_precision = 1e-12,
                    double gradient_precision = 1e-12);

        /// Choose update method (currently only "BFGS" or "SR1")
        void set_update(const std::string& method);

        /// Runs the optimizer

        /// @return True if converged
        bool optimize(Tensor<double>& x);

        /// After running the optimizer returns true if converged

        /// @return True if converged
        bool converged() const;

        /// Value of objective function

        /// @return Value of objective function
        double value() const;

        /// Value of gradient norm

        /// @return Norm of gradient of objective function
        double gradient_norm() const;

	virtual ~QuasiNewton() {}
    };
}

#endif // MADNESS_LINALG_SOLVERS_H__INCLUDED
