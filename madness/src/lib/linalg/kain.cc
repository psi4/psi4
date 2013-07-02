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

  $Id$
*/
#include <iostream>
#include <tensor/tensor.h>
#include <linalg/tensor_lapack.h>
#include <world/print.h>

using namespace madness;
using namespace std;



/// Solves the KAIN equations for coefficients to compute the next vector

/// \verbatim
///   Q(i,j) = <xi|fj>
///   A(i,j) = <xi-xm | fj-fm> = Q(i,j) - Q(m,j) - Q(i,m) + Q(m,m)
///   b(i,j) =-<xi-xm | fm> = -Q(i,m) + Q(m,m)
///   A c = b
///
///   . Correction to vector m
///   .   interior = sum(i<m)[ c(i)*(x(i)-x(m)) ] = sum(i<m)[c(i)*x(i)] - x[m]*sum(i<m)[c(i)]
///   .   exterior = -f(m) - sum(i<m)[ c(i)*(f(i)-f(m)) ] = -f(m) - sum(i<m)[c(i)*f(i)] + f(m)*sum(i<m)[c(i)]
///   . New vector
///   .   define C = sum(i<m)(c(i))  (note less than)
///   .   define c(m) = 1.0 - C
///   .   xnew = sum(i<=m) [ c(i)*(x(i) - f(i)) ]
/// \endverbatim
template <typename T>
Tensor<T> KAIN(const Tensor<T>& Q) {
    const int nvec = Q.dim[0];
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

    double rcond = 1e-12;
    Tensor<T> x;
    Tensor<double> s;
    long rank;
    gelss(A, b, rcond, &x, &s, &rank);
    print("singular values", s);
    print("rank", rank);
    print("solution", x);

    Tensor<T> c(nvec);
    T sumC = 0.0;
    for (long i=0; i<m; ++i) sumC += x(i);
    c(Slice(0,m-1)) = x;
    print("SUMC", nvec, m, sumC);
    c(m) = 1.0 - sumC;

    print("returned C", c);

    return c;
}


/// The interface to be provided by
struct SolverTargetInterface {
    virtual bool provides_jacobian() const = 0;
    virtual Tensor<double> residual(const Tensor<double>& x) = 0;
    virtual Tensor<double> jacobian(const Tensor<double>& x) {
        throw "not implemented";
    }
    virtual void residual_and_jacobian(const Tensor<double>& x,
                                       Tensor<double>& residual, Tensor<double>& jacobian) {
        residual = this->residual(x);
        jacobian = this->jacobian(x);
    }
    virtual ~SolverTargetInterface() {};
};


/// The interface to be provided by functions to be optimized
struct OptimizationTargetInterface {

    /// Should return true if the gradient is implemented
    virtual bool provides_gradient() const = 0;

    /// Should return the value of the objective function
    virtual double value(const Tensor<double>& x) = 0;

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
    double test_gradient(Tensor<double>& x, double value_precision, bool doprint=true) {
        const double eps = pow(value_precision,0.3333);
        if (doprint) {
            printf("\n");
            printf("Testing gradient eps=%.1e\n----------------\n", eps);
            printf("  i            f-                 f+               ganalytic           gnumeric         err\n");
            printf(" ----  ------------------  ------------------  ------------------  ------------------  -------\n");
        }
        Tensor<double> tt = gradient(x);
        int n = int(tt.dim[0]);
        double maxerr = 0.0;
        for (int i=0; i<n; ++i) {
            x[i] += eps;
            double fp = value(x);
            x[i] -= 2.0*eps;
            double fm = value(x);
            x[i] += eps;

            double gg = 0.5*(fp-fm)/eps;
            if (doprint)
                printf("% 5d%20.12e%20.12e%20.12e%20.12e  %.1e\n", i, fm, fp, tt(i), gg, abs(tt(i)-gg));
            maxerr = max(abs(gg-tt(i)),maxerr);
        }
        if (doprint) printf("\n");

        return maxerr;
    }
};


/// The interface to be provided by solvers
struct SolverInterface {
    virtual bool solve(Tensor<double>& x) = 0;
    virtual bool converged() const = 0;
    virtual double residual_norm() const = 0;
};

/// The interface to be provided by optimizers
struct OptimizerInterface {
    virtual bool optimize(Tensor<double>& x) = 0;
    virtual bool converged() const = 0;
    virtual double value() const = 0;
    virtual double gradient_norm() const = 0;
};


/// Optimization via steepest descent
class SteepestDescent : public OptimizerInterface {
    std::shared_ptr<OptimizationTargetInterface> target;
    const double tol;
    const double value_precision;  // Numerical precision of value
    const double gradient_precision; // Numerical precision of each element of residual
    double f;
    double gnorm;

public:
    SteepestDescent(const std::shared_ptr<OptimizationTargetInterface>& target,
                    double tol = 1e-6,
                    double value_precision = 1e-12,
                    double gradient_precision = 1e-12)
        : target(target)
        , tol(tol)
        , value_precision(value_precision)
        , gradient_precision(gradient_precision)
        , gnorm(tol*1e16)
    {
        if (!target->provides_gradient()) throw "Steepest descent requires the gradient";
    }

    bool optimize(Tensor<double>& x) {
        double step = 10.0;
        double  fnew;
        Tensor<double> g;
        target->value_and_gradient(x,f,g);
        gnorm = g.normf();
        for (int i=0; i<100; ++i) {
            while (1) {
                Tensor<double> gnew;
                x.gaxpy(1.0, g, -step);
                target->value_and_gradient(x,fnew,gnew);
                if (fnew < f) {
                    f = fnew;
                    g = gnew;
                    break;
                }
                x.gaxpy(1.0, g, step);
                step *= 0.5;
                print("reducing step size", f, fnew, step);
            }
            Tensor<double> g = target->gradient(x);
            gnorm = g.normf();
            print("iteration",i,"value",f,"gradient",gnorm);
            if (converged()) break;
        }
        return converged();
    }

    bool converged() const {return gnorm < tol;}

    double gradient_norm() const {return gnorm;}

    double value() const {return f;}
};


/// Optimization via quasi-Newton (BFGS or SR1 update)
class QuasiNewton : public OptimizerInterface {
private:
    string update;              // One of BFGS or SR1
    std::shared_ptr<OptimizationTargetInterface> target;
    const double tol;
    const double value_precision;  // Numerical precision of value
    const double gradient_precision; // Numerical precision of each element of residual
    double f;
    double gnorm;
    Tensor<double> h;
    int n;


    double line_search(double a1, double f0, double dxgrad, const Tensor<double>& x, const Tensor<double>& dx) {
        double f1, f2p;
        double hess, a2;
        const char* lsmode = "";

        if (dxgrad*a1 > 0.0) {
            print("    line search gradient +ve ", a1, dxgrad);
            a1 = -a1;
        }

        f1 = target->value(x + a1 * dx);

        // Fit to a parabola using f0, g0, f1
        hess = 2.0*(f1-f0-a1*dxgrad)/(a1*a1);
        a2 = -dxgrad/hess;

        if (abs(f1-f0) < value_precision) { // Insufficient precision
            a2 = a1;
            lsmode = "fixed";
        }
        else if (hess > 0.0) { // Positive curvature
            if ((f1 - f0) <= -value_precision) { // Downhill
                lsmode = "downhill";
                if (abs(a2) > 4.0*abs(a1)) {
                    lsmode = "restrict";
                    a2 = 4.0*a1;
                }
            }
            else { // Uphill
                lsmode = "bracket";
            }
        }
        else { // Negative curvature
            if ((f1 - f0) < value_precision) { // Downhill
                lsmode = "negative";
                a2 = 2e0*a1;
            }
            else {
                lsmode = "punt";
                a2 = a1;
            }
        }

        f2p = f0 + dxgrad*a2 + 0.5*hess*a2*a2;
        printf("   line search grad=%.2e hess=%.2e mode=%s newstep=%.3f\n", dxgrad, hess, lsmode, a2);
        printf("                      predicted %.12f\n", f2p);

        return a2;
    }

    void hessian_update_sr1(const Tensor<double>& s, const Tensor<double>& y) {
        Tensor<double> q = y - inner(h,s);
        double qds = q.trace(s);
        if (abs(qds) > 1e-8 * s.normf() * q.normf()) {
            h += outer(q,q).scale(1.0/qds);
        }
        else {
            printf("   SR1 not updating\n");
        }
    }


    void hessian_update_bfgs(const Tensor<double>& dx,
                             const Tensor<double>& dg)
    {
        /*
          Apply the BFGS update to the approximate Hessian h[][].

          h[][] = Hessian matrix from previous iteration
          dx[]  = Step from previous iteration
          .       (dx[] = x[] - xp[] where xp[] is the previous point)
          dg[]  = gradient difference (dg = g - gp)
        */

        Tensor<double> hdx  = inner(h,dx);

        double dxhdx = dx.trace(hdx);
        double dxdx  = dx.trace(dx);
        double dxdg  = dx.trace(dg);
        double dgdg  = dg.trace(dg);

        if ( (dxdx > 0.0) && (dgdg > 0.0) && (abs(dxdg/sqrt(dxdx*dgdg)) > 1.e-8) ) {
            for (int i=0; i<n; ++i) {
                for (int j=0; j<n; ++j) {
                    h(i,j) += dg[i]*dg[j]/dxdg - hdx[i]*hdx[j]/dxhdx;
                }
            }
        }
        else {
            printf("   BFGS not updating dxdg (%e), dgdg (%e), dxhdx (%f), dxdx(%e)\n" , dxdg, dgdg, dxhdx, dxdx);
        }
    }

    Tensor<double> new_search_direction(const Tensor<double>& g) {
        Tensor<double> dx, s;
        double tol = gradient_precision;
        double trust = 1.0; // This applied in spectral basis

        Tensor<double> v, e;
        syev(h, &v, &e);

        // Transform gradient into spectral basis
        Tensor<double> gv = inner(g,v);

        // Take step applying restriction
        for (int i=0; i<n; ++i) {
            if (e[i] < -tol) {
                printf("   forcing negative eigenvalue to be positive %d %.1e\n", i, e[i]);
                e[i] = -2.0*e[i]; // Enforce positive search direction
            }
            else if (e[i] < -tol) {
                printf("   forcing small eigenvalue to be positive %d %.1e\n", i, e[i]);
                e[i] = tol;
            }

            gv[i] = -gv[i] / e[i];
            if (abs(gv[i]) > trust) { // Step restriction
                double gvnew = trust*abs(gv(i))/gv[i];
                printf("   restricting step in spectral direction %d %.1e --> %.1e\n", i, gv[i], gvnew);
                gv[i] = gvnew;
            }
        }

        // Transform back from spectral basis
        return inner(v,gv);
    }

public:
    QuasiNewton(const std::shared_ptr<OptimizationTargetInterface>& target,
         double tol = 1e-6,
         double value_precision = 1e-12,
         double gradient_precision = 1e-12)
        : update("BFGS")
        , target(target)
        , tol(tol)
        , value_precision(value_precision)
        , gradient_precision(gradient_precision)
        , gnorm(tol*1e16)
        , n(0)
    {
        if (!target->provides_gradient()) throw "QuasiNewton requires the gradient";
    }

    void set_update(const string& method) {
        if (method == "BFGS" || method == "SR1") update=method;
        else throw "QuasiNewton: unknown update mthod";
    }

    bool optimize(Tensor<double>& x) {
        if (n != x.dim[0]) {
            n = x.dim[0];
            h = Tensor<double>();
        }

        target->test_gradient(x, value_precision);

        bool h_is_identity = (h.size == 0);
        if (h_is_identity) {
            h = Tensor<double>(n,n);
            for (int i=0; i<n; ++i) h(i,i) = 1.0;
        }

        Tensor<double> gp, dx;
        double fp;
        for (int iter=0; iter<20; ++iter) {
            Tensor<double> g;
            target->value_and_gradient(x, f, g);
            gnorm = g.normf();
            printf(" QuasiNewton iteration %2d value %.12f gradient %.2e\n",iter,f,gnorm);
            if (converged()) break;

            if (iter == 1 && h_is_identity) {
                // Default initial Hessian is scaled identity but
                // prefer to reuse any existing approximation.
                h.scale(g.trace(gp)/gp.trace(dx));
            }

            if (iter > 0) {
                if (update == "BFGS") hessian_update_bfgs(dx, g-gp);
                else hessian_update_sr1(dx, g-gp);
            }

            dx = new_search_direction(g);

            double step = line_search(1.0, f, dx.trace(g), x, dx);

            dx.scale(step);
            x += dx;
            gp = g;
            fp = f;

        }
        print("final hessian");
        print(h);
        return converged();
    }

    bool converged() const {return gnorm < tol;}

    double value() const {return f;}

    double gradient_norm() const {return gnorm;}
};


struct Test : public OptimizationTargetInterface {
    bool provides_gradient() const {return true;}

    double value(const Tensor<double>& x) {
        return 0.5*3.0*x.sumsq();
    }

    Tensor<double> gradient(const Tensor<double>& x) {
        return x*3.0;
    }
};


struct Test2 : public OptimizationTargetInterface {
    bool provides_gradient() const {return true;}

    double value(const Tensor<double>& x) {
        double v = 1.0;
        for (int i=0; i<x.dim[0]; ++i) {
            v *= cos((i+1)*x[i]);
        }
        return v;
    }

    Tensor<double> gradient(const Tensor<double>& x) {
        double v = value(x);
        Tensor<double> g(x.dim[0]);
        for (int i=0; i<x.dim[0]; ++i) {
            g[i]= -v*(i+1)*sin((i+1)*x[i])/cos((i+1)*x[i]);
        }
        return g;
    }
};



Tensor<double> op(const Tensor<double>& x) {
    const long n = x.dim[0];
    Tensor<double> f(n);
    for (long i=0; i<n; ++i) {
        f(i) = (i + 1)*x[i]; // + 0.01*i*x[i]*x[i]*x[i];
        for (long j=0; j<n; ++j)
            f(i) += 0.0001*i*j*x[i]*x[i]*x[j]*x[j]/((i+1)*(j+1));
    }
    return f;
}

double dot_product(const Tensor<double>& a, const Tensor<double>& b) {
    return a.trace(b);
}

int main() {

    Tensor<double> x(5);
    x.fillrandom();
    QuasiNewton solver(std::shared_ptr<OptimizationTargetInterface>(new Test2));
    solver.set_update("SR1");
    solver.optimize(x);
    return 0;


//     int n= 40;
//     int maxiter = 100;
//     int maxnvec = 20;
//     Tensor<double> f(maxnvec,n), x(maxnvec,n), Q(maxnvec,maxnvec);

//     int m = 0;
//     x(0,_).fillrandom();
//     for (int iter=0; iter<maxiter; ++iter) {
//         print("\nITERATION", iter, m);
//         f(m,_) = op(x(m,_));
//         print("x");
//         print(x(m,_));
//         print(f(m,_));

//         for (int j=0; j<=m; ++j) {
//             Q(j,m) = dot_product(x(j,_),f(m,_));
//             Q(m,j) = dot_product(x(m,_),f(j,_));
//         }

//         Tensor<double> c = KAIN(Q(Slice(0,m),Slice(0,m)));
//         print("THIS IS THE C I GOT");
//         print(c);

//         {
//             ++m;

//             Tensor<double> xnew(n);
//             for (int j=0; j<m; ++j) {
//                 xnew += c(j)*(x(j,_) - f(j,_));
//             }

//             double steplen = (xnew-x(m-1,_)).normf();
//             double xnorm = xnew.normf();
//             if (steplen > 0.5) {
//                 double damp = 0.5/steplen;
//                 if (damp < 0.1) damp = 0.1;
//                 print("restrictING", steplen, xnorm, damp);
//                 xnew = damp*xnew + (1.0-damp)*x(m-1,_);
//             }

//             if (m == maxnvec) {
//                 for (int i=1; i<m; ++i) {
//                     f(i-1,_) = f(i,_);
//                     x(i-1,_) = f(i,_);
//                 }
//                 Q(Slice(0,-2),Slice(0,-2)) = copy(Q(Slice(1,-1),Slice(1,-1)));

//                 m--;
//             }
//             x(m,_) = xnew;
//         }
//     }
//     return 0;

}
