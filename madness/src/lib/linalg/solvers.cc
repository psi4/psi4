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
#include <linalg/solvers.h>

namespace madness {


    double OptimizationTargetInterface::test_gradient(Tensor<double>& x, double value_precision, bool doprint) {
        const double eps = pow(value_precision,0.3333);
        if (doprint) {
            printf("\n");
            printf("Testing gradient eps=%.1e\n----------------\n", eps);
            printf("  i            f-                 f+               ganalytic           gnumeric         err\n");
            printf(" ----  ------------------  ------------------  ------------------  ------------------  -------\n");
        }
        Tensor<double> tt = gradient(x);
        int n = int(tt.dim(0));
        double maxerr = 0.0;
        for (int i=0; i<n; ++i) {
            x[i] += eps;
            double fp = value(x);
            x[i] -= 2.0*eps;
            double fm = value(x);
            x[i] += eps;

            double gg = 0.5*(fp-fm)/eps;
            if (doprint)
                printf("% 5d%20.12e%20.12e%20.12e%20.12e  %.1e\n", i, fm, fp, tt(i), gg, std::abs(tt(i)-gg));
            maxerr = std::max(std::abs(gg-tt(i)),maxerr);
        }
        if (doprint) printf("\n");

        return maxerr;
    }


    SteepestDescent::SteepestDescent(const std::shared_ptr<OptimizationTargetInterface>& tar,
                                     double tol,
                                     double value_precision,
                                     double gradient_precision)
        : target(tar)
        , tol(tol)
        , value_precision(value_precision)
        , gradient_precision(gradient_precision)
        , gnorm(tol*1e16)
    {
        if (!target->provides_gradient()) throw "Steepest descent requires the gradient";
    }

    bool SteepestDescent::optimize(Tensor<double>& x) {
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

    bool SteepestDescent::converged() const {return gnorm < tol;}

    double SteepestDescent::gradient_norm() const {return gnorm;}

    double SteepestDescent::value() const {return f;}

    double QuasiNewton::line_search(double a1, double f0, double dxgrad, const Tensor<double>& x, const Tensor<double>& dx) {
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

        if (std::abs(f1-f0) < value_precision) { // Insufficient precision
            a2 = a1;
            lsmode = "fixed";
        }
        else if (hess > 0.0) { // Positive curvature
            if ((f1 - f0) <= -value_precision) { // Downhill
                lsmode = "downhill";
                if (std::abs(a2) > 4.0*std::abs(a1)) {
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

    void QuasiNewton::hessian_update_sr1(const Tensor<double>& s, const Tensor<double>& y) {
        Tensor<double> q = y - inner(h,s);
        double qds = q.trace(s);
        if (std::abs(qds) > 1e-8 * s.normf() * q.normf()) {
            h += outer(q,q).scale(1.0/qds);
        }
        else {
            printf("   SR1 not updating\n");
        }
    }


    void QuasiNewton::hessian_update_bfgs(const Tensor<double>& dx,
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

        if ( (dxdx > 0.0) && (dgdg > 0.0) && (std::abs(dxdg/std::sqrt(dxdx*dgdg)) > 1.e-8) ) {
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

    Tensor<double> QuasiNewton::new_search_direction(const Tensor<double>& g) {
        Tensor<double> dx, s;
        double tol = gradient_precision;
        double trust = 1.0; // This applied in spectral basis

        Tensor<double> v, e;
        syev(h, v, e);

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
            if (std::abs(gv[i]) > trust) { // Step restriction
                double gvnew = trust*std::abs(gv(i))/gv[i];
                printf("   restricting step in spectral direction %d %.1e --> %.1e\n", i, gv[i], gvnew);
                gv[i] = gvnew;
            }
        }

        // Transform back from spectral basis
        return inner(v,gv);
    }

    QuasiNewton::QuasiNewton(const std::shared_ptr<OptimizationTargetInterface>& tar,
                             double tol,
                             double value_precision,
                             double gradient_precision)
        : update("BFGS")
        , target(tar)
        , tol(tol)
        , value_precision(value_precision)
        , gradient_precision(gradient_precision)
        , gnorm(tol*1e16)
        , n(0)
    {
        if (!target->provides_gradient()) throw "QuasiNewton requires the gradient";
    }

    void QuasiNewton::set_update(const std::string& method) {
        if (method == "BFGS" || method == "SR1") update=method;
        else throw "QuasiNewton: unknown update mthod";
    }

    bool QuasiNewton::optimize(Tensor<double>& x) {
        int maxiter = 20; // !!!!!!!!! dumb
        if (n != x.dim(0)) {
            n = x.dim(0);
            h = Tensor<double>();
        }

        target->test_gradient(x, value_precision);

        bool h_is_identity = (h.size() == 0);
        if (h_is_identity) {
            h = Tensor<double>(n,n);
            for (int i=0; i<n; ++i) h(i,i) = 1.0;
        }

        Tensor<double> gp, dx;
        double fp;
        for (int iter=0; iter<maxiter; ++iter) {
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

        bool QuasiNewton::converged() const {return gnorm < tol;}

        double QuasiNewton::value() const {return f;}

        double QuasiNewton::gradient_norm() const {return gnorm;}
}
