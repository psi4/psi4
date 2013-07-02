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


  $Id: legendre.cc 2203 2011-03-04 15:09:41Z edward.valeev $
*/


#include <world/world.h>
#include <iostream>
using std::cout;
using std::endl;

#include <cmath>
#include <mra/legendre.h>
#include <tensor/tensor.h>

/// \file legendre.cc
/// \brief Legendre quadrature, polynomials and scaling functions.

// static void pn(int n, double x, double *p) {
// // old slow but simple recursion for polynomials
//   p[0] = 1;
//   if (n == 0) return;
//   p[1] = x;
//   for (int i=0; i<n; ++i) {
//     p[i+1] = i*(x*p[i] - p[i-1])/(i+1) + x*p[i];
//   }
// }


namespace madness {

    /// Evaluate the Legendre polynomials up to the given order at x in [-1,1].

    /// p should be an array of order+1 elements.
    void legendre_polynomials(double x, long order, double *p) {

        static double nn1[100];
        static long firstcall=1;
        long n;

        p[0] = 1.0;
        if (order == 0) return;

        if (firstcall) {
            for (n=0; n<100; ++n) nn1[n] = n/((double)(n+1));
            firstcall=0;
        }

        p[1] = x;
        for (n=1; n<order; ++n)
            p[n+1] = (x*p[n] - p[n-1])*nn1[n] + x*p[n];
    }

    /// Evaluate the first k Legendre scaling functions.

    /// p should be an array of k elements.
    void legendre_scaling_functions(double x, long k, double *p) {

        static double phi_norms[100];
        static long first_call = 1;
        long n;

        if (first_call) {
            for (n=0; n<100; ++n) phi_norms[n] = sqrt(2.0*n+1.0);
            first_call = 0;
        }

        legendre_polynomials(2.*x-1,k-1,p);
        for (n=0; n<k; ++n) {
            p[n] = p[n]*phi_norms[n];
        }
    }

    static bool data_is_read = false;
    static const int max_npt = 64;

    static const char *filename = "gaussleg";   // Is overridden by
    // These are the points and weights on [-1,1]
    static Tensor<double> points[max_npt+1];
    static Tensor<double> weights[max_npt+1];

    /// read_data loads the precomputed Gauss-Legendre data
    static bool read_data() {
        if (data_is_read) return true;
        FILE *f = fopen(filename,"r");
        if (!f) {
            cout << "legendre: read_data: could not find file " << filename << endl;
            return false;
        }
        for (int npt=0; npt<=max_npt; ++npt) {
            points[npt] = Tensor<double>(npt);
            weights[npt] = Tensor<double>(npt);

            int nnpt;
            if (fscanf(f,"%d",&nnpt) != 1) {
                cout << "legendre: read_data: failed reading " << npt << endl;
                fclose(f);
                return false;
            }
            if (nnpt != npt) {
                cout << "legendre: read_data: npt did not match " << npt << endl;
                fclose(f);
                return false;
            }
            for (int i=0; i<npt; ++i) {
                int ii;
                if (fscanf(f,"%d %lf %lf",&ii,&points[npt][i],&weights[npt][i]) != 3) {
                    cout << "legendre: read_data: failed reading data " << npt << " " << i << endl;
                    fclose(f);
                    return false;
                }
            }
        }
        fclose(f);
        data_is_read = true;
        return true;
    }

    /// Collective routine to pre-load and cache the quadrature points and weights

    /// Only process rank 0 will access the file.
    void load_quadrature(World& world, const char* dir) {
        if (data_is_read) return;
        if (world.rank() == 0) {
            char buf[32768];
            buf[0] = 0;
            strcat(buf,dir);
            strcat(buf,"/");
            strcat(buf,filename);
            filename = strdup(buf);
            if (!read_data()) throw "load_quadrature: failed reading quadrature coefficients";
        }
        else {
            for (int npt=0; npt<=max_npt; ++npt) {
                points[npt] = Tensor<double>(npt);
                weights[npt] = Tensor<double>(npt);
            }
        }
        for (int npt=1; npt<=max_npt; ++npt) {
            world.gop.broadcast(points[npt].ptr(), npt, 0);
            if (world.rank() == 99999999) world.gop.fence(); // Work aroung g++ 4.2.* bug on x86-64
            world.gop.broadcast(weights[npt].ptr(), npt, 0);
            if (world.rank() == 99999999) world.gop.fence(); // Work aroung g++ 4.2.* bug on x86-64
        }
        data_is_read = true;
    }

    /// Compute the Gauss-Legendre quadrature points and weights

    /// Return in x and w, which should be arrays of n elements, the
    /// points and weights for the n-point Gauss Legendre quadrature rule
    /// in [xlo,xhi].
    ///
    /// !!!! This routine delivers accurate points, but some weights are accurate
    /// only to about 1e-13 for higher order rules.  Need a more intelligent way
    /// to compute them.  See internal comments for more info.
    bool gauss_legendre_numeric(int n, double xlo, double xhi, double *x, double *w) {

        throw "gauss_legendre_numeric: why are we in here?";

#if 0
        double midpoint = (xhi + xlo)*0.5;
        double scale    = (xhi - xlo)*0.5;
        double acc = 1e-16;
        double p[100];

        double pi = atan(1.0)*4.0;

        // References made to the equation numbers in Davis & Rabinowitz 2nd ed.
        for (int k=0; k<n; ++k) {
            // Initial guess for the root using 2.7.3.3b
            // x = (1.0 - 1.0/(8*n*n) + 1.0/(8*n*n*n))*cos(((4*k+3.0)/(4*n+2))*pi)
            double r = (1.0 - 1.0/(8.0*n*n) + 1.0/(8.0*n*n*n))*cos(((4*k+3.0)/(4*n+2))*pi);
            // Refine the initial guess using a 3rd-order Newton correction 2.7.3.9
            for (int iter=0; iter<3; ++iter) {
                legendre_polynomials(r,n,p);
                double p0 = p[n];		// Value
                double p1 = n*(p[n-1] - r*p[n])/(1.0-r*r); // First derivative ... 2.7.3.5
                double p2 = (2*r*p1 - n*(n+1)*p0)/(1-r*r); // Second derivative
                double delta = - (p0/p1)*(1 + (p0*p2)/(2*p1*p1));
                r += delta;
                if (std::abs(delta) < acc) break;
                if (iter >= 2) return false;
            }
            legendre_polynomials(r,n,p);

            // This expression seems to have some cancellation of
            // significant digits for higher n.  Probably in the
            // 1-r*r for r close to 1 ... (1-r)*(1+r) gives us just
            // one more bit.  There are other expressions ...
            w[k] = 2.0*(1.0-r)*(1.0+r)/(n*n*p[n-1]*p[n-1]);
            x[k] = r;
        }

        for (int i=0; i<n; ++i) {
            x[i] = x[i]*scale + midpoint;
            w[i] = w[i]*scale;
        }
#endif

        return true;
    }

    /// Return precomputed (most accurate) or if not available computed
    /// (not quite as accurate) points and weights for Gauss Legendre
    /// quadrature rule on the specified interval.
    bool gauss_legendre(int n, double xlo, double xhi, double *x, double *w) {
        if (!read_data())
            return false;

        if (n < 1)
            return false;

        if (n > max_npt)
            return gauss_legendre_numeric(n, xlo, xhi, x, w);

        // Cached are the points and weights for the interval [0,1]

        // int(f(x),x=xlo..xhi) = int(f(x(z)),z=0..1)*(xhi-xlo)
        // z = (x-xlo)/(xhi-xlo)  ->  x = xlo + z*(xhi-xlo)

        double range = xhi - xlo;
        for (int i=0; i<n; ++i) {
            x[i] = xlo + points[n][i] * range;
            w[i] = weights[n][i] * range;
        }

        return true;
    }

    static double testf(int n, double x) {
        /// test function for gauss_legendre_test
        double sum = 0.0;
        double xn = 1.0;
        for (int i=0; i<=(2*n-1); ++i) {
            sum += xn;
            xn *= x;
        }
        return sum;
    }

    bool gauss_legendre_test(bool print) {
        /// Check error in numerical integ(x+x^2+...+x^(2n-1),x=0..1)
        /// using n-pt rule.

        const int maxnpt = 64;
        double x[maxnpt], w[maxnpt];

        for (int npt=1; npt<maxnpt; ++npt) {
            double sum = 0.0;
            gauss_legendre(npt,0,1,x,w);
            for (int i=0; i<npt; ++i) sum += testf(npt, x[i])*w[i];
            for (int i=0; i<=(2*npt-1); ++i) sum -= 1.0/(i+1);
            bool err = (std::abs(sum/npt) > 1.3e-14);
            if (err || print)
                cout << "gauss_leg_test: " << npt << " " << sum << " " << sum/npt << endl;
            if (err) return false;
        }

        return true;
    }

    // int main() {
    //   if (gauss_legendre_test())
    //     cout << " gauss_legendre seems OK\n";
    //   else
    //     cout << " gauss_legendre failed !\n";
    //   return 0;
    // }

}
