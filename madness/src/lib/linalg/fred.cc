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
#include <cmath>
#include <tensor/tensor.h>
#include <linalg/tensor_lapack.h>
#include <world/print.h>

using namespace madness;
using namespace std;

void drot(long n, double* restrict a, double* restrict b, double s, double c) {
    for (long i=0; i<n; ++i) {
        double aa = a[i]*c + b[i]*s;
        double bb = b[i]*c - a[i]*s;
        a[i] = aa;
        b[i] = bb;
    }
}

int main() {
    const double pi = 3.14; // No need for precision
    const double thetamax = pi*0.2;

    const long n = 400;
    Tensor<double> d(n);
    d.fillrandom();
    ITERATOR(d, d(IND) -= 0.5);
    Tensor<double> dorig = copy(d);

    Tensor<double> U(n,n);
    for (long i=0; i<n; ++i) U(i,i) = 1.0;

    cout.precision(8);

    const double thresh = 1e-9; // Final convergence test
    double tol = 0.1;           // Current test
    long ndone=0;               // Total no. of rotations performed

    for (int iter=0; iter<300; ++iter) {
        double sum=d.sum();  // Current sum
        print("iter", iter, sum);

        // Screening not only reduces work (in same manner as it does
        // for Jacobi diagonalization) it also dramatically increases
        // the rate of convergence
        tol *= 0.25;
        if (tol < thresh) tol = thresh;

        long ndone_iter=0; // No. of rotations done this iteration
        double screen = tol*sum;
        for (long i=0; i<n; ++i) {
            for (long j=0; j<i; ++j) {

                // construct rotation
                // |inew>  =  cos theta |i>  +  sin theta |j>
                // |jnew>  =  cos theta |j>  -  sin theta |i>
                //
                double numerator = d[j] - d[i];
                double denominator = d[j] + d[i];
                double theta;

                if (fabs(denominator) < 1e-2*fabs(numerator)) {
                    // Need a nearly pi/2 rotation in the direction of increase
                    // Just take a step hopefully in the right direction.
                    //print("     small", i, j);
                    theta = thetamax*0.7;
                    if (numerator < 0.0) theta = -theta;
                }
                else {
                    theta = atan(numerator/fabs(denominator));
                }

                if (fabs(theta) > screen) {
                    if (fabs(theta) > thetamax) {
                        //print("     restricting", i, j, theta);
                        if (theta > 0.0) theta = thetamax;
                        else theta = -thetamax;
                    }
                    //print("    rotn", i, j, theta);

                    double c = cos(theta);
                    double s = sin(theta);

                    double dbefore = d.sum();
                    double di = d[i]*c + d[j]*s;
                    double dj = d[j]*c - d[i]*s;
                    d[i] = di;
                    d[j] = dj;
                    double dafter = d.sum();
                    if (dbefore-dafter > 1e-12*fabs(dbefore))
                        print("BAD BAD", i, j, dbefore, dafter);

                    // Accumulate the transpose of U for efficient inner loop
                    //                 for (long k=0; k<n; ++k) {
                    //                     double uik = U(i,k)*c + U(j,k)*s;
                    //                     double ujk = U(j,k)*c - U(i,k)*s;
                    //                     U(i,k) = uik;
                    //                     U(j,k) = ujk;
                    //                 }

                    ++ndone_iter;
                    drot(n, &U(i,0), &U(j,0), s, c);
                }
            }
        }
        print(ndone_iter);
        ndone += ndone_iter;
        if (ndone_iter == 0 && tol==thresh) {
            print("CONVERGED");
            break;
        }
    }

    U = transpose(U);
    double err = (inner(dorig,U)-d).normf();
    print("Error ", err);
    print("Ndone ", ndone);

    return 0;
}
