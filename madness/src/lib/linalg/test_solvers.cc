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

using namespace madness;
using namespace std;

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
        for (int i=0; i<x.dim(0); ++i) {
            v *= cos((i+1)*x[i]);
        }
        return v;
    }

    Tensor<double> gradient(const Tensor<double>& x) {
        double v = value(x);
        Tensor<double> g(x.dim(0));
        for (int i=0; i<x.dim(0); ++i) {
            g[i]= -v*(i+1)*sin((i+1)*x[i])/cos((i+1)*x[i]);
        }
        return g;
    }
};



Tensor<double> op(const Tensor<double>& x) {
    const long n = x.dim(0);
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
