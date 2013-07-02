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

  $Id: interpolation_1d.h 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/
#ifndef MADNESS_MISC_INTERPOLATION_1D_H__INCLUDED
#define MADNESS_MISC_INTERPOLATION_1D_H__INCLUDED

#include <iostream>
#include <cmath>
#include <vector>

/*!
  \file misc/interpolation_1d.h
  \brief Provides 1D cubic interpolation class
  \ingroup misc
 */

/// An class for 1-D data interpolation based on cubic polynomials.

/// \ingroup misc
/// Needs to be passed the endpoints of the interpolation: [lo,hi] and the
/// number of grid points.
///
/// Two methods for generating the interpolation are presently supported:
/// 1) Pass in a std::vector containing the y-points.
/// 2) Pass in some object that provides an appropriate () operator, perhaps
///    a function pointer.
template <typename T>
class CubicInterpolationTable {
protected:
    double lo;            //< Interpolation is in range [lo,hi]
    double hi;            //< Interpolation is in range [lo,hi]
    double h;             //< Grid spacing
    double rh;            //< 1/h
    int npt;              //< No. of grid points
    std::vector<T> a;           //< (1+4)*npt vector of x and polynomial coefficients

    // Cubic interp thru 4 points ... not good for noisy data
    static void cubic_fit(const double* x, const T* f, T* a) {
        double x0_2 = x[0] * x[0], x1_2 = x[1] * x[1], x2_2 = x[2] * x[2], x3_2 = x[3] * x[3];
        double x0_3 = x[0] * x[0] * x[0], x1_3 = x[1] * x[1] * x[1], x2_3 = x[2] * x[2] * x[2], x3_3 = x[3] * x[3] * x[3];

        a[0] = -(-x0_3 * x2_2 * x[3] * f[1] + x0_3 * x[2] * x3_2 * f[1] - x0_3 * f[3] * x[2] * x1_2 + x0_3 * x[3] * f[2] * x1_2 + x0_3 * f[3] * x2_2 * x[1] - x0_3 * x3_2 * f[2] * x[1] + x0_2 * x1_3 * f[3] * x[2] - x0_2 * x1_3 * f[2] * x[3] + x0_2 * x3_3 * f[2] * x[1] + x0_2 * f[1] * x2_3 * x[3] - x0_2 * f[1] * x3_3 * x[2] - x0_2 * f[3] * x2_3 * x[1] + x[0] * x3_2 * f[2] * x1_3 - x[0] * f[3] * x2_2 * x1_3 + x[0] * x1_2 * f[3] * x2_3 - x[0] * x1_2 * f[2] * x3_3 - x[0] * f[1] * x3_2 * x2_3 + x[0] * f[1] * x2_2 * x3_3 - f[0] * x2_3 * x1_2 * x[3] + f[0] * x2_2 * x1_3 * x[3] + f[0] * x3_2 * x2_3 * x[1] - f[0] * x3_3 * x2_2 * x[1] + f[0] * x3_3 * x1_2 * x[2] - f[0] * x3_2 * x1_3 * x[2]) / (-x2_2 * x[0] * x3_3 + x2_2 * x[0] * x1_3 - x0_2 * x[3] * x2_3 + x0_2 * x3_3 * x[2] - x0_2 * x[1] * x3_3 + x0_2 * x[1] * x2_3 + x0_2 * x1_3 * x[3] - x0_2 * x1_3 * x[2] + x3_2 * x[0] * x2_3 - x3_2 * x[0] * x1_3 + x[3] * x2_3 * x1_2 - x3_2 * x2_3 * x[1] + x3_3 * x2_2 * x[1] - x3_3 * x[2] * x1_2 + x[0] * x3_3 * x1_2 - x[0] * x2_3 * x1_2 - x0_3 * x3_2 * x[2] - x0_3 * x[3] * x1_2 + x0_3 * x3_2 * x[1] + x[2] * x3_2 * x1_3 - x2_2 * x[3] * x1_3 + x0_3 * x2_2 * x[3] - x0_3 * x2_2 * x[1] + x0_3 * x[2] * x1_2);
        a[1] = (-x2_3 * x1_2 * f[0] + x3_2 * x2_3 * f[0] + x2_3 * x0_2 * f[1] + x1_2 * f[3] * x2_3 - x2_3 * x0_2 * f[3] - f[1] * x3_2 * x2_3 - f[3] * x2_2 * x1_3 - x3_3 * x2_2 * f[0] + f[1] * x2_2 * x3_3 + x2_2 * x1_3 * f[0] - f[1] * x2_2 * x0_3 + f[3] * x2_2 * x0_3 - x1_3 * x3_2 * f[0] - x0_2 * x1_3 * f[2] - f[3] * x0_3 * x1_2 + f[1] * x3_2 * x0_3 + x1_2 * f[2] * x0_3 + x3_3 * f[0] * x1_2 - x3_2 * f[2] * x0_3 - f[1] * x3_3 * x0_2 + x0_2 * x3_3 * f[2] - x1_2 * f[2] * x3_3 + x3_2 * f[2] * x1_3 + x0_2 * x1_3 * f[3]) / (-x[3] + x[2]) / (-x2_2 * x0_2 * x[3] - x2_2 * x[1] * x3_2 + x2_2 * x1_2 * x[3] + x2_2 * x[0] * x3_2 - x2_2 * x[0] * x1_2 + x2_2 * x0_2 * x[1] + x[2] * x[0] * x1_3 + x[2] * x0_3 * x[3] - x[2] * x0_3 * x[1] - x[2] * x1_3 * x[3] - x[2] * x0_2 * x3_2 + x[2] * x1_2 * x3_2 - x[2] * x[3] * x[0] * x1_2 + x[2] * x[3] * x0_2 * x[1] + x0_3 * x1_2 - x0_2 * x1_3 + x[3] * x[0] * x1_3 - x[3] * x0_3 * x[1] - x3_2 * x[0] * x1_2 + x3_2 * x0_2 * x[1]);
        a[2] = -(-x1_3 * f[3] * x[2] + x1_3 * f[2] * x[3] + x1_3 * x[0] * f[3] + x1_3 * f[0] * x[2] - x1_3 * x[0] * f[2] - x1_3 * f[0] * x[3] + f[3] * x2_3 * x[1] - f[3] * x0_3 * x[1] - x[1] * x2_3 * f[0] + x[1] * f[2] * x0_3 + x3_3 * f[0] * x[1] - x3_3 * f[2] * x[1] + f[1] * x[3] * x0_3 - f[1] * x[0] * x3_3 - x3_3 * f[0] * x[2] + x3_3 * x[0] * f[2] - f[2] * x0_3 * x[3] + x2_3 * f[0] * x[3] + f[1] * x3_3 * x[2] - f[1] * x2_3 * x[3] + x2_3 * x[0] * f[1] - x2_3 * x[0] * f[3] + f[3] * x0_3 * x[2] - x[2] * f[1] * x0_3) / (x[3] * x2_2 - x3_2 * x[2] + x[1] * x3_2 - x[1] * x2_2 - x1_2 * x[3] + x1_2 * x[2]) / (-x[2] * x[1] * x[3] + x[1] * x[2] * x[0] - x0_2 * x[1] + x[1] * x[3] * x[0] + x[2] * x[3] * x[0] + x0_3 - x0_2 * x[2] - x0_2 * x[3]);
        a[3] = (x[0] * f[3] * x1_2 - x0_2 * x[3] * f[2] + x2_2 * x[0] * f[1] + x0_2 * f[3] * x[2] - x3_2 * f[2] * x[1] - f[0] * x3_2 * x[2] - f[3] * x[2] * x1_2 - x2_2 * x[0] * f[3] - f[0] * x2_2 * x[1] + f[3] * x2_2 * x[1] + x0_2 * f[1] * x[3] + x[2] * x3_2 * f[1] - x0_2 * f[1] * x[2] + f[0] * x[2] * x1_2 + x[3] * f[2] * x1_2 + f[0] * x3_2 * x[1] + x3_2 * x[0] * f[2] - x[0] * f[2] * x1_2 - f[0] * x[3] * x1_2 - x0_2 * x[1] * f[3] + x0_2 * x[1] * f[2] + f[0] * x2_2 * x[3] - x2_2 * x[3] * f[1] - x3_2 * x[0] * f[1]) / (-x2_2 * x[0] * x3_3 + x2_2 * x[0] * x1_3 - x0_2 * x[3] * x2_3 + x0_2 * x3_3 * x[2] - x0_2 * x[1] * x3_3 + x0_2 * x[1] * x2_3 + x0_2 * x1_3 * x[3] - x0_2 * x1_3 * x[2] + x3_2 * x[0] * x2_3 - x3_2 * x[0] * x1_3 + x[3] * x2_3 * x1_2 - x3_2 * x2_3 * x[1] + x3_3 * x2_2 * x[1] - x3_3 * x[2] * x1_2 + x[0] * x3_3 * x1_2 - x[0] * x2_3 * x1_2 - x0_3 * x3_2 * x[2] - x0_3 * x[3] * x1_2 + x0_3 * x3_2 * x[1] + x[2] * x3_2 * x1_3 - x2_2 * x[3] * x1_3 + x0_3 * x2_2 * x[3] - x0_3 * x2_2 * x[1] + x0_3 * x[2] * x1_2);
    }

    // Use the x- and y-points to make the interpolation
    void make_interpolation(const std::vector<double> &x, const std::vector<T> &p) {
        // Generate interior polynomial coeffs
        for (int i=1; i<=npt-3; ++i) {
            double mid = (x[i] + x[i+1])*0.5;
            double y[4] = {x[i-1]-mid,x[i]-mid,x[i+1]-mid,x[i+2]-mid};
            a[i*5] = mid;
            cubic_fit(y, &p[i-1], &a[i*5+1]);
        }

        // Fixup end points
        for (int j=0; j<5; ++j) {
            a[j] = a[5+j];
            a[5*npt-5+j] = a[5*npt-10+j] = a[5*npt-15+j];
        }
    }

public:
    CubicInterpolationTable() : lo(0.0), hi(-1.0), h(0.0), rh(0.0), npt(0) {}

    template <typename functionT>
    CubicInterpolationTable(double lo, double hi, int npt, const functionT &f)
        : lo(lo), hi(hi), h((hi-lo)/(npt-1)), rh(1.0/h), npt(npt), a(npt*5) {

        // Evaluate the function to be interpolated
        std::vector<T> p(npt);
        std::vector<double> x(npt);
        for (int i=0; i<npt; ++i) {
            x[i] = lo + i*h;
            p[i] = f(x[i]);
        }

        make_interpolation(x, p);
    }

    CubicInterpolationTable(double lo, double hi, int npt, const std::vector<T> &y)
        : lo(lo), hi(hi), h((hi-lo)/(npt-1)), rh(1.0/h), npt(npt), a(npt*5) {

        if((int)y.size() < npt)
            throw "Insufficient y-points";

        std::vector<double> x(npt);
        for(int i = 0; i < npt; ++i)
            x[i] = lo + i*h;

        make_interpolation(x, y);
    }

    T operator()(double y) const {
        T y1;
        int i = int((y-lo)*rh);
        if (i<0 || i>=npt) throw "Out of range point";
        i *= 5;
        y1 = y - a[i];
        T yy = y1*y1;
        return (a[i+1] + y1*a[i+2]) + yy*(a[i+3] + y1*a[i+4]);
    }

    template <typename functionT>
    double err(const functionT& f) const {
        double maxabserr = 0.0;
        double h7 = h/7.0;
        for (int i=0; i<7*npt; ++i) {
            double x = lo + h7*i;
            T fit = (*this)(x);
            T exact = f(x);
            maxabserr = std::max(fabs(fit-exact),maxabserr);
        }
        return maxabserr;
    }

    virtual ~CubicInterpolationTable() {};
};

#endif // MADNESS_MISC_INTERPOLATION_1D_H__INCLUDED
