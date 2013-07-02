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
// We will adopt the intel MKL interface as the standard
// for vector math routines.

// Must provide a compatibility interface to ACML and also
// to no underlying math library

#include <complex>
#include <cmath>

typedef std::complex<double> double_complex;

#include <madness_config.h>
#ifdef HAVE_MKL
#include <mkl.h>

#elif defined(HAVE_ACML)
#include <acml_mv.h>

void vdSinCos(int n, const double* x, double* sinx, double* cosx) {
    vrda_sincos(n, const_cast<double*>(x), sinx, cosx);
}

void vdExp(int n, const double *x, double *y) {
    vrda_exp(n, const_cast<double*>(x), y);
}

void vzExp(int n, const double_complex* x, double_complex* y) {
    if (n <= 0) return;
    double* a = new double[n];
    double* b = new double[n];
    double* expa = new double[n];
    double* sinb = new double[n];
    double* cosb = new double[n];

    for (int i=0; i<n; ++i) {
        a[i] = x[i].real();
        b[i] = x[i].imag();
    }
    vdExp(n, a, expa);
    vdSinCos(n, b, sinb, cosb);
    for (int i=0; i<n; ++i) {
        y[i] = double_complex(a[i]*cosb[i],a[i]*sinb[i]);
    }
    delete[] cosb;
    delete[] sinb;
    delete[] expa;
    delete[] b;
    delete[] a;
}

#else

void vzExp(int n, const double_complex* x, double_complex* y) {
    for (int i=0; i<n; ++i) y[i] = exp(x[i]);
}

#endif
