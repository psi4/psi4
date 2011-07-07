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

  
  $Id: cblas.cc 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/

  
#include <cstdio>
#include <complex>

#include <linalg/cblas.h>

typedef std::complex<float> complex_float;
typedef std::complex<double> complex_double;

extern "C" void xerbla_(char *message, integer *info, int length) {
    std::fprintf(stderr,
                 " ** On entry to  %6s, parameter number %2ld had an illegal value\n", message,
                 (long) *info);
    throw "XERBLA";
}

/// \file cblas.cc
/// \brief This file provides gemm template BLAS.

namespace madness {


    template <> void gemm<double> (bool transa, bool transb,
                                   integer m, integer n, integer k,
                                   double alpha, const double* a, integer lda,
                                   const double* b, integer ldb,
                                   double beta, double* c, integer ldc) {
        const char *op[] = {"n","t"
                           };
        dgemm_(op[transa], op[transb], &m, &n, &k,
               &alpha, a, &lda, b, &ldb,
               &beta,  c, &ldc,
               1, 1);
    }

    template <> void gemm<float> (bool transa, bool transb,
                                  integer m, integer n, integer k,
                                  float alpha, const float* a, integer lda,
                                  const float* b, integer ldb,
                                  float beta, float* c, integer ldc) {
        const char *op[] = {"n","t"
                           };
        sgemm_(op[transa], op[transb], &m, &n, &k,
               &alpha, a, &lda, b, &ldb,
               &beta,  c, &ldc,
               1, 1);
    }

    // complex_double == complex_real8
    template <> void gemm<complex_double> (bool transa, bool transb,
                                           integer m, integer n, integer k,
                                           complex_double alpha, const complex_double* a, integer lda,
                                           const complex_double* b, integer ldb,
                                           complex_double beta, complex_double* c, integer ldc) {
        const char *op[] = {"n","t"
                           };
        zgemm_(op[transa], op[transb], &m, &n, &k,
               &alpha, a, &lda, b, &ldb,
               &beta,  c, &ldc,
               1, 1);
    }

    // complex_float == complex_real4
    template <> void gemm<complex_float> (bool transa, bool transb,
                                          integer m, integer n, integer k,
                                          complex_float alpha, const complex_float* a, integer lda,
                                          const complex_float* b, integer ldb,
                                          complex_float beta, complex_float* c, integer ldc) {
        const char *op[] = {"n","t"
                           };
        cgemm_(op[transa], op[transb], &m, &n, &k,
               &alpha, a, &lda, b, &ldb,
               &beta,  c, &ldc,
               1, 1);
    }

}


