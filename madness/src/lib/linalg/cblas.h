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

  
  $Id: cblas.h 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/

  
#ifndef MADNESS_LINALG_CBLAS_H__INCLUDED
#define MADNESS_LINALG_CBLAS_H__INCLUDED

#include <fortran_ctypes.h>

#ifdef FORTRAN_LINKAGE_LC
#  define sgemm_ sgemm
#  define dgemm_ dgemm
#  define cgemm_ cgemm
#  define zgemm_ zgemm
#else
  // only lowercase with zero and one underscores are handled -- if detected another convention complain loudly
#  ifndef FORTRAN_LINKAGE_LCU
#    error "cblas.h does not support the current Fortran symbol convention -- please, edit and check in the changes."
#  endif
#endif

extern "C" void dgemm_(const char *opa, const char *opb, const integer *m, const integer *n, const integer *k,
                       const real8 *alpha, const real8 *a, const integer *lda, const real8 *b, const integer *ldb,
                       const real8 *beta, real8 *c, const integer *ldc, char_len opalen, char_len opblen);

extern "C" void sgemm_(const char *opa, const char *opb, const integer *m, const integer *n, const integer *k,
                       const real4 *alpha, const real4 *a, const integer *lda, const real4 *b, const integer *ldb,
                       const real4 *beta, real4 *c, const integer *ldc, char_len opalen, char_len opblen);

extern "C" void zgemm_(const char *opa, const char *opb, const integer *m, const integer *n, const integer *k,
                       const complex_real8 *alpha,
                       const complex_real8 *a, const integer *lda, const complex_real8 *b, const integer *ldb,
                       const complex_real8 *beta, complex_real8 *c, const integer *ldc,  char_len opalen, char_len opblen);

extern "C" void cgemm_(const char *opa, const char *opb, const integer *m, const integer *n, const integer *k,
                       const complex_real4 *alpha,
                       const complex_real4 *a, const integer *lda, const complex_real4 *b, const integer *ldb,
                       const complex_real4 *beta, complex_real4 *c, const integer *ldc, char_len opalen, char_len opblen);

namespace madness {

    /// BLAS gemm exists only for float, double, complex<float> and complex<double>
    /// these types are defined in fortran_ctypes.h

    template <typename T> void gemm(bool transa, bool transb,
                                    integer m, integer n, integer k,
                                    T alpha, const T* a, integer lda, const T*b, integer ldb,
                                    T beta, T* c, integer ldc);
}

#endif // MADNESS_LINALG_CBLAS_H__INCLUDED

