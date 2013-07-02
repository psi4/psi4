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


  $Id: mxm.h 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/

#ifndef MADNESS_TENSOR_MXM_H__INCLUDED
#define MADNESS_TENSOR_MXM_H__INCLUDED

/// \file tensor/mxm.h
/// \brief Internal use only

// This file is ONLY included into tensor.cc ... separated here just
// to shrink file size.  Don't try to include anywhere else

// Due to both flakey compilers and performance concerns,
// we use a simple reference implementation of the mxm
// routines for all except T=double.


// Here undef _CRAY since no longer care about Cray X1
// and want to ignore it for the XT3/4
#ifdef _CRAY
#undef _CRAY
#endif


/// Matrix * matrix reference implementation (slow but correct)
template <typename T, typename Q, typename S>
static inline void mxm(long dimi, long dimj, long dimk,
                       T* restrict c, const Q* restrict a,
                       const S* restrict b) {
    /*
      c(i,j) = c(i,j) + sum(k) a(i,k)*b(k,j)

      where it is assumed that the last index in each array is has unit
      stride and the dimensions are as provided.
    */

    for (long i=0; i<dimi; ++i) {
        for (long k=0; k<dimk; ++k) {
#ifdef _CRAY
#pragma _CRI prefervector
#endif
            for (long j=0; j<dimj; ++j) {
                c[i*dimj+j] += a[i*dimk+k]*b[k*dimj+j];
            }
        }
    }
}


/// Matrix transpose * matrix ... reference implementation (slow but correct)
template <typename T, typename Q, typename S>
static inline
void mTxm(long dimi, long dimj, long dimk,
          T* restrict c, const Q* restrict a,
          const S* restrict b) {
    /*
      c(i,j) = c(i,j) + sum(k) a(k,i)*b(k,j)

      where it is assumed that the last index in each array is has unit
      stride and the dimensions are as provided.

      i loop might be long in anticipated application
    */

    for (long k=0; k<dimk; ++k) {
        for (long j=0; j<dimj; ++j) {
#ifdef _CRAY
#pragma _CRI prefervector
#endif
            for (long i=0; i<dimi; ++i) {
                c[i*dimj+j] += a[k*dimi+i]*b[k*dimj+j];
            }
        }
    }
}

/// Matrix * matrix transpose ... reference implementation (slow but correct)
template <typename T, typename Q, typename S>
static inline void mxmT(long dimi, long dimj, long dimk,
                        T* restrict c, const Q* restrict a,
                        const S* restrict b) {
    /*
      c(i,j) = c(i,j) + sum(k) a(i,k)*b(j,k)

      where it is assumed that the last index in each array is has unit
      stride and the dimensions are as provided.

      i loop might be long in anticipated application
    */

    for (long i=0; i<dimi; ++i) {
        for (long j=0; j<dimj; ++j) {
            T sum = 0;
            for (long k=0; k<dimk; ++k) {
                sum += a[i*dimk+k]*b[j*dimk+k];
            }
            c[i*dimj+j] += sum;
        }
    }
}

/// Matrix transpose * matrix transpose reference implementation (slow but correct)
template <typename T, typename Q, typename S>
static inline void mTxmT(long dimi, long dimj, long dimk,
                         T* restrict c, const Q* restrict a,
                         const S* restrict b) {
    /*
      c(i,j) = c(i,j) + sum(k) a(k,i)*b(j,k)

      where it is assumed that the last index in each array is has unit
      stride and the dimensions are as provided.
    */

    for (long i=0; i<dimi; ++i) {
#ifdef _CRAY
#pragma _CRI prefervector
#endif
        for (long j=0; j<dimj; ++j) {
            for (long k=0; k<dimk; ++k) {
                c[i*dimj+j] += a[k*dimi+i]*b[j*dimk+k];
            }
        }
    }
}

#ifdef _CRAY
// Simple loop structure best on the Cray X1
template
void mTxm(long dimi, long dimj, long dimk, double* restrict c,
          const double* restrict a, const double* restrict b);
#else
// The following are restricted to double only

/// Matrix transpose * matrix (hand unrolled version)

template <>
inline void mTxm(long dimi, long dimj, long dimk,
                 double* restrict c, const double* restrict a,
                 const double* restrict b) {
    /*
    c(i,j) = c(i,j) + sum(k) a(k,i)*b(k,j)  <--- NOTE ACCUMULATION INTO C

    where it is assumed that the last index in each array is has unit
    stride and the dimensions are as provided.

    i loop might be long in anticipated application

    4-way unrolled k loop ... empirically fastest on PIII
    compared to 2/3 way unrolling (though not by much).
    */

    long dimk4 = (dimk/4)*4;
    for (long i=0; i<dimi; ++i,c+=dimj) {
        const double* ai = a+i;
        const double* p = b;
        for (long k=0; k<dimk4; k+=4,ai+=4*dimi,p+=4*dimj) {
            double ak0i = ai[0   ];
            double ak1i = ai[dimi];
            double ak2i = ai[dimi+dimi];
            double ak3i = ai[dimi+dimi+dimi];
            const double* bk0 = p;
            const double* bk1 = p+dimj;
            const double* bk2 = p+dimj+dimj;
            const double* bk3 = p+dimj+dimj+dimj;
            for (long j=0; j<dimj; ++j) {
                c[j] += ak0i*bk0[j] + ak1i*bk1[j] + ak2i*bk2[j] + ak3i*bk3[j];
            }
        }
        for (long k=dimk4; k<dimk; ++k) {
            double aki = a[k*dimi+i];
            const double* bk = b+k*dimj;
            for (long j=0; j<dimj; ++j) {
                c[j] += aki*bk[j];
            }
        }
    }
}



/// Matrix * matrix transpose (hand unrolled version)

template <>
inline void mxmT(long dimi, long dimj, long dimk,
                 double* restrict c,
                 const double* restrict a, const double* restrict b) {
    /*
    c(i,j) = c(i,j) + sum(k) a(i,k)*b(j,k)

    where it is assumed that the last index in each array is has unit
    stride and the dimensions are as provided.

    j loop might be long in anticipated application

    Unrolled i loop.  Empirically fastest on PIII compared
    to unrolling j, or both i&j.
    */

    long dimi2 = (dimi/2)*2;
    for (long i=0; i<dimi2; i+=2) {
        const double* ai0 = a+i*dimk;
        const double* ai1 = a+i*dimk+dimk;
        double* restrict ci0 = c+i*dimj;
        double* restrict ci1 = c+i*dimj+dimj;
        for (long j=0; j<dimj; ++j) {
            double sum0 = 0;
            double sum1 = 0;
            const double* bj = b + j*dimk;
            for (long k=0; k<dimk; ++k) {
                sum0 += ai0[k]*bj[k];
                sum1 += ai1[k]*bj[k];
            }
            ci0[j] += sum0;
            ci1[j] += sum1;
        }
    }
    for (long i=dimi2; i<dimi; ++i) {
        const double* ai = a+i*dimk;
        double* restrict ci = c+i*dimj;
        for (long j=0; j<dimj; ++j) {
            double sum = 0;
            const double* bj = b+j*dimk;
            for (long k=0; k<dimk; ++k) {
                sum += ai[k]*bj[k];
            }
            ci[j] += sum;
        }
    }
}

/// Matrix * matrix (hand unrolled version)
template <>
inline void mxm(long dimi, long dimj, long dimk,
                double* restrict c, const double* restrict a, const double* restrict b) {
    /*
    c(i,j) = c(i,j) + sum(k) a(i,k)*b(k,j)

    where it is assumed that the last index in each array is has unit
    stride and the dimensions are as provided.

    4-way unrolled k loop ... empirically fastest on PIII
    compared to 2/3 way unrolling (though not by much).
    */

    long dimk4 = (dimk/4)*4;
    for (long i=0; i<dimi; ++i, c+=dimj,a+=dimk) {
        const double* p = b;
        for (long k=0; k<dimk4; k+=4,p+=4*dimj) {
            double aik0 = a[k  ];
            double aik1 = a[k+1];
            double aik2 = a[k+2];
            double aik3 = a[k+3];
            const double* bk0 = p;
            const double* bk1 = bk0+dimj;
            const double* bk2 = bk1+dimj;
            const double* bk3 = bk2+dimj;
            for (long j=0; j<dimj; ++j) {
                c[j] += aik0*bk0[j] + aik1*bk1[j] + aik2*bk2[j] + aik3*bk3[j];
            }
        }
        for (long k=dimk4; k<dimk; ++k) {
            double aik = a[k];
            for (long j=0; j<dimj; ++j) {
                c[j] += aik*b[k*dimj+j];
            }
        }
    }
}

/// Matrix transpose * matrix transpose (hand tiled and unrolled)
template <>
inline void mTxmT(long dimi, long dimj, long dimk,
                  double* restrict c, const double* restrict a, const double* restrict b) {
    /*
    c(i,j) = c(i,j) + sum(k) a(k,i)*b(j,k)

    where it is assumed that the last index in each array is has unit
    stride and the dimensions are as provided.

    Tiled k, copy row of a into temporary, and unroll j once.
    */

    const int ktile=32;
    double ai[ktile];
    long dimj2 = (dimj/2)*2;

    double* restrict csave = c;
    const double* asave = a;
    for (long klo=0; klo<dimk; klo+=ktile, asave+=ktile*dimi, b+=ktile) {
        long khi = klo+ktile;
        if (khi > dimk) khi = dimk;
        long nk = khi-klo;

        a = asave;
        c = csave;
        for (long i=0; i<dimi; ++i,c+=dimj,++a) {
            const double* q = a;
            for (long k=0; k<nk; ++k,q+=dimi) ai[k] = *q;

            const double* bj0 = b;
            for (long j=0; j<dimj2; j+=2,bj0+=2*dimk) {
                const double* bj1 = bj0+dimk;
                double sum0 = 0;
                double sum1 = 0;
                for (long k=0; k<nk; ++k) {
                    sum0 += ai[k]*bj0[k];
                    sum1 += ai[k]*bj1[k];
                }
                c[j  ] += sum0;
                c[j+1] += sum1;
            }

            for (long j=dimj2; j<dimj; ++j,bj0+=dimk) {
                double sum = 0;
                for (long k=0; k<nk; ++k) {
                    sum += ai[k]*bj0[k];
                }
                c[j] += sum;
            }
        }
    }
}
#endif
#endif // MADNESS_TENSOR_MXM_H__INCLUDED
