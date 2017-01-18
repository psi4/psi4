/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \defgroup QT libqt: The Quantum-Trio Miscellaneous Library */

/*!
 \file
 \ingroup QT
 \brief The PSI3 BLAS1 interface routines

 Interface to the BLAS routines

 C. David Sherrill
 Anna I. Krylov

 May 1998

 NOTE: Refactored by Rob Parrish on 1/24/2010
 This file now contains all relevant BLAS1
 routines, with provisions made for >2^31
 elements (ULI sizes).

 All BLAS2 and BLAS3 routines are now wrapped
 and are in blas_inftc23.cc

*/

#include <cstdio>
#include <limits.h>
#include <cmath>

#include "psi4/libqt/blas_intfc_mangle.h"

extern "C" {

extern void F_DSWAP(int *length, double *x, int *incx, double *y, int *inc_y);
extern void F_DAXPY(int *length, double *a, double *x, int *inc_x,
  double *y, int *inc_y);
extern void F_DCOPY(int *length, double *x, int *inc_x,
  double *y, int *inc_y);
extern void F_DGEMM(char *transa, char *transb, int *m, int *n, int *k,
  double *alpha, double *A, int *lda, double *B, int *ldb,
  double *beta, double *C, int *ldc);
extern void F_DSYMM(char* side, char *uplo, int *m, int *n,
  double *alpha, double *A, int *lda, double *B, int *ldb,
  double *beta, double *C, int *ldc);
extern void F_DROT(int *ntot, double *x, int *incx, double *y, int *incy,
  double *cotheta, double *sintheta);
extern void F_DSCAL(int *n, double *alpha, double *vec, int *inc);
extern void F_DGEMV(char *transa, int *m, int *n, double *alpha, double *A,
  int *lda, double *X, int *inc_x, double *beta, double *Y, int *inc_y);
extern void F_DSYMV(char *uplo, int *n, double *alpha, double *A,
  int *lda, double *X, int *inc_x, double *beta, double *Y, int *inc_y);
extern void F_DSPMV(char *uplo, int *n, double *alpha, double *A, double *X,
  int *inc_x, double *beta, double *Y, int *inc_y);
extern double F_DDOT(int *n, double *x, int *incx, double *y, int *incy);
extern double F_DNRM2(int *n, double *x, int *incx);
extern double F_DASUM(int *n, double *x, int *incx);
extern int F_IDAMAX(int *n, double *x, int *incx);

}

namespace psi {

/**
 * Swaps a vector with another vector.
 *
 * @param length Specifies the number of elements in vectors x and y.
 * @param x Array, DIMENSION at least (1 + (n-1)*abs(incx)).
 * @param inc_x Specifies the increment for the elements of x.
 * @param y Array, DIMENSION at least (1 + (n-1)*abs(incy)).
 * @param inc_y Specifies the increment for the elements of y.
 *
 * @ingroup QT
 */
void C_DSWAP(unsigned long int length, double *x, int inc_x, double *y, int inc_y)
{
    int big_blocks = (int)(length / INT_MAX);
    int small_size = (int)(length % INT_MAX);
    for (int block = 0; block <= big_blocks; block++) {
        double* x_s = &x[block*inc_x*(unsigned long int)INT_MAX];
        double* y_s = &y[block*inc_y*(unsigned long int)INT_MAX];
        signed int length_s = (block == big_blocks) ? small_size : INT_MAX;
        ::F_DSWAP(&length_s, x_s, &inc_x, y_s, &inc_y);
    }
}

/*!
 * This function performs y = a * x + y.
 *
 * Steps every inc_x in x and every inc_y in y (normally both 1).
 *
 * \param length   length of arrays
 * \param a        scalar a to multiply vector x
 * \param x        vector x
 * \param inc_x    how many places to skip to get to next element in x
 * \param y        vector y
 * \param inc_y    how many places to skip to get to next element in y
 *
 * \ingroup QT
 */
void C_DAXPY(unsigned long int length, double a, double *x, int inc_x,
             double *y, int inc_y)
{
    int big_blocks = (int)(length / INT_MAX);
    int small_size = (int)(length % INT_MAX);
    for (int block = 0; block <= big_blocks; block++) {
        double* x_s = &x[block*inc_x*(unsigned long int)INT_MAX];
        double* y_s = &y[block*inc_y*(unsigned long int)INT_MAX];
        signed int length_s = (block == big_blocks) ? small_size : INT_MAX;
        ::F_DAXPY(&length_s, &a, x_s, &inc_x, y_s, &inc_y);
    }
}

/*!
 * This function copies x into y.
 *
 * Steps every inc_x in x and every inc_y in y (normally both 1).
 *
 * \param length  = length of array
 * \param x       = vector x
 * \param inc_x   = how many places to skip to get to next element in x
 * \param y       = vector y
 * \param inc_y   = how many places to skip to get to next element in y
 *
 * \ingroup QT
 */
void C_DCOPY(unsigned long int length, double *x, int inc_x,
             double *y, int inc_y)
{
    int big_blocks = (int)(length / INT_MAX);
    int small_size = (int)(length % INT_MAX);
    for (int block = 0; block <= big_blocks; block++) {
        double* x_s = &x[block*inc_x*(unsigned long int)INT_MAX];
        double* y_s = &y[block*inc_y*(unsigned long int)INT_MAX];
        signed int length_s = (block == big_blocks) ? small_size : INT_MAX;
        ::F_DCOPY(&length_s, x_s, &inc_x, y_s, &inc_y);
    }
}


/*!
 * This function scales a vector by a real scalar.
 *
 * \param length length of array
 * \param alpha  scale factor
 * \param vec    vector to scale
 * \param inc    how many places to skip to get to next element in vec
 *
 * \ingroup QT
 */
void C_DSCAL(unsigned long int length, double alpha, double *vec, int inc)
{
    int big_blocks = (int)(length / INT_MAX);
    int small_size = (int)(length % INT_MAX);
    for (int block = 0; block <= big_blocks; block++) {
        double* vec_s = &vec[block*inc*(unsigned long int)INT_MAX];
        signed int length_s = (block == big_blocks) ? small_size : INT_MAX;
        ::F_DSCAL(&length_s, &alpha, vec_s, &inc);
    }
}


/*!
 *Calculates a plane Givens rotation for vectors x, y and
 * angle theta.  x = x*cos + y*sin, y = -x*sin + y*cos.
 *
 * \param x      vector x
 * \param y      vector Y
 * \param length length of x,y
 * \param inc_x  how many places to skip to get to the next element of x
 * \param inc_y  how many places to skip to get to the next element of y
 *
 * \ingroup QT
 */
void C_DROT(unsigned long int length, double *x, int inc_x, double *y, int inc_y,
            double costheta, double sintheta)
{
    int big_blocks = (int)(length / INT_MAX);
    int small_size = (int)(length % INT_MAX);
    for (int block = 0; block <= big_blocks; block++) {
        double* x_s = &x[block*inc_x*(unsigned long int)INT_MAX];
        double* y_s = &y[block*inc_y*(unsigned long int)INT_MAX];
        signed int length_s = (block == big_blocks) ? small_size : INT_MAX;
        ::F_DROT(&length_s, x_s, &inc_x, y_s, &inc_y, &costheta, &sintheta);
    }

}

/*!
 * This function returns the dot product of two vectors, x and y.
 *
 * \param length Number of elements in x and y.
 * \param x      A pointer to the beginning of the data in x.
 *               Must be of at least length (1+(N-1)*abs(inc_x).
 * \param inc_x  how many places to skip to get to next element in x
 * \param y      A pointer to the beginning of the data in y.
 * \param inc_y  how many places to skip to get to next element in y
 *
 * @returns the dot product
 *
 * \ingroup QT
 */

double C_DDOT(unsigned long int length, double *x, int inc_x, double *y, int inc_y)
{
    if(length == 0) return 0.0;

    double reg = 0.0;

    int big_blocks = (int)(length / INT_MAX);
    int small_size = (int)(length % INT_MAX);
    for (int block = 0; block <= big_blocks; block++) {
        double* x_s = &x[block*inc_x*(unsigned long int)INT_MAX];
        double* y_s = &y[block*inc_y*(unsigned long int)INT_MAX];
        signed int length_s = (block == big_blocks) ? small_size : INT_MAX;
        reg += ::F_DDOT(&length_s, x_s, &inc_x, y_s, &inc_y);
    }

    return reg;
}
/*!
 * This function returns the square of the norm of this vector.
 *
 * \param length Number of elements in x.
 * \param x      A pointer to the beginning of the data in x.
 *               Must be of at least length (1+(N-1)*abs(inc_x).
 * \param inc_x  how many places to skip to get to next element in x
 *
 * @returns the norm squared product
 *
 * \ingroup QT
 */

double C_DNRM2(unsigned long int length, double *x, int inc_x)
{
    if(length == 0) return 0.0;

    double reg = 0.0;

    int big_blocks = (int)(length / INT_MAX);
    int small_size = (int)(length % INT_MAX);
    for (int block = 0; block <= big_blocks; block++) {
        double* x_s = &x[block*inc_x*(unsigned long int)INT_MAX];
        signed int length_s = (block == big_blocks) ? small_size : INT_MAX;
        reg += ::F_DNRM2(&length_s, x_s, &inc_x);
    }

    return reg;
}
/*!
 * This function returns the sum of the absolute value of this vector.
 *
 * \param length Number of elements in x.
 * \param x      A pointer to the beginning of the data in x.
 *               Must be of at least length (1+(N-1)*abs(inc_x).
 * \param inc_x  how many places to skip to get to next element in x
 *
 * @returns the sum of the absolute value
 *
 * \ingroup QT
 */

double C_DASUM(unsigned long int length, double *x, int inc_x)
{
    if(length == 0) return 0.0;

    double reg = 0.0;

    int big_blocks = (int)(length / INT_MAX);
    int small_size = (int)(length % INT_MAX);
    for (int block = 0; block <= big_blocks; block++) {
        double* x_s = &x[block*inc_x*(unsigned long int)INT_MAX];
        signed int length_s = (block == big_blocks) ? small_size : INT_MAX;
        reg += ::F_DASUM(&length_s, x_s, &inc_x);
    }

    return reg;
}
/*!
 * This function returns the index of the largest absolute value compoment of this vector.
 *
 * \param length Number of elements in x.
 * \param x      A pointer to the beginning of the data in x.
 *               Must be of at least length (1+(N-1)*abs(inc_x).
 * \param inc_x  how many places to skip to get to next element in x
 *
 * @returns the index of the largest absolute value
 *
 * \ingroup QT
 */

unsigned long int C_IDAMAX(unsigned long int length, double *x, int inc_x)
{
    if(length == 0) return 0L;

    unsigned long int reg = 0L;
    unsigned long int reg2 = 0L;

    int big_blocks = (int)(length / INT_MAX);
    int small_size = (int)(length % INT_MAX);
    for (int block = 0; block <= big_blocks; block++) {
        double* x_s = &x[block*inc_x*(unsigned long int)INT_MAX];
        signed int length_s = (block == big_blocks) ? small_size : INT_MAX;
        reg2 = ::F_IDAMAX(&length_s, x_s, &inc_x) + block*inc_x*(unsigned long int)INT_MAX;
        if (fabs(x[reg]) > fabs(x[reg2]))
            reg = reg2;
    }

    return reg;
}


}
