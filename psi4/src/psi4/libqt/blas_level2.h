/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*!
 * \file
 * \brief Function signatures for wrappers to BLAS level 2 subroutines
 * \ingroup QT
 */

#pragma once

#include <string>

#ifdef USING_LAPACK_MKL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

#include "psi4/pragma.h"

namespace psi {
namespace detail {
CBLAS_TRANSPOSE dispatch_trans(char trans, std::string func, const char * file, int line);

CBLAS_UPLO dispatch_uplo(char uplo, std::string func, const char * file, int line);

CBLAS_DIAG dispatch_diag(char diag, std::string func, const char * file, int line);

CBLAS_SIDE dispatch_side(char side, std::string func, const char * file, int line);
}  // namespace detail

/**
 *  Purpose
 *  =======
 *
 *  DGBMV  performs one of the matrix-vector operations
 *
 *     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
 *
 *  where alpha and beta are scalars, x and y are vectors and A is an
 *  m by n band matrix, with kl sub-diagonals and ku super-diagonals.
 *
 *  Arguments
 *  ==========
 *
 *  TRANS  - CHARACTER*1.
 *           On entry, TRANS specifies the operation to be performed as
 *           follows:
 *
 *              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
 *
 *              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
 *
 *              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
 *
 *           Unchanged on exit.
 *
 *  M      - INTEGER.
 *           On entry, M specifies the number of rows of the matrix A.
 *           M must be at least zero.
 *           Unchanged on exit.
 *
 *  N      - INTEGER.
 *           On entry, N specifies the number of columns of the matrix A.
 *           N must be at least zero.
 *           Unchanged on exit.
 *
 *  KL     - INTEGER.
 *           On entry, KL specifies the number of sub-diagonals of the
 *           matrix A. KL must satisfy  0 .le. KL.
 *           Unchanged on exit.
 *
 *  KU     - INTEGER.
 *           On entry, KU specifies the number of super-diagonals of the
 *           matrix A. KU must satisfy  0 .le. KU.
 *           Unchanged on exit.
 *
 *  ALPHA  - DOUBLE PRECISION.
 *           On entry, ALPHA specifies the scalar alpha.
 *           Unchanged on exit.
 *
 *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
 *           Before entry, the leading ( kl + ku + 1 ) by n part of the
 *           array A must contain the matrix of coefficients, supplied
 *           column by column, with the leading diagonal of the matrix in
 *           row ( ku + 1 ) of the array, the first super-diagonal
 *           starting at position 2 in row ku, the first sub-diagonal
 *           starting at position 1 in row ( ku + 2 ), and so on.
 *           Elements in the array A that do not correspond to elements
 *           in the band matrix (such as the top left ku by ku triangle)
 *           are not referenced.
 *           The following program segment will transfer a band matrix
 *           from conventional full matrix storage to band storage:
 *
 *                 DO 20, J = 1, N
 *                    K = KU + 1 - J
 *                    DO 10, I = MAX( 1, J - KU ), MIN( M, J + KL )
 *                       A( K + I, J ) = matrix( I, J )
 *              10    CONTINUE
 *              20 CONTINUE
 *
 *           Unchanged on exit.
 *
 *  LDA    - INTEGER.
 *           On entry, LDA specifies the first dimension of A as declared
 *           in the calling (sub) program. LDA must be at least
 *           ( kl + ku + 1 ).
 *           Unchanged on exit.
 *
 *  X      - DOUBLE PRECISION array of DIMENSION at least
 *           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
 *           and at least
 *           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
 *           Before entry, the incremented array X must contain the
 *           vector x.
 *           Unchanged on exit.
 *
 *  INCX   - INTEGER.
 *           On entry, INCX specifies the increment for the elements of
 *           X. INCX must not be zero.
 *           Unchanged on exit.
 *
 *  BETA   - DOUBLE PRECISION.
 *           On entry, BETA specifies the scalar beta. When BETA is
 *           supplied as zero then Y need not be set on input.
 *           Unchanged on exit.
 *
 *  Y      - DOUBLE PRECISION array of DIMENSION at least
 *           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
 *           and at least
 *           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
 *           Before entry, the incremented array Y must contain the
 *           vector y. On exit, Y is overwritten by the updated vector y.
 *
 *  INCY   - INTEGER.
 *           On entry, INCY specifies the increment for the elements of
 *           Y. INCY must not be zero.
 *           Unchanged on exit.
 *
 *
 *  Level 2 Blas routine.
 *
 *  -- Written on 22-October-1986.
 *     Jack Dongarra, Argonne National Lab.
 *     Jeremy Du Croz, Nag Central Office.
 *     Sven Hammarling, Nag Central Office.
 *     Richard Hanson, Sandia National Labs.
 *
 *     .. Parameters ..
 */
PSI_API
void C_DGBMV(char trans, int m, int n, int kl, int ku, double alpha, double* a, int lda, double* x, int incx,
             double beta, double* y, int incy);

/**
 *  Purpose
 *  =======
 *
 *  DGEMV  performs one of the matrix-vector operations
 *
 *     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
 *
 *  where alpha and beta are scalars, x and y are vectors and A is an
 *  m by n matrix.
 *
 *  Arguments
 *  ==========
 *
 *  TRANS  - CHARACTER*1.
 *           On entry, TRANS specifies the operation to be performed as
 *           follows:
 *
 *              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
 *
 *              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
 *
 *              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
 *
 *           Unchanged on exit.
 *
 *  M      - INTEGER.
 *           On entry, M specifies the number of rows of the matrix A.
 *           M must be at least zero.
 *           Unchanged on exit.
 *
 *  N      - INTEGER.
 *           On entry, N specifies the number of columns of the matrix A.
 *           N must be at least zero.
 *           Unchanged on exit.
 *
 *  ALPHA  - DOUBLE PRECISION.
 *           On entry, ALPHA specifies the scalar alpha.
 *           Unchanged on exit.
 *
 *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
 *           Before entry, the leading m by n part of the array A must
 *           contain the matrix of coefficients.
 *           Unchanged on exit.
 *
 *  LDA    - INTEGER.
 *           On entry, LDA specifies the first dimension of A as declared
 *           in the calling (sub) program. LDA must be at least
 *           max( 1, m ).
 *           Unchanged on exit.
 *
 *  X      - DOUBLE PRECISION array of DIMENSION at least
 *           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
 *           and at least
 *           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
 *           Before entry, the incremented array X must contain the
 *           vector x.
 *           Unchanged on exit.
 *
 *  INCX   - INTEGER.
 *           On entry, INCX specifies the increment for the elements of
 *           X. INCX must not be zero.
 *           Unchanged on exit.
 *
 *  BETA   - DOUBLE PRECISION.
 *           On entry, BETA specifies the scalar beta. When BETA is
 *           supplied as zero then Y need not be set on input.
 *           Unchanged on exit.
 *
 *  Y      - DOUBLE PRECISION array of DIMENSION at least
 *           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
 *           and at least
 *           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
 *           Before entry with BETA non-zero, the incremented array Y
 *           must contain the vector y. On exit, Y is overwritten by the
 *           updated vector y.
 *
 *  INCY   - INTEGER.
 *           On entry, INCY specifies the increment for the elements of
 *           Y. INCY must not be zero.
 *           Unchanged on exit.
 *
 *
 *  Level 2 Blas routine.
 *
 *  -- Written on 22-October-1986.
 *     Jack Dongarra, Argonne National Lab.
 *     Jeremy Du Croz, Nag Central Office.
 *     Sven Hammarling, Nag Central Office.
 *     Richard Hanson, Sandia National Labs.
 *
 *
 *     .. Parameters ..
 */
PSI_API
void C_DGEMV(char trans, int m, int n, double alpha, double* a, int lda, double* x, int incx, double beta, double* y,
             int incy);

/**
 *  Purpose
 *  =======
 *
 *  DGER   performs the rank 1 operation
 *
 *     A := alpha*x*y' + A,
 *
 *  where alpha is a scalar, x is an m element vector, y is an n element
 *  vector and A is an m by n matrix.
 *
 *  Arguments
 *  ==========
 *
 *  M      - INTEGER.
 *           On entry, M specifies the number of rows of the matrix A.
 *           M must be at least zero.
 *           Unchanged on exit.
 *
 *  N      - INTEGER.
 *           On entry, N specifies the number of columns of the matrix A.
 *           N must be at least zero.
 *           Unchanged on exit.
 *
 *  ALPHA  - DOUBLE PRECISION.
 *           On entry, ALPHA specifies the scalar alpha.
 *           Unchanged on exit.
 *
 *  X      - DOUBLE PRECISION array of dimension at least
 *           ( 1 + ( m - 1 )*abs( INCX ) ).
 *           Before entry, the incremented array X must contain the m
 *           element vector x.
 *           Unchanged on exit.
 *
 *  INCX   - INTEGER.
 *           On entry, INCX specifies the increment for the elements of
 *           X. INCX must not be zero.
 *           Unchanged on exit.
 *
 *  Y      - DOUBLE PRECISION array of dimension at least
 *           ( 1 + ( n - 1 )*abs( INCY ) ).
 *           Before entry, the incremented array Y must contain the n
 *           element vector y.
 *           Unchanged on exit.
 *
 *  INCY   - INTEGER.
 *           On entry, INCY specifies the increment for the elements of
 *           Y. INCY must not be zero.
 *           Unchanged on exit.
 *
 *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
 *           Before entry, the leading m by n part of the array A must
 *           contain the matrix of coefficients. On exit, A is
 *           overwritten by the updated matrix.
 *
 *  LDA    - INTEGER.
 *           On entry, LDA specifies the first dimension of A as declared
 *           in the calling (sub) program. LDA must be at least
 *           max( 1, m ).
 *           Unchanged on exit.
 *
 *
 *  Level 2 Blas routine.
 *
 *  -- Written on 22-October-1986.
 *     Jack Dongarra, Argonne National Lab.
 *     Jeremy Du Croz, Nag Central Office.
 *     Sven Hammarling, Nag Central Office.
 *     Richard Hanson, Sandia National Labs.
 *
 *
 *     .. Parameters ..
 */
PSI_API
void C_DGER(int m, int n, double alpha, double* x, int incx, double* y, int incy, double* a, int lda);

/**
 *  Purpose
 *  =======
 *
 *  DSBMV  performs the matrix-vector  operation
 *
 *     y := alpha*A*x + beta*y,
 *
 *  where alpha and beta are scalars, x and y are n element vectors and
 *  A is an n by n symmetric band matrix, with k super-diagonals.
 *
 *  Arguments
 *  ==========
 *
 *  UPLO   - CHARACTER*1.
 *           On entry, UPLO specifies whether the upper or lower
 *           triangular part of the band matrix A is being supplied as
 *           follows:
 *
 *              UPLO = 'U' or 'u'   The upper triangular part of A is
 *                                  being supplied.
 *
 *              UPLO = 'L' or 'l'   The lower triangular part of A is
 *                                  being supplied.
 *
 *           Unchanged on exit.
 *
 *  N      - INTEGER.
 *           On entry, N specifies the order of the matrix A.
 *           N must be at least zero.
 *           Unchanged on exit.
 *
 *  K      - INTEGER.
 *           On entry, K specifies the number of super-diagonals of the
 *           matrix A. K must satisfy  0 .le. K.
 *           Unchanged on exit.
 *
 *  ALPHA  - DOUBLE PRECISION.
 *           On entry, ALPHA specifies the scalar alpha.
 *           Unchanged on exit.
 *
 *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
 *           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
 *           by n part of the array A must contain the upper triangular
 *           band part of the symmetric matrix, supplied column by
 *           column, with the leading diagonal of the matrix in row
 *           ( k + 1 ) of the array, the first super-diagonal starting at
 *           position 2 in row k, and so on. The top left k by k triangle
 *           of the array A is not referenced.
 *           The following program segment will transfer the upper
 *           triangular part of a symmetric band matrix from conventional
 *           full matrix storage to band storage:
 *
 *                 DO 20, J = 1, N
 *                    M = K + 1 - J
 *                    DO 10, I = MAX( 1, J - K ), J
 *                       A( M + I, J ) = matrix( I, J )
 *              10    CONTINUE
 *              20 CONTINUE
 *
 *           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
 *           by n part of the array A must contain the lower triangular
 *           band part of the symmetric matrix, supplied column by
 *           column, with the leading diagonal of the matrix in row 1 of
 *           the array, the first sub-diagonal starting at position 1 in
 *           row 2, and so on. The bottom right k by k triangle of the
 *           array A is not referenced.
 *           The following program segment will transfer the lower
 *           triangular part of a symmetric band matrix from conventional
 *           full matrix storage to band storage:
 *
 *                 DO 20, J = 1, N
 *                    M = 1 - J
 *                    DO 10, I = J, MIN( N, J + K )
 *                       A( M + I, J ) = matrix( I, J )
 *              10    CONTINUE
 *              20 CONTINUE
 *
 *           Unchanged on exit.
 *
 *  LDA    - INTEGER.
 *           On entry, LDA specifies the first dimension of A as declared
 *           in the calling (sub) program. LDA must be at least
 *           ( k + 1 ).
 *           Unchanged on exit.
 *
 *  X      - DOUBLE PRECISION array of DIMENSION at least
 *           ( 1 + ( n - 1 )*abs( INCX ) ).
 *           Before entry, the incremented array X must contain the
 *           vector x.
 *           Unchanged on exit.
 *
 *  INCX   - INTEGER.
 *           On entry, INCX specifies the increment for the elements of
 *           X. INCX must not be zero.
 *           Unchanged on exit.
 *
 *  BETA   - DOUBLE PRECISION.
 *           On entry, BETA specifies the scalar beta.
 *           Unchanged on exit.
 *
 *  Y      - DOUBLE PRECISION array of DIMENSION at least
 *           ( 1 + ( n - 1 )*abs( INCY ) ).
 *           Before entry, the incremented array Y must contain the
 *           vector y. On exit, Y is overwritten by the updated vector y.
 *
 *  INCY   - INTEGER.
 *           On entry, INCY specifies the increment for the elements of
 *           Y. INCY must not be zero.
 *           Unchanged on exit.
 *
 *
 *  Level 2 Blas routine.
 *
 *  -- Written on 22-October-1986.
 *     Jack Dongarra, Argonne National Lab.
 *     Jeremy Du Croz, Nag Central Office.
 *     Sven Hammarling, Nag Central Office.
 *     Richard Hanson, Sandia National Labs.
 *
 *
 *     .. Parameters ..
 */
PSI_API
void C_DSBMV(char uplo, int n, int k, double alpha, double* a, int lda, double* x, int incx, double beta, double* y,
             int incy);

/**
 *  Purpose
 *  =======
 *
 *  DSPMV  performs the matrix-vector operation
 *
 *     y := alpha*A*x + beta*y,
 *
 *  where alpha and beta are scalars, x and y are n element vectors and
 *  A is an n by n symmetric matrix, supplied in packed form.
 *
 *  Arguments
 *  ==========
 *
 *  UPLO   - CHARACTER*1.
 *           On entry, UPLO specifies whether the upper or lower
 *           triangular part of the matrix A is supplied in the packed
 *           array AP as follows:
 *
 *              UPLO = 'U' or 'u'   The upper triangular part of A is
 *                                  supplied in AP.
 *
 *              UPLO = 'L' or 'l'   The lower triangular part of A is
 *                                  supplied in AP.
 *
 *           Unchanged on exit.
 *
 *  N      - INTEGER.
 *           On entry, N specifies the order of the matrix A.
 *           N must be at least zero.
 *           Unchanged on exit.
 *
 *  ALPHA  - DOUBLE PRECISION.
 *           On entry, ALPHA specifies the scalar alpha.
 *           Unchanged on exit.
 *
 *  AP     - DOUBLE PRECISION array of DIMENSION at least
 *           ( ( n*( n + 1 ) )/2 ).
 *           Before entry with UPLO = 'U' or 'u', the array AP must
 *           contain the upper triangular part of the symmetric matrix
 *           packed sequentially, column by column, so that AP( 1 )
 *           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
 *           and a( 2, 2 ) respectively, and so on.
 *           Before entry with UPLO = 'L' or 'l', the array AP must
 *           contain the lower triangular part of the symmetric matrix
 *           packed sequentially, column by column, so that AP( 1 )
 *           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
 *           and a( 3, 1 ) respectively, and so on.
 *           Unchanged on exit.
 *
 *  X      - DOUBLE PRECISION array of dimension at least
 *           ( 1 + ( n - 1 )*abs( INCX ) ).
 *           Before entry, the incremented array X must contain the n
 *           element vector x.
 *           Unchanged on exit.
 *
 *  INCX   - INTEGER.
 *           On entry, INCX specifies the increment for the elements of
 *           X. INCX must not be zero.
 *           Unchanged on exit.
 *
 *  BETA   - DOUBLE PRECISION.
 *           On entry, BETA specifies the scalar beta. When BETA is
 *           supplied as zero then Y need not be set on input.
 *           Unchanged on exit.
 *
 *  Y      - DOUBLE PRECISION array of dimension at least
 *           ( 1 + ( n - 1 )*abs( INCY ) ).
 *           Before entry, the incremented array Y must contain the n
 *           element vector y. On exit, Y is overwritten by the updated
 *           vector y.
 *
 *  INCY   - INTEGER.
 *           On entry, INCY specifies the increment for the elements of
 *           Y. INCY must not be zero.
 *           Unchanged on exit.
 *
 *
 *  Level 2 Blas routine.
 *
 *  -- Written on 22-October-1986.
 *     Jack Dongarra, Argonne National Lab.
 *     Jeremy Du Croz, Nag Central Office.
 *     Sven Hammarling, Nag Central Office.
 *     Richard Hanson, Sandia National Labs.
 *
 *
 *     .. Parameters ..
 */
PSI_API
void C_DSPMV(char uplo, int n, double alpha, double* ap, double* x, int incx, double beta, double* y, int incy);

/**
 *  Purpose
 *  =======
 *
 *  DSPR    performs the symmetric rank 1 operation
 *
 *     A := alpha*x*x' + A,
 *
 *  where alpha is a real scalar, x is an n element vector and A is an
 *  n by n symmetric matrix, supplied in packed form.
 *
 *  Arguments
 *  ==========
 *
 *  UPLO   - CHARACTER*1.
 *           On entry, UPLO specifies whether the upper or lower
 *           triangular part of the matrix A is supplied in the packed
 *           array AP as follows:
 *
 *              UPLO = 'U' or 'u'   The upper triangular part of A is
 *                                  supplied in AP.
 *
 *              UPLO = 'L' or 'l'   The lower triangular part of A is
 *                                  supplied in AP.
 *
 *           Unchanged on exit.
 *
 *  N      - INTEGER.
 *           On entry, N specifies the order of the matrix A.
 *           N must be at least zero.
 *           Unchanged on exit.
 *
 *  ALPHA  - DOUBLE PRECISION.
 *           On entry, ALPHA specifies the scalar alpha.
 *           Unchanged on exit.
 *
 *  X      - DOUBLE PRECISION array of dimension at least
 *           ( 1 + ( n - 1 )*abs( INCX ) ).
 *           Before entry, the incremented array X must contain the n
 *           element vector x.
 *           Unchanged on exit.
 *
 *  INCX   - INTEGER.
 *           On entry, INCX specifies the increment for the elements of
 *           X. INCX must not be zero.
 *           Unchanged on exit.
 *
 *  AP     - DOUBLE PRECISION array of DIMENSION at least
 *           ( ( n*( n + 1 ) )/2 ).
 *           Before entry with  UPLO = 'U' or 'u', the array AP must
 *           contain the upper triangular part of the symmetric matrix
 *           packed sequentially, column by column, so that AP( 1 )
 *           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
 *           and a( 2, 2 ) respectively, and so on. On exit, the array
 *           AP is overwritten by the upper triangular part of the
 *           updated matrix.
 *           Before entry with UPLO = 'L' or 'l', the array AP must
 *           contain the lower triangular part of the symmetric matrix
 *           packed sequentially, column by column, so that AP( 1 )
 *           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
 *           and a( 3, 1 ) respectively, and so on. On exit, the array
 *           AP is overwritten by the lower triangular part of the
 *           updated matrix.
 *
 *
 *  Level 2 Blas routine.
 *
 *  -- Written on 22-October-1986.
 *     Jack Dongarra, Argonne National Lab.
 *     Jeremy Du Croz, Nag Central Office.
 *     Sven Hammarling, Nag Central Office.
 *     Richard Hanson, Sandia National Labs.
 *
 *
 *     .. Parameters ..
 */
PSI_API
void C_DSPR(char uplo, int n, double alpha, double* x, int incx, double* ap);

/**
 *  Purpose
 *  =======
 *
 *  DSPR2  performs the symmetric rank 2 operation
 *
 *     A := alpha*x*y' + alpha*y*x' + A,
 *
 *  where alpha is a scalar, x and y are n element vectors and A is an
 *  n by n symmetric matrix, supplied in packed form.
 *
 *  Arguments
 *  ==========
 *
 *  UPLO   - CHARACTER*1.
 *           On entry, UPLO specifies whether the upper or lower
 *           triangular part of the matrix A is supplied in the packed
 *           array AP as follows:
 *
 *              UPLO = 'U' or 'u'   The upper triangular part of A is
 *                                  supplied in AP.
 *
 *              UPLO = 'L' or 'l'   The lower triangular part of A is
 *                                  supplied in AP.
 *
 *           Unchanged on exit.
 *
 *  N      - INTEGER.
 *           On entry, N specifies the order of the matrix A.
 *           N must be at least zero.
 *           Unchanged on exit.
 *
 *  ALPHA  - DOUBLE PRECISION.
 *           On entry, ALPHA specifies the scalar alpha.
 *           Unchanged on exit.
 *
 *  X      - DOUBLE PRECISION array of dimension at least
 *           ( 1 + ( n - 1 )*abs( INCX ) ).
 *           Before entry, the incremented array X must contain the n
 *           element vector x.
 *           Unchanged on exit.
 *
 *  INCX   - INTEGER.
 *           On entry, INCX specifies the increment for the elements of
 *           X. INCX must not be zero.
 *           Unchanged on exit.
 *
 *  Y      - DOUBLE PRECISION array of dimension at least
 *           ( 1 + ( n - 1 )*abs( INCY ) ).
 *           Before entry, the incremented array Y must contain the n
 *           element vector y.
 *           Unchanged on exit.
 *
 *  INCY   - INTEGER.
 *           On entry, INCY specifies the increment for the elements of
 *           Y. INCY must not be zero.
 *           Unchanged on exit.
 *
 *  AP     - DOUBLE PRECISION array of DIMENSION at least
 *           ( ( n*( n + 1 ) )/2 ).
 *           Before entry with  UPLO = 'U' or 'u', the array AP must
 *           contain the upper triangular part of the symmetric matrix
 *           packed sequentially, column by column, so that AP( 1 )
 *           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
 *           and a( 2, 2 ) respectively, and so on. On exit, the array
 *           AP is overwritten by the upper triangular part of the
 *           updated matrix.
 *           Before entry with UPLO = 'L' or 'l', the array AP must
 *           contain the lower triangular part of the symmetric matrix
 *           packed sequentially, column by column, so that AP( 1 )
 *           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
 *           and a( 3, 1 ) respectively, and so on. On exit, the array
 *           AP is overwritten by the lower triangular part of the
 *           updated matrix.
 *
 *
 *  Level 2 Blas routine.
 *
 *  -- Written on 22-October-1986.
 *     Jack Dongarra, Argonne National Lab.
 *     Jeremy Du Croz, Nag Central Office.
 *     Sven Hammarling, Nag Central Office.
 *     Richard Hanson, Sandia National Labs.
 *
 *
 *     .. Parameters ..
 */
PSI_API
void C_DSPR2(char uplo, int n, double alpha, double* x, int incx, double* y, int incy, double* ap);

/**
 *  Purpose
 *  =======
 *
 *  DSYMV  performs the matrix-vector  operation
 *
 *     y := alpha*A*x + beta*y,
 *
 *  where alpha and beta are scalars, x and y are n element vectors and
 *  A is an n by n symmetric matrix.
 *
 *  Arguments
 *  ==========
 *
 *  UPLO   - CHARACTER*1.
 *           On entry, UPLO specifies whether the upper or lower
 *           triangular part of the array A is to be referenced as
 *           follows:
 *
 *              UPLO = 'U' or 'u'   Only the upper triangular part of A
 *                                  is to be referenced.
 *
 *              UPLO = 'L' or 'l'   Only the lower triangular part of A
 *                                  is to be referenced.
 *
 *           Unchanged on exit.
 *
 *  N      - INTEGER.
 *           On entry, N specifies the order of the matrix A.
 *           N must be at least zero.
 *           Unchanged on exit.
 *
 *  ALPHA  - DOUBLE PRECISION.
 *           On entry, ALPHA specifies the scalar alpha.
 *           Unchanged on exit.
 *
 *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
 *           Before entry with  UPLO = 'U' or 'u', the leading n by n
 *           upper triangular part of the array A must contain the upper
 *           triangular part of the symmetric matrix and the strictly
 *           lower triangular part of A is not referenced.
 *           Before entry with UPLO = 'L' or 'l', the leading n by n
 *           lower triangular part of the array A must contain the lower
 *           triangular part of the symmetric matrix and the strictly
 *           upper triangular part of A is not referenced.
 *           Unchanged on exit.
 *
 *  LDA    - INTEGER.
 *           On entry, LDA specifies the first dimension of A as declared
 *           in the calling (sub) program. LDA must be at least
 *           max( 1, n ).
 *           Unchanged on exit.
 *
 *  X      - DOUBLE PRECISION array of dimension at least
 *           ( 1 + ( n - 1 )*abs( INCX ) ).
 *           Before entry, the incremented array X must contain the n
 *           element vector x.
 *           Unchanged on exit.
 *
 *  INCX   - INTEGER.
 *           On entry, INCX specifies the increment for the elements of
 *           X. INCX must not be zero.
 *           Unchanged on exit.
 *
 *  BETA   - DOUBLE PRECISION.
 *           On entry, BETA specifies the scalar beta. When BETA is
 *           supplied as zero then Y need not be set on input.
 *           Unchanged on exit.
 *
 *  Y      - DOUBLE PRECISION array of dimension at least
 *           ( 1 + ( n - 1 )*abs( INCY ) ).
 *           Before entry, the incremented array Y must contain the n
 *           element vector y. On exit, Y is overwritten by the updated
 *           vector y.
 *
 *  INCY   - INTEGER.
 *           On entry, INCY specifies the increment for the elements of
 *           Y. INCY must not be zero.
 *           Unchanged on exit.
 *
 *
 *  Level 2 Blas routine.
 *
 *  -- Written on 22-October-1986.
 *     Jack Dongarra, Argonne National Lab.
 *     Jeremy Du Croz, Nag Central Office.
 *     Sven Hammarling, Nag Central Office.
 *     Richard Hanson, Sandia National Labs.
 *
 *
 *     .. Parameters ..
 */
PSI_API
void C_DSYMV(char uplo, int n, double alpha, double* a, int lda, double* x, int incx, double beta, double* y, int incy);

/**
 *  Purpose
 *  =======
 *
 *  DSYR   performs the symmetric rank 1 operation
 *
 *     A := alpha*x*x' + A,
 *
 *  where alpha is a real scalar, x is an n element vector and A is an
 *  n by n symmetric matrix.
 *
 *  Arguments
 *  ==========
 *
 *  UPLO   - CHARACTER*1.
 *           On entry, UPLO specifies whether the upper or lower
 *           triangular part of the array A is to be referenced as
 *           follows:
 *
 *              UPLO = 'U' or 'u'   Only the upper triangular part of A
 *                                  is to be referenced.
 *
 *              UPLO = 'L' or 'l'   Only the lower triangular part of A
 *                                  is to be referenced.
 *
 *           Unchanged on exit.
 *
 *  N      - INTEGER.
 *           On entry, N specifies the order of the matrix A.
 *           N must be at least zero.
 *           Unchanged on exit.
 *
 *  ALPHA  - DOUBLE PRECISION.
 *           On entry, ALPHA specifies the scalar alpha.
 *           Unchanged on exit.
 *
 *  X      - DOUBLE PRECISION array of dimension at least
 *           ( 1 + ( n - 1 )*abs( INCX ) ).
 *           Before entry, the incremented array X must contain the n
 *           element vector x.
 *           Unchanged on exit.
 *
 *  INCX   - INTEGER.
 *           On entry, INCX specifies the increment for the elements of
 *           X. INCX must not be zero.
 *           Unchanged on exit.
 *
 *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
 *           Before entry with  UPLO = 'U' or 'u', the leading n by n
 *           upper triangular part of the array A must contain the upper
 *           triangular part of the symmetric matrix and the strictly
 *           lower triangular part of A is not referenced. On exit, the
 *           upper triangular part of the array A is overwritten by the
 *           upper triangular part of the updated matrix.
 *           Before entry with UPLO = 'L' or 'l', the leading n by n
 *           lower triangular part of the array A must contain the lower
 *           triangular part of the symmetric matrix and the strictly
 *           upper triangular part of A is not referenced. On exit, the
 *           lower triangular part of the array A is overwritten by the
 *           lower triangular part of the updated matrix.
 *
 *  LDA    - INTEGER.
 *           On entry, LDA specifies the first dimension of A as declared
 *           in the calling (sub) program. LDA must be at least
 *           max( 1, n ).
 *           Unchanged on exit.
 *
 *
 *  Level 2 Blas routine.
 *
 *  -- Written on 22-October-1986.
 *     Jack Dongarra, Argonne National Lab.
 *     Jeremy Du Croz, Nag Central Office.
 *     Sven Hammarling, Nag Central Office.
 *     Richard Hanson, Sandia National Labs.
 *
 *
 *     .. Parameters ..
 */
PSI_API
void C_DSYR(char uplo, int n, double alpha, double* x, int incx, double* a, int lda);

/**
 *  Purpose
 *  =======
 *
 *  DSYR2  performs the symmetric rank 2 operation
 *
 *     A := alpha*x*y' + alpha*y*x' + A,
 *
 *  where alpha is a scalar, x and y are n element vectors and A is an n
 *  by n symmetric matrix.
 *
 *  Arguments
 *  ==========
 *
 *  UPLO   - CHARACTER*1.
 *           On entry, UPLO specifies whether the upper or lower
 *           triangular part of the array A is to be referenced as
 *           follows:
 *
 *              UPLO = 'U' or 'u'   Only the upper triangular part of A
 *                                  is to be referenced.
 *
 *              UPLO = 'L' or 'l'   Only the lower triangular part of A
 *                                  is to be referenced.
 *
 *           Unchanged on exit.
 *
 *  N      - INTEGER.
 *           On entry, N specifies the order of the matrix A.
 *           N must be at least zero.
 *           Unchanged on exit.
 *
 *  ALPHA  - DOUBLE PRECISION.
 *           On entry, ALPHA specifies the scalar alpha.
 *           Unchanged on exit.
 *
 *  X      - DOUBLE PRECISION array of dimension at least
 *           ( 1 + ( n - 1 )*abs( INCX ) ).
 *           Before entry, the incremented array X must contain the n
 *           element vector x.
 *           Unchanged on exit.
 *
 *  INCX   - INTEGER.
 *           On entry, INCX specifies the increment for the elements of
 *           X. INCX must not be zero.
 *           Unchanged on exit.
 *
 *  Y      - DOUBLE PRECISION array of dimension at least
 *           ( 1 + ( n - 1 )*abs( INCY ) ).
 *           Before entry, the incremented array Y must contain the n
 *           element vector y.
 *           Unchanged on exit.
 *
 *  INCY   - INTEGER.
 *           On entry, INCY specifies the increment for the elements of
 *           Y. INCY must not be zero.
 *           Unchanged on exit.
 *
 *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
 *           Before entry with  UPLO = 'U' or 'u', the leading n by n
 *           upper triangular part of the array A must contain the upper
 *           triangular part of the symmetric matrix and the strictly
 *           lower triangular part of A is not referenced. On exit, the
 *           upper triangular part of the array A is overwritten by the
 *           upper triangular part of the updated matrix.
 *           Before entry with UPLO = 'L' or 'l', the leading n by n
 *           lower triangular part of the array A must contain the lower
 *           triangular part of the symmetric matrix and the strictly
 *           upper triangular part of A is not referenced. On exit, the
 *           lower triangular part of the array A is overwritten by the
 *           lower triangular part of the updated matrix.
 *
 *  LDA    - INTEGER.
 *           On entry, LDA specifies the first dimension of A as declared
 *           in the calling (sub) program. LDA must be at least
 *           max( 1, n ).
 *           Unchanged on exit.
 *
 *
 *  Level 2 Blas routine.
 *
 *  -- Written on 22-October-1986.
 *     Jack Dongarra, Argonne National Lab.
 *     Jeremy Du Croz, Nag Central Office.
 *     Sven Hammarling, Nag Central Office.
 *     Richard Hanson, Sandia National Labs.
 *
 *
 *     .. Parameters ..
 */
PSI_API
void C_DSYR2(char uplo, int n, double alpha, double* x, int incx, double* y, int incy, double* a, int lda);

/**
 *  Purpose
 *  =======
 *
 *  DTBMV  performs one of the matrix-vector operations
 *
 *     x := A*x,   or   x := A'*x,
 *
 *  where x is an n element vector and  A is an n by n unit, or non-unit,
 *  upper or lower triangular band matrix, with ( k + 1 ) diagonals.
 *
 *  Arguments
 *  ==========
 *
 *  UPLO   - CHARACTER*1.
 *           On entry, UPLO specifies whether the matrix is an upper or
 *           lower triangular matrix as follows:
 *
 *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
 *
 *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
 *
 *           Unchanged on exit.
 *
 *  TRANS  - CHARACTER*1.
 *           On entry, TRANS specifies the operation to be performed as
 *           follows:
 *
 *              TRANS = 'N' or 'n'   x := A*x.
 *
 *              TRANS = 'T' or 't'   x := A'*x.
 *
 *              TRANS = 'C' or 'c'   x := A'*x.
 *
 *           Unchanged on exit.
 *
 *  DIAG   - CHARACTER*1.
 *           On entry, DIAG specifies whether or not A is unit
 *           triangular as follows:
 *
 *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
 *
 *              DIAG = 'N' or 'n'   A is not assumed to be unit
 *                                  triangular.
 *
 *           Unchanged on exit.
 *
 *  N      - INTEGER.
 *           On entry, N specifies the order of the matrix A.
 *           N must be at least zero.
 *           Unchanged on exit.
 *
 *  K      - INTEGER.
 *           On entry with UPLO = 'U' or 'u', K specifies the number of
 *           super-diagonals of the matrix A.
 *           On entry with UPLO = 'L' or 'l', K specifies the number of
 *           sub-diagonals of the matrix A.
 *           K must satisfy  0 .le. K.
 *           Unchanged on exit.
 *
 *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
 *           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
 *           by n part of the array A must contain the upper triangular
 *           band part of the matrix of coefficients, supplied column by
 *           column, with the leading diagonal of the matrix in row
 *           ( k + 1 ) of the array, the first super-diagonal starting at
 *           position 2 in row k, and so on. The top left k by k triangle
 *           of the array A is not referenced.
 *           The following program segment will transfer an upper
 *           triangular band matrix from conventional full matrix storage
 *           to band storage:
 *
 *                 DO 20, J = 1, N
 *                    M = K + 1 - J
 *                    DO 10, I = MAX( 1, J - K ), J
 *                       A( M + I, J ) = matrix( I, J )
 *              10    CONTINUE
 *              20 CONTINUE
 *
 *           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
 *           by n part of the array A must contain the lower triangular
 *           band part of the matrix of coefficients, supplied column by
 *           column, with the leading diagonal of the matrix in row 1 of
 *           the array, the first sub-diagonal starting at position 1 in
 *           row 2, and so on. The bottom right k by k triangle of the
 *           array A is not referenced.
 *           The following program segment will transfer a lower
 *           triangular band matrix from conventional full matrix storage
 *           to band storage:
 *
 *                 DO 20, J = 1, N
 *                    M = 1 - J
 *                    DO 10, I = J, MIN( N, J + K )
 *                       A( M + I, J ) = matrix( I, J )
 *              10    CONTINUE
 *              20 CONTINUE
 *
 *           Note that when DIAG = 'U' or 'u' the elements of the array A
 *           corresponding to the diagonal elements of the matrix are not
 *           referenced, but are assumed to be unity.
 *           Unchanged on exit.
 *
 *  LDA    - INTEGER.
 *           On entry, LDA specifies the first dimension of A as declared
 *           in the calling (sub) program. LDA must be at least
 *           ( k + 1 ).
 *           Unchanged on exit.
 *
 *  X      - DOUBLE PRECISION array of dimension at least
 *           ( 1 + ( n - 1 )*abs( INCX ) ).
 *           Before entry, the incremented array X must contain the n
 *           element vector x. On exit, X is overwritten with the
 *           tranformed vector x.
 *
 *  INCX   - INTEGER.
 *           On entry, INCX specifies the increment for the elements of
 *           X. INCX must not be zero.
 *           Unchanged on exit.
 *
 *
 *  Level 2 Blas routine.
 *
 *  -- Written on 22-October-1986.
 *     Jack Dongarra, Argonne National Lab.
 *     Jeremy Du Croz, Nag Central Office.
 *     Sven Hammarling, Nag Central Office.
 *     Richard Hanson, Sandia National Labs.
 *
 *
 *     .. Parameters ..
 */
PSI_API
void C_DTBMV(char uplo, char trans, char diag, int n, int k, double* a, int lda, double* x, int incx);

/**
 *  Purpose
 *  =======
 *
 *  DTBSV  solves one of the systems of equations
 *
 *     A*x = b,   or   A'*x = b,
 *
 *  where b and x are n element vectors and A is an n by n unit, or
 *  non-unit, upper or lower triangular band matrix, with ( k + 1 )
 *  diagonals.
 *
 *  No test for singularity or near-singularity is included in this
 *  routine. Such tests must be performed before calling this routine.
 *
 *  Arguments
 *  ==========
 *
 *  UPLO   - CHARACTER*1.
 *           On entry, UPLO specifies whether the matrix is an upper or
 *           lower triangular matrix as follows:
 *
 *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
 *
 *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
 *
 *           Unchanged on exit.
 *
 *  TRANS  - CHARACTER*1.
 *           On entry, TRANS specifies the equations to be solved as
 *           follows:
 *
 *              TRANS = 'N' or 'n'   A*x = b.
 *
 *              TRANS = 'T' or 't'   A'*x = b.
 *
 *              TRANS = 'C' or 'c'   A'*x = b.
 *
 *           Unchanged on exit.
 *
 *  DIAG   - CHARACTER*1.
 *           On entry, DIAG specifies whether or not A is unit
 *           triangular as follows:
 *
 *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
 *
 *              DIAG = 'N' or 'n'   A is not assumed to be unit
 *                                  triangular.
 *
 *           Unchanged on exit.
 *
 *  N      - INTEGER.
 *           On entry, N specifies the order of the matrix A.
 *           N must be at least zero.
 *           Unchanged on exit.
 *
 *  K      - INTEGER.
 *           On entry with UPLO = 'U' or 'u', K specifies the number of
 *           super-diagonals of the matrix A.
 *           On entry with UPLO = 'L' or 'l', K specifies the number of
 *           sub-diagonals of the matrix A.
 *           K must satisfy  0 .le. K.
 *           Unchanged on exit.
 *
 *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
 *           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
 *           by n part of the array A must contain the upper triangular
 *           band part of the matrix of coefficients, supplied column by
 *           column, with the leading diagonal of the matrix in row
 *           ( k + 1 ) of the array, the first super-diagonal starting at
 *           position 2 in row k, and so on. The top left k by k triangle
 *           of the array A is not referenced.
 *           The following program segment will transfer an upper
 *           triangular band matrix from conventional full matrix storage
 *           to band storage:
 *
 *                 DO 20, J = 1, N
 *                    M = K + 1 - J
 *                    DO 10, I = MAX( 1, J - K ), J
 *                       A( M + I, J ) = matrix( I, J )
 *              10    CONTINUE
 *              20 CONTINUE
 *
 *           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
 *           by n part of the array A must contain the lower triangular
 *           band part of the matrix of coefficients, supplied column by
 *           column, with the leading diagonal of the matrix in row 1 of
 *           the array, the first sub-diagonal starting at position 1 in
 *           row 2, and so on. The bottom right k by k triangle of the
 *           array A is not referenced.
 *           The following program segment will transfer a lower
 *           triangular band matrix from conventional full matrix storage
 *           to band storage:
 *
 *                 DO 20, J = 1, N
 *                    M = 1 - J
 *                    DO 10, I = J, MIN( N, J + K )
 *                       A( M + I, J ) = matrix( I, J )
 *              10    CONTINUE
 *              20 CONTINUE
 *
 *           Note that when DIAG = 'U' or 'u' the elements of the array A
 *           corresponding to the diagonal elements of the matrix are not
 *           referenced, but are assumed to be unity.
 *           Unchanged on exit.
 *
 *  LDA    - INTEGER.
 *           On entry, LDA specifies the first dimension of A as declared
 *           in the calling (sub) program. LDA must be at least
 *           ( k + 1 ).
 *           Unchanged on exit.
 *
 *  X      - DOUBLE PRECISION array of dimension at least
 *           ( 1 + ( n - 1 )*abs( INCX ) ).
 *           Before entry, the incremented array X must contain the n
 *           element right-hand side vector b. On exit, X is overwritten
 *           with the solution vector x.
 *
 *  INCX   - INTEGER.
 *           On entry, INCX specifies the increment for the elements of
 *           X. INCX must not be zero.
 *           Unchanged on exit.
 *
 *
 *  Level 2 Blas routine.
 *
 *  -- Written on 22-October-1986.
 *     Jack Dongarra, Argonne National Lab.
 *     Jeremy Du Croz, Nag Central Office.
 *     Sven Hammarling, Nag Central Office.
 *     Richard Hanson, Sandia National Labs.
 *
 *
 *     .. Parameters ..
 */
PSI_API
void C_DTBSV(char uplo, char trans, char diag, int n, int k, double* a, int lda, double* x, int incx);

/**
 *  Purpose
 *  =======
 *
 *  DTPMV  performs one of the matrix-vector operations
 *
 *     x := A*x,   or   x := A'*x,
 *
 *  where x is an n element vector and  A is an n by n unit, or non-unit,
 *  upper or lower triangular matrix, supplied in packed form.
 *
 *  Arguments
 *  ==========
 *
 *  UPLO   - CHARACTER*1.
 *           On entry, UPLO specifies whether the matrix is an upper or
 *           lower triangular matrix as follows:
 *
 *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
 *
 *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
 *
 *           Unchanged on exit.
 *
 *  TRANS  - CHARACTER*1.
 *           On entry, TRANS specifies the operation to be performed as
 *           follows:
 *
 *              TRANS = 'N' or 'n'   x := A*x.
 *
 *              TRANS = 'T' or 't'   x := A'*x.
 *
 *              TRANS = 'C' or 'c'   x := A'*x.
 *
 *           Unchanged on exit.
 *
 *  DIAG   - CHARACTER*1.
 *           On entry, DIAG specifies whether or not A is unit
 *           triangular as follows:
 *
 *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
 *
 *              DIAG = 'N' or 'n'   A is not assumed to be unit
 *                                  triangular.
 *
 *           Unchanged on exit.
 *
 *  N      - INTEGER.
 *           On entry, N specifies the order of the matrix A.
 *           N must be at least zero.
 *           Unchanged on exit.
 *
 *  AP     - DOUBLE PRECISION array of DIMENSION at least
 *           ( ( n*( n + 1 ) )/2 ).
 *           Before entry with  UPLO = 'U' or 'u', the array AP must
 *           contain the upper triangular matrix packed sequentially,
 *           column by column, so that AP( 1 ) contains a( 1, 1 ),
 *           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )
 *           respectively, and so on.
 *           Before entry with UPLO = 'L' or 'l', the array AP must
 *           contain the lower triangular matrix packed sequentially,
 *           column by column, so that AP( 1 ) contains a( 1, 1 ),
 *           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )
 *           respectively, and so on.
 *           Note that when  DIAG = 'U' or 'u', the diagonal elements of
 *           A are not referenced, but are assumed to be unity.
 *           Unchanged on exit.
 *
 *  X      - DOUBLE PRECISION array of dimension at least
 *           ( 1 + ( n - 1 )*abs( INCX ) ).
 *           Before entry, the incremented array X must contain the n
 *           element vector x. On exit, X is overwritten with the
 *           tranformed vector x.
 *
 *  INCX   - INTEGER.
 *           On entry, INCX specifies the increment for the elements of
 *           X. INCX must not be zero.
 *           Unchanged on exit.
 *
 *
 *  Level 2 Blas routine.
 *
 *  -- Written on 22-October-1986.
 *     Jack Dongarra, Argonne National Lab.
 *     Jeremy Du Croz, Nag Central Office.
 *     Sven Hammarling, Nag Central Office.
 *     Richard Hanson, Sandia National Labs.
 *
 *
 *     .. Parameters ..
 */
PSI_API
void C_DTPMV(char uplo, char trans, char diag, int n, double* ap, double* x, int incx);

/**
 *  Purpose
 *  =======
 *
 *  DTPSV  solves one of the systems of equations
 *
 *     A*x = b,   or   A'*x = b,
 *
 *  where b and x are n element vectors and A is an n by n unit, or
 *  non-unit, upper or lower triangular matrix, supplied in packed form.
 *
 *  No test for singularity or near-singularity is included in this
 *  routine. Such tests must be performed before calling this routine.
 *
 *  Arguments
 *  ==========
 *
 *  UPLO   - CHARACTER*1.
 *           On entry, UPLO specifies whether the matrix is an upper or
 *           lower triangular matrix as follows:
 *
 *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
 *
 *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
 *
 *           Unchanged on exit.
 *
 *  TRANS  - CHARACTER*1.
 *           On entry, TRANS specifies the equations to be solved as
 *           follows:
 *
 *              TRANS = 'N' or 'n'   A*x = b.
 *
 *              TRANS = 'T' or 't'   A'*x = b.
 *
 *              TRANS = 'C' or 'c'   A'*x = b.
 *
 *           Unchanged on exit.
 *
 *  DIAG   - CHARACTER*1.
 *           On entry, DIAG specifies whether or not A is unit
 *           triangular as follows:
 *
 *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
 *
 *              DIAG = 'N' or 'n'   A is not assumed to be unit
 *                                  triangular.
 *
 *           Unchanged on exit.
 *
 *  N      - INTEGER.
 *           On entry, N specifies the order of the matrix A.
 *           N must be at least zero.
 *           Unchanged on exit.
 *
 *  AP     - DOUBLE PRECISION array of DIMENSION at least
 *           ( ( n*( n + 1 ) )/2 ).
 *           Before entry with  UPLO = 'U' or 'u', the array AP must
 *           contain the upper triangular matrix packed sequentially,
 *           column by column, so that AP( 1 ) contains a( 1, 1 ),
 *           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )
 *           respectively, and so on.
 *           Before entry with UPLO = 'L' or 'l', the array AP must
 *           contain the lower triangular matrix packed sequentially,
 *           column by column, so that AP( 1 ) contains a( 1, 1 ),
 *           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )
 *           respectively, and so on.
 *           Note that when  DIAG = 'U' or 'u', the diagonal elements of
 *           A are not referenced, but are assumed to be unity.
 *           Unchanged on exit.
 *
 *  X      - DOUBLE PRECISION array of dimension at least
 *           ( 1 + ( n - 1 )*abs( INCX ) ).
 *           Before entry, the incremented array X must contain the n
 *           element right-hand side vector b. On exit, X is overwritten
 *           with the solution vector x.
 *
 *  INCX   - INTEGER.
 *           On entry, INCX specifies the increment for the elements of
 *           X. INCX must not be zero.
 *           Unchanged on exit.
 *
 *
 *  Level 2 Blas routine.
 *
 *  -- Written on 22-October-1986.
 *     Jack Dongarra, Argonne National Lab.
 *     Jeremy Du Croz, Nag Central Office.
 *     Sven Hammarling, Nag Central Office.
 *     Richard Hanson, Sandia National Labs.
 *
 *
 *     .. Parameters ..
 */
PSI_API
void C_DTPSV(char uplo, char trans, char diag, int n, double* ap, double* x, int incx);

/**
 *  Purpose
 *  =======
 *
 *  DTRMV  performs one of the matrix-vector operations
 *
 *     x := A*x,   or   x := A'*x,
 *
 *  where x is an n element vector and  A is an n by n unit, or non-unit,
 *  upper or lower triangular matrix.
 *
 *  Arguments
 *  ==========
 *
 *  UPLO   - CHARACTER*1.
 *           On entry, UPLO specifies whether the matrix is an upper or
 *           lower triangular matrix as follows:
 *
 *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
 *
 *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
 *
 *           Unchanged on exit.
 *
 *  TRANS  - CHARACTER*1.
 *           On entry, TRANS specifies the operation to be performed as
 *           follows:
 *
 *              TRANS = 'N' or 'n'   x := A*x.
 *
 *              TRANS = 'T' or 't'   x := A'*x.
 *
 *              TRANS = 'C' or 'c'   x := A'*x.
 *
 *           Unchanged on exit.
 *
 *  DIAG   - CHARACTER*1.
 *           On entry, DIAG specifies whether or not A is unit
 *           triangular as follows:
 *
 *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
 *
 *              DIAG = 'N' or 'n'   A is not assumed to be unit
 *                                  triangular.
 *
 *           Unchanged on exit.
 *
 *  N      - INTEGER.
 *           On entry, N specifies the order of the matrix A.
 *           N must be at least zero.
 *           Unchanged on exit.
 *
 *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
 *           Before entry with  UPLO = 'U' or 'u', the leading n by n
 *           upper triangular part of the array A must contain the upper
 *           triangular matrix and the strictly lower triangular part of
 *           A is not referenced.
 *           Before entry with UPLO = 'L' or 'l', the leading n by n
 *           lower triangular part of the array A must contain the lower
 *           triangular matrix and the strictly upper triangular part of
 *           A is not referenced.
 *           Note that when  DIAG = 'U' or 'u', the diagonal elements of
 *           A are not referenced either, but are assumed to be unity.
 *           Unchanged on exit.
 *
 *  LDA    - INTEGER.
 *           On entry, LDA specifies the first dimension of A as declared
 *           in the calling (sub) program. LDA must be at least
 *           max( 1, n ).
 *           Unchanged on exit.
 *
 *  X      - DOUBLE PRECISION array of dimension at least
 *           ( 1 + ( n - 1 )*abs( INCX ) ).
 *           Before entry, the incremented array X must contain the n
 *           element vector x. On exit, X is overwritten with the
 *           tranformed vector x.
 *
 *  INCX   - INTEGER.
 *           On entry, INCX specifies the increment for the elements of
 *           X. INCX must not be zero.
 *           Unchanged on exit.
 *
 *
 *  Level 2 Blas routine.
 *
 *  -- Written on 22-October-1986.
 *     Jack Dongarra, Argonne National Lab.
 *     Jeremy Du Croz, Nag Central Office.
 *     Sven Hammarling, Nag Central Office.
 *     Richard Hanson, Sandia National Labs.
 *
 *
 *     .. Parameters ..
 */
PSI_API
void C_DTRMV(char uplo, char trans, char diag, int n, double* a, int lda, double* x, int incx);

/**
 *  Purpose
 *  =======
 *
 *  DTRSM  solves one of the matrix equations
 *
 *     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
 *
 *  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
 *  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
 *
 *     op( A ) = A   or   op( A ) = A'.
 *
 *  The matrix X is overwritten on B.
 *
 *  Arguments
 *  ==========
 *
 *  SIDE   - CHARACTER*1.
 *           On entry, SIDE specifies whether op( A ) appears on the left
 *           or right of X as follows:
 *
 *              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
 *
 *              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
 *
 *           Unchanged on exit.
 *
 *  UPLO   - CHARACTER*1.
 *           On entry, UPLO specifies whether the matrix A is an upper or
 *           lower triangular matrix as follows:
 *
 *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
 *
 *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
 *
 *           Unchanged on exit.
 *
 *  TRANSA - CHARACTER*1.
 *           On entry, TRANSA specifies the form of op( A ) to be used in
 *           the matrix multiplication as follows:
 *
 *              TRANSA = 'N' or 'n'   op( A ) = A.
 *
 *              TRANSA = 'T' or 't'   op( A ) = A'.
 *
 *              TRANSA = 'C' or 'c'   op( A ) = A'.
 *
 *           Unchanged on exit.
 *
 *  DIAG   - CHARACTER*1.
 *           On entry, DIAG specifies whether or not A is unit triangular
 *           as follows:
 *
 *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
 *
 *              DIAG = 'N' or 'n'   A is not assumed to be unit
 *                                  triangular.
 *
 *           Unchanged on exit.
 *
 *  M      - INTEGER.
 *           On entry, M specifies the number of rows of B. M must be at
 *           least zero.
 *           Unchanged on exit.
 *
 *  N      - INTEGER.
 *           On entry, N specifies the number of columns of B.  N must be
 *           at least zero.
 *           Unchanged on exit.
 *
 *  ALPHA  - DOUBLE PRECISION.
 *           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
 *           zero then  A is not referenced and  B need not be set before
 *           entry.
 *           Unchanged on exit.
 *
 *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
 *           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
 *           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
 *           upper triangular part of the array  A must contain the upper
 *           triangular matrix  and the strictly lower triangular part of
 *           A is not referenced.
 *           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
 *           lower triangular part of the array  A must contain the lower
 *           triangular matrix  and the strictly upper triangular part of
 *           A is not referenced.
 *           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
 *           A  are not referenced either,  but are assumed to be  unity.
 *           Unchanged on exit.
 *
 *  LDA    - INTEGER.
 *           On entry, LDA specifies the first dimension of A as declared
 *           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
 *           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
 *           then LDA must be at least max( 1, n ).
 *           Unchanged on exit.
 *
 *  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
 *           Before entry,  the leading  m by n part of the array  B must
 *           contain  the  right-hand  side  matrix  B,  and  on exit  is
 *           overwritten by the solution matrix  X.
 *
 *  LDB    - INTEGER.
 *           On entry, LDB specifies the first dimension of B as declared
 *           in  the  calling  (sub)  program.   LDB  must  be  at  least
 *           max( 1, m ).
 *           Unchanged on exit.
 *
 *
 *  Level 3 Blas routine.
 *
 *
 *  -- Written on 8-February-1989.
 *     Jack Dongarra, Argonne National Laboratory.
 *     Iain Duff, AERE Harwell.
 *     Jeremy Du Croz, Numerical Algorithms Group Ltd.
 *     Sven Hammarling, Numerical Algorithms Group Ltd.
 *
 *
 *     .. External Functions ..
 */
PSI_API
void C_DTRSM(char side, char uplo, char transa, char diag, int m, int n, double alpha, double* a, int lda, double* b,
             int ldb);

/**
 *  Purpose
 *  =======
 *
 *  DTRSV  solves one of the systems of equations
 *
 *     A*x = b,   or   A'*x = b,
 *
 *  where b and x are n element vectors and A is an n by n unit, or
 *  non-unit, upper or lower triangular matrix.
 *
 *  No test for singularity or near-singularity is included in this
 *  routine. Such tests must be performed before calling this routine.
 *
 *  Arguments
 *  ==========
 *
 *  UPLO   - CHARACTER*1.
 *           On entry, UPLO specifies whether the matrix is an upper or
 *           lower triangular matrix as follows:
 *
 *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
 *
 *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
 *
 *           Unchanged on exit.
 *
 *  TRANS  - CHARACTER*1.
 *           On entry, TRANS specifies the equations to be solved as
 *           follows:
 *
 *              TRANS = 'N' or 'n'   A*x = b.
 *
 *              TRANS = 'T' or 't'   A'*x = b.
 *
 *              TRANS = 'C' or 'c'   A'*x = b.
 *
 *           Unchanged on exit.
 *
 *  DIAG   - CHARACTER*1.
 *           On entry, DIAG specifies whether or not A is unit
 *           triangular as follows:
 *
 *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
 *
 *              DIAG = 'N' or 'n'   A is not assumed to be unit
 *                                  triangular.
 *
 *           Unchanged on exit.
 *
 *  N      - INTEGER.
 *           On entry, N specifies the order of the matrix A.
 *           N must be at least zero.
 *           Unchanged on exit.
 *
 *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
 *           Before entry with  UPLO = 'U' or 'u', the leading n by n
 *           upper triangular part of the array A must contain the upper
 *           triangular matrix and the strictly lower triangular part of
 *           A is not referenced.
 *           Before entry with UPLO = 'L' or 'l', the leading n by n
 *           lower triangular part of the array A must contain the lower
 *           triangular matrix and the strictly upper triangular part of
 *           A is not referenced.
 *           Note that when  DIAG = 'U' or 'u', the diagonal elements of
 *           A are not referenced either, but are assumed to be unity.
 *           Unchanged on exit.
 *
 *  LDA    - INTEGER.
 *           On entry, LDA specifies the first dimension of A as declared
 *           in the calling (sub) program. LDA must be at least
 *           max( 1, n ).
 *           Unchanged on exit.
 *
 *  X      - DOUBLE PRECISION array of dimension at least
 *           ( 1 + ( n - 1 )*abs( INCX ) ).
 *           Before entry, the incremented array X must contain the n
 *           element right-hand side vector b. On exit, X is overwritten
 *           with the solution vector x.
 *
 *  INCX   - INTEGER.
 *           On entry, INCX specifies the increment for the elements of
 *           X. INCX must not be zero.
 *           Unchanged on exit.
 *
 *
 *  Level 2 Blas routine.
 *
 *  -- Written on 22-October-1986.
 *     Jack Dongarra, Argonne National Lab.
 *     Jeremy Du Croz, Nag Central Office.
 *     Sven Hammarling, Nag Central Office.
 *     Richard Hanson, Sandia National Labs.
 *
 *
 *     .. Parameters ..
 */
PSI_API
void C_DTRSV(char uplo, char trans, char diag, int n, double* a, int lda, double* x, int incx);
}  // namespace psi
