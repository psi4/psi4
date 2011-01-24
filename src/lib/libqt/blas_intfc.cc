/*! \defgroup QT libqt: The Quantum-Trio Miscellaneous Library */

/*!
 \file
 \ingroup QT
 \brief The PSI3 BLAS interface routines

 Interface to the BLAS routines

 C. David Sherrill
 Anna I. Krylov

 May 1998

 Additions by TD Crawford and EF Valeev, June 1999.

 Modifications to support BLASI calls with > 2B elements
 Rob Parrish, August 2010. The future has finally arrived.
*/

#include <cstdio>
#include <limits.h>
#include <cmath>

#if FC_SYMBOL==2
#define F_DSWAP dswap_
#define F_DAXPY daxpy_
#define F_DCOPY dcopy_
#define F_DGEMM dgemm_
#define F_DSYMM dsymm_
#define F_DROT drot_
#define F_DSCAL dscal_
#define F_DGEMV dgemv_
#define F_DSYMV dsymv_
#define F_DSPMV dspmv_
#define F_DDOT  ddot_
#define F_DASUM  dasum_ 
#define F_DNRM2  dnrm2_
#define F_IDAMAX  idamax_
#elif FC_SYMBOL==1
#define F_DSWAP dswap
#define F_DAXPY daxpy
#define F_DCOPY dcopy
#define F_DGEMM dgemm
#define F_DSYMM dsymm
#define F_DROT drot
#define F_DSCAL dscal
#define F_DGEMV dgemv
#define F_DSYMV dsymv
#define F_DSPMV dspmv
#define F_DDOT  ddot
#define F_DASUM  dasum 
#define F_DNRM2  dnrm2
#define F_IDAMAX  idamax
#elif FC_SYMBOL==3
#define F_DSWAP DSWAP
#define F_DAXPY DAXPY
#define F_DCOPY DCOPY
#define F_DGEMM DGEMM
#define F_DSYMM DSYMM
#define F_DROT DROT
#define F_DSCAL DSCAL
#define F_DGEMV DGEMV
#define F_DSYMV DSYMV
#define F_DSPMV DSPMV
#define F_DDOT  DDOT
#define F_DASUM  DASUM
#define F_DNRM2  DNRM2
#define F_IDAMAX  IDAMAX
#elif FC_SYMBOL==4
#define F_DSWAP DSWAP_
#define F_DAXPY DAXPY_
#define F_DCOPY DCOPY_
#define F_DGEMM DGEMM_
#define F_DSYMM DSYMM_
#define F_DROT DROT_
#define F_DSCAL DSCAL_
#define F_DGEMV DGEMV_
#define F_DSYMV DSYMV_
#define F_DSPMV DSPMV_
#define F_DDOT  DDOT_
#define F_DASUM  DASUM_
#define F_DNRM2  DNRM2_
#define F_IDAMAX  IDAMAX_
#endif

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

/*!
** This function calculates C(m,n)=alpha*(opT)A(m,k)*(opT)B(k,n)+
** beta*C(m,n)
**
** These arguments mimic their Fortran conterparts; parameters have been
** reversed (nca, ncb, ncc, A, B, C), to make it correct for C.
**
** \param transa On entry, specifies the form of (op)A used in the
**               matrix multiplication:
**               If transa = 'N' or 'n', (op)A = A.
**               If transa = 'T' or 't', (op)A = transp(A).
**               If transa = 'R' or 'r', (op)A = conjugate(A).
**               If transa = 'C' or 'c', (op)A = conjug_transp(A).
**               On exit, transa is unchanged.
**
** \param transb On entry, specifies the form of (op)B used in the
**               matrix multiplication:
**               If transb = 'N' or 'n', (op)B = B.
**               If transb = 'T' or 't', (op)B = transp(B).
**               If transb = 'R' or 'r', (op)B = conjugate(B)
**
** \param m      On entry, the number of rows of the matrix (op)A and of
**               the matrix C; m >= 0. On exit, m is unchanged.
**
** \param n      On entry, the number of columns of the matrix (op)B and
**               of the matrix C; n >= 0. On exit, n is unchanged.
**
** \param k      On entry, the number of columns of the matrix (op)A and
**               the number of rows of the matrix (op)B; k >= 0. On exit,
**               k is unchanged.
**
** \param alpha  On entry, specifies the scalar alpha. On exit, alpha is
**               unchanged.
**
** \param A      On entry, a two-dimensional array A with dimensions ka
**               by nca. For (op)A = A  or  conjugate(A), nca >= k and the
**               leading m by k portion of the array A contains the matrix
**               A. For (op)A = transp(A) or conjug_transp(A), nca >= m
**               and the leading k by m part of the array A contains the
**               matrix A. On exit, a is unchanged.
**
** \param nca    On entry, the second dimension of array A.
**               For (op)A = A  or conjugate(A), nca >= MAX(1,k).
**               For (op)A=transp(A) or conjug_transp(A), nca >= MAX(1,m).
**               On exit, nca is unchanged.
**
** \param B      On entry, a two-dimensional array B with dimensions kb
**               by ncb. For (op)B = B or conjugate(B), kb >= k and the
**               leading k by n portion of the array contains the matrix
**               B. For (op)B = transp(B) or conjug_transp(B), ncb >= k and
**               the leading n by k part of the array contains the matrix
**               B. On exit, B is unchanged.
**
** \param ncb    On entry, the second dimension of array B.
**               For (op)B = B or <conjugate(B), ncb >= MAX(1,n).
**               For (op)B = transp(B) or conjug_transp(B), ncb >=
**               MAX(1,k). On exit, ncb is unchanged.
**
** \param beta   On entry, specifies the scalar beta. On exit, beta is
**               unchanged.
**
** \param C      On entry, a two-dimensional array with the dimension
**               at least m by ncc. On exit,  the leading  m by n part of
**               array C is overwritten by the matrix alpha*(op)A*(op)B +
**               beta*C.
**
** \param ncc    On entry, the second dimension  of array C;
**               ncc >=MAX(1,n).  On exit, ncc is unchanged.
**
** \ingroup QT
*/
void C_DGEMM(char transa, char transb, int m, int n, int k, double alpha,
           double *A, int nca, double *B, int ncb, double beta, double *C,
           int ncc)
{

  /* the only strange thing we need to do is reverse everything
     since the stride runs differently in C vs. Fortran
   */

  /* also, do nothing if a dimension is 0 */
  if (m == 0 || n == 0 || k == 0) return;

  ::F_DGEMM(&transb,&transa,&n,&m,&k,&alpha,B,&ncb,A,&nca,&beta,C,&ncc);

}
/*
*  Purpose
*  =======
*  Note: C-style ordering implemented via wrapper
*  DSYMM  performs one of the matrix-matrix operations
*
*     C := alpha*A*B + beta*C,
*
*  or
*
*     C := alpha*B*A + beta*C,
*
*  where alpha and beta are scalars,  A is a symmetric matrix and  B and
*  C are  m by n matrices.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry,  SIDE  specifies whether  the  symmetric matrix  A
*           appears on the  left or right  in the  operation as follows:
*
*              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
*
*              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On  entry,   UPLO  specifies  whether  the  upper  or  lower
*           triangular  part  of  the  symmetric  matrix   A  is  to  be
*           referenced as follows:
*
*              UPLO = 'U' or 'u'   Only the upper triangular part of the
*                                  symmetric matrix is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the lower triangular part of the
*                                  symmetric matrix is to be referenced.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies the number of rows of the matrix  C.
*           M  must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix C.
*           N  must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           m  when  SIDE = 'L' or 'l'  and is  n otherwise.
*           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
*           the array  A  must contain the  symmetric matrix,  such that
*           when  UPLO = 'U' or 'u', the leading m by m upper triangular
*           part of the array  A  must contain the upper triangular part
*           of the  symmetric matrix and the  strictly  lower triangular
*           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
*           the leading  m by m  lower triangular part  of the  array  A
*           must  contain  the  lower triangular part  of the  symmetric
*           matrix and the  strictly upper triangular part of  A  is not
*           referenced.
*           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
*           the array  A  must contain the  symmetric matrix,  such that
*           when  UPLO = 'U' or 'u', the leading n by n upper triangular
*           part of the array  A  must contain the upper triangular part
*           of the  symmetric matrix and the  strictly  lower triangular
*           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
*           the leading  n by n  lower triangular part  of the  array  A
*           must  contain  the  lower triangular part  of the  symmetric
*           matrix and the  strictly upper triangular part of  A  is not
*           referenced.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
*           Before entry, the leading  m by n part of the array  B  must
*           contain the matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n updated
*           matrix.
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*/
void C_DSYMM(char side, char uplo, int m, int n, double alpha,
           double *A, int nca, double *B, int ncb, double beta, double *C,
           int ncc)
{

  /* the only strange thing we need to do is reverse everything
     since the stride runs differently in C vs. Fortran
   */

  /* also, do nothing if a dimension is 0 */
  if (m == 0 || n == 0 ) return;

  if (side == 'R' || side == 'r')
    side = 'L';
  else 
    side = 'R';
  if (uplo == 'U' || uplo == 'u')
    uplo = 'L';
  else 
    uplo = 'U';

  ::F_DSYMM(&side,&uplo,&n,&m,&alpha,A,&nca,B,&ncb,&beta,C,&ncc);

}
/*
*  Purpose
*  =======
*  Note: C-style indexing implemented in wrapper
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
*/
void C_DSYMV(char uplo, int n, double alpha, double *A,
             int nca, double *X, int inc_x, double beta, double *Y,
             int inc_y)
{
  if (n == 0) return;

  if(uplo == 'L' || uplo == 'l') uplo = 'U';
  else uplo = 'L';

  ::F_DSYMV(&uplo,&n,&alpha,A,&nca,X,&inc_x,&beta,Y,&inc_y);

}

/*!
** This function calculates the matrix-vector product.
**
** Y = alpha * A * X + beta * Y
**
** where X and Y are vectors, A is a matrix, and alpha and beta are
** constants.
**
** \param transa     =     Indicates whether the matrix A should be
**                         transposed ('t') or left alone ('n').
** \param m          =     The row dimension of A (regardless of transa).
** \param n          =     The column dimension of A (regardless of transa).
** \param alpha      =     The scalar alpha.
** \param A          =     A pointer to the beginning of the data in A.
** \param nca        =     The number of columns *actually* in A.  This is
**                         useful if one only wishes to multiply the first
**                         n columns of A times X even though A
**                         contains nca columns.
** \param X          =     A pointer to the beginning of the data in X.
** \param inc_x      =     The desired stride for X.  Useful for skipping
**                         sections of data to treat only one column of a
**                         complete matrix.  Usually 1, though.
** \param beta       =     The scalar beta.
** \param Y          =     A pointer to the beginning of the data in Y.
** \param inc_y      =     The desired stride for Y.
**
** \ingroup QT
*/
void C_DGEMV(char transa, int m, int n, double alpha, double *A,
             int nca, double *X, int inc_x, double beta, double *Y,
             int inc_y)
{
  if (m == 0 || n == 0) return;

  if(transa == 'n' || transa == 'N') transa = 't';
  else transa = 'n';

  ::F_DGEMV(&transa,&n,&m,&alpha,A,&nca,X,&inc_x,&beta,Y,&inc_y);

}


/*!
** This function calculates the matrix-vector product
**
**  Y = alpha * A * X + beta * Y
**
** where X and Y are vectors, A is a matrix, and alpha and beta are
** constants.
**
** \param uplo  Indicates whether the matrix A is packed in
**              upper ('U' or 'u') or lower ('L' or 'l')
**              triangular form.  We reverse what is passed
**              before sending it on to Fortran because of
**              the different Fortran/C conventions
** \param n     The order of the matrix A (number of rows/columns)
** \param alpha The scalar alpha.
** \param A     A pointer to the beginning of the data in A.
** \param X     A pointer to the beginning of the data in X.
** \param inc_x The desired stride for X.  Useful for skipping
**              sections of data to treat only one column of a
**              complete matrix.  Usually 1, though.
** \param beta  The scalar beta.
** \param Y     A pointer to the beginning of the data in Y.
** \param inc_y The desired stride for Y.
**
** \ingroup QT
*/
void C_DSPMV(char uplo, int n, double alpha, double *A,
             double *X, int inc_x, double beta, double *Y,
             int inc_y)
{
  if (n == 0) return;

  if (uplo != 'U' && uplo != 'u' && uplo != 'L' && uplo != 'l')
    fprintf(stderr, "C_DSPMV: called with unrecognized option for uplo!\n");

  if (uplo == 'U' || uplo == 'u') uplo = 'L';
  else uplo = 'U';

  ::F_DSPMV(&uplo,&n,&alpha,A,X,&inc_x,&beta,Y,&inc_y);

}


}

