/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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
** \file
** \brief Interface to LAPACK routines
** \ingroup QT
**
** Rollin A. King and T. Daniel Crawford
** August 2001 - January 2002
**
** 03/08/2002 EFV Added DGETRF since DGETRI isn't useful without it
**
** 05/29/2010 RMP added DPOTRF, DPOTRS and DPOTRI for density fitting applications
**
** Written to work similarly to the BLAS C interface in blas_intfc.c
*/

#include "lapack_intfc_mangle.h"

extern "C" {
extern int F_DGEEV(char *, char *, int *, double *, int *, double *, double *, double *, int *, double *, int *,
                   double *, int *, int *);

extern int F_DGESV(int *, int *, double *, int *, int *, double *, int *, int *);

extern int F_DGETRF(int *, int *, double *, int *, int *, int *);

extern int F_DGETRI(int *, double *, int *, int *, double *, int *, int *);

extern int F_DPOTRF(char *, int *, double *, int *, int *);

extern int F_DPOTRS(char *, int *, int *, double *, int *, double *, int *, int *);

extern int F_DPOTRI(char *, int *, double *, int *, int *);

extern int F_DGESVD(char *, char *, int *, int *, double *, int *, double *, double *, int *, double *, int *, double *,
                    int *, int *);

extern int F_DSYEV(char *, char *, int *, double *, int *, double *, double *, int *, int *);
}

namespace psi {
/*!
** C_DGEEV()
**
** This function computes the eigenvalues and the left and right
** right eigenvectors of a real, nonsymmetric matrix A.  For
** symmetric matrices, refer to C_DSYEV().
**
** \param n = The order of the matrix A.  n >= 0.
**
** \param a =  The n-by-n matrix A.  As with all other lapack
** routines, must be allocated by contiguous memory (i.e., block_matrix()).
**
** \param lda = The leading dimension of matrix A.  lda >= max(1,n).
**
** \param wr = array of length n containing the real parts of
** computed eigenvalues.
**
** \param wi = array of length n containing the imaginary parts of
** computed eigenvalues.
**
** \param vl = matrix of dimensions ldvl*n.  The columns store
** the left eigenvectors u(j).  If the j-th eigenvalues is real, then
** u(j) = vl(:,j), the j-th column of vl.  If the j-th and (j+1)-st
** eigenvalues form a complex conjugate pair, then u(j) = vl(:,j) +
** i*vl(:,j+1) and u(j+1) = vl(:,j) - i*vl(:,j+1).  Note: this is
** the Fortran documentation, may need to change cols <-> rows.
**
** \param ldvl = The leading dimension of matrix vl.
**
** \param vr = matrix of dimensions ldvr*n.  The columns store
** the right eigenvectors v(j).  If the j-th eigenvalues is real, then
** v(j) = vr(:,j), the j-th column of vr.  If the j-th and (j+1)-st
** eigenvalues form a complex conjugate pair, then v(j) = vr(:,j) +
** i*vr(:,j+1) and v(j+1) = vr(:,j) - i*vr(:,j+1).  Note: this is
** the Fortran documentation, may need to change cols <-> rows.
**
** \param ldvr = The leading dimension of matrix vr.
**
** \param work =  Array for scratch computations, of dimension
** lwork.  On successful exit, work[0] returns the optimal value of lwork.
**
** \param lwork = The dimension of the array work.  lwork >= max(1,3*n),
** and if eigenvectors are required (default for this wrapper at present)
** then actually lwork >= 4*n.  For good performance, lwork must generally
** be larger.  If lwork = -1, then a workspace query is assumed.  The
** routine only calculates the optimal size of the work array, returns
** this value ans the first entry of the work array, and no error
** message related to lword is issued by xerbla.
**
** \param info = On output (returned by C_DGEEV), a status flag.
** info = 0 for successful exit, if info = -i, the ith argument had
** an illegal value.  If info = i, the QR algorithm failed to
** compute all the eigenvalues, and no eigenvectors have been computed.
** Elements i+1:n of wr and wi contain eigenvalues which have converged.
**
** \ingroup QT
*/
int C_DGEEV(int n, double **a, int lda, double *wr, double *wi, double **vl, int ldvl, double **vr, int ldvr,
            double *work, int lwork, int info) {
    char jobvl, jobvr;
    jobvl = 'V';
    jobvr = 'V';
    ::F_DGEEV(&jobvl, &jobvr, &n, &(a[0][0]), &lda, &(wr[0]), &(wi[0]), &(vl[0][0]), &ldvl, &(vr[0][0]), &ldvr,
              &(work[0]), &lwork, &info);

    return info;
}

/*!
** C_DGESV()
**
** This function solves a system of linear equations A * X = B, where
** A is an n x n matrix and X and B are n x nrhs matrices.
**
** \param n = The number of linear equations, i.e., the order of
** the matrix A.  n >= 0.
**
** \param nrhs = The number of right hand sides, i.e., the number
** of columns of the matrix B.  nrhs >= 0.
**
** \param A = On entry, the n-by-n coefficient matrix A.  On
** exit, the factors L and U from the factorization A = P*L*U; the
** unit diagonal elements of L are not stored.
**
** \params lda = The leading dimension of the array A. lda >=
** max(1,n).
**
** \param int ipiv = An integer array of length n.  The pivot indices
** that define the permutation matrix P; row i of the matrix was
** interchanged with row ipiv(i).
**
** \param B = On entry, the n-by-nrhs matrix of right hand side
** matrix B.  On exit, if info = 0, the n-by-nrhs solution matrix X.
**
** \param ldb = The leading dimension of the array B.  ldb >=
** max(1,n).
**
** Returns: int info. info = 0: successful exit.  info < 0: if info = -i,
** the i-th argument had an illegal value.  info > 0: if info = i,
** U(i,i) is exactly zero.  The factorization has been completed, but the
** factor U is exactly singular, so the solution could not be computed.
**
** \ingroup QT
*/
int C_DGESV(int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb) {
    int info;

    ::F_DGESV(&n, &nrhs, &(a[0]), &lda, &(ipiv[0]), &(b[0]), &ldb, &info);

    return info;
}

/*!
** C_DGETRF(): Compute an LU factorization of a general M-by-N matrix a
** using partial pivoting with row interchanges.
**
** The factorization has the form
**   A = P * L * U
** where P is a permutation matrix, L is lower triangular with unit
** diagonal elements (lower trapezoidal if nrow > ncol), and U is upper
** triangular (upper trapezoidal if nrow < ncol).
**
** \param nrow = number of rows
** \param ncol = number of colums
** \param a    = matrix to factorize
** \param lda  = leading dimension of a, lda >= max(1,ncol)
** \param ipiv = output integer array of the pivot indices; for 1 <= i <=
**               min(nrow, ncol), col i of the matrix was interchanged with
**               col ipiv(i).
**
** Returns: int info.  info = 0: successful exit.  info < 0: if info=-i,
**   the ith-argument had an illegal vale.  info>0: if info=i, U(i,i) is
**   exactly zero.  The factorization has been completed, but the factor
**   U is exactly singular, and division by zero will occur if it is used
**   to solve a system of equations.
**
** \ingroup QT
*/
int C_DGETRF(int nrow, int ncol, double *a, int lda, int *ipiv) {
    int info;

    ::F_DGETRF(&ncol, &nrow, &(a[0]), &lda, &(ipiv[0]), &info);

    return info;
}

/*!
** C_DPOTRF(): Compute an Cholesky factorization of a Hermitian positive definite N-by-N matrix A.
**
** The factorization has the form
**   A = U**T* U if uplo = 'U'
**   or  A = L*L**T if oplo = 'L'
**
** \param uplo: uplo = 'U' upper triangle of A is stored.
**              uplo = 'L' lower triangle of A is stored
** \param n = number of rows and columns
** \param A = matrix to factorize
** \param lda  = leading dimension of a, lda >= max(1,ncol)
**
** Returns: int info.  info = 0: successful exit.  info < 0: if info=-i,
**   the ith-argument had an illegal vale.  info>0: if info=i, the i-th
**   minor is not positive definite
**
** \ingroup QT
*/
int C_DPOTRF(char uplo, int n, double *A, int lda) {
    int info;

    ::F_DPOTRF(&uplo, &n, &(A[0]), &lda, &info);

    return info;
}
/*!
** C_DGETRI(): computes the inverse of a matrix using the LU factorization
** computed by DGETRF
**
** This method inverts U and then computes inv(A) by solving the system
** inv(A)*L = inv(U) for inv(A).
**
** \param n     = order of the matrix a.  n >= 0.
** \param a     = array of dimension (n,lda).  On entry, the factors L and U
**                from the factorization A = P*L*U as computed by DGETRF.
**                On exit, if info=0, the inverse of hte original matrix A.
** \param lda   = the leading dimension of the array a.  lda >= max(1,n).
** \param ipiv  = The pivot indices from DGETRF.  For 1<=i<=n, row i of the
**                matrix was interchanged with row ipiv(i)
** \param work  = workspace of dimension max(1,lwork).  On exit, if info=0,
**                then work(0) returns the optimal lwork.
** \param lwork = the dimension of the array work.  lwork >= max(1,n).
**                For optimal performance, lwork > n*nb, where nb is the
**                optimal blocksize returned by ilaenv.  If lwork=-1, then
**                a workspace query is assumed; the routine only calculates
**                the optimal size of the work array, returns this value
**                as the first entry of the work array, and no error
**                message related to lwork is issued by xerbla.
**
** Returns: int info.  info=0: successful exit.  info<0: if info=-i, the
**   i-th argument had an illegal value.  info>0: if info=i, U(i,i) is
**   exactly zero; the matrix is singular and its inverse could not be
**   computed.
**
** \ingroup QT
*/
int C_DGETRI(int n, double *a, int lda, int *ipiv, double *work, int lwork) {
    int info;

    ::F_DGETRI(&n, &(a[0]), &lda, &(ipiv[0]), &(work[0]), &lwork, &info);

    return info;
}
/*!
** C_DPOTRI(): Compute inverse of a Hermitian positive deifinite matrix A AFTER a Cholesky
** decomposition is computed by C_DPORF.
**
**
** \param uplo: uplo = 'U' upper triangle of A is stored.
**              uplo = 'L' lower triangle of A is stored
** \param n = number of rows and columns
** \param A = matrix to factorize
** \param lda  = leading dimension of a, lda >= max(1,ncol)
**
** Returns: int info.  info = 0: successful exit.  info < 0: if info=-i,
**   the ith-argument had an illegal vale.  info>0: if info=i, the i-th
**   minor is no positive definite
**
** \ingroup QT
*/
int C_DPOTRI(char uplo, int n, double *A, int lda) {
    int info;

    ::F_DPOTRI(&uplo, &n, &(A[0]), &lda, &info);

    return info;
}
/*!
** C_DPOTRS(): Solve a Hermitian positive definite matrix Ax = B AFTER a Cholesky
** decomposition is computed for A by C_DPORF.
**
**
** \param uplo: uplo = 'U' upper triangle of A is stored.
**              uplo = 'L' lower triangle of A is stored
** \param n = number of rows and columns of A
** \param nrhs = number of vectors to solve in B
** \param A = matrix to solve
** \param lda  = leading dimension of a, lda >= max(1,ncol)
** \param B = (input)forcing vector/(output)solution vector
** \param ldb  = leading dimension of b, ldb >= max(1,ncol)
**
** Returns: int info.  info = 0: successful exit.  info < 0: if info=-i,
**   the ith-argument had an illegal value.
**
** \ingroup QT
*/
int C_DPOTRS(char uplo, int n, int nrhs, double *A, int lda, double *B, int ldb) {
    int info;

    ::F_DPOTRS(&uplo, &n, &nrhs, &(A[0]), &lda, &(B[0]), &ldb, &info);

    return info;
}

/*!
** C_DGESVD()
** This function computes the singular value decomposition (SVD) of a
** real mxn matrix A, ** optionally computing the left and/or right
** singular vectors.  The SVD is written
**
**  A = U * S * transpose(V)
**
** where S is an mxn matrix which is zero except for its min(m,n)
** diagonal elements, U is an M-by-M orthogonal matrix, and V is an
** N-by-N orthogonal matrix.  The diagonal elements of SIGMA are the
** singular values of A; they are real and non-negative, and are returned
** in descending order.  The first min(m,n) columns of U and V are the
** left and right singular vectors of A.
**
** Note that the routine returns V^t, not V;
**
** These arguments mimic their Fortran counterparts.  See the LAPACK
** manual for additional information.
**
** \param jobu   = 'A' = all m columns of U are returned
**                 'S' = the first min(m,n) columns of U are returned
**                 'O' = the first min(m,n) columns of U are
**                       returned in the input matrix A
**                 'N' = no columns of U are returned
**
** \param jobvt  = 'A' = all n rows of VT are returned
**                 'S' = the first min(m,n) rows of VT are
**                       returned
**                 'O' = the first min(m,n) rows of VT are
**                       returned in the input matrix A
**                 'N' = no rows of VT are returned
**
** Obviously jobu and jobvt cannot both be 'O'.
**
** \param m      = The row dimension of the matrix A.
**
** \param n      = The column dimension of the matrix A.
**
** \param A      = On entry, the mxn matrix A with dimensions m
**                 by lda.  On exit, if jobu='O' the first
**                 min(m,n) columns of A are overwritten with the
**                 left singular vectors; if jobvt='O' the first
**                 min(m,n) rows of A are overwritten with the
**                 right singular vectors; otherwise, the
**                 contents of A are destroyed.
**
** \param lda    = The number of columns allocated for A.
**
** \param s      = The singular values of A.
**
** \param u      = The right singular vectors of A, returned as
**                 column of the matrix u if jobu='A'.  If
**                 jobu='N' or 'O', u is not referenced.
**
** \param ldu    = The number of columns allocated for u.
**
** \param vt     = The left singular vectors of A, returned as
**                 rows of the matrix VT is jobvt='A'.  If
**                 jobvt='N' or 'O', vt is not referenced.
**
** \param ldvt   = The number of columns allocated for vt.
**
** \param work   = Workspace array of length lwork.
**
** \param lwork  = The length of the workspace array work, which
**                 should be at least as large as
**                 max(3*min(m,n)+max(m,n),5*min(m,n)).  For good
**                 performance, lwork should generally be
**                 larger.  If lwork=-1, a workspace query is
**                 assumed, and the value of work[0] upon return
**                 is the optimal size of the work array.
**
** Return value: int info.
**                 If info=0, successful exit.
**                 If info = -i, the i-th argument had an illegal value.
**                 If info > 0,  Related to failure in the convergence of
**                 the upper bidiagonal matrix B.  See the LAPACK
**                 manual for additiona information.
**
** Interface written by TDC, July 2001, updated April 2004
** \ingroup QT
*/
int C_DGESVD(char jobu, char jobvt, int m, int n, double *A, int lda, double *s, double *u, int ldu, double *vt,
             int ldvt, double *work, int lwork) {
    int info;

    ::F_DGESVD(&jobvt, &jobu, &n, &m, A, &lda, s, vt, &ldvt, u, &ldu, work, &lwork, &info);

    return info;
}

/*!
** C_DSYEV(): Computes all eigenvalues and, optionally, eigenvectors of
** a real symmetric matrix A.
**
** These arguments mimic their Fortran counterparts.
**
** \param jobz =   'N' or 'n' = compute eigenvalues only;
**                 'V' or 'v' = compute both eigenvalues and eigenvectors.
**
** \param uplo =   'U' or 'u' = A contains the upper triangular part
**                 of the matrix;
**                 'L' or 'l' = A contains the lower triangular part
**                 of the matrix.
**
** \param n    =   The order of the matrix A.
**
** \param A    =   On entry, the two-dimensional array with dimensions
**                 n by lda.
**                 On exit, if jobz = 'V', the columns of the matrix
**                 contain the eigenvectors of A, but if jobz = 'N',
**                 the contents of the matrix are destroyed.
**
** \param lda   =  The second dimension of A (i.e., the number of columns
**                 allocated for A).
**
** \param w     =  The computed eigenvalues in ascending order.
**
** \param work  =  An array of length lwork.  On exit, if the return
**                 value is 0, work[0] contains the optimal value of lwork.
**
** \param lwork =  The length of the array work.  A useful value of
**                 lwork seems to be 3*N.
**
** Returns:  0 = successful exit
**          <0 = the value of the i-th argument to the function was illegal
**          >0 = the algorithm failed to converge.
**
** Interface written by TDC, 10/2002
**
** \ingroup QT
*/
int C_DSYEV(char jobz, char uplo, int n, double *A, int lda, double *w, double *work, int lwork) {
    int info;

    ::F_DSYEV(&jobz, &uplo, &n, A, &lda, w, work, &lwork, &info);

    return info;
}

}  // namespace psi
