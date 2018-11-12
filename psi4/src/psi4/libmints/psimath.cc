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

#include "psimath.h"

#include <memory>

#include "psi4/libqt/blas_level1.h"
#include "psi4/libqt/blas_level2.h"
#include "psi4/libqt/blas_level3.h"
#include "psi4/libqt/lapack.h"

#include "matrix.h"
#include "vector.h"

namespace psi {
/// PSI_DGBMV, a wrapper to C_DGBMV using objects
void PSI_DGBMV(int irrep, char trans, int m, int n, int kl, int ku, double alpha, SharedMatrix a, int lda,
               std::shared_ptr<Vector> x, int incx, double beta, std::shared_ptr<Vector> y, int incy) {
    C_DGBMV(trans, m, n, kl, ku, alpha, a->pointer(irrep)[0], lda, x->pointer(irrep), incx, beta, y->pointer(irrep),
            incy);
}
/// PSI_DGEMM, a wrapper to C_DGEMM using objects
void PSI_DGEMM(int irrep, char transa, char transb, int m, int n, int k, double alpha, SharedMatrix a, int lda,
               SharedMatrix b, int ldb, double beta, SharedMatrix c, int ldc) {
    C_DGEMM(transa, transb, m, n, k, alpha, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, beta,
            c->pointer(irrep)[0], ldc);
}
/// PSI_DGEMV, a wrapper to C_DGEMV using objects
void PSI_DGEMV(int irrep, char trans, int m, int n, double alpha, SharedMatrix a, int lda, std::shared_ptr<Vector> x,
               int incx, double beta, std::shared_ptr<Vector> y, int incy) {
    C_DGEMV(trans, m, n, alpha, a->pointer(irrep)[0], lda, x->pointer(irrep), incx, beta, y->pointer(irrep), incy);
}
/// PSI_DGER, a wrapper to C_DGER using objects
void PSI_DGER(int irrep, int m, int n, double alpha, std::shared_ptr<Vector> x, int incx, std::shared_ptr<Vector> y,
              int incy, SharedMatrix a, int lda) {
    C_DGER(m, n, alpha, x->pointer(irrep), incx, y->pointer(irrep), incy, a->pointer(irrep)[0], lda);
}
/// PSI_DSBMV, a wrapper to C_DSBMV using objects
void PSI_DSBMV(int irrep, char uplo, int n, int k, double alpha, SharedMatrix a, int lda, std::shared_ptr<Vector> x,
               int incx, double beta, std::shared_ptr<Vector> y, int incy) {
    C_DSBMV(uplo, n, k, alpha, a->pointer(irrep)[0], lda, x->pointer(irrep), incx, beta, y->pointer(irrep), incy);
}
/// PSI_DSYMM, a wrapper to C_DSYMM using objects
void PSI_DSYMM(int irrep, char side, char uplo, int m, int n, double alpha, SharedMatrix a, int lda, SharedMatrix b,
               int ldb, double beta, SharedMatrix c, int ldc) {
    C_DSYMM(side, uplo, m, n, alpha, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, beta, c->pointer(irrep)[0],
            ldc);
}
/// PSI_DSYMV, a wrapper to C_DSYMV using objects
void PSI_DSYMV(int irrep, char uplo, int n, double alpha, SharedMatrix a, int lda, std::shared_ptr<Vector> x, int incx,
               double beta, std::shared_ptr<Vector> y, int incy) {
    C_DSYMV(uplo, n, alpha, a->pointer(irrep)[0], lda, x->pointer(irrep), incx, beta, y->pointer(irrep), incy);
}
/// PSI_DSYR, a wrapper to C_DSYR using objects
void PSI_DSYR(int irrep, char uplo, int n, double alpha, std::shared_ptr<Vector> x, int incx, SharedMatrix a, int lda) {
    C_DSYR(uplo, n, alpha, x->pointer(irrep), incx, a->pointer(irrep)[0], lda);
}
/// PSI_DSYR2, a wrapper to C_DSYR2 using objects
void PSI_DSYR2(int irrep, char uplo, int n, double alpha, std::shared_ptr<Vector> x, int incx,
               std::shared_ptr<Vector> y, int incy, SharedMatrix a, int lda) {
    C_DSYR2(uplo, n, alpha, x->pointer(irrep), incx, y->pointer(irrep), incy, a->pointer(irrep)[0], lda);
}
/// PSI_DSYR2K, a wrapper to C_DSYR2K using objects
void PSI_DSYR2K(int irrep, char uplo, char trans, int n, int k, double alpha, SharedMatrix a, int lda, SharedMatrix b,
                int ldb, double beta, SharedMatrix c, int ldc) {
    C_DSYR2K(uplo, trans, n, k, alpha, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, beta, c->pointer(irrep)[0],
             ldc);
}
/// PSI_DSYRK, a wrapper to C_DSYRK using objects
void PSI_DSYRK(int irrep, char uplo, char trans, int n, int k, double alpha, SharedMatrix a, int lda, double beta,
               SharedMatrix c, int ldc) {
    C_DSYRK(uplo, trans, n, k, alpha, a->pointer(irrep)[0], lda, beta, c->pointer(irrep)[0], ldc);
}
/// PSI_DTBMV, a wrapper to C_DTBMV using objects
void PSI_DTBMV(int irrep, char uplo, char trans, char diag, int n, int k, SharedMatrix a, int lda,
               std::shared_ptr<Vector> x, int incx) {
    C_DTBMV(uplo, trans, diag, n, k, a->pointer(irrep)[0], lda, x->pointer(irrep), incx);
}
/// PSI_DTBSV, a wrapper to C_DTBSV using objects
void PSI_DTBSV(int irrep, char uplo, char trans, char diag, int n, int k, SharedMatrix a, int lda,
               std::shared_ptr<Vector> x, int incx) {
    C_DTBSV(uplo, trans, diag, n, k, a->pointer(irrep)[0], lda, x->pointer(irrep), incx);
}
/// PSI_DTRMM, a wrapper to C_DTRMM using objects
void PSI_DTRMM(int irrep, char side, char uplo, char transa, char diag, int m, int n, double alpha, SharedMatrix a,
               int lda, SharedMatrix b, int ldb) {
    C_DTRMM(side, uplo, transa, diag, m, n, alpha, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb);
}
/// PSI_DTRMV, a wrapper to C_DTRMV using objects
void PSI_DTRMV(int irrep, char uplo, char trans, char diag, int n, SharedMatrix a, int lda, std::shared_ptr<Vector> x,
               int incx) {
    C_DTRMV(uplo, trans, diag, n, a->pointer(irrep)[0], lda, x->pointer(irrep), incx);
}
/// PSI_DTRSM, a wrapper to C_DTRSM using objects
void PSI_DTRSM(int irrep, char side, char uplo, char transa, char diag, int m, int n, double alpha, SharedMatrix a,
               int lda, SharedMatrix b, int ldb) {
    C_DTRSM(side, uplo, transa, diag, m, n, alpha, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb);
}
/// PSI_DTRSV, a wrapper to C_DTRSV using objects
void PSI_DTRSV(int irrep, char uplo, char trans, char diag, int n, SharedMatrix a, int lda, std::shared_ptr<Vector> x,
               int incx) {
    C_DTRSV(uplo, trans, diag, n, a->pointer(irrep)[0], lda, x->pointer(irrep), incx);
}

/// PSI_DROT, a wrapper to C_DROT using objects
void PSI_DROT(int irrep, size_t n, std::shared_ptr<Vector> x, int incx, std::shared_ptr<Vector> y, int incy, double c,
              double s) {
    C_DROT(n, x->pointer(irrep), incx, y->pointer(irrep), incy, c, s);
}
/// PSI_DSWAP, a wrapper to C_DSWAP using objects
void PSI_DSWAP(int irrep, size_t n, std::shared_ptr<Vector> x, int incx, std::shared_ptr<Vector> y, int incy) {
    C_DSWAP(n, x->pointer(irrep), incx, y->pointer(irrep), incy);
}
/// PSI_DCOPY, a wrapper to C_DCOPY using objects
void PSI_DCOPY(int irrep, size_t n, std::shared_ptr<Vector> x, int incx, std::shared_ptr<Vector> y, int incy) {
    C_DCOPY(n, x->pointer(irrep), incx, y->pointer(irrep), incy);
}
/// PSI_DSCAL, a wrapper to C_DSCAL using objects
void PSI_DSCAL(int irrep, size_t n, double alpha, std::shared_ptr<Vector> x, int incx) {
    C_DSCAL(n, alpha, x->pointer(irrep), incx);
}
/// PSI_DAXPY, a wrapper to C_DAXPY using objects
void PSI_DAXPY(int irrep, size_t n, double alpha, std::shared_ptr<Vector> x, int incx, std::shared_ptr<Vector> y,
               int incy) {
    C_DAXPY(n, alpha, x->pointer(irrep), incx, y->pointer(irrep), incy);
}
/// PSI_DDOT, a wrapper to C_DDOT using objects
double PSI_DDOT(int irrep, size_t n, std::shared_ptr<Vector> x, int incx, std::shared_ptr<Vector> y, int incy) {
    return C_DDOT(n, x->pointer(irrep), incx, y->pointer(irrep), incy);
}
/// PSI_DNRM2, a wrapper to C_DNRM2 using objects
double PSI_DNRM2(int irrep, size_t n, std::shared_ptr<Vector> x, int incx) {
    return C_DNRM2(n, x->pointer(irrep), incx);
}
/// PSI_DASUM, a wrapper to C_DASUM using objects
double PSI_DASUM(int irrep, size_t n, std::shared_ptr<Vector> x, int incx) {
    return C_DASUM(n, x->pointer(irrep), incx);
}
/// PSI_IDAMAX, a wrapper to C_IDAMAX using objects
size_t PSI_IDAMAX(int irrep, size_t n, std::shared_ptr<Vector> x, int incx) {
    return C_IDAMAX(n, x->pointer(irrep), incx);
}

/// Selection of LAPACK subroutines
/// PSI_DGEEV, a wrapper to return C_DGEEV using objects
int PSI_DGEEV(int irrep, char jobvl, char jobvr, int n, SharedMatrix a, int lda, std::shared_ptr<Vector> wr,
              std::shared_ptr<Vector> wi, SharedMatrix vl, int ldvl, SharedMatrix vr, int ldvr) {
    return C_DGEEV(jobvl, jobvr, n, a->pointer(irrep)[0], lda, wr->pointer(irrep), wi->pointer(irrep),
                   vl->pointer(irrep)[0], ldvl, vr->pointer(irrep)[0], ldvr);
}
/// PSI_DSYEV, a wrapper to return C_DSYEV using objects
int PSI_DSYEV(int irrep, char jobz, char uplo, int n, SharedMatrix a, int lda, std::shared_ptr<Vector> w) {
    return C_DSYEV(jobz, uplo, n, a->pointer(irrep)[0], lda, w->pointer(irrep));
}
/// PSI_DSYSV, a wrapper to return C_DSYSV using objects
int PSI_DSYSV(int irrep, char uplo, int n, int nrhs, SharedMatrix a, int lda, std::shared_ptr<IntVector> ipiv,
              SharedMatrix b, int ldb) {
    return C_DSYSV(uplo, n, nrhs, a->pointer(irrep)[0], lda, ipiv->pointer(irrep), b->pointer(irrep)[0], ldb);
}
/// PSI_DGETRF, a wrapper to return C_DGETRF using objects
int PSI_DGETRF(int irrep, int m, int n, SharedMatrix a, int lda, std::shared_ptr<IntVector> ipiv) {
    return C_DGETRF(m, n, a->pointer(irrep)[0], lda, ipiv->pointer(irrep));
}
/// PSI_DGETRI, a wrapper to return C_DGETRI using objects
int PSI_DGETRI(int irrep, int n, SharedMatrix a, int lda, std::shared_ptr<IntVector> ipiv, std::shared_ptr<Vector> work,
               int lwork) {
    return C_DGETRI(n, a->pointer(irrep)[0], lda, ipiv->pointer(irrep));
}
/// PSI_DGETRS, a wrapper to return C_DGETRS using objects
int PSI_DGETRS(int irrep, char trans, int n, int nrhs, SharedMatrix a, int lda, std::shared_ptr<IntVector> ipiv,
               SharedMatrix b, int ldb) {
    return C_DGETRS(trans, n, nrhs, a->pointer(irrep)[0], lda, ipiv->pointer(irrep), b->pointer(irrep)[0], ldb);
}
/// PSI_DPOTRF, a wrapper to return C_DPOTRF using objects
int PSI_DPOTRF(int irrep, char uplo, int n, SharedMatrix a, int lda) {
    return C_DPOTRF(uplo, n, a->pointer(irrep)[0], lda);
}
/// PSI_DPOTRI, a wrapper to return C_DPOTRI using objects
int PSI_DPOTRI(int irrep, char uplo, int n, SharedMatrix a, int lda) {
    return C_DPOTRI(uplo, n, a->pointer(irrep)[0], lda);
}
/// PSI_DPOTRS, a wrapper to return C_DPOTRS using objects
int PSI_DPOTRS(int irrep, char uplo, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix b, int ldb) {
    return C_DPOTRS(uplo, n, nrhs, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb);
}
}  // namespace psi
