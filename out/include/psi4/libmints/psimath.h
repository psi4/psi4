/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
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

#ifndef LIBMINTS_PSIMATH
#define LIBMINTS_PSIMATH

#include "typedefs.h"

/**
 * psimath.h
 * A set of wrappers to BLAS and LAPACK for use with Matrix, Vector, and IntVector
 * Rob Parrish 1/25/2010
 **/

namespace psi {

class Vector;
class IntVector;

/// PSI_DGBMV, a wrapper to C_DGBMV using objects
void PSI_DGBMV(int irrep, char trans, int m, int n, int kl, int ku, double alpha, SharedMatrix a, int lda,
               std::shared_ptr<Vector> x, int incx, double beta, std::shared_ptr<Vector> y, int incy);
/// PSI_DGEMM, a wrapper to C_DGEMM using objects
void PSI_DGEMM(int irrep, char transa, char transb, int m, int n, int k, double alpha, SharedMatrix a, int lda,
               SharedMatrix b, int ldb, double beta, SharedMatrix c, int ldc);
/// PSI_DGEMV, a wrapper to C_DGEMV using objects
void PSI_DGEMV(int irrep, char trans, int m, int n, double alpha, SharedMatrix a, int lda, std::shared_ptr<Vector> x,
               int incx, double beta, std::shared_ptr<Vector> y, int incy);
/// PSI_DGER, a wrapper to C_DGER using objects
void PSI_DGER(int irrep, int m, int n, double alpha, std::shared_ptr<Vector> x, int incx, std::shared_ptr<Vector> y,
              int incy, SharedMatrix a, int lda);
/// PSI_DSBMV, a wrapper to C_DSBMV using objects
void PSI_DSBMV(int irrep, char uplo, int n, int k, double alpha, SharedMatrix a, int lda, std::shared_ptr<Vector> x,
               int incx, double beta, std::shared_ptr<Vector> y, int incy);
/// PSI_DSYMM, a wrapper to C_DSYMM using objects
void PSI_DSYMM(int irrep, char side, char uplo, int m, int n, double alpha, SharedMatrix a, int lda, SharedMatrix b,
               int ldb, double beta, SharedMatrix c, int ldc);
/// PSI_DSYMV, a wrapper to C_DSYMV using objects
void PSI_DSYMV(int irrep, char uplo, int n, double alpha, SharedMatrix a, int lda, std::shared_ptr<Vector> x, int incx,
               double beta, std::shared_ptr<Vector> y, int incy);
/// PSI_DSYR, a wrapper to C_DSYR using objects
void PSI_DSYR(int irrep, char uplo, int n, double alpha, std::shared_ptr<Vector> x, int incx, SharedMatrix a, int lda);
/// PSI_DSYR2, a wrapper to C_DSYR2 using objects
void PSI_DSYR2(int irrep, char uplo, int n, double alpha, std::shared_ptr<Vector> x, int incx,
               std::shared_ptr<Vector> y, int incy, SharedMatrix a, int lda);
/// PSI_DSYR2K, a wrapper to C_DSYR2K using objects
void PSI_DSYR2K(int irrep, char uplo, char trans, int n, int k, double alpha, SharedMatrix a, int lda, SharedMatrix b,
                int ldb, double beta, SharedMatrix c, int ldc);
/// PSI_DSYRK, a wrapper to C_DSYRK using objects
void PSI_DSYRK(int irrep, char uplo, char trans, int n, int k, double alpha, SharedMatrix a, int lda, double beta,
               SharedMatrix c, int ldc);
/// PSI_DTBMV, a wrapper to C_DTBMV using objects
void PSI_DTBMV(int irrep, char uplo, char trans, char diag, int n, int k, SharedMatrix a, int lda,
               std::shared_ptr<Vector> x, int incx);
/// PSI_DTBSV, a wrapper to C_DTBSV using objects
void PSI_DTBSV(int irrep, char uplo, char trans, char diag, int n, int k, SharedMatrix a, int lda,
               std::shared_ptr<Vector> x, int incx);
/// PSI_DTRMM, a wrapper to C_DTRMM using objects
void PSI_DTRMM(int irrep, char side, char uplo, char transa, char diag, int m, int n, double alpha, SharedMatrix a,
               int lda, SharedMatrix b, int ldb);
/// PSI_DTRMV, a wrapper to C_DTRMV using objects
void PSI_DTRMV(int irrep, char uplo, char trans, char diag, int n, SharedMatrix a, int lda, std::shared_ptr<Vector> x,
               int incx);
/// PSI_DTRSM, a wrapper to C_DTRSM using objects
void PSI_DTRSM(int irrep, char side, char uplo, char transa, char diag, int m, int n, double alpha, SharedMatrix a,
               int lda, SharedMatrix b, int ldb);
/// PSI_DTRSV, a wrapper to C_DTRSV using objects
void PSI_DTRSV(int irrep, char uplo, char trans, char diag, int n, SharedMatrix a, int lda, std::shared_ptr<Vector> x,
               int incx);
/// PSI_DROT, a wrapper to C_DROT using objects
void PSI_DROT(int irrep, size_t n, std::shared_ptr<Vector> x, int incx, std::shared_ptr<Vector> y, int incy, double c,
              double s);
/// PSI_DSWAP, a wrapper to C_DSWAP using objects
void PSI_DSWAP(int irrep, size_t n, std::shared_ptr<Vector> x, int incx, std::shared_ptr<Vector> y, int incy);
/// PSI_DCOPY, a wrapper to C_DCOPY using objects
void PSI_DCOPY(int irrep, size_t n, std::shared_ptr<Vector> x, int incx, std::shared_ptr<Vector> y, int incy);
/// PSI_DSCAL, a wrapper to C_DSCAL using objects
void PSI_DSCAL(int irrep, size_t n, double alpha, std::shared_ptr<Vector> x, int incx);
/// PSI_DAXPY, a wrapper to C_DAXPY using objects
void PSI_DAXPY(int irrep, size_t n, double alpha, std::shared_ptr<Vector> x, int incx, std::shared_ptr<Vector> y,
               int incy);
/// PSI_DDOT, a wrapper to C_DDOT using objects
double PSI_DDOT(int irrep, size_t n, std::shared_ptr<Vector> x, int incx, std::shared_ptr<Vector> y, int incy);
/// PSI_DNRM2, a wrapper to C_DNRM2 using objects
double PSI_DNRM2(int irrep, size_t n, std::shared_ptr<Vector> x, int incx);
/// PSI_DASUM, a wrapper to C_DASUM using objects
double PSI_DASUM(int irrep, size_t n, std::shared_ptr<Vector> x, int incx);
/// PSI_IDAMAX, a wrapper to C_IDAMAX using objects
size_t PSI_IDAMAX(int irrep, size_t n, std::shared_ptr<Vector> x, int incx);
/// LAPACK
/// PSI_DBDSDC, a wrapper to return C_DBDSDC using objects
int PSI_DBDSDC(int irrep, char uplo, char compq, int n, std::shared_ptr<Vector> d, std::shared_ptr<Vector> e,
               SharedMatrix u, int ldu, SharedMatrix vt, int ldvt, std::shared_ptr<Vector> q,
               std::shared_ptr<IntVector> iq, std::shared_ptr<Vector> work, std::shared_ptr<IntVector> iwork);
/// PSI_DBDSQR, a wrapper to return C_DBDSQR using objects
int PSI_DBDSQR(int irrep, char uplo, int n, int ncvt, int nru, int ncc, std::shared_ptr<Vector> d,
               std::shared_ptr<Vector> e, SharedMatrix vt, int ldvt, SharedMatrix u, int ldu, SharedMatrix c, int ldc,
               std::shared_ptr<Vector> work);
/// PSI_DDISNA, a wrapper to return C_DDISNA using objects
int PSI_DDISNA(int irrep, char job, int m, int n, std::shared_ptr<Vector> d, std::shared_ptr<Vector> sep);
/// PSI_DGBBRD, a wrapper to return C_DGBBRD using objects
int PSI_DGBBRD(int irrep, char vect, int m, int n, int ncc, int kl, int ku, SharedMatrix ab, int ldab,
               std::shared_ptr<Vector> d, std::shared_ptr<Vector> e, SharedMatrix q, int ldq, SharedMatrix pt, int ldpt,
               SharedMatrix c, int ldc, std::shared_ptr<Vector> work);
/// PSI_DGBCON, a wrapper to return C_DGBCON using objects
int PSI_DGBCON(int irrep, char norm, int n, int kl, int ku, SharedMatrix ab, int ldab, std::shared_ptr<IntVector> ipiv,
               double anorm, std::shared_ptr<Vector> rcond, std::shared_ptr<Vector> work,
               std::shared_ptr<IntVector> iwork);
/// PSI_DGBEQU, a wrapper to return C_DGBEQU using objects
int PSI_DGBEQU(int irrep, int m, int n, int kl, int ku, SharedMatrix ab, int ldab, std::shared_ptr<Vector> r,
               std::shared_ptr<Vector> c, std::shared_ptr<Vector> rowcnd, std::shared_ptr<Vector> colcnd,
               std::shared_ptr<Vector> amax);
/// PSI_DGBRFS, a wrapper to return C_DGBRFS using objects
int PSI_DGBRFS(int irrep, char trans, int n, int kl, int ku, int nrhs, SharedMatrix ab, int ldab, SharedMatrix afb,
               int ldafb, std::shared_ptr<IntVector> ipiv, SharedMatrix b, int ldb, SharedMatrix x, int ldx,
               std::shared_ptr<Vector> ferr, std::shared_ptr<Vector> berr, std::shared_ptr<Vector> work,
               std::shared_ptr<IntVector> iwork);
/// PSI_DGBSV, a wrapper to return C_DGBSV using objects
int PSI_DGBSV(int irrep, int n, int kl, int ku, int nrhs, SharedMatrix ab, int ldab, std::shared_ptr<IntVector> ipiv,
              SharedMatrix b, int ldb);
/// PSI_DGBSVX, a wrapper to return C_DGBSVX using objects
int PSI_DGBSVX(int irrep, char fact, char trans, int n, int kl, int ku, int nrhs, SharedMatrix ab, int ldab,
               SharedMatrix afb, int ldafb, std::shared_ptr<IntVector> ipiv, char equed, std::shared_ptr<Vector> r,
               std::shared_ptr<Vector> c, SharedMatrix b, int ldb, SharedMatrix x, int ldx,
               std::shared_ptr<Vector> rcond, std::shared_ptr<Vector> ferr, std::shared_ptr<Vector> berr,
               std::shared_ptr<Vector> work, std::shared_ptr<IntVector> iwork);
/// PSI_DGBTRF, a wrapper to return C_DGBTRF using objects
int PSI_DGBTRF(int irrep, int m, int n, int kl, int ku, SharedMatrix ab, int ldab, std::shared_ptr<IntVector> ipiv);
/// PSI_DGBTRS, a wrapper to return C_DGBTRS using objects
int PSI_DGBTRS(int irrep, char trans, int n, int kl, int ku, int nrhs, SharedMatrix ab, int ldab,
               std::shared_ptr<IntVector> ipiv, SharedMatrix b, int ldb);
/// PSI_DGEBAK, a wrapper to return C_DGEBAK using objects
int PSI_DGEBAK(int irrep, char job, char side, int n, int ilo, int ihi, std::shared_ptr<Vector> scale, int m,
               SharedMatrix v, int ldv);
/// PSI_DGEBAL, a wrapper to return C_DGEBAL using objects
int PSI_DGEBAL(int irrep, char job, int n, SharedMatrix a, int lda, std::shared_ptr<IntVector> ilo,
               std::shared_ptr<IntVector> ihi, std::shared_ptr<Vector> scale);
/// PSI_DGEBRD, a wrapper to return C_DGEBRD using objects
int PSI_DGEBRD(int irrep, int m, int n, SharedMatrix a, int lda, std::shared_ptr<Vector> d, std::shared_ptr<Vector> e,
               std::shared_ptr<Vector> tauq, std::shared_ptr<Vector> taup, std::shared_ptr<Vector> work, int lwork);
/// PSI_DGECON, a wrapper to return C_DGECON using objects
int PSI_DGECON(int irrep, char norm, int n, SharedMatrix a, int lda, double anorm, std::shared_ptr<Vector> rcond,
               std::shared_ptr<Vector> work, std::shared_ptr<IntVector> iwork);
/// PSI_DGEEQU, a wrapper to return C_DGEEQU using objects
int PSI_DGEEQU(int irrep, int m, int n, SharedMatrix a, int lda, std::shared_ptr<Vector> r, std::shared_ptr<Vector> c,
               std::shared_ptr<Vector> rowcnd, std::shared_ptr<Vector> colcnd, std::shared_ptr<Vector> amax);
/// PSI_DGEES, a wrapper to return C_DGEES using objects
int PSI_DGEES(int irrep, char jobvs, char sort, int n, SharedMatrix a, int lda, std::shared_ptr<IntVector> sdim,
              std::shared_ptr<Vector> wr, std::shared_ptr<Vector> wi, SharedMatrix vs, int ldvs,
              std::shared_ptr<Vector> work, int lwork);
/// PSI_DGEESX, a wrapper to return C_DGEESX using objects
int PSI_DGEESX(int irrep, char jobvs, char sort, char sense, int n, SharedMatrix a, int lda,
               std::shared_ptr<IntVector> sdim, std::shared_ptr<Vector> wr, std::shared_ptr<Vector> wi, SharedMatrix vs,
               int ldvs, std::shared_ptr<Vector> rconde, std::shared_ptr<Vector> rcondv, std::shared_ptr<Vector> work,
               int lwork, std::shared_ptr<IntVector> iwork, int liwork);
/// PSI_DGEEV, a wrapper to return C_DGEEV using objects
int PSI_DGEEV(int irrep, char jobvl, char jobvr, int n, SharedMatrix a, int lda, std::shared_ptr<Vector> wr,
              std::shared_ptr<Vector> wi, SharedMatrix vl, int ldvl, SharedMatrix vr, int ldvr,
              std::shared_ptr<Vector> work, int lwork);
/// PSI_DGEEVX, a wrapper to return C_DGEEVX using objects
int PSI_DGEEVX(int irrep, char balanc, char jobvl, char jobvr, char sense, int n, SharedMatrix a, int lda,
               std::shared_ptr<Vector> wr, std::shared_ptr<Vector> wi, SharedMatrix vl, int ldvl, SharedMatrix vr,
               int ldvr, std::shared_ptr<IntVector> ilo, std::shared_ptr<IntVector> ihi, std::shared_ptr<Vector> scale,
               std::shared_ptr<Vector> abnrm, std::shared_ptr<Vector> rconde, std::shared_ptr<Vector> rcondv,
               std::shared_ptr<Vector> work, int lwork, std::shared_ptr<IntVector> iwork);
/// PSI_DGEGS, a wrapper to return C_DGEGS using objects
int PSI_DGEGS(int irrep, char jobvsl, char jobvsr, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb,
              std::shared_ptr<Vector> alphar, std::shared_ptr<Vector> alphai, std::shared_ptr<Vector> beta,
              SharedMatrix vsl, int ldvsl, SharedMatrix vsr, int ldvsr, std::shared_ptr<Vector> work, int lwork);
/// PSI_DGEGV, a wrapper to return C_DGEGV using objects
int PSI_DGEGV(int irrep, char jobvl, char jobvr, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb,
              std::shared_ptr<Vector> alphar, std::shared_ptr<Vector> alphai, std::shared_ptr<Vector> beta,
              SharedMatrix vl, int ldvl, SharedMatrix vr, int ldvr, std::shared_ptr<Vector> work, int lwork);
/// PSI_DGEHRD, a wrapper to return C_DGEHRD using objects
int PSI_DGEHRD(int irrep, int n, int ilo, int ihi, SharedMatrix a, int lda, std::shared_ptr<Vector> tau,
               std::shared_ptr<Vector> work, int lwork);
/// PSI_DGELQF, a wrapper to return C_DGELQF using objects
int PSI_DGELQF(int irrep, int m, int n, SharedMatrix a, int lda, std::shared_ptr<Vector> tau,
               std::shared_ptr<Vector> work, int lwork);
/// PSI_DGELS, a wrapper to return C_DGELS using objects
int PSI_DGELS(int irrep, char trans, int m, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix b, int ldb,
              std::shared_ptr<Vector> work, int lwork);
/// PSI_DGELSD, a wrapper to return C_DGELSD using objects
int PSI_DGELSD(int irrep, int m, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix b, int ldb,
               std::shared_ptr<Vector> s, double rcond, std::shared_ptr<IntVector> rank, std::shared_ptr<Vector> work,
               int lwork, std::shared_ptr<IntVector> iwork);
/// PSI_DGELSS, a wrapper to return C_DGELSS using objects
int PSI_DGELSS(int irrep, int m, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix b, int ldb,
               std::shared_ptr<Vector> s, double rcond, std::shared_ptr<IntVector> rank, std::shared_ptr<Vector> work,
               int lwork);
/// PSI_DGELSX, a wrapper to return C_DGELSX using objects
int PSI_DGELSX(int irrep, int m, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix b, int ldb,
               std::shared_ptr<IntVector> jpvt, double rcond, std::shared_ptr<IntVector> rank,
               std::shared_ptr<Vector> work);
/// PSI_DGELSY, a wrapper to return C_DGELSY using objects
int PSI_DGELSY(int irrep, int m, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix b, int ldb,
               std::shared_ptr<IntVector> jpvt, double rcond, std::shared_ptr<IntVector> rank,
               std::shared_ptr<Vector> work, int lwork);
/// PSI_DGEQLF, a wrapper to return C_DGEQLF using objects
int PSI_DGEQLF(int irrep, int m, int n, SharedMatrix a, int lda, std::shared_ptr<Vector> tau,
               std::shared_ptr<Vector> work, int lwork);
/// PSI_DGEQP3, a wrapper to return C_DGEQP3 using objects
int PSI_DGEQP3(int irrep, int m, int n, SharedMatrix a, int lda, std::shared_ptr<IntVector> jpvt,
               std::shared_ptr<Vector> tau, std::shared_ptr<Vector> work, int lwork);
/// PSI_DGEQPF, a wrapper to return C_DGEQPF using objects
int PSI_DGEQPF(int irrep, int m, int n, SharedMatrix a, int lda, std::shared_ptr<IntVector> jpvt,
               std::shared_ptr<Vector> tau, std::shared_ptr<Vector> work);
/// PSI_DGERFS, a wrapper to return C_DGERFS using objects
int PSI_DGERFS(int irrep, char trans, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix af, int ldaf,
               std::shared_ptr<IntVector> ipiv, SharedMatrix b, int ldb, SharedMatrix x, int ldx,
               std::shared_ptr<Vector> ferr, std::shared_ptr<Vector> berr, std::shared_ptr<Vector> work,
               std::shared_ptr<IntVector> iwork);
/// PSI_DGERQF, a wrapper to return C_DGERQF using objects
int PSI_DGERQF(int irrep, int m, int n, SharedMatrix a, int lda, std::shared_ptr<Vector> tau,
               std::shared_ptr<Vector> work, int lwork);
/// PSI_DGESDD, a wrapper to return C_DGESDD using objects
int PSI_DGESDD(int irrep, char jobz, int m, int n, SharedMatrix a, int lda, std::shared_ptr<Vector> s, SharedMatrix u,
               int ldu, SharedMatrix vt, int ldvt, std::shared_ptr<Vector> work, int lwork,
               std::shared_ptr<IntVector> iwork);
/// PSI_DGESV, a wrapper to return C_DGESV using objects
int PSI_DGESV(int irrep, int n, int nrhs, SharedMatrix a, int lda, std::shared_ptr<IntVector> ipiv, SharedMatrix b,
              int ldb);
/// PSI_DGESVX, a wrapper to return C_DGESVX using objects
int PSI_DGESVX(int irrep, char fact, char trans, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix af, int ldaf,
               std::shared_ptr<IntVector> ipiv, char equed, std::shared_ptr<Vector> r, std::shared_ptr<Vector> c,
               SharedMatrix b, int ldb, SharedMatrix x, int ldx, std::shared_ptr<Vector> rcond,
               std::shared_ptr<Vector> ferr, std::shared_ptr<Vector> berr, std::shared_ptr<Vector> work,
               std::shared_ptr<IntVector> iwork);
/// PSI_DGETRF, a wrapper to return C_DGETRF using objects
int PSI_DGETRF(int irrep, int m, int n, SharedMatrix a, int lda, std::shared_ptr<IntVector> ipiv);
/// PSI_DGETRI, a wrapper to return C_DGETRI using objects
int PSI_DGETRI(int irrep, int n, SharedMatrix a, int lda, std::shared_ptr<IntVector> ipiv, std::shared_ptr<Vector> work,
               int lwork);
/// PSI_DGETRS, a wrapper to return C_DGETRS using objects
int PSI_DGETRS(int irrep, char trans, int n, int nrhs, SharedMatrix a, int lda, std::shared_ptr<IntVector> ipiv,
               SharedMatrix b, int ldb);
/// PSI_DGGBAK, a wrapper to return C_DGGBAK using objects
int PSI_DGGBAK(int irrep, char job, char side, int n, int ilo, int ihi, std::shared_ptr<Vector> lscale,
               std::shared_ptr<Vector> rscale, int m, SharedMatrix v, int ldv);
/// PSI_DGGBAL, a wrapper to return C_DGGBAL using objects
int PSI_DGGBAL(int irrep, char job, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb,
               std::shared_ptr<IntVector> ilo, std::shared_ptr<IntVector> ihi, std::shared_ptr<Vector> lscale,
               std::shared_ptr<Vector> rscale, std::shared_ptr<Vector> work);
/// PSI_DGGES, a wrapper to return C_DGGES using objects
int PSI_DGGES(int irrep, char jobvsl, char jobvsr, char sort, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb,
              std::shared_ptr<IntVector> sdim, std::shared_ptr<Vector> alphar, std::shared_ptr<Vector> alphai,
              std::shared_ptr<Vector> beta, SharedMatrix vsl, int ldvsl, SharedMatrix vsr, int ldvsr,
              std::shared_ptr<Vector> work, int lwork);
/// PSI_DGGESX, a wrapper to return C_DGGESX using objects
int PSI_DGGESX(int irrep, char jobvsl, char jobvsr, char sort, char sense, int n, SharedMatrix a, int lda,
               SharedMatrix b, int ldb, std::shared_ptr<IntVector> sdim, std::shared_ptr<Vector> alphar,
               std::shared_ptr<Vector> alphai, std::shared_ptr<Vector> beta, SharedMatrix vsl, int ldvsl,
               SharedMatrix vsr, int ldvsr, std::shared_ptr<Vector> rconde, std::shared_ptr<Vector> rcondv,
               std::shared_ptr<Vector> work, int lwork, std::shared_ptr<IntVector> iwork, int liwork);
/// PSI_DGGEV, a wrapper to return C_DGGEV using objects
int PSI_DGGEV(int irrep, char jobvl, char jobvr, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb,
              std::shared_ptr<Vector> alphar, std::shared_ptr<Vector> alphai, std::shared_ptr<Vector> beta,
              SharedMatrix vl, int ldvl, SharedMatrix vr, int ldvr, std::shared_ptr<Vector> work, int lwork);
/// PSI_DGGEVX, a wrapper to return C_DGGEVX using objects
int PSI_DGGEVX(int irrep, char balanc, char jobvl, char jobvr, char sense, int n, SharedMatrix a, int lda,
               SharedMatrix b, int ldb, std::shared_ptr<Vector> alphar, std::shared_ptr<Vector> alphai,
               std::shared_ptr<Vector> beta, SharedMatrix vl, int ldvl, SharedMatrix vr, int ldvr,
               std::shared_ptr<IntVector> ilo, std::shared_ptr<IntVector> ihi, std::shared_ptr<Vector> lscale,
               std::shared_ptr<Vector> rscale, std::shared_ptr<Vector> abnrm, std::shared_ptr<Vector> bbnrm,
               std::shared_ptr<Vector> rconde, std::shared_ptr<Vector> rcondv, std::shared_ptr<Vector> work, int lwork,
               std::shared_ptr<IntVector> iwork);
/// PSI_DGGGLM, a wrapper to return C_DGGGLM using objects
int PSI_DGGGLM(int irrep, int n, int m, int p, SharedMatrix a, int lda, SharedMatrix b, int ldb,
               std::shared_ptr<Vector> d, std::shared_ptr<Vector> x, std::shared_ptr<Vector> y,
               std::shared_ptr<Vector> work, int lwork);
/// PSI_DGGHRD, a wrapper to return C_DGGHRD using objects
int PSI_DGGHRD(int irrep, char compq, char compz, int n, int ilo, int ihi, SharedMatrix a, int lda, SharedMatrix b,
               int ldb, SharedMatrix q, int ldq, SharedMatrix z, int ldz);
/// PSI_DGGLSE, a wrapper to return C_DGGLSE using objects
int PSI_DGGLSE(int irrep, int m, int n, int p, SharedMatrix a, int lda, SharedMatrix b, int ldb,
               std::shared_ptr<Vector> c, std::shared_ptr<Vector> d, std::shared_ptr<Vector> x,
               std::shared_ptr<Vector> work, int lwork);
/// PSI_DGGQRF, a wrapper to return C_DGGQRF using objects
int PSI_DGGQRF(int irrep, int n, int m, int p, SharedMatrix a, int lda, std::shared_ptr<Vector> taua, SharedMatrix b,
               int ldb, std::shared_ptr<Vector> taub, std::shared_ptr<Vector> work, int lwork);
/// PSI_DGGRQF, a wrapper to return C_DGGRQF using objects
int PSI_DGGRQF(int irrep, int m, int p, int n, SharedMatrix a, int lda, std::shared_ptr<Vector> taua, SharedMatrix b,
               int ldb, std::shared_ptr<Vector> taub, std::shared_ptr<Vector> work, int lwork);
/// PSI_DGGSVD, a wrapper to return C_DGGSVD using objects
int PSI_DGGSVD(int irrep, char jobu, char jobv, char jobq, int m, int n, int p, std::shared_ptr<IntVector> k,
               std::shared_ptr<IntVector> l, SharedMatrix a, int lda, SharedMatrix b, int ldb,
               std::shared_ptr<Vector> alpha, std::shared_ptr<Vector> beta, SharedMatrix u, int ldu, SharedMatrix v,
               int ldv, SharedMatrix q, int ldq, std::shared_ptr<Vector> work, std::shared_ptr<IntVector> iwork);
/// PSI_DGGSVP, a wrapper to return C_DGGSVP using objects
int PSI_DGGSVP(int irrep, char jobu, char jobv, char jobq, int m, int p, int n, SharedMatrix a, int lda, SharedMatrix b,
               int ldb, double tola, double tolb, std::shared_ptr<IntVector> k, std::shared_ptr<IntVector> l,
               SharedMatrix u, int ldu, SharedMatrix v, int ldv, SharedMatrix q, int ldq,
               std::shared_ptr<IntVector> iwork, std::shared_ptr<Vector> tau, std::shared_ptr<Vector> work);
/// PSI_DGTCON, a wrapper to return C_DGTCON using objects
int PSI_DGTCON(int irrep, char norm, int n, std::shared_ptr<Vector> dl, std::shared_ptr<Vector> d,
               std::shared_ptr<Vector> du, std::shared_ptr<Vector> du2, std::shared_ptr<IntVector> ipiv, double anorm,
               std::shared_ptr<Vector> rcond, std::shared_ptr<Vector> work, std::shared_ptr<IntVector> iwork);
/// PSI_DGTRFS, a wrapper to return C_DGTRFS using objects
int PSI_DGTRFS(int irrep, char trans, int n, int nrhs, std::shared_ptr<Vector> dl, std::shared_ptr<Vector> d,
               std::shared_ptr<Vector> du, std::shared_ptr<Vector> dlf, std::shared_ptr<Vector> df,
               std::shared_ptr<Vector> duf, std::shared_ptr<Vector> du2, std::shared_ptr<IntVector> ipiv,
               SharedMatrix b, int ldb, SharedMatrix x, int ldx, std::shared_ptr<Vector> ferr,
               std::shared_ptr<Vector> berr, std::shared_ptr<Vector> work, std::shared_ptr<IntVector> iwork);
/// PSI_DGTSV, a wrapper to return C_DGTSV using objects
int PSI_DGTSV(int irrep, int n, int nrhs, std::shared_ptr<Vector> dl, std::shared_ptr<Vector> d,
              std::shared_ptr<Vector> du, SharedMatrix b, int ldb);
/// PSI_DGTSVX, a wrapper to return C_DGTSVX using objects
int PSI_DGTSVX(int irrep, char fact, char trans, int n, int nrhs, std::shared_ptr<Vector> dl, std::shared_ptr<Vector> d,
               std::shared_ptr<Vector> du, std::shared_ptr<Vector> dlf, std::shared_ptr<Vector> df,
               std::shared_ptr<Vector> duf, std::shared_ptr<Vector> du2, std::shared_ptr<IntVector> ipiv,
               SharedMatrix b, int ldb, SharedMatrix x, int ldx, std::shared_ptr<Vector> rcond);
/// PSI_DGTTRF, a wrapper to return C_DGTTRF using objects
int PSI_DGTTRF(int irrep, int n, std::shared_ptr<Vector> dl, std::shared_ptr<Vector> d, std::shared_ptr<Vector> du,
               std::shared_ptr<Vector> du2, std::shared_ptr<IntVector> ipiv);
/// PSI_DGTTRS, a wrapper to return C_DGTTRS using objects
int PSI_DGTTRS(int irrep, char trans, int n, int nrhs, std::shared_ptr<Vector> dl, std::shared_ptr<Vector> d,
               std::shared_ptr<Vector> du, std::shared_ptr<Vector> du2, std::shared_ptr<IntVector> ipiv, SharedMatrix b,
               int ldb);
/// PSI_DHGEQZ, a wrapper to return C_DHGEQZ using objects
int PSI_DHGEQZ(int irrep, char job, char compq, char compz, int n, int ilo, int ihi, SharedMatrix h, int ldh,
               SharedMatrix t, int ldt, std::shared_ptr<Vector> alphar, std::shared_ptr<Vector> alphai,
               std::shared_ptr<Vector> beta, SharedMatrix q, int ldq, SharedMatrix z, int ldz,
               std::shared_ptr<Vector> work, int lwork);
/// PSI_DHSEIN, a wrapper to return C_DHSEIN using objects
int PSI_DHSEIN(int irrep, char side, char eigsrc, char initv, int n, SharedMatrix h, int ldh,
               std::shared_ptr<Vector> wr, std::shared_ptr<Vector> wi, SharedMatrix vl, int ldvl, SharedMatrix vr,
               int ldvr, int mm, std::shared_ptr<IntVector> m, std::shared_ptr<Vector> work,
               std::shared_ptr<IntVector> ifaill, std::shared_ptr<IntVector> ifailr);
/// PSI_DHSEQR, a wrapper to return C_DHSEQR using objects
int PSI_DHSEQR(int irrep, char job, char compz, int n, int ilo, int ihi, SharedMatrix h, int ldh,
               std::shared_ptr<Vector> wr, std::shared_ptr<Vector> wi, SharedMatrix z, int ldz,
               std::shared_ptr<Vector> work, int lwork);
/// PSI_DORGBR, a wrapper to return C_DORGBR using objects
int PSI_DORGBR(int irrep, char vect, int m, int n, int k, SharedMatrix a, int lda, std::shared_ptr<Vector> tau,
               std::shared_ptr<Vector> work, int lwork);
/// PSI_DORGHR, a wrapper to return C_DORGHR using objects
int PSI_DORGHR(int irrep, int n, int ilo, int ihi, SharedMatrix a, int lda, std::shared_ptr<Vector> tau,
               std::shared_ptr<Vector> work, int lwork);
/// PSI_DORGLQ, a wrapper to return C_DORGLQ using objects
int PSI_DORGLQ(int irrep, int m, int n, int k, SharedMatrix a, int lda, std::shared_ptr<Vector> tau,
               std::shared_ptr<Vector> work, int lwork);
/// PSI_DORGQL, a wrapper to return C_DORGQL using objects
int PSI_DORGQL(int irrep, int m, int n, int k, SharedMatrix a, int lda, std::shared_ptr<Vector> tau,
               std::shared_ptr<Vector> work, int lwork);
/// PSI_DORGQR, a wrapper to return C_DORGQR using objects
int PSI_DORGQR(int irrep, int m, int n, int k, SharedMatrix a, int lda, std::shared_ptr<Vector> tau,
               std::shared_ptr<Vector> work, int lwork);
/// PSI_DORGRQ, a wrapper to return C_DORGRQ using objects
int PSI_DORGRQ(int irrep, int m, int n, int k, SharedMatrix a, int lda, std::shared_ptr<Vector> tau,
               std::shared_ptr<Vector> work, int lwork);
/// PSI_DORGTR, a wrapper to return C_DORGTR using objects
int PSI_DORGTR(int irrep, char uplo, int n, SharedMatrix a, int lda, std::shared_ptr<Vector> tau,
               std::shared_ptr<Vector> work, int lwork);
/// PSI_DORMBR, a wrapper to return C_DORMBR using objects
int PSI_DORMBR(int irrep, char vect, char side, char trans, int m, int n, int k, std::shared_ptr<Vector> a, int lda,
               std::shared_ptr<Vector> tau, SharedMatrix c, int ldc, std::shared_ptr<Vector> work, int lwork);
/// PSI_DORMHR, a wrapper to return C_DORMHR using objects
int PSI_DORMHR(int irrep, char side, char trans, int m, int n, int ilo, int ihi, std::shared_ptr<Vector> a, int lda,
               std::shared_ptr<Vector> tau, SharedMatrix c, int ldc, std::shared_ptr<Vector> work, int lwork);
/// PSI_DORMLQ, a wrapper to return C_DORMLQ using objects
int PSI_DORMLQ(int irrep, char side, char trans, int m, int n, int k, std::shared_ptr<Vector> a, int lda,
               std::shared_ptr<Vector> tau, SharedMatrix c, int ldc, std::shared_ptr<Vector> work, int lwork);
/// PSI_DORMQL, a wrapper to return C_DORMQL using objects
int PSI_DORMQL(int irrep, char side, char trans, int m, int n, int k, SharedMatrix a, int lda,
               std::shared_ptr<Vector> tau, SharedMatrix c, int ldc, std::shared_ptr<Vector> work, int lwork);
/// PSI_DORMQR, a wrapper to return C_DORMQR using objects
int PSI_DORMQR(int irrep, char side, char trans, int m, int n, int k, SharedMatrix a, int lda,
               std::shared_ptr<Vector> tau, SharedMatrix c, int ldc, std::shared_ptr<Vector> work, int lwork);
/// PSI_DORMR3, a wrapper to return C_DORMR3 using objects
int PSI_DORMR3(int irrep, char side, char trans, int m, int n, int k, int l, std::shared_ptr<Vector> a, int lda,
               std::shared_ptr<Vector> tau, SharedMatrix c, int ldc, std::shared_ptr<Vector> work);
/// PSI_DORMRQ, a wrapper to return C_DORMRQ using objects
int PSI_DORMRQ(int irrep, char side, char trans, int m, int n, int k, std::shared_ptr<Vector> a, int lda,
               std::shared_ptr<Vector> tau, SharedMatrix c, int ldc, std::shared_ptr<Vector> work, int lwork);
/// PSI_DORMRZ, a wrapper to return C_DORMRZ using objects
int PSI_DORMRZ(int irrep, char side, char trans, int m, int n, int k, int l, std::shared_ptr<Vector> a, int lda,
               std::shared_ptr<Vector> tau, SharedMatrix c, int ldc, std::shared_ptr<Vector> work, int lwork);
/// PSI_DORMTR, a wrapper to return C_DORMTR using objects
int PSI_DORMTR(int irrep, char side, char uplo, char trans, int m, int n, std::shared_ptr<Vector> a, int lda,
               std::shared_ptr<Vector> tau, SharedMatrix c, int ldc, std::shared_ptr<Vector> work, int lwork);
/// PSI_DPBCON, a wrapper to return C_DPBCON using objects
int PSI_DPBCON(int irrep, char uplo, int n, int kd, SharedMatrix ab, int ldab, double anorm,
               std::shared_ptr<Vector> rcond, std::shared_ptr<Vector> work, std::shared_ptr<IntVector> iwork);
/// PSI_DPBEQU, a wrapper to return C_DPBEQU using objects
int PSI_DPBEQU(int irrep, char uplo, int n, int kd, SharedMatrix ab, int ldab, std::shared_ptr<Vector> s,
               std::shared_ptr<Vector> scond, std::shared_ptr<Vector> amax);
/// PSI_DPBRFS, a wrapper to return C_DPBRFS using objects
int PSI_DPBRFS(int irrep, char uplo, int n, int kd, int nrhs, SharedMatrix ab, int ldab, SharedMatrix afb, int ldafb,
               SharedMatrix b, int ldb, SharedMatrix x, int ldx, std::shared_ptr<Vector> ferr,
               std::shared_ptr<Vector> berr, std::shared_ptr<Vector> work, std::shared_ptr<IntVector> iwork);
/// PSI_DPBSTF, a wrapper to return C_DPBSTF using objects
int PSI_DPBSTF(int irrep, char uplo, int n, int kd, SharedMatrix ab, int ldab);
/// PSI_DPBSV, a wrapper to return C_DPBSV using objects
int PSI_DPBSV(int irrep, char uplo, int n, int kd, int nrhs, SharedMatrix ab, int ldab, SharedMatrix b, int ldb);
/// PSI_DPBSVX, a wrapper to return C_DPBSVX using objects
int PSI_DPBSVX(int irrep, char fact, char uplo, int n, int kd, int nrhs, SharedMatrix ab, int ldab, SharedMatrix afb,
               int ldafb, char equed, std::shared_ptr<Vector> s, SharedMatrix b, int ldb, SharedMatrix x, int ldx,
               std::shared_ptr<Vector> rcond, std::shared_ptr<Vector> ferr, std::shared_ptr<Vector> berr,
               std::shared_ptr<Vector> work, std::shared_ptr<IntVector> iwork);
/// PSI_DPBTRF, a wrapper to return C_DPBTRF using objects
int PSI_DPBTRF(int irrep, char uplo, int n, int kd, SharedMatrix ab, int ldab);
/// PSI_DPBTRS, a wrapper to return C_DPBTRS using objects
int PSI_DPBTRS(int irrep, char uplo, int n, int kd, int nrhs, SharedMatrix ab, int ldab, SharedMatrix b, int ldb);
/// PSI_DPOCON, a wrapper to return C_DPOCON using objects
int PSI_DPOCON(int irrep, char uplo, int n, SharedMatrix a, int lda, double anorm, std::shared_ptr<Vector> rcond,
               std::shared_ptr<Vector> work, std::shared_ptr<IntVector> iwork);
/// PSI_DPOEQU, a wrapper to return C_DPOEQU using objects
int PSI_DPOEQU(int irrep, int n, SharedMatrix a, int lda, std::shared_ptr<Vector> s, std::shared_ptr<Vector> scond,
               std::shared_ptr<Vector> amax);
/// PSI_DPORFS, a wrapper to return C_DPORFS using objects
int PSI_DPORFS(int irrep, char uplo, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix af, int ldaf,
               SharedMatrix b, int ldb, SharedMatrix x, int ldx, std::shared_ptr<Vector> ferr,
               std::shared_ptr<Vector> berr, std::shared_ptr<Vector> work, std::shared_ptr<IntVector> iwork);
/// PSI_DPOSV, a wrapper to return C_DPOSV using objects
int PSI_DPOSV(int irrep, char uplo, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix b, int ldb);
/// PSI_DPOSVX, a wrapper to return C_DPOSVX using objects
int PSI_DPOSVX(int irrep, char fact, char uplo, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix af, int ldaf,
               char equed, std::shared_ptr<Vector> s, SharedMatrix b, int ldb, SharedMatrix x, int ldx,
               std::shared_ptr<Vector> rcond, std::shared_ptr<Vector> ferr, std::shared_ptr<Vector> berr,
               std::shared_ptr<Vector> work, std::shared_ptr<IntVector> iwork);
/// PSI_DPOTRF, a wrapper to return C_DPOTRF using objects
int PSI_DPOTRF(int irrep, char uplo, int n, SharedMatrix a, int lda);
/// PSI_DPOTRI, a wrapper to return C_DPOTRI using objects
int PSI_DPOTRI(int irrep, char uplo, int n, SharedMatrix a, int lda);
/// PSI_DPOTRS, a wrapper to return C_DPOTRS using objects
int PSI_DPOTRS(int irrep, char uplo, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix b, int ldb);
/// PSI_DPTCON, a wrapper to return C_DPTCON using objects
int PSI_DPTCON(int irrep, int n, std::shared_ptr<Vector> d, std::shared_ptr<Vector> e, double anorm,
               std::shared_ptr<Vector> rcond, std::shared_ptr<Vector> work);
/// PSI_DPTEQR, a wrapper to return C_DPTEQR using objects
int PSI_DPTEQR(int irrep, char compz, int n, std::shared_ptr<Vector> d, std::shared_ptr<Vector> e, SharedMatrix z,
               int ldz, std::shared_ptr<Vector> work);
/// PSI_DPTRFS, a wrapper to return C_DPTRFS using objects
int PSI_DPTRFS(int irrep, int n, int nrhs, std::shared_ptr<Vector> d, std::shared_ptr<Vector> e,
               std::shared_ptr<Vector> df, std::shared_ptr<Vector> ef, SharedMatrix b, int ldb, SharedMatrix x, int ldx,
               std::shared_ptr<Vector> ferr, std::shared_ptr<Vector> berr, std::shared_ptr<Vector> work);
/// PSI_DPTSV, a wrapper to return C_DPTSV using objects
int PSI_DPTSV(int irrep, int n, int nrhs, std::shared_ptr<Vector> d, std::shared_ptr<Vector> e, SharedMatrix b,
              int ldb);
/// PSI_DPTSVX, a wrapper to return C_DPTSVX using objects
int PSI_DPTSVX(int irrep, char fact, int n, int nrhs, std::shared_ptr<Vector> d, std::shared_ptr<Vector> e,
               std::shared_ptr<Vector> df, std::shared_ptr<Vector> ef, SharedMatrix b, int ldb, SharedMatrix x, int ldx,
               std::shared_ptr<Vector> rcond, std::shared_ptr<Vector> ferr, std::shared_ptr<Vector> berr,
               std::shared_ptr<Vector> work);
/// PSI_DPTTRF, a wrapper to return C_DPTTRF using objects
int PSI_DPTTRF(int irrep, int n, std::shared_ptr<Vector> d, std::shared_ptr<Vector> e);
/// PSI_DPTTRS, a wrapper to return C_DPTTRS using objects
int PSI_DPTTRS(int irrep, int n, int nrhs, std::shared_ptr<Vector> d, std::shared_ptr<Vector> e, SharedMatrix b,
               int ldb);
/// PSI_DSBEV, a wrapper to return C_DSBEV using objects
int PSI_DSBEV(int irrep, char jobz, char uplo, int n, int kd, SharedMatrix ab, int ldab, std::shared_ptr<Vector> w,
              SharedMatrix z, int ldz, std::shared_ptr<Vector> work);
/// PSI_DSBEVD, a wrapper to return C_DSBEVD using objects
int PSI_DSBEVD(int irrep, char jobz, char uplo, int n, int kd, SharedMatrix ab, int ldab, std::shared_ptr<Vector> w,
               SharedMatrix z, int ldz, std::shared_ptr<Vector> work, int lwork, std::shared_ptr<IntVector> iwork,
               int liwork);
/// PSI_DSBEVX, a wrapper to return C_DSBEVX using objects
int PSI_DSBEVX(int irrep, char jobz, char range, char uplo, int n, int kd, SharedMatrix ab, int ldab, SharedMatrix q,
               int ldq, double vl, double vu, int il, int iu, double abstol, std::shared_ptr<IntVector> m,
               std::shared_ptr<Vector> w, SharedMatrix z, int ldz, std::shared_ptr<Vector> work,
               std::shared_ptr<IntVector> iwork, std::shared_ptr<IntVector> ifail);
/// PSI_DSBGST, a wrapper to return C_DSBGST using objects
int PSI_DSBGST(int irrep, char vect, char uplo, int n, int ka, int kb, SharedMatrix ab, int ldab, SharedMatrix bb,
               int ldbb, SharedMatrix x, int ldx, std::shared_ptr<Vector> work);
/// PSI_DSBGV, a wrapper to return C_DSBGV using objects
int PSI_DSBGV(int irrep, char jobz, char uplo, int n, int ka, int kb, SharedMatrix ab, int ldab, SharedMatrix bb,
              int ldbb, std::shared_ptr<Vector> w, SharedMatrix z, int ldz, std::shared_ptr<Vector> work);
/// PSI_DSBGVD, a wrapper to return C_DSBGVD using objects
int PSI_DSBGVD(int irrep, char jobz, char uplo, int n, int ka, int kb, SharedMatrix ab, int ldab, SharedMatrix bb,
               int ldbb, std::shared_ptr<Vector> w, SharedMatrix z, int ldz, std::shared_ptr<Vector> work, int lwork,
               std::shared_ptr<IntVector> iwork, int liwork);
/// PSI_DSBGVX, a wrapper to return C_DSBGVX using objects
int PSI_DSBGVX(int irrep, char jobz, char range, char uplo, int n, int ka, int kb, SharedMatrix ab, int ldab,
               SharedMatrix bb, int ldbb, SharedMatrix q, int ldq, double vl, double vu, int il, int iu, double abstol,
               std::shared_ptr<IntVector> m, std::shared_ptr<Vector> w, SharedMatrix z, int ldz,
               std::shared_ptr<Vector> work, std::shared_ptr<IntVector> iwork, std::shared_ptr<IntVector> ifail);
/// PSI_DSBTRD, a wrapper to return C_DSBTRD using objects
int PSI_DSBTRD(int irrep, char vect, char uplo, int n, int kd, SharedMatrix ab, int ldab, std::shared_ptr<Vector> d,
               std::shared_ptr<Vector> e, SharedMatrix q, int ldq, std::shared_ptr<Vector> work);
/// PSI_DSGESV, a wrapper to return C_DSGESV using objects
int PSI_DSGESV(int irrep, int n, int nrhs, std::shared_ptr<Vector> a, int lda, std::shared_ptr<IntVector> ipiv,
               SharedMatrix b, int ldb, SharedMatrix x, int ldx, SharedMatrix work, std::shared_ptr<IntVector> iter);
/// PSI_DSTEBZ, a wrapper to return C_DSTEBZ using objects
int PSI_DSTEBZ(int irrep, char range, char order, int n, double vl, double vu, int il, int iu, double abstol,
               std::shared_ptr<Vector> d, std::shared_ptr<Vector> e, std::shared_ptr<IntVector> m,
               std::shared_ptr<IntVector> nsplit, std::shared_ptr<Vector> w, std::shared_ptr<IntVector> iblock,
               std::shared_ptr<IntVector> isplit, std::shared_ptr<Vector> work, std::shared_ptr<IntVector> iwork);
/// PSI_DSTEDC, a wrapper to return C_DSTEDC using objects
int PSI_DSTEDC(int irrep, char compz, int n, std::shared_ptr<Vector> d, std::shared_ptr<Vector> e, SharedMatrix z,
               int ldz, std::shared_ptr<Vector> work, int lwork, std::shared_ptr<IntVector> iwork, int liwork);
/// PSI_DSTEGR, a wrapper to return C_DSTEGR using objects
int PSI_DSTEGR(int irrep, char jobz, char range, int n, std::shared_ptr<Vector> d, std::shared_ptr<Vector> e, double vl,
               double vu, int il, int iu, double abstol, std::shared_ptr<IntVector> m, std::shared_ptr<Vector> w,
               SharedMatrix z, int ldz, std::shared_ptr<IntVector> isuppz, std::shared_ptr<Vector> work, int lwork,
               std::shared_ptr<IntVector> iwork, int liwork);
/// PSI_DSTEIN, a wrapper to return C_DSTEIN using objects
int PSI_DSTEIN(int irrep, int n, std::shared_ptr<Vector> d, std::shared_ptr<Vector> e, int m, std::shared_ptr<Vector> w,
               std::shared_ptr<IntVector> iblock, std::shared_ptr<IntVector> isplit, SharedMatrix z, int ldz,
               std::shared_ptr<Vector> work, std::shared_ptr<IntVector> iwork, std::shared_ptr<IntVector> ifail);
/// PSI_DSTEQR, a wrapper to return C_DSTEQR using objects
int PSI_DSTEQR(int irrep, char compz, int n, std::shared_ptr<Vector> d, std::shared_ptr<Vector> e, SharedMatrix z,
               int ldz, std::shared_ptr<Vector> work);
/// PSI_DSTERF, a wrapper to return C_DSTERF using objects
int PSI_DSTERF(int irrep, int n, std::shared_ptr<Vector> d, std::shared_ptr<Vector> e);
/// PSI_DSTEV, a wrapper to return C_DSTEV using objects
int PSI_DSTEV(int irrep, char jobz, int n, std::shared_ptr<Vector> d, std::shared_ptr<Vector> e, SharedMatrix z,
              int ldz, std::shared_ptr<Vector> work);
/// PSI_DSTEVD, a wrapper to return C_DSTEVD using objects
int PSI_DSTEVD(int irrep, char jobz, int n, std::shared_ptr<Vector> d, std::shared_ptr<Vector> e, SharedMatrix z,
               int ldz, std::shared_ptr<Vector> work, int lwork, std::shared_ptr<IntVector> iwork, int liwork);
/// PSI_DSTEVR, a wrapper to return C_DSTEVR using objects
int PSI_DSTEVR(int irrep, char jobz, char range, int n, std::shared_ptr<Vector> d, std::shared_ptr<Vector> e, double vl,
               double vu, int il, int iu, double abstol, std::shared_ptr<IntVector> m, std::shared_ptr<Vector> w,
               SharedMatrix z, int ldz, std::shared_ptr<IntVector> isuppz, std::shared_ptr<Vector> work, int lwork,
               std::shared_ptr<IntVector> iwork, int liwork);
/// PSI_DSTEVX, a wrapper to return C_DSTEVX using objects
int PSI_DSTEVX(int irrep, char jobz, char range, int n, std::shared_ptr<Vector> d, std::shared_ptr<Vector> e, double vl,
               double vu, int il, int iu, double abstol, std::shared_ptr<IntVector> m, std::shared_ptr<Vector> w,
               SharedMatrix z, int ldz, std::shared_ptr<Vector> work, std::shared_ptr<IntVector> iwork,
               std::shared_ptr<IntVector> ifail);
/// PSI_DSYCON, a wrapper to return C_DSYCON using objects
int PSI_DSYCON(int irrep, char uplo, int n, SharedMatrix a, int lda, std::shared_ptr<IntVector> ipiv, double anorm,
               std::shared_ptr<Vector> rcond, std::shared_ptr<Vector> work, std::shared_ptr<IntVector> iwork);
/// PSI_DSYEV, a wrapper to return C_DSYEV using objects
int PSI_DSYEV(int irrep, char jobz, char uplo, int n, SharedMatrix a, int lda, std::shared_ptr<Vector> w,
              std::shared_ptr<Vector> work, int lwork);
/// PSI_DSYEVD, a wrapper to return C_DSYEVD using objects
int PSI_DSYEVD(int irrep, char jobz, char uplo, int n, SharedMatrix a, int lda, std::shared_ptr<Vector> w,
               std::shared_ptr<Vector> work, int lwork, std::shared_ptr<IntVector> iwork, int liwork);
/// PSI_DSYEVR, a wrapper to return C_DSYEVR using objects
int PSI_DSYEVR(int irrep, char jobz, char range, char uplo, int n, SharedMatrix a, int lda, double vl, double vu,
               int il, int iu, double abstol, std::shared_ptr<IntVector> m, std::shared_ptr<Vector> w, SharedMatrix z,
               int ldz, std::shared_ptr<IntVector> isuppz, std::shared_ptr<Vector> work, int lwork,
               std::shared_ptr<IntVector> iwork, int liwork);
/// PSI_DSYEVX, a wrapper to return C_DSYEVX using objects
int PSI_DSYEVX(int irrep, char jobz, char range, char uplo, int n, SharedMatrix a, int lda, double vl, double vu,
               int il, int iu, double abstol, std::shared_ptr<IntVector> m, std::shared_ptr<Vector> w, SharedMatrix z,
               int ldz, std::shared_ptr<Vector> work, int lwork, std::shared_ptr<IntVector> iwork,
               std::shared_ptr<IntVector> ifail);
/// PSI_DSYGST, a wrapper to return C_DSYGST using objects
int PSI_DSYGST(int irrep, int itype, char uplo, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb);
/// PSI_DSYGV, a wrapper to return C_DSYGV using objects
int PSI_DSYGV(int irrep, int itype, char jobz, char uplo, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb,
              std::shared_ptr<Vector> w, std::shared_ptr<Vector> work, int lwork);
/// PSI_DSYGVD, a wrapper to return C_DSYGVD using objects
int PSI_DSYGVD(int irrep, int itype, char jobz, char uplo, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb,
               std::shared_ptr<Vector> w, std::shared_ptr<Vector> work, int lwork, std::shared_ptr<IntVector> iwork,
               int liwork);
/// PSI_DSYGVX, a wrapper to return C_DSYGVX using objects
int PSI_DSYGVX(int irrep, int itype, char jobz, char range, char uplo, int n, SharedMatrix a, int lda, SharedMatrix b,
               int ldb, double vl, double vu, int il, int iu, double abstol, std::shared_ptr<IntVector> m,
               std::shared_ptr<Vector> w, SharedMatrix z, int ldz, std::shared_ptr<Vector> work, int lwork,
               std::shared_ptr<IntVector> iwork, std::shared_ptr<IntVector> ifail);
/// PSI_DSYRFS, a wrapper to return C_DSYRFS using objects
int PSI_DSYRFS(int irrep, char uplo, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix af, int ldaf,
               std::shared_ptr<IntVector> ipiv, SharedMatrix b, int ldb, SharedMatrix x, int ldx,
               std::shared_ptr<Vector> ferr, std::shared_ptr<Vector> berr, std::shared_ptr<Vector> work,
               std::shared_ptr<IntVector> iwork);
/// PSI_DSYSV, a wrapper to return C_DSYSV using objects
int PSI_DSYSV(int irrep, char uplo, int n, int nrhs, SharedMatrix a, int lda, std::shared_ptr<IntVector> ipiv,
              SharedMatrix b, int ldb, std::shared_ptr<Vector> work, int lwork);
/// PSI_DSYSVX, a wrapper to return C_DSYSVX using objects
int PSI_DSYSVX(int irrep, char fact, char uplo, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix af, int ldaf,
               std::shared_ptr<IntVector> ipiv, SharedMatrix b, int ldb, SharedMatrix x, int ldx,
               std::shared_ptr<Vector> rcond);
/// PSI_DSYTRD, a wrapper to return C_DSYTRD using objects
int PSI_DSYTRD(int irrep, char uplo, int n, SharedMatrix a, int lda, std::shared_ptr<Vector> d,
               std::shared_ptr<Vector> e, std::shared_ptr<Vector> tau, std::shared_ptr<Vector> work, int lwork);
/// PSI_DSYTRF, a wrapper to return C_DSYTRF using objects
int PSI_DSYTRF(int irrep, char uplo, int n, SharedMatrix a, int lda, std::shared_ptr<IntVector> ipiv,
               std::shared_ptr<Vector> work, int lwork);
/// PSI_DSYTRI, a wrapper to return C_DSYTRI using objects
int PSI_DSYTRI(int irrep, char uplo, int n, SharedMatrix a, int lda, std::shared_ptr<IntVector> ipiv,
               std::shared_ptr<Vector> work);
/// PSI_DSYTRS, a wrapper to return C_DSYTRS using objects
int PSI_DSYTRS(int irrep, char uplo, int n, int nrhs, SharedMatrix a, int lda, std::shared_ptr<IntVector> ipiv,
               SharedMatrix b, int ldb);
/// PSI_DTBCON, a wrapper to return C_DTBCON using objects
int PSI_DTBCON(int irrep, char norm, char uplo, char diag, int n, int kd, SharedMatrix ab, int ldab,
               std::shared_ptr<Vector> rcond, std::shared_ptr<Vector> work, std::shared_ptr<IntVector> iwork);
/// PSI_DTBRFS, a wrapper to return C_DTBRFS using objects
int PSI_DTBRFS(int irrep, char uplo, char trans, char diag, int n, int kd, int nrhs, SharedMatrix ab, int ldab,
               SharedMatrix b, int ldb, SharedMatrix x, int ldx, std::shared_ptr<Vector> ferr,
               std::shared_ptr<Vector> berr, std::shared_ptr<Vector> work, std::shared_ptr<IntVector> iwork);
/// PSI_DTBTRS, a wrapper to return C_DTBTRS using objects
int PSI_DTBTRS(int irrep, char uplo, char trans, char diag, int n, int kd, int nrhs, SharedMatrix ab, int ldab,
               SharedMatrix b, int ldb);
/// PSI_DTGEVC, a wrapper to return C_DTGEVC using objects
int PSI_DTGEVC(int irrep, char side, char howmny, int n, SharedMatrix s, int lds, SharedMatrix p, int ldp,
               SharedMatrix vl, int ldvl, SharedMatrix vr, int ldvr, int mm, std::shared_ptr<IntVector> m,
               std::shared_ptr<Vector> work);
/// PSI_DTGEXC, a wrapper to return C_DTGEXC using objects
int PSI_DTGEXC(int irrep, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb, SharedMatrix q, int ldq,
               SharedMatrix z, int ldz, std::shared_ptr<IntVector> ifst, std::shared_ptr<IntVector> ilst,
               std::shared_ptr<Vector> work, int lwork);
/// PSI_DTGSEN, a wrapper to return C_DTGSEN using objects
int PSI_DTGSEN(int irrep, int ijob, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb,
               std::shared_ptr<Vector> alphar, std::shared_ptr<Vector> alphai, std::shared_ptr<Vector> beta,
               SharedMatrix q, int ldq, SharedMatrix z, int ldz, std::shared_ptr<IntVector> m,
               std::shared_ptr<Vector> pl, std::shared_ptr<Vector> pr, std::shared_ptr<Vector> dif,
               std::shared_ptr<Vector> work, int lwork, std::shared_ptr<IntVector> iwork, int liwork);
/// PSI_DTGSJA, a wrapper to return C_DTGSJA using objects
int PSI_DTGSJA(int irrep, char jobu, char jobv, char jobq, int m, int p, int n, int k, int l, SharedMatrix a, int lda,
               SharedMatrix b, int ldb, double tola, double tolb, std::shared_ptr<Vector> alpha,
               std::shared_ptr<Vector> beta, SharedMatrix u, int ldu, SharedMatrix v, int ldv, SharedMatrix q, int ldq,
               std::shared_ptr<Vector> work, std::shared_ptr<IntVector> ncycle);
/// PSI_DTGSNA, a wrapper to return C_DTGSNA using objects
int PSI_DTGSNA(int irrep, char job, char howmny, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb,
               SharedMatrix vl, int ldvl, SharedMatrix vr, int ldvr, std::shared_ptr<Vector> s,
               std::shared_ptr<Vector> dif, int mm, std::shared_ptr<IntVector> m, std::shared_ptr<Vector> work,
               int lwork, std::shared_ptr<IntVector> iwork);
/// PSI_DTGSYL, a wrapper to return C_DTGSYL using objects
int PSI_DTGSYL(int irrep, char trans, int ijob, int m, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb,
               SharedMatrix c, int ldc, SharedMatrix d, int ldd, SharedMatrix e, int lde, SharedMatrix f, int ldf,
               std::shared_ptr<Vector> dif, std::shared_ptr<Vector> scale, std::shared_ptr<Vector> work, int lwork,
               std::shared_ptr<IntVector> iwork);
/// PSI_DTRCON, a wrapper to return C_DTRCON using objects
int PSI_DTRCON(int irrep, char norm, char uplo, char diag, int n, SharedMatrix a, int lda,
               std::shared_ptr<Vector> rcond, std::shared_ptr<Vector> work, std::shared_ptr<IntVector> iwork);
/// PSI_DTREVC, a wrapper to return C_DTREVC using objects
int PSI_DTREVC(int irrep, char side, char howmny, int n, SharedMatrix t, int ldt, SharedMatrix vl, int ldvl,
               SharedMatrix vr, int ldvr, int mm, std::shared_ptr<IntVector> m, std::shared_ptr<Vector> work);
/// PSI_DTREXC, a wrapper to return C_DTREXC using objects
int PSI_DTREXC(int irrep, char compq, int n, SharedMatrix t, int ldt, SharedMatrix q, int ldq,
               std::shared_ptr<IntVector> ifst, std::shared_ptr<IntVector> ilst, std::shared_ptr<Vector> work);
/// PSI_DTRRFS, a wrapper to return C_DTRRFS using objects
int PSI_DTRRFS(int irrep, char uplo, char trans, char diag, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix b,
               int ldb, SharedMatrix x, int ldx, std::shared_ptr<Vector> ferr, std::shared_ptr<Vector> berr,
               std::shared_ptr<Vector> work, std::shared_ptr<IntVector> iwork);
/// PSI_DTRSEN, a wrapper to return C_DTRSEN using objects
int PSI_DTRSEN(int irrep, char job, char compq, int n, SharedMatrix t, int ldt, SharedMatrix q, int ldq,
               std::shared_ptr<Vector> wr, std::shared_ptr<Vector> wi, std::shared_ptr<IntVector> m,
               std::shared_ptr<Vector> s, std::shared_ptr<Vector> sep, std::shared_ptr<Vector> work, int lwork,
               std::shared_ptr<IntVector> iwork, int liwork);
/// PSI_DTRSNA, a wrapper to return C_DTRSNA using objects
int PSI_DTRSNA(int irrep, char job, char howmny, int n, SharedMatrix t, int ldt, SharedMatrix vl, int ldvl,
               SharedMatrix vr, int ldvr, std::shared_ptr<Vector> s, std::shared_ptr<Vector> sep, int mm,
               std::shared_ptr<IntVector> m, SharedMatrix work, int ldwork, std::shared_ptr<IntVector> iwork);
/// PSI_DTRSYL, a wrapper to return C_DTRSYL using objects
int PSI_DTRSYL(int irrep, char trana, char tranb, int isgn, int m, int n, SharedMatrix a, int lda, SharedMatrix b,
               int ldb, SharedMatrix c, int ldc, std::shared_ptr<Vector> scale);
/// PSI_DTRTRI, a wrapper to return C_DTRTRI using objects
int PSI_DTRTRI(int irrep, char uplo, char diag, int n, SharedMatrix a, int lda);
/// PSI_DTRTRS, a wrapper to return C_DTRTRS using objects
int PSI_DTRTRS(int irrep, char uplo, char trans, char diag, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix b,
               int ldb);
/// PSI_DTZRQF, a wrapper to return C_DTZRQF using objects
int PSI_DTZRQF(int irrep, int m, int n, SharedMatrix a, int lda, std::shared_ptr<Vector> tau);
/// PSI_DTZRZF, a wrapper to return C_DTZRZF using objects
int PSI_DTZRZF(int irrep, int m, int n, SharedMatrix a, int lda, std::shared_ptr<Vector> tau,
               std::shared_ptr<Vector> work, int lwork);
}  // namespace psi
#endif
