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

#include "psi4/pragma.h"

namespace psi {
PSI_API
void C_DGBMV(char trans, int m, int n, int kl, int ku, double alpha, double* a, int lda, double* x, int incx,
             double beta, double* y, int incy);
PSI_API
void C_DGEMV(char trans, int m, int n, double alpha, double* a, int lda, double* x, int incx, double beta, double* y,
             int incy);
PSI_API
void C_DGER(int m, int n, double alpha, double* x, int incx, double* y, int incy, double* a, int lda);
PSI_API
void C_DSBMV(char uplo, int n, int k, double alpha, double* a, int lda, double* x, int incx, double beta, double* y,
             int incy);
PSI_API
void C_DSPMV(char uplo, int n, double alpha, double* ap, double* x, int incx, double beta, double* y, int incy);
PSI_API
void C_DSPR(char uplo, int n, double alpha, double* x, int incx, double* ap);
PSI_API
void C_DSPR2(char uplo, int n, double alpha, double* x, int incx, double* y, int incy, double* ap);
PSI_API
void C_DSYMV(char uplo, int n, double alpha, double* a, int lda, double* x, int incx, double beta, double* y, int incy);
PSI_API
void C_DSYR(char uplo, int n, double alpha, double* x, int incx, double* a, int lda);
PSI_API
void C_DSYR2(char uplo, int n, double alpha, double* x, int incx, double* y, int incy, double* a, int lda);
PSI_API
void C_DTBMV(char uplo, char trans, char diag, int n, int k, double* a, int lda, double* x, int incx);
PSI_API
void C_DTBSV(char uplo, char trans, char diag, int n, int k, double* a, int lda, double* x, int incx);
PSI_API
void C_DTPMV(char uplo, char trans, char diag, int n, double* ap, double* x, int incx);
PSI_API
void C_DTPSV(char uplo, char trans, char diag, int n, double* ap, double* x, int incx);
PSI_API
void C_DTRMV(char uplo, char trans, char diag, int n, double* a, int lda, double* x, int incx);
PSI_API
void C_DTRSM(char side, char uplo, char transa, char diag, int m, int n, double alpha, double* a, int lda, double* b,
             int ldb);
}  // namespace psi
