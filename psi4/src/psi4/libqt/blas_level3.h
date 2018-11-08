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
 * \brief Function signatures for wrappers to BLAS level 3 subroutines
 * \ingroup QT
 */

#pragma once

#include "psi4/pragma.h"

namespace psi {
PSI_API
void C_DGEMM(char transa, char transb, int m, int n, int k, double alpha, double* a, int lda, double* b, int ldb,
             double beta, double* c, int ldc);
PSI_API
void C_DSYMM(char side, char uplo, int m, int n, double alpha, double* a, int lda, double* b, int ldb, double beta,
             double* c, int ldc);
PSI_API
void C_DTRMM(char side, char uplo, char transa, char diag, int m, int n, double alpha, double* a, int lda, double* b,
             int ldb);
PSI_API
void C_DSYRK(char uplo, char trans, int n, int k, double alpha, double* a, int lda, double beta, double* c, int ldc);
PSI_API
void C_DSYR2K(char uplo, char trans, int n, int k, double alpha, double* a, int lda, double* b, int ldb, double beta,
              double* c, int ldc);
PSI_API
void C_DTRSV(char uplo, char trans, char diag, int n, double* a, int lda, double* x, int incx);
}  // namespace psi
