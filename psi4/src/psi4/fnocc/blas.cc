/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#include "blas.h"
#include <cstdlib>

#include "psi4/pragma.h"

namespace psi {
namespace fnocc {

// position in a symmetric packed matrix
long int Position(long int i, long int j) {
    if (i < j) {
        return ((j * (j + 1)) >> 1) + i;
    }
    return ((i * (i + 1)) >> 1) + j;
}

/**
 * fortran-ordered dgemv
 */
void PSI_API F_DGEMV(char trans, integer m, integer n, doublereal alpha, doublereal* A, integer lda, doublereal* X,
                     integer incx, doublereal beta, doublereal* Y, integer incy) {
    DGEMV(trans, m, n, alpha, A, lda, X, incx, beta, Y, incy);
}
/**
 * fortran-ordered dgemm
 */
void PSI_API F_DGEMM(char transa, char transb, integer m, integer n, integer k, doublereal alpha, doublereal* A,
                     integer lda, doublereal* B, integer ldb, doublereal beta, doublereal* C, integer ldc) {
    DGEMM(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

/**
 *  Diagonalize a real symmetric matrix
 */
void Diagonalize(integer N, doublereal* A, doublereal* W) {
    char JOBZ = 'V';
    char UPLO = 'U';
    integer LDA = N;
    integer LWORK = 3 * N - 1;
    doublereal* WORK = (doublereal*)malloc(LWORK * sizeof(doublereal));
    integer INFO = 0;
    DSYEV(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO);
    free(WORK);
}
void Diagonalize2(integer N, doublereal* AP, doublereal* W, doublereal* Z) {
    char JOBZ = 'V';
    char UPLO = 'U';
    integer LDZ = N;
    doublereal* WORK = (doublereal*)malloc(3 * N * sizeof(doublereal));
    integer INFO = 0;
    DSPEV(JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO);
    free(WORK);
}

/**
 *  General SVD
 */
void SVD(integer M, integer N, doublereal* A, doublereal* U, doublereal* VT, doublereal* S) {
    char JOBU = 'S';   // all M columns of U are returned in array U
    char JOBVT = 'A';  // all N rows of V**T are returned in the array VT
    integer LDA = M;
    integer LDU = M;
    integer LDVT = N;

    integer min = (M < N) ? M : N;
    integer max = (M < N) ? N : M;
    integer LWORK = (3 * min + max > 5 * min) ? (3 * min + max) : (5 * min);
    doublereal* WORK = (doublereal*)malloc(LWORK * sizeof(doublereal));

    integer INFO = 0;
    DGESVD(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO);
    free(WORK);
}
}
}
