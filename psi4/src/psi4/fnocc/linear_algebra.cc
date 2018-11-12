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

#include "psi4/libqt/lapack.h"

namespace psi {
namespace fnocc {
long int Position(long int i, long int j) {
    if (i < j) {
        return ((j * (j + 1)) >> 1) + i;
    }
    return ((i * (i + 1)) >> 1) + j;
}

void Diagonalize(long int N, double* A, double* W) {
    char JOBZ = 'V';
    char UPLO = 'U';
    long int LDA = N;
    auto info = C_DSYEV(JOBZ, UPLO, N, A, LDA, W);
}

void Diagonalize2(long int N, double* AP, double* W, double* Z) {
    char JOBZ = 'V';
    char UPLO = 'U';
    long int LDZ = N;
    auto info = C_DSPEV(JOBZ, UPLO, N, AP, W, Z, LDZ);
}

void SVD(long int M, long int N, double* A, double* U, double* VT, double* S) {
    char JOBU = 'S';   // all M columns of U are returned in array U
    char JOBVT = 'A';  // all N rows of V**T are returned in the array VT
    long int LDA = M;
    long int LDU = M;
    long int LDVT = N;

    auto info = C_DGESVD(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT);
}
}  // namespace fnocc
}  // namespace psi
