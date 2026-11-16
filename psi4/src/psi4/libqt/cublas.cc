/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#include "psi4/libqt/qt.h"

#ifdef USING_cuEST
#include <cublas_v2.h>
#include "psi4/libfock/cuESTCommon.h"

namespace psi {

// TODO: add device pointer mode to each function
void C_DGEMM_CUBLAS(
    cublasHandle_t handle,
    bool transa, bool transb,
    int m, int n, int k,
    double alpha, const double *A, int offa, int lda,
    const double *B, int offb, int ldb, double beta,
    double *C, int ldc
) {

    double alpha_arg = alpha;
    double beta_arg = beta;
    
    // swap argument order as cublas expects column-major buffers
    CHECK_CUBLAS(cublasDgemm(handle,
                (transb ? CUBLAS_OP_T : CUBLAS_OP_N), 
                (transa ? CUBLAS_OP_T : CUBLAS_OP_N),
                n, m, k,
                &alpha_arg,
                (B + offb), ldb,
                (A + offa), lda,
                &beta_arg,
                C, ldc));
}

void C_DGEAM_CUBLAS(
    cublasHandle_t handle,
    bool transa, bool transb,
    int m, int n,
    double alpha, const double *A, int lda,
    double beta, const double *B, int ldb,
    double *C, int ldc
) {

    double alpha_arg = alpha;
    double beta_arg = beta;

    CHECK_CUBLAS(cublasDgeam(handle,
                (transa ? CUBLAS_OP_T : CUBLAS_OP_N),
                (transb ? CUBLAS_OP_T : CUBLAS_OP_N), 
                m, n,
                &alpha_arg, A, lda,
                &beta_arg, B, ldb,
                C, ldc));
}

void C_DDOT_CUBLAS(
    cublasHandle_t handle, int n,
    const double *x, int offx, int incx,
    const double *y, int offy, int incy,
    double *result
) {

    CHECK_CUBLAS(cublasDdot(handle, n,
                (x + offx), incx,
                (y + offy), incy,
                result));
}

} // namespace psi
#endif
