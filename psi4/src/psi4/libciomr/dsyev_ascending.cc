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

/*!
** \file
** \brief Diagnoalize a symmetrix square matrix
** \ingroup CIOMR
*/

#include "psi4/libqt/qt.h"
#include "libciomr.h"
#include <vector>

#ifdef USING_cuEST
#include <cublas_v2.h>
#include <cusolverDn.h>
extern cusolverDnHandle_t cusolver_handle;
#endif

// Helper macro for brevity (you may want a nicer error handler)
#define CHECK_CUDA(call) \
    do { \
        cudaError_t err__ = (call); \
        if (err__ != cudaSuccess) return static_cast<int>(err__); \
    } while (0)

#define CHECK_CUSOLVER(call) \
    do { \
        cusolverStatus_t st__ = (call); \
        if (st__ != CUSOLVER_STATUS_SUCCESS) return static_cast<int>(st__); \
    } while (0)

namespace psi {

namespace {

[[nodiscard]] int DSYEV_ascending_lapack(const int N, const double* const* const array, double* e_vals,
                                         double* const* const e_vecs) {
    // We need to make a copy of the matrix before diagonalization, because LAPACK overwrites it.
    // LAPACK also needs the mtx to be flattened to a 1D array, so a copy is inevitable.
    // The new 1D array will correspond to a column-major array, suitable for LAPACK.
    std::vector<double> tmp_matrix(N * N);
    for (int64_t i = 0, ij = 0; i < N; i++) {
        for (int64_t j = 0; j < N; j++, ij++) {
            tmp_matrix[ij] = array[j][i];
        }
    }
    // LAPACK also needs some extra memory to store temporaries in
    // TODO: query C_DSYEV for optimal workspace size
    const int64_t workspace_size = 3 * N;
    std::vector<double> tmp_work(workspace_size);
    const char jobtype = (e_vecs != nullptr) ? 'V' : 'N';
    const int info = C_DSYEV(jobtype, 'U', N, tmp_matrix.data(), N, e_vals, tmp_work.data(), workspace_size);
    if ((info == 0) && (e_vecs != nullptr)) {
        // tmp_matrix has now been overwritten with the eigenvecs as the columns, flattened as column-major
        // Copy them to the columns of a row-major 2D array
        for (int64_t j = 0, ij = 0; j < N; j++) {
            for (int64_t i = 0; i < N; i++, ij++) {
                e_vecs[i][j] = tmp_matrix[ij];
            }
        }
    }
    return info;
}

}  // namespace

#ifdef USING_cuEST
[[nodiscard]] int DSYEV_ascending(
    const int N,
    const double* const* const array,   // row-major 2D, CPU
    double* e_vals,                     // length N, CPU
    double* const* const e_vecs         // row-major 2D, CPU (optional)
) {
    if (cusolver_handle == nullptr) {
        return DSYEV_ascending_lapack(N, array, e_vals, e_vecs);
    }

    // 1. Flatten input matrix from row-major array[i][j] to column-major tmp_A(j,i)
    std::vector<double> h_A(N * N);
    for (int i = 0, ij = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j, ++ij) {
            h_A[ij] = array[j][i];  // column-major: A(j,i)
        }
    }

    // 3. Allocate device memory
    double* d_A = nullptr;  // N x N, column-major
    double* d_W = nullptr;  // eigenvalues
    int*    d_info = nullptr;

    CHECK_CUDA(cudaMalloc((void**)&d_A, sizeof(double) * N * N));
    CHECK_CUDA(cudaMalloc((void**)&d_W, sizeof(double) * N));
    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));

    // 4. Copy matrix to device
    CHECK_CUDA(cudaMemcpy(d_A, h_A.data(),
                          sizeof(double) * N * N,
                          cudaMemcpyHostToDevice));

    // 5. Query workspace size for Dsyevd
    int lwork = 0;
    CHECK_CUSOLVER(cusolverDnDsyevd_bufferSize(
        cusolver_handle,
        CUSOLVER_EIG_MODE_VECTOR,      // we want eigenvectors always; we can discard later
        CUBLAS_FILL_MODE_UPPER,
        N,
        d_A,
        N,
        d_W,
        &lwork));

    double* d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&d_work, sizeof(double) * lwork));

    // 6. Compute eigen-decomposition on device
    CHECK_CUSOLVER(cusolverDnDsyevd(
        cusolver_handle,
        (e_vecs ? CUSOLVER_EIG_MODE_VECTOR : CUSOLVER_EIG_MODE_NOVECTOR),
        CUBLAS_FILL_MODE_UPPER,
        N,
        d_A,
        N,
        d_W,
        d_work,
        lwork,
        d_info));

    // 7. Copy info and eigenvalues back
    int info = 0;
    CHECK_CUDA(cudaMemcpy(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (info != 0) {
        // Clean up and return LAPACK-style error code
        cudaFree(d_A);
        cudaFree(d_W);
        cudaFree(d_work);
        cudaFree(d_info);
        return info;
    }

    CHECK_CUDA(cudaMemcpy(e_vals, d_W,
                          sizeof(double) * N,
                          cudaMemcpyDeviceToHost));

    // 8. If eigenvectors requested, copy them back and reorder to row-major
    if (e_vecs) {
        std::vector<double> h_V(N * N);
        CHECK_CUDA(cudaMemcpy(h_V.data(), d_A,
                              sizeof(double) * N * N,
                              cudaMemcpyDeviceToHost));
        // h_V is column-major: columns = eigenvectors
        // Copy to row-major e_vecs[i][j]
        for (int j = 0, ij = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i, ++ij) {
                e_vecs[i][j] = h_V[ij];
            }
        }
    }

    // 9. Clean up
    cudaFree(d_A);
    cudaFree(d_W);
    cudaFree(d_work);
    cudaFree(d_info);

    return info;
}
#else
/*!
** DSYEV_ascending(): diagonalize a symmetric square matrix ('array') using LAPACK DSYEV
**
** \param n      = number of rows (and columns)
** \param array  = matrix to diagonalize (2D row major array)
** \param e_vals = array to hold eigenvalues (returned in ascending order)
** \param e_vecs = (optional) matrix of eigenvectors (2D row major array, one column for each eigvector)
**
** \ingroup CIOMR
*/
[[nodiscard]] int DSYEV_ascending(const int N, const double* const* const array, double* e_vals,
                                  double* const* const e_vecs /* = nullptr*/) {
    return DSYEV_ascending_lapack(N, array, e_vals, e_vecs);
}
#endif

}  // namespace psi
