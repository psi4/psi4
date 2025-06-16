/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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

#include <iostream>
#include <cmath>

#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"

#include "algebra_interface.h"
#include "blas.h"
#include "index.h"
#include "matrix.h"
#include "operation.h"

namespace psi {

namespace psimrcc {

void CCOperation::setup_contractions() {
    Timer PartA;

    CCMatTmp AMatTmp = wfn_->blas()->get_MatTmp(A_Matrix, none);
    check_and_zero_target();

    bool need_sort = false;
    if (reindexing.size() > 0) need_sort = true;
    CCIndex* T_left = A_Matrix->get_left();
    CCIndex* T_right = A_Matrix->get_right();
    auto T_matrix = A_Matrix->get_matrix();

    // Determine the target indexing if sorting is needed
    if (need_sort) {
        if (operation[0] == '1') {
            T_left = B_Matrix->get_right();
        } else {
            T_left = B_Matrix->get_left();
        }
        if (operation[2] == '1') {
            T_right = C_Matrix->get_right();
        } else {
            T_right = C_Matrix->get_left();
        }
        size_t T_matrix_offset = 0;
        T_matrix = std::vector<double**>(wfn_->nirrep(), nullptr);
        for (int irrep = 0; irrep < wfn_->nirrep(); irrep++) {
            T_matrix[irrep] = new double*[T_left->get_pairpi(irrep)];
            size_t block_size = T_left->get_pairpi(irrep) * T_right->get_pairpi(irrep);
            for (size_t i = 0; i < T_left->get_pairpi(irrep); i++) {
                T_matrix[irrep][i] = &(local_work[T_matrix_offset + i * T_right->get_pairpi(irrep)]);
            }
            T_matrix_offset += block_size;
        }
        if (T_matrix_offset > 0) zero_arr(&(local_work[0]), T_matrix_offset);
    }

    PartA_timing += PartA.get();
    Timer PartB;

    double** A_matrix;
    double** B_matrix;
    double** C_matrix;

    for (int h = 0; h < wfn_->moinfo()->get_nirreps(); h++) {
        bool B_on_disk = false;
        bool C_on_disk = false;
        if (!B_Matrix->is_block_allocated(h)) {
            if (B_Matrix->is_integral()) {
                B_on_disk = true;
            } else {
                B_Matrix->load_irrep(h);
            }
        }
        if (!C_Matrix->is_block_allocated(h)) {
            if (C_Matrix->is_integral()) {
                C_on_disk = true;
            } else {
                C_Matrix->load_irrep(h);
            }
        }
        if (B_on_disk && C_on_disk) throw PSIEXCEPTION("BOTH ON DISK MULTIPLY");

        //////////////////////////////////////////////////////////
        // Case I. A,B,C are in core. Perform direct a BLAS call
        // in one pass.
        //////////////////////////////////////////////////////////
        if (!B_on_disk && !C_on_disk) {
            size_t offset = 0;
            A_matrix = T_matrix[h];
            B_matrix = B_Matrix->get_matrix()[h];
            C_matrix = C_Matrix->get_matrix()[h];
            size_t rows_A = T_left->get_pairpi(h);
            size_t cols_A = T_right->get_pairpi(h);
            size_t rows_B = B_Matrix->get_left_pairpi(h);
            size_t cols_B = B_Matrix->get_right_pairpi(h);
            size_t rows_C = C_Matrix->get_left_pairpi(h);
            size_t cols_C = C_Matrix->get_right_pairpi(h);
            // Now call BLAS
            // Start a timer
            Timer timer;
            // Do the job
            contract_in_core(A_matrix, B_matrix, C_matrix, B_on_disk, C_on_disk, rows_A, rows_B, rows_C, cols_A, cols_B,
                             cols_C, offset);
            // Store the timing in moinfo
            wfn_->moinfo()->add_dgemm_timing(timer.get());
        }

        //////////////////////////////////////////////////////////
        // Case II. A,C are in core. B is on disk. Perform several
        // BLAS calls.
        //////////////////////////////////////////////////////////
        if (B_on_disk && !C_on_disk) {
            // Assign pointers to in core matrices
            double** A_matrix = T_matrix[h];
            double** B_matrix = new double*[1];
            B_matrix[0] = &out_of_core_buffer[0];
            double** C_matrix = C_Matrix->get_matrix()[h];

            int strip = 0;
            size_t offset = 0;
            bool done = false;
            while (!done) {
                size_t strip_length = B_Matrix->read_strip_from_disk(h, strip, out_of_core_buffer);
                if (strip_length == 0) {
                    done = true;
                } else {
                    // Compute the size of the sub-matrix B
                    size_t rows_A = T_left->get_pairpi(h);
                    size_t cols_A = T_right->get_pairpi(h);
                    size_t rows_B = strip_length;
                    size_t cols_B = B_Matrix->get_right_pairpi(h);
                    size_t rows_C = C_Matrix->get_left_pairpi(h);
                    size_t cols_C = C_Matrix->get_right_pairpi(h);

                    // Now call BLAS
                    // Start a timer
                    Timer timer;
                    // Do the job
                    contract_in_core(A_matrix, B_matrix, C_matrix, B_on_disk, C_on_disk, rows_A, rows_B, rows_C, cols_A,
                                     cols_B, cols_C, offset);
                    // Store the timing in moinfo
                    wfn_->moinfo()->add_dgemm_timing(timer.get());
                    offset += strip_length;
                }
                strip++;
            }
            delete[] B_matrix;
        }
        //////////////////////////////////////////////////////////
        // Case III. A,B are in core. C is on disk. Perform several
        // BLAS calls.
        //////////////////////////////////////////////////////////
        if (!B_on_disk && C_on_disk) {
            // Assign pointers to in core matrices
            double** A_matrix = T_matrix[h];
            double** B_matrix = B_Matrix->get_matrix()[h];
            double** C_matrix = new double*[1];
            C_matrix[0] = &out_of_core_buffer[0];

            int strip = 0;
            size_t offset = 0;
            bool done = false;
            while (!done) {
                size_t strip_length = C_Matrix->read_strip_from_disk(h, strip, out_of_core_buffer);
                if (strip_length == 0) {
                    done = true;
                } else {
                    // Compute the size of the sub-matrix B
                    size_t rows_A = T_left->get_pairpi(h);
                    size_t cols_A = T_right->get_pairpi(h);
                    size_t rows_B = B_Matrix->get_left_pairpi(h);
                    size_t cols_B = B_Matrix->get_right_pairpi(h);
                    size_t rows_C = strip_length;
                    size_t cols_C = C_Matrix->get_right_pairpi(h);

                    // Now call BLAS
                    // Start a timer
                    Timer timer;
                    // Do the job
                    contract_in_core(A_matrix, B_matrix, C_matrix, B_on_disk, C_on_disk, rows_A, rows_B, rows_C, cols_A,
                                     cols_B, cols_C, offset);
                    // Store the timing in moinfo
                    wfn_->moinfo()->add_dgemm_timing(timer.get());
                    offset += strip_length;
                }
                strip++;
            }
            delete[] C_matrix;
        }
    }  // end of for loop over irreps

    PartB_timing += PartB.get();
    Timer PartC;
    if (need_sort) {
        sort(T_left, T_right, T_matrix, 1.0);
        for (int h = 0; h < wfn_->moinfo()->get_nirreps(); h++)
            //       if(T_left->get_pairpi(h)*T_right->get_pairpi(h)>0)
            delete[] T_matrix[h];
    }
    PartC_timing += PartC.get();
}

void CCOperation::contract_in_core(double** A_matrix, double** B_matrix, double** C_matrix, bool B_on_disk,
                                   bool C_on_disk, int rows_A, int rows_B, int rows_C, int cols_A, int cols_B,
                                   int cols_C, int offset) {
    // CASE A, worst case scenario
    //
    //            -----------------------------
    //           |         ---------------     |
    //           |        |               |    |
    //   A[A_l][A_r] = B[B_i][B_c] 1@1 C[C_c][C_i]
    //      |                  |
    //       ------------------
    //  Solution:
    //  1) Transpose B,C
    //  2) Perform contraction of the indices
    //  3) Transpose the result and store in A ???
    if (operation == "1@1") {
        double beta = 1.0;
        int m = rows_A;
        int n = cols_A;
        int k;
        if (m * n == 0) return;
        if (!B_on_disk && !C_on_disk) {
            k = rows_B;
        }
        if (B_on_disk && !C_on_disk) {
            k = rows_B;
        }
        if (!B_on_disk && C_on_disk) {
            k = rows_C;
        }
        if (k != 0) {
            F_DGEMM("n", "t", &n, &m, &k, &factor, &(C_matrix[(B_on_disk ? offset : 0)][0]), &n,
                    &(B_matrix[(C_on_disk ? offset : 0)][0]), &m, &beta, &(A_matrix[0][0]), &n);
        }
    }
    // CASE B
    //
    //            ------------------------
    //           |         ---------------|----
    //           |        |               |    |
    //   A[A_l][A_r] = B[B_i][B_c] 1@2 C[C_c][C_i]
    //      |                  |
    //       ------------------
    //  Solution:
    //  1) Transpose B
    //  2) Perform contraction of the indices ???
    //  3) Transpose the result and store in A ??
    if (operation == "1@2") {
        double beta = 1.0;
        int m = rows_A;
        int n = rows_C;
        int k;
        if (m * n == 0) return;
        k = rows_B;
        if (k != 0) {
            F_DGEMM("t", "t", &n, &m, &k, &factor, &(C_matrix[0][(B_on_disk ? offset : 0)]), &cols_C, &(B_matrix[0][0]),
                    &m, &beta, &(A_matrix[0][(C_on_disk ? offset : 0)]), &cols_A);
        }
    }
    // CASE C,
    //
    //       -------------       ---------
    //      |             |     |         |
    //   A[A_l][A_r] = B[B_i][B_c] 2@1 C[C_c][C_i]
    //           |                             |
    //            -----------------------------
    //  Solution:
    //  Perform contraction of the indices by calling BLAS with swapped C and B
    if (operation == "2@1") {
        int m = rows_B;
        int n = cols_A;
        int k = rows_C;
        if (m * n == 0) return;
        double beta = 1.0;
        if (k != 0) {
            F_DGEMM("n", "n", &n, &m, &k, &factor, &(C_matrix[0][0]), &n, &(B_matrix[0][(C_on_disk ? offset : 0)]),
                    &cols_B, &beta, &(A_matrix[(B_on_disk ? offset : 0)][0]), &n);
        }
    }
    // CASE D, best case scenario
    //
    //       -------------       ---------------
    //      |             |     |               |
    //   A[A_l][A_r] = B[B_i][B_c] 2@2 C[C_c][C_i]
    //           |                        |
    //            ------------------------
    //  Solution:
    //  1) Perform contraction of the indices
    //  2) Store the result in A
    if (operation == "2@2") {
        int m = rows_B;
        int n = rows_C;
        int k = cols_C;
        if (m * n == 0) return;
        double beta = 1.0;
        if (k != 0) {
            F_DGEMM("t", "n", &n, &m, &k, &factor, &(C_matrix[0][0]), &k, &(B_matrix[0][0]), &k, &beta,
                    &(A_matrix[(B_on_disk ? offset : 0)][(C_on_disk ? offset : 0)]), &cols_A);
        }
    }
}

}  // namespace psimrcc
}  // namespace psi
