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

#include <iostream>
#include <cmath>
#include <algorithm>

#include "psi4/libciomr/libciomr.h"
#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libpsi4util/libpsi4util.h"

#include "index.h"
#include "matrix.h"
#include "psi4/libpsi4util/PsiOutStream.h"
namespace psi {

namespace psimrcc {

double CCMatrix::fraction_of_memory_for_buffer = 0.05;

CCMatrix::CCMatrix(std::string& str, CCIndex* left_index, CCIndex* right_index)
    : label(str),
      memory2(0),
      naccess(0),
      reference(-1),
      symmetry(-1),
      fock(false),
      integral(false),
      chemist_notation(false),
      antisymmetric(false),
      out_of_core(false),
      right(right_index),
      left(left_index) {
    if (left_index->wfn() != right_index->wfn()) {
        throw PSIEXCEPTION("\n Error in CCMatrix construction. The CCIndexes must belong to the same wavefunction.");
    } else {
        wfn_ = left_index->wfn();
    }
    nirreps = wfn_->moinfo()->get_nirreps();

    if (str.find("(") != std::string::npos || str.find("<") != std::string::npos) integral = true;
    if (str.find("(") != std::string::npos) chemist_notation = true;
    if (label.find(":") != std::string::npos) antisymmetric = true;
    if (str.find("fock") != std::string::npos) fock = true;

    // Copy the pairpi arrays from the CCIndex object
    // Compute the memory required to store the matrix in core

    matrix = std::vector<double**>(nirreps, nullptr);
    left_pairpi = std::vector<size_t>(nirreps, 0);
    right_pairpi = std::vector<size_t>(nirreps, 0);
    block_sizepi = std::vector<size_t>(nirreps, 0);

    for (int h = 0; h < nirreps; h++) {
        left_pairpi[h] = left->get_pairpi(h);
        right_pairpi[h] = right->get_pairpi(h);
        block_sizepi[h] = left_pairpi[h] * right_pairpi[h];
        memorypi2.push_back(static_cast<size_t>(sizeof(double)) * block_sizepi[h]);
        memory2 += memorypi2[h];
        out_of_core.push_back(false);
    }

    // Get a copy of the indices
    index_label = compute_index_label();

    // Parse the curly braces to get the reference
    std::string::size_type left_curly = str.find("{");
    std::string::size_type right_curly = str.find("}");
    if (left_curly != std::string::npos &&
        right_curly != std::string::npos)  // TODO add check on the size of the string
        reference = to_integer(str.substr(left_curly + 1, right_curly - left_curly - 1));
}

CCMatrix::~CCMatrix() { free_memory(); }

/*********************************************************
  Printing Routines
*********************************************************/

/*!
    \fn CCMatrix::print()
 */
void CCMatrix::print() {
    outfile->Printf("\n\n\t\t\t\t\t%s Matrix\n", label.c_str());
    for (int i = 0; i < nirreps; i++) {
        if (left->get_pairpi(i) * right->get_pairpi(i)) {
            outfile->Printf("\nBlock %d (%s,%s)", i, wfn_->moinfo()->irr_lab_ref(i).c_str(),
                            wfn_->moinfo()->irr_lab_ref(i).c_str());
            print_dpdmatrix(i, "outfile");
        }
    }
}

void CCMatrix::add_numerical_factor(double factor) {
    for (int h = 0; h < nirreps; h++) add_numerical_factor(factor, h);
}

void CCMatrix::add_numerical_factor(double factor, int h) {
    if (block_sizepi[h] > 0) {
        double* matrix_block = &matrix[h][0][0];
        for (size_t i = 0; i != block_sizepi[h]; ++i) matrix_block[i] += factor;
    }
}

void CCMatrix::scale(double factor) {
    for (int h = 0; h < nirreps; h++) scale(factor, h);
}

void CCMatrix::scale(double factor, int h) {
    if (block_sizepi[h] > 0) {
        double* matrix_block = &matrix[h][0][0];
        for (size_t i = 0; i != block_sizepi[h]; ++i) matrix_block[i] *= factor;
    }
}

void CCMatrix::zero_matrix() {
    for (int h = 0; h < nirreps; h++) zero_matrix_block(h);
}

void CCMatrix::zero_matrix_block(int h) {
    if (block_sizepi[h] > 0) zero_arr(&(matrix[h][0][0]), block_sizepi[h]);
}

void CCMatrix::zero_two_diagonal() {
    auto* pq = new short[2];
    for (int h = 0; h < nirreps; h++)
        for (int i = 0; i < left->get_pairpi(h); i++)
            for (int j = 0; j < right->get_pairpi(h); j++) {
                get_two_indices(pq, h, i, j);
                if (pq[0] == pq[1]) matrix[h][i][j] = 0.0;
            }
    delete[] pq;
}

void CCMatrix::zero_non_doubly_occupied() {
    const boolvec& is_act_in_occ = wfn_->moinfo()->get_is_actv_in_occ();
    auto* pq = new short[2];
    for (int h = 0; h < nirreps; h++)
        for (int i = 0; i < left->get_pairpi(h); i++)
            for (int j = 0; j < right->get_pairpi(h); j++) {
                get_two_indices(pq, h, i, j);
                if (is_act_in_occ[pq[0]] && !is_act_in_occ[pq[1]]) matrix[h][i][j] = 0.0;
                if (!is_act_in_occ[pq[0]] && is_act_in_occ[pq[1]]) matrix[h][i][j] = 0.0;
            }
    delete[] pq;
}

void CCMatrix::zero_non_external() {
    const boolvec& is_act_in_vir = wfn_->moinfo()->get_is_actv_in_vir();
    auto* pq = new short[2];
    for (int h = 0; h < nirreps; h++)
        for (int i = 0; i < left->get_pairpi(h); i++)
            for (int j = 0; j < right->get_pairpi(h); j++) {
                get_two_indices(pq, h, i, j);
                if (is_act_in_vir[pq[0]] && !is_act_in_vir[pq[1]]) matrix[h][i][j] = 0.0;
                if (!is_act_in_vir[pq[0]] && is_act_in_vir[pq[1]]) matrix[h][i][j] = 0.0;
            }
    delete[] pq;
}

void CCMatrix::zero_right_four_diagonal() {
    auto* pqrs = new short[4];
    for (int h = 0; h < nirreps; h++)
        for (int j = 0; j < right->get_pairpi(h); j++) {
            if (left->get_pairpi(h) > 0) {
                get_four_indices(pqrs, h, 0, j);
                if (pqrs[2] == pqrs[3])
                    for (int i = 0; i < left->get_pairpi(h); i++) matrix[h][i][j] = 0.0;
            }
        }
    delete[] pqrs;
}

void CCMatrix::zero_left_four_diagonal() {
    auto* pqrs = new short[4];
    for (int h = 0; h < nirreps; h++)
        for (int i = 0; i < left->get_pairpi(h); i++) {
            if (right->get_pairpi(h) > 0) {
                get_four_indices(pqrs, h, i, 0);
                if (pqrs[0] == pqrs[1])
                    for (int j = 0; j < right->get_pairpi(h); j++) matrix[h][i][j] = 0.0;
            }
        }
    delete[] pqrs;
}

void CCMatrix::element_by_element_product(double factor, CCMatrix* B_Matrix, CCMatrix* C_Matrix, int h) {
    if (block_sizepi[h] > 0) {
        double* A_matrix = &(matrix[h][0][0]);
        double* B_matrix = &(B_Matrix->get_matrix()[h][0][0]);
        double* C_matrix = &(C_Matrix->get_matrix()[h][0][0]);
        for (size_t i = 0; i < block_sizepi[h]; i++) A_matrix[i] += factor * B_matrix[i] * C_matrix[i];
    }
}

void CCMatrix::element_by_element_division(double factor, CCMatrix* B_Matrix, CCMatrix* C_Matrix, int h) {
    if (block_sizepi[h] > 0) {
        double* A_matrix = &(matrix[h][0][0]);
        double* B_matrix = &(B_Matrix->get_matrix()[h][0][0]);
        double* C_matrix = &(C_Matrix->get_matrix()[h][0][0]);
        for (size_t i = 0; i < block_sizepi[h]; i++) A_matrix[i] += factor * B_matrix[i] / C_matrix[i];
    }
}

void CCMatrix::element_by_element_addition(double factor, CCMatrix* B_Matrix, int h) {
    if (block_sizepi[h] > 0) {
        double* A_matrix = &(matrix[h][0][0]);
        double* B_matrix = &(B_Matrix->get_matrix()[h][0][0]);
        for (size_t i = 0; i < block_sizepi[h]; i++) A_matrix[i] += factor * B_matrix[i];
    }
}

void CCMatrix::tensor_product(std::string& reindexing, double factor, CCMatrix* B_Matrix, CCMatrix* C_Matrix) {
    auto* reindexing_array = new short[4];

    intpairvec pairs;
    for (int i = 0; i < reindexing.size(); i++) pairs.push_back(std::make_pair(to_integer(reindexing.substr(i, 1)), i));
    sort(pairs.begin(), pairs.end());
    for (int i = 0; i < reindexing.size(); i++) reindexing_array[i] = pairs[i].second;

    // This assumes that the reindexing starts from 1 !!! This can cost you an headache
    auto* pqrs = new short[4];
    auto* pq = new short[2];
    auto* rs = new short[2];
    auto B_matrix = B_Matrix->get_matrix();
    auto C_matrix = C_Matrix->get_matrix();
    double value;
    for (int b_n = 0; b_n < wfn_->moinfo()->get_nirreps(); b_n++) {
        for (int c_n = 0; c_n < wfn_->moinfo()->get_nirreps(); c_n++) {
            for (int b_i = 0; b_i < B_Matrix->get_left_pairpi(b_n); b_i++) {
                for (int b_j = 0; b_j < B_Matrix->get_right_pairpi(b_n); b_j++) {
                    for (int c_i = 0; c_i < C_Matrix->get_left_pairpi(c_n); c_i++) {
                        for (int c_j = 0; c_j < C_Matrix->get_right_pairpi(c_n); c_j++) {
                            value = factor * B_matrix[b_n][b_i][b_j] * C_matrix[c_n][c_i][c_j];
                            B_Matrix->get_two_indices(pq, b_n, b_i, b_j);
                            C_Matrix->get_two_indices(rs, c_n, c_i, c_j);
                            pqrs[0] = pq[0];
                            pqrs[1] = pq[1];
                            pqrs[2] = rs[0];
                            pqrs[3] = rs[1];
                            add_four_address_element(pqrs[reindexing_array[0]], pqrs[reindexing_array[1]],
                                                     pqrs[reindexing_array[2]], pqrs[reindexing_array[3]], value);
                        }
                    }
                }
            }
        }
    }
    delete[] pqrs;
    delete[] pq;
    delete[] rs;
    delete[] reindexing_array;
}

double CCMatrix::dot_product(CCMatrix* B_Matrix, CCMatrix* C_Matrix, int h) {
    double value = 0.0;
    size_t block_size = B_Matrix->get_block_sizepi(h);
    if (block_size > 0) {
        size_t i;
        double* B_matrix = &(B_Matrix->get_matrix()[h][0][0]);
        double* C_matrix = &(C_Matrix->get_matrix()[h][0][0]);
        for (size_t i = 0; i < block_size; i++) value += B_matrix[i] * C_matrix[i];
    }
    return (value);
}

void CCMatrix::print_dpdmatrix(int irrep, std::string out_fname) {
    int ii, jj, kk, nn, ll;
    int i, j;
    auto mode = std::ostream::app;
    auto printer = out_fname == "outfile" ? psi::outfile : std::make_shared<psi::PsiOutStream>(out_fname, mode);
    double** mat = matrix[irrep];
    int left_offset = left->get_first(irrep);
    int right_offset = right->get_first(irrep);

    int m = left->get_pairpi(irrep);
    int n = right->get_pairpi(irrep);

    ii = 0;
    jj = 0;
L200:
    ii++;
    jj++;
    kk = 10 * jj;
    nn = n;
    if (nn > kk) nn = kk;
    ll = 2 * (nn - ii + 1) + 1;
    printer->Printf("\n            ");
    for (i = ii; i <= nn; i++) {
        auto& right_indices = right->get_tuple(i + right_offset - 1);
        printer->Printf("(");
        for (int p = 0; p < right->get_nelements(); p++) printer->Printf("%3d", right_indices[p]);
        printer->Printf(")");
        int nspaces = 10 - 3 * right->get_nelements();
        for (int p = 0; p < nspaces; p++) printer->Printf(" ");
    }
    printer->Printf("\n");
    for (i = 0; i < m; i++) {
        auto& left_indices = left->get_tuple(i + left_offset);
        printer->Printf("\n(");
        for (int p = 0; p < left->get_nelements(); p++) printer->Printf("%3d", left_indices[p]);
        printer->Printf(")  ");

        for (j = ii - 1; j < nn; j++) {
            if (std::fabs(mat[i][j]) < 100.0)
                printer->Printf("%12.7f", mat[i][j]);
            else
                printer->Printf("    infinity");
        }
    }
    printer->Printf("\n");
    if (n <= kk) {
        return;
    }
    ii = kk;
    goto L200;
}

void CCMatrix::set_scalar(double val) { matrix[0][0][0] = val; }

void CCMatrix::add_scalar(double val) { matrix[0][0][0] += val; }

double CCMatrix::get_scalar() { return (matrix[0][0][0]); }

std::string CCMatrix::compute_index_label() {
    std::string label;
    int left_indices = left->get_label().size();
    if (left_indices > 2) label += left->get_label().substr(1, left_indices - 2);
    int right_indices = right->get_label().size();
    if (right_indices > 2) label += right->get_label().substr(1, right_indices - 2);
    return (label);
}

}  // namespace psimrcc
}  // namespace psi
