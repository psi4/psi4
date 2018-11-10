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

// Latest revision on April 38, 2013.
#include <algorithm>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include "psi4/libqt/qt.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/libiwl/iwl.hpp"
#include "tensors.h"
#include "tensors_float.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"

namespace psi {
namespace dfoccwave {

/********************************************************************************************/
/************************** 1d array ********************************************************/
/********************************************************************************************/
Tensor1f::Tensor1f(int d1) {
    A1d_ = NULL;
    dim1_ = d1;
    memalloc();
}  //

Tensor1f::Tensor1f(std::string name, int d1) {
    A1d_ = NULL;
    dim1_ = d1;
    name_ = name;
    memalloc();
}  //

Tensor1f::Tensor1f() {
    A1d_ = NULL;
    dim1_ = 0;

}  //

Tensor1f::~Tensor1f() { release(); }  //

void Tensor1f::memalloc() {
    if (A1d_) release();
    A1d_ = new float[dim1_];
    zero();
}  //

void Tensor1f::release() {
    if (!A1d_) return;
    delete[] A1d_;
    A1d_ = NULL;
}  //

void Tensor1f::init(int d1) {
    dim1_ = d1;
    if (A1d_) release();
    A1d_ = new float[dim1_];
}  //

void Tensor1f::init(std::string name, int d1) {
    dim1_ = d1;
    name_ = name;
    if (A1d_) release();
    A1d_ = new float[dim1_];
}  //

void Tensor1f::zero() { memset(A1d_, 0, sizeof(float) * dim1_); }  //

void Tensor1f::print() {
    if (name_.length()) outfile->Printf("\n ## %s ##\n", name_.c_str());
    for (int p = 0; p < dim1_; p++) {
        outfile->Printf(" %3d %10.7f \n", p, A1d_[p]);
    }

}  //

void Tensor1f::print(std::string out_fname) {
    std::shared_ptr<psi::PsiOutStream> printer =
        (out_fname == "outfile" ? outfile : std::shared_ptr<PsiOutStream>(new PsiOutStream(out_fname)));
    if (name_.length()) printer->Printf("\n ## %s ##\n", name_.c_str());
    for (int p = 0; p < dim1_; p++) {
        printer->Printf(" %3d %10.7f \n", p, A1d_[p]);
    }
}  //

void Tensor1f::print(FILE *out) {
    if (name_.length()) fprintf(out, "\n ## %s ##\n", name_.c_str());
    for (int p = 0; p < dim1_; p++) {
        fprintf(out, " %3d %10.7f \n", p, A1d_[p]);
    }
    fflush(out);
}  //

void Tensor1f::print(const char *outfile) {
    // Open the file
    std::ofstream out(outfile, std::ios::app);
    out.precision(6);

    if (name_.length()) out << "\n ## %s ##\n" << name_.c_str();
    for (int p = 0; p < dim1_; p++) {
        out << " %3d %10.7f \n" << p << A1d_[p];
    }

    // Close output file
    out.close();
}  //

void Tensor1f::set(int i, float value) { A1d_[i] = value; }  //

void Tensor1f::set(float *vec) {
    for (int i = 0; i < dim1_; ++i) A1d_[i] = vec[i];
}  //

void Tensor1f::set(const SharedTensor1f &vec) {
    for (int i = 0; i < dim1_; ++i) A1d_[i] = vec->A1d_[i];
}  //

float Tensor1f::get(int i) { return A1d_[i]; }  //

void Tensor1f::add(const SharedTensor1f &a) {
/*
float *lhs, *rhs;
size_t size = dim1_;
if (size) {
    lhs = A1d_;
    rhs = Adum->A1d_;
    for (size_t ij=0; ij<size; ++ij) {
        *lhs += *rhs;
        lhs++; rhs++;
    }
}
*/
#pragma omp parallel for
    for (int i = 0; i < dim1_; ++i) A1d_[i] += a->A1d_[i];

}  //

void Tensor1f::add(int i, float value) { A1d_[i] += value; }  //

void Tensor1f::subtract(const SharedTensor1f &a) {
/*
float *lhs, *rhs;
size_t size = dim1_;
if (size) {
    lhs = A1d_;
    rhs = Adum->A1d_;
    for (size_t ij=0; ij<size; ++ij) {
        *lhs -= *rhs;
        lhs++; rhs++;
    }
}
*/
#pragma omp parallel for
    for (int i = 0; i < dim1_; ++i) A1d_[i] -= a->A1d_[i];
}  //

void Tensor1f::subtract(int i, float value) { A1d_[i] -= value; }  //

void Tensor1f::to_shared_vector(SharedVector A) {
#pragma omp parallel for
    for (int i = 0; i < dim1_; ++i) {
        A->set(0, i, A1d_[i]);
    }
}  //

float Tensor1f::rms() {
    float summ = 0.0;
    for (int i = 0; i < dim1_; ++i) summ += A1d_[i] * A1d_[i];
    summ = std::sqrt(summ / dim1_);

    return summ;
}  //

float Tensor1f::rms(const SharedTensor1f &Atemp) {
    float summ = 0.0;
    for (int i = 0; i < dim1_; ++i) summ += (A1d_[i] - Atemp->A1d_[i]) * (A1d_[i] - Atemp->A1d_[i]);
    summ = std::sqrt(summ / dim1_);

    return summ;
}  //

float Tensor1f::dot(const SharedTensor1f &y) {
    float value = 0.0;
    int incx = 1;
    int incy = 1;
    if (dim1_ == y->dim1_) value = C_SDOT((size_t)dim1_, A1d_, incx, y->A1d_, incy);
    return value;
}  //

// void Tensor1f::gbmv(bool transa, const SharedTensor2f &a, const SharedTensor1f &b, float alpha, float beta) {
//     char ta = transa ? 't' : 'n';
//     int m, n, k, kl, ku, incx, incy, lda;

//     m = a->dim1_;
//     n = a->dim2_;
//     k = b->dim1_;
//     kl = m - 1;  // # of subdiagonal of matrix A, at most kl = m - 1
//     ku = n - 1;  // # of superdiagonal of matrix A, at most ku = n - 1
//     lda = kl + ku + 1;
//     incx = 1;  // increments in elements of b vector
//     incy = 1;  // increments in elements of A1d_

//     // A1d_ = alpha * A * b + beta, where A is a band matrix
//     if (m && n) {
//         C_SGBMV(ta, m, n, kl, ku, alpha, &(a->A2d_[0][0]), lda, b->A1d_, incx, beta, A1d_, incy);
//     }
// }  //

void Tensor1f::gemv(bool transa, const SharedTensor2f &a, const SharedTensor1f &b, float alpha, float beta) {
    char ta = transa ? 't' : 'n';
    int m, n, k, incx, incy, lda;

    m = a->dim1();
    n = a->dim2();
    lda = n;
    incx = 1;  // increments in elements of b vector
    incy = 1;  // increments in elements of A1d_

    // A1d_ = alpha * A * b + beta, where A is a general matrix
    if (m && n) {
        C_SGEMV(ta, m, n, alpha, &(a->A2d_[0][0]), lda, b->A1d_, incx, beta, A1d_, incy);
    }
}  //

void Tensor1f::gemv(bool transa, int m, int n, const SharedTensor2f &a, const SharedTensor2f &b, float alpha,
                    float beta) {
    char ta = transa ? 't' : 'n';
    int incx, incy, lda;

    lda = n;
    incx = 1;  // increments in elements of b vector
    incy = 1;  // increments in elements of A1d_

    // A1d_ = alpha * A * b + beta, where A is a general matrix
    if (m && n) {
        C_SGEMV(ta, m, n, alpha, &(a->A2d_[0][0]), lda, &(b->A2d_[0][0]), incx, beta, A1d_, incy);
    }
}  //

void Tensor1f::gemv(bool transa, const SharedTensor2f &a, const SharedTensor2f &b, float alpha, float beta) {
    char ta = transa ? 't' : 'n';
    int incx, incy, lda, m, n;

    m = a->dim1();
    n = a->dim2();

    /*
    if (transa) {
        m = a->dim1();
        n = a->dim2();
    }

    else {
        n = a->dim1();
        m = a->dim2();
    }
    */

    lda = n;
    incx = 1;  // increments in elements of b vector
    incy = 1;  // increments in elements of A1d_

    // A1d_ = alpha * A * b + beta, where A is a general matrix
    if (m && n) {
        C_SGEMV(ta, m, n, alpha, &(a->A2d_[0][0]), lda, &(b->A2d_[0][0]), incx, beta, A1d_, incy);
    }
}  //

void Tensor1f::gemv(bool transa, int m, int n, const SharedTensor2f &a, const SharedTensor2f &b, int start_a,
                    int start_b, float alpha, float beta) {
    char ta = transa ? 't' : 'n';
    int incx, incy, lda;

    lda = n;
    incx = 1;  // increments in elements of b vector
    incy = 1;  // increments in elements of A1d_

    // A1d_ = alpha * A * b + beta, where A is a general matrix
    if (m && n) {
        C_SGEMV(ta, m, n, alpha, a->A2d_[0] + start_a, lda, b->A2d_[0] + start_b, incx, beta, A1d_, incy);
    }
}  //

void Tensor1f::gemv(bool transa, int m, int n, const SharedTensor2f &a, const SharedTensor2f &b, int start_a,
                    int start_b, int start_c, float alpha, float beta) {
    char ta = transa ? 't' : 'n';
    int incx, incy, lda;

    lda = n;
    incx = 1;  // increments in elements of b vector
    incy = 1;  // increments in elements of A1d_

    // A1d_ = alpha * A * b + beta, where A is a general matrix
    if (m && n) {
        C_SGEMV(ta, m, n, alpha, a->A2d_[0] + start_a, lda, b->A2d_[0] + start_b, incx, beta, A1d_ + start_c, incy);
    }
}  //

float Tensor1f::xay(const SharedTensor2f &a, const SharedTensor1f &y) {
    float value = 0.0;
    SharedTensor1f ay = SharedTensor1f(new Tensor1f(a->dim1_));
    ay->gemv(false, a, y, 1.0, 0.0);
    value = dot(ay);
    return value;
}  //

void Tensor1f::axpy(const SharedTensor1f &a, float alpha) {
    size_t length = (size_t)dim1_;
    C_SAXPY(length, alpha, a->A1d_, 1, A1d_, 1);
}

void Tensor1f::scale(float a) {
    // size_t size = dim1_ ;
    size_t size = (size_t)dim1_;
    if (size) C_SSCAL(size, a, A1d_, 1);
}  //

void Tensor1f::copy(float *a) {
    // size_t size;
    // size = dim1_ * sizeof(float);
    // if (size) memcpy(&(A1d_[0]), &(x[0]), size);
    size_t size = (size_t)dim1_;
    C_SCOPY(size, a, 1, A1d_, 1);
}  //

void Tensor1f::copy(const SharedTensor1f &a) {
    // size_t size;
    // size = dim1_ * sizeof(float);
    // if (size) memcpy(&(A1d_[0]), &(x->A1d_[0]), size);
    size_t size = (size_t)dim1_;
    C_SCOPY(size, a->A1d_, 1, A1d_, 1);
}  //

void Tensor1f::row_vector(SharedTensor2f &A, int n) {
    int dim = A->dim2();
    for (int i = 0; i < dim; i++) A1d_[i] = A->get(n, i);
}  //

void Tensor1f::column_vector(SharedTensor2f &A, int n) {
    int dim = A->dim1();
    for (int i = 0; i < dim; i++) A1d_[i] = A->get(i, n);
}  //

void Tensor1f::dirprd(SharedTensor1f &a, SharedTensor1f &b) {
    int dima = a->dim1();
    int dimb = b->dim1();

    if (dima == dimb && dima == dim1_) {
        for (int i = 0; i < dim1_; i++) A1d_[i] = a->get(i) * b->get(i);
    } else
        throw SanityCheckError("Vector dimensions do NOT match!", __FILE__, __LINE__);
}  //

void Tensor1f::symm_packed(const SharedTensor2f &A) {
// Form Lower triangular part
#pragma omp parallel for
    for (int p = 0; p < A->dim1(); p++) {
        for (int q = 0; q <= p; q++) {
            int pq = index2(p, q);
            float perm = (p == q ? 1.0 : 2.0);
            A1d_[pq] = perm * A->get(p, q);
        }
    }

}  //

void Tensor1f::ltm(const SharedTensor2f &A) {
// Form Lower triangular part
#pragma omp parallel for
    for (int p = 0; p < A->dim1(); p++) {
        for (int q = 0; q <= p; q++) {
            int pq = index2(p, q);
            A1d_[pq] = A->get(p, q);
        }
    }

}  //

/********************************************************************************************/
/************************** 2d array ********************************************************/
/********************************************************************************************/
Tensor2f::Tensor2f(int d1, int d2) {
    A2d_ = NULL;
    row_idx_ = NULL;
    col_idx_ = NULL;
    row2d1_ = NULL;
    row2d2_ = NULL;
    col2d1_ = NULL;
    col2d2_ = NULL;
    d1_ = 0;
    d2_ = 0;
    dim1_ = d1;
    dim2_ = d2;
    memalloc();
}  //

Tensor2f::Tensor2f(std::string name, int d1, int d2) {
    A2d_ = NULL;
    row_idx_ = NULL;
    col_idx_ = NULL;
    row2d1_ = NULL;
    row2d2_ = NULL;
    col2d1_ = NULL;
    col2d2_ = NULL;
    d1_ = 0;
    d2_ = 0;
    dim1_ = d1;
    dim2_ = d2;
    name_ = name;
    memalloc();
}  //

Tensor2f::Tensor2f() {
    A2d_ = NULL;
    row_idx_ = NULL;
    col_idx_ = NULL;
    row2d1_ = NULL;
    row2d2_ = NULL;
    col2d1_ = NULL;
    col2d2_ = NULL;
    d1_ = 0;
    d2_ = 0;
    dim1_ = 0;
    dim2_ = 0;

}  //

Tensor2f::Tensor2f(psi::PSIO *psio, size_t fileno, std::string name, int d1, int d2) {
    A2d_ = NULL;
    row_idx_ = NULL;
    col_idx_ = NULL;
    row2d1_ = NULL;
    row2d2_ = NULL;
    col2d1_ = NULL;
    col2d2_ = NULL;
    d1_ = 0;
    d2_ = 0;
    dim1_ = d1;
    dim2_ = d2;
    name_ = name;
    memalloc();
    read(psio, fileno);
}

Tensor2f::Tensor2f(std::shared_ptr<psi::PSIO> psio, size_t fileno, std::string name, int d1, int d2) {
    A2d_ = NULL;
    row_idx_ = NULL;
    col_idx_ = NULL;
    row2d1_ = NULL;
    row2d2_ = NULL;
    col2d1_ = NULL;
    col2d2_ = NULL;
    d1_ = 0;
    d2_ = 0;
    dim1_ = d1;
    dim2_ = d2;
    name_ = name;
    memalloc();
    read(psio, fileno);
}

Tensor2f::Tensor2f(psi::PSIO &psio, size_t fileno, std::string name, int d1, int d2) {
    A2d_ = NULL;
    row_idx_ = NULL;
    col_idx_ = NULL;
    row2d1_ = NULL;
    row2d2_ = NULL;
    col2d1_ = NULL;
    col2d2_ = NULL;
    d1_ = 0;
    d2_ = 0;
    dim1_ = d1;
    dim2_ = d2;
    name_ = name;
    memalloc();
    read(&psio, fileno);
}  //

Tensor2f::Tensor2f(std::string name, int d1, int d2, int d3, int d4) {
    A2d_ = NULL;
    row_idx_ = NULL;
    col_idx_ = NULL;
    row2d1_ = NULL;
    row2d2_ = NULL;
    col2d1_ = NULL;
    col2d2_ = NULL;
    d1_ = d1;
    d2_ = d2;
    d3_ = d3;
    d4_ = d4;
    dim1_ = d1 * d2;
    dim2_ = d3 * d4;
    name_ = name;

    // memalloc
    if (A2d_) release();
    A2d_ = block_matrix_float   (dim1_, dim2_);
    zero();

    // row idx
    row_idx_ = init_int_matrix(d1_, d2_);
    memset(row_idx_[0], 0, sizeof(int) * d1_ * d2_);
    row2d1_ = new int[dim1_];
    row2d2_ = new int[dim1_];
    memset(row2d1_, 0, sizeof(int) * dim1_);
    memset(row2d2_, 0, sizeof(int) * dim1_);
    for (int i = 0; i < d1_; i++) {
        for (int a = 0; a < d2_; a++) {
            int ia = a + (i * d2_);
            row_idx_[i][a] = ia;
            row2d1_[ia] = i;
            row2d2_[ia] = a;
        }
    }

    // col idx
    col_idx_ = init_int_matrix(d3_, d4_);
    memset(col_idx_[0], 0, sizeof(int) * d3_ * d4_);
    col2d1_ = new int[dim2_];
    col2d2_ = new int[dim2_];
    memset(col2d1_, 0, sizeof(int) * dim2_);
    memset(col2d2_, 0, sizeof(int) * dim2_);
    for (int i = 0; i < d3_; i++) {
        for (int a = 0; a < d4_; a++) {
            int ia = a + (i * d4_);
            col_idx_[i][a] = ia;
            col2d1_[ia] = i;
            col2d2_[ia] = a;
        }
    }

}  //

Tensor2f::Tensor2f(std::string name, int d1, int d2, int d3) {
    A2d_ = NULL;
    row_idx_ = NULL;
    col_idx_ = NULL;
    row2d1_ = NULL;
    row2d2_ = NULL;
    col2d1_ = NULL;
    col2d2_ = NULL;
    d1_ = d1;
    d2_ = d2;
    d3_ = d3;
    d4_ = 0;
    dim1_ = d1;
    dim2_ = d2 * d3;
    name_ = name;

    // memalloc
    if (A2d_) release();
    A2d_ = block_matrix_float   (dim1_, dim2_);
    zero();

    // col idx
    col_idx_ = init_int_matrix(d2_, d3_);
    memset(col_idx_[0], 0, sizeof(int) * d2_ * d3_);
    col2d1_ = new int[dim2_];
    col2d2_ = new int[dim2_];
    memset(col2d1_, 0, sizeof(int) * dim2_);
    memset(col2d2_, 0, sizeof(int) * dim2_);
    for (int i = 0; i < d2_; i++) {
        for (int a = 0; a < d3_; a++) {
            int ia = a + (i * d3_);
            col_idx_[i][a] = ia;
            col2d1_[ia] = i;
            col2d2_[ia] = a;
        }
    }

}  //

Tensor2f::~Tensor2f() { release(); }  //

void Tensor2f::memalloc() {
    if (A2d_) release();
    A2d_ = block_matrix_float   (dim1_, dim2_);
    zero();
}  //

void Tensor2f::release() {
    // if (!A2d_) return;
    // free_block_float(A2d_);
    if (A2d_) free_block_float(A2d_);
    if (row_idx_) free_int_matrix(row_idx_);
    if (col_idx_) free_int_matrix(col_idx_);
    if (row2d1_) delete[] row2d1_;
    if (row2d2_) delete[] row2d2_;
    if (col2d1_) delete[] col2d1_;
    if (col2d2_) delete[] col2d2_;

    A2d_ = NULL;
    row_idx_ = NULL;
    col_idx_ = NULL;
    row2d1_ = NULL;
    row2d2_ = NULL;
    col2d1_ = NULL;
    col2d2_ = NULL;
}  //

void Tensor2f::init(int d1, int d2) {
    dim1_ = d1;
    dim2_ = d2;
    if (A2d_) release();
    A2d_ = block_matrix_float   (dim1_, dim2_);
}  //

void Tensor2f::init(std::string name, int d1, int d2) {
    dim1_ = d1;
    dim2_ = d2;
    name_ = name;
    if (A2d_) release();
    A2d_ = block_matrix_float   (dim1_, dim2_);
}  //

void Tensor2f::zero() { memset(A2d_[0], 0, sizeof(float) * dim1_ * dim2_); }  //

void Tensor2f::zero_diagonal() {
    if (dim1_ == dim2_) {
        for (int i = 0; i < dim1_; i++) A2d_[i][i] = 0.0;
    }
}  //

void Tensor2f::zero_off_diagonal() {
    for (int i = 0; i < dim1_; i++) {
        for (int j = 0; j < dim2_; j++) {
            if (i != j) A2d_[i][j] = 0.0;
        }
    }
}  //

void Tensor2f::print() {
    if (A2d_) {
        if (name_.length()) outfile->Printf("\n ## %s ##\n", name_.c_str());
        print_mat(A2d_, dim1_, dim2_, "outfile");
    }
}  //

void Tensor2f::print(std::string out_fname) {
    std::shared_ptr<psi::PsiOutStream> printer =
        (out_fname == "outfile" ? outfile : std::shared_ptr<PsiOutStream>(new PsiOutStream(out_fname)));
    if (A2d_) {
        if (name_.length()) printer->Printf("\n ## %s ##\n", name_.c_str());
        print_mat(A2d_, dim1_, dim2_, out_fname);
    }
}  //

/*
void Tensor2f::print(FILE *out)
{
  if (A2d_) {
      if (name_.length()) fprintf(out, "\n ## %s ##\n", name_.c_str());
      print_mat(A2d_,dim1_,dim2_,out);
      fflush(out);
  }
}//
*/

void Tensor2f::print(const char *outfile) {
    // Open the file
    std::ofstream out(outfile, std::ios::app);
    out.precision(6);
    if (name_.length()) out << "\n ## %s ##\n" << name_.c_str();

    int m = dim1_;
    int n = dim2_;

    int num_frames = int(n / 10);
    int num_frames_rem = n % 10;  // adding one for changing 0->1 start
    int num_frame_counter = 0;
    // for each frame
    for (num_frame_counter = 0; num_frame_counter < num_frames; num_frame_counter++) {
        out << "\n";
        for (int j = 10 * num_frame_counter + 1; j < 10 * num_frame_counter + 11; j++) {
            if (j == 10 * num_frame_counter + 1) {
                out << "%18d" << j;
            } else {
                out << "        %5d" << j;
            }
        }
        out << "\n\n";

        for (int k = 1; k <= m; ++k) {
            for (int j = 10 * num_frame_counter + 1; j < 10 * num_frame_counter + 12; j++) {
                if (j == 10 * num_frame_counter + 1) {
                    printf("%5d", k);
                }  // printf left here
                else {
                    out << " %12.7f" << A2d_[k - 1][j - 2];
                }
            }
            out << "\n";
        }
    }

    // ALREADY DID THE FULL FRAMES BY THIS POINT
    // NEED TO TAKE CARE OF THE REMAINDER
    if (num_frames_rem != 0) {
        out << "\n";
        for (int j = 10 * num_frame_counter + 1; j <= n; j++) {
            if (j == 10 * num_frame_counter + 1) {
                out << "%18d" << j;
            } else {
                out << "        %5d" << j;
            }
        }
        out << "\n\n";

        for (int k = 1; k <= m; ++k) {
            for (int j = 10 * num_frame_counter + 1; j < n + 2; j++) {
                if (j == 10 * num_frame_counter + 1) {
                    out << "%5d" << k;
                } else {
                    out << " %12.7f" << A2d_[k - 1][j - 2];
                }
            }
            out << "\n";
        }
    }
    out << "\n\n";

    // Close output file
    out.close();

}  //

void Tensor2f::set(int i, int j, float value) { A2d_[i][j] = value; }  //

void Tensor2f::set(float **A) {
    if (A == NULL) return;
#pragma omp parallel for
    for (int i = 0; i < dim1_; ++i) {
        for (int j = 0; j < dim2_; ++j) {
            A2d_[i][j] = A[i][j];
        }
    }
}  //

void Tensor2f::set(SharedTensor2f &A) {
    if (A == NULL) return;
#pragma omp parallel for
    for (int i = 0; i < dim1_; ++i) {
        for (int j = 0; j < dim2_; ++j) {
            A2d_[i][j] = A->A2d_[i][j];
        }
    }
}  //

void Tensor2f::set(SharedMatrix A) {
#pragma omp parallel for
    for (int i = 0; i < dim1_; ++i) {
        for (int j = 0; j < dim2_; ++j) {
            A2d_[i][j] = A->get(0, i, j);
        }
    }
}  //

void Tensor2f::set2(SharedMatrix A) {
#pragma omp parallel for
    for (int i = 0; i < dim1_; ++i) {
        for (int j = 0; j < dim2_; ++j) {
            A2d_[i][j] = A->get(i, j);
        }
    }
}  //

void Tensor2f::set(SharedTensor1f &A) {
#pragma omp parallel for
    for (int i = 0; i < dim1_; i++) {
        for (int j = 0; j < dim2_; j++) {
            int ij = j + (i * dim2_);
            A2d_[i][j] = A->get(ij);
        }
    }
}  //

void Tensor2f::set(float *A) {
#pragma omp parallel for
    for (int i = 0; i < dim1_; i++) {
        for (int j = 0; j < dim2_; j++) {
            int ij = j + (i * dim2_);
            A2d_[i][j] = A[ij];
        }
    }
}  //

float Tensor2f::get(int i, int j) { return A2d_[i][j]; }  //

void Tensor2f::gemm(bool transa, bool transb, const SharedTensor2f &a, const SharedTensor2f &b, float alpha,
                    float beta) {
    char ta = transa ? 't' : 'n';
    char tb = transb ? 't' : 'n';
    int m, n, k, nca, ncb, ncc;

    m = dim1_;
    n = dim2_;
    k = transa ? a->dim1_ : a->dim2_;
    nca = transa ? m : k;  // lda
    ncb = transb ? k : n;  // ldb
    ncc = n;               // ldc

    if (m && n && k) {
        C_SGEMM(ta, tb, m, n, k, alpha, &(a->A2d_[0][0]), nca, &(b->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
    }
}  //

void Tensor2f::contract(bool transa, bool transb, int m, int n, int k, const SharedTensor2f &a, const SharedTensor2f &b,
                        float alpha, float beta) {
    char ta = transa ? 't' : 'n';
    char tb = transb ? 't' : 'n';
    int lda, ldb, ldc;

    lda = transa ? m : k;
    ldb = transb ? k : n;
    ldc = n;

    if (m && n && k) {
        C_SGEMM(ta, tb, m, n, k, alpha, &(a->A2d_[0][0]), lda, &(b->A2d_[0][0]), ldb, beta, &(A2d_[0][0]), ldc);
    }
}  //

void Tensor2f::contract(bool transa, bool transb, int m, int n, int k, const SharedTensor2f &a, const SharedTensor2f &b,
                        int start_a, int start_b, float alpha, float beta) {
    char ta = transa ? 't' : 'n';
    char tb = transb ? 't' : 'n';
    int lda, ldb, ldc;

    lda = transa ? m : k;
    ldb = transb ? k : n;
    ldc = n;

    if (m && n && k) {
        C_SGEMM(ta, tb, m, n, k, alpha, a->A2d_[0] + start_a, lda, b->A2d_[0] + start_b, ldb, beta, A2d_[0], ldc);
    }
}  //

void Tensor2f::contract(bool transa, bool transb, int m, int n, int k, const SharedTensor2f &a, const SharedTensor2f &b,
                        int start_a, int start_b, int start_c, float alpha, float beta) {
    char ta = transa ? 't' : 'n';
    char tb = transb ? 't' : 'n';
    int lda, ldb, ldc;

    lda = transa ? m : k;
    ldb = transb ? k : n;
    ldc = n;

    if (m && n && k) {
        C_SGEMM(ta, tb, m, n, k, alpha, a->A2d_[0] + start_a, lda, b->A2d_[0] + start_b, ldb, beta, A2d_[0] + start_c,
                ldc);
    }
}  //

void Tensor2f::contract323(bool transa, bool transb, int m, int n, const SharedTensor2f &a, const SharedTensor2f &b,
                           float alpha, float beta) {
    char ta = transa ? 't' : 'n';
    char tb = transb ? 't' : 'n';
    int k, nca, ncb, ncc;

    k = transb ? b->dim2_ : b->dim1_;
    nca = transa ? m : k;
    ncb = transb ? k : n;
    ncc = n;

    if (m && n && k) {
#pragma omp parallel for
        for (int Q = 0; Q < dim1_; Q++) {
            C_SGEMM(ta, tb, m, n, k, alpha, a->A2d_[Q], nca, b->A2d_[0], ncb, beta, A2d_[Q], ncc);
        }
    }
}  //

void Tensor2f::contract233(bool transa, bool transb, int m, int n, const SharedTensor2f &a, const SharedTensor2f &b,
                           float alpha, float beta) {
    char ta = transa ? 't' : 'n';
    char tb = transb ? 't' : 'n';
    int k, lda, ldb, ldc;

    k = transa ? a->dim1_ : a->dim2_;
    lda = transa ? m : k;
    ldb = transb ? k : n;
    ldc = n;

    if (m && n && k) {
#pragma omp parallel for
        for (int Q = 0; Q < dim1_; Q++) {
            C_SGEMM(ta, tb, m, n, k, alpha, a->A2d_[0], lda, b->A2d_[Q], ldb, beta, A2d_[Q], ldc);
        }
    }
}  //

void Tensor2f::contract332(bool transa, bool transb, int k, const SharedTensor2f &a, const SharedTensor2f &b,
                           float alpha, float beta) {
    char ta = transa ? 't' : 'n';
    char tb = transb ? 't' : 'n';
    int m, n, nca, ncb, ncc;

    m = dim1_;
    n = dim2_;
    nca = transa ? m : k;
    ncb = transb ? k : n;
    ncc = n;

    if (m && n && k) {
        //#pragma omp parallel for
        for (int Q = 0; Q < a->dim1(); Q++) {
            C_SGEMM(ta, tb, m, n, k, alpha, a->A2d_[Q], nca, b->A2d_[Q], ncb, beta, A2d_[0], ncc);
        }
    }
}  //

void Tensor2f::contract424(int target_x, int target_y, const SharedTensor2f &a, const SharedTensor2f &b, float alpha,
                           float beta) {
    char ta;
    char tb;
    int lda, ldb, ldc;
    int m, n, k;

    // C(pq,rs) = \sum_{o} A(oq,rs) B(o,p)
    if (target_x == 1 && target_y == 1) {
        ta = 't';
        tb = 'n';
        m = d1_;
        n = d2_ * d3_ * d4_;
        k = b->dim1();
        lda = m;
        ldb = n;
        ldc = n;

        if (m && n && k) {
            C_SGEMM(ta, tb, m, n, k, alpha, b->A2d_[0], lda, a->A2d_[0], ldb, beta, A2d_[0], ldc);
        }
    }

    // C(pq,rs) = \sum_{o} A(oq,rs) B(p,o)
    else if (target_x == 1 && target_y == 2) {
        ta = 'n';
        tb = 'n';
        m = d1_;
        n = d2_ * d3_ * d4_;
        k = b->dim2();
        lda = k;
        ldb = n;
        ldc = n;

        if (m && n && k) {
            C_SGEMM(ta, tb, m, n, k, alpha, b->A2d_[0], lda, a->A2d_[0], ldb, beta, A2d_[0], ldc);
        }
    }

    // C(pq,rs) = \sum_{o} A(po,rs) B(o,q)
    else if (target_x == 2 && target_y == 1) {
        ta = 'n';
        tb = 'n';
        m = d1_ * d3_ * d4_;
        n = d2_;
        k = b->dim1();
        lda = k;
        ldb = n;
        ldc = n;

        SharedTensor2f temp1 = SharedTensor2f(new Tensor2f("temp1", a->d1_, a->d3_, a->d4_, a->d2_));
        SharedTensor2f temp2 = SharedTensor2f(new Tensor2f("temp2", d1_, d3_, d4_, d2_));
        temp1->sort(1342, a, 1.0, 0.0);

        if (m && n && k) {
            C_SGEMM(ta, tb, m, n, k, alpha, temp1->A2d_[0], lda, b->A2d_[0], ldb, 0.0, temp2->A2d_[0], ldc);
        }
        temp1.reset();
        SharedTensor2f temp3 = SharedTensor2f(new Tensor2f("temp3", d1_, d2_, d3_, d4_));
        temp3->sort(1423, temp2, 1.0, 0.0);
        temp2.reset();
        scale(beta);
        add(temp3);
        temp3.reset();
    }

    // C(pq,rs) = \sum_{o} A(po,rs) B(q,o)
    else if (target_x == 2 && target_y == 2) {
        ta = 'n';
        tb = 't';
        m = d1_ * d3_ * d4_;
        n = d2_;
        k = b->dim2();
        lda = k;
        ldb = k;
        ldc = n;

        SharedTensor2f temp1 = SharedTensor2f(new Tensor2f("temp1", a->d1_, a->d3_, a->d4_, a->d2_));
        SharedTensor2f temp2 = SharedTensor2f(new Tensor2f("temp2", d1_, d3_, d4_, d2_));
        temp1->sort(1342, a, 1.0, 0.0);

        if (m && n && k) {
            C_SGEMM(ta, tb, m, n, k, alpha, temp1->A2d_[0], lda, b->A2d_[0], ldb, 0.0, temp2->A2d_[0], ldc);
        }
        temp1.reset();
        SharedTensor2f temp3 = SharedTensor2f(new Tensor2f("temp3", d1_, d2_, d3_, d4_));
        temp3->sort(1423, temp2, 1.0, 0.0);
        temp2.reset();
        scale(beta);
        add(temp3);
        temp3.reset();
    }

    // C(pq,rs) = \sum_{o} A(pq,os) B(o,r)
    else if (target_x == 3 && target_y == 1) {
        ta = 'n';
        tb = 'n';
        m = d1_ * d2_ * d4_;
        n = d3_;
        k = b->dim1();
        lda = k;
        ldb = n;
        ldc = n;

        SharedTensor2f temp1 = SharedTensor2f(new Tensor2f("temp1", a->d1_, a->d2_, a->d4_, a->d3_));
        SharedTensor2f temp2 = SharedTensor2f(new Tensor2f("temp2", d1_, d2_, d4_, d3_));
        temp1->sort(1243, a, 1.0, 0.0);

        if (m && n && k) {
            C_SGEMM(ta, tb, m, n, k, alpha, temp1->A2d_[0], lda, b->A2d_[0], ldb, 0.0, temp2->A2d_[0], ldc);
        }
        temp1.reset();
        SharedTensor2f temp3 = SharedTensor2f(new Tensor2f("temp3", d1_, d2_, d3_, d4_));
        temp3->sort(1243, temp2, 1.0, 0.0);
        temp2.reset();
        scale(beta);
        add(temp3);
        temp3.reset();
    }

    // C(pq,rs) = \sum_{o} A(pq,os) B(r,o)
    else if (target_x == 3 && target_y == 2) {
        ta = 'n';
        tb = 't';
        m = d1_ * d2_ * d4_;
        n = d3_;
        k = b->dim2();
        lda = k;
        ldb = k;
        ldc = n;

        SharedTensor2f temp1 = SharedTensor2f(new Tensor2f("temp1", a->d1_, a->d2_, a->d4_, a->d3_));
        SharedTensor2f temp2 = SharedTensor2f(new Tensor2f("temp2", d1_, d2_, d4_, d3_));
        temp1->sort(1243, a, 1.0, 0.0);

        if (m && n && k) {
            C_SGEMM(ta, tb, m, n, k, alpha, temp1->A2d_[0], lda, b->A2d_[0], ldb, 0.0, temp2->A2d_[0], ldc);
        }
        temp1.reset();
        SharedTensor2f temp3 = SharedTensor2f(new Tensor2f("temp3", d1_, d2_, d3_, d4_));
        temp3->sort(1243, temp2, 1.0, 0.0);
        temp2.reset();
        scale(beta);
        add(temp3);
        temp3.reset();
    }

    // C(pq,rs) = \sum_{o} A(pq,ro) B(o,s)
    else if (target_x == 4 && target_y == 1) {
        ta = 'n';
        tb = 'n';
        m = d1_ * d2_ * d3_;
        n = d4_;
        k = b->dim1();
        lda = k;
        ldb = n;
        ldc = n;
        if (m && n && k) {
            C_SGEMM(ta, tb, m, n, k, alpha, &(a->A2d_[0][0]), lda, &(b->A2d_[0][0]), ldb, beta, &(A2d_[0][0]), ldc);
        }
    }

    // C(pq,rs) = \sum_{o} A(pq,ro) B(s,o)
    else if (target_x == 4 && target_y == 2) {
        ta = 'n';
        tb = 't';
        m = d1_ * d2_ * d3_;
        n = d4_;
        k = b->dim2();
        lda = k;
        ldb = k;
        ldc = n;
        if (m && n && k) {
            C_SGEMM(ta, tb, m, n, k, alpha, &(a->A2d_[0][0]), lda, &(b->A2d_[0][0]), ldb, beta, &(A2d_[0][0]), ldc);
        }
    }

    /*
    // C(pq,rs) = \sum_{o} A(pq,ro) B(o,s)
    else if (target_x == 4 && target_y == 1) {
    #pragma omp parallel for
    for (int p = 0; p < d1_; p++) {
         for (int q = 0; q < d2_; q++) {
              int pq = row_idx_[p][q];
              for (int r = 0; r < d3_; r++) {
                   for (int s = 0; s < d4_; s++) {
                        int rs = col_idx_[r][s];
                        float sum = 0.0;
                        for (int o = 0; o < b->dim1(); o++) {
                             int ro = a->col_idx_[r][o];
                             sum += a->get(pq,ro) * b->get(o,s);
                        }
                        A2d_[pq][rs] = (alpha * sum) + (beta * A2d_[pq][rs]);
                   }
              }
         }
    }
    }

    // C(pq,rs) = \sum_{o} A(pq,ro) B(s,o)
    else if (target_x == 4 && target_y == 2) {
    #pragma omp parallel for
    for (int p = 0; p < d1_; p++) {
         for (int q = 0; q < d2_; q++) {
              int pq = row_idx_[p][q];
              for (int r = 0; r < d3_; r++) {
                   for (int s = 0; s < d4_; s++) {
                        int rs = col_idx_[r][s];
                        float sum = 0.0;
                        for (int o = 0; o < b->dim2(); o++) {
                             int ro = a->col_idx_[r][o];
                             sum += a->get(pq,ro) * b->get(s,o);
                        }
                        A2d_[pq][rs] = (alpha * sum) + (beta * A2d_[pq][rs]);
                   }
              }
         }
    }
    }

    // C(pq,rs) = \sum_{o} A(oq,rs) B(o,p)
    if (target_x == 1 && target_y == 1) {
    #pragma omp parallel for
    for (int p = 0; p < d1_; p++) {
         for (int q = 0; q < d2_; q++) {
              int pq = row_idx_[p][q];
              for (int r = 0; r < d3_; r++) {
                   for (int s = 0; s < d4_; s++) {
                        int rs = col_idx_[r][s];
                        float sum = 0.0;
                        for (int o = 0; o < b->dim1(); o++) {
                             int oq = a->row_idx_[o][q];
                             sum += a->get(oq,rs) * b->get(o,p);
                        }
                        A2d_[pq][rs] = (alpha * sum) + (beta * A2d_[pq][rs]);
                   }
              }
         }
    }
    }

    // C(pq,rs) = \sum_{o} A(oq,rs) B(p,o)
    else if (target_x == 1 && target_y == 2) {
    #pragma omp parallel for
    for (int p = 0; p < d1_; p++) {
         for (int q = 0; q < d2_; q++) {
              int pq = row_idx_[p][q];
              for (int r = 0; r < d3_; r++) {
                   for (int s = 0; s < d4_; s++) {
                        int rs = col_idx_[r][s];
                        float sum = 0.0;
                        for (int o = 0; o < b->dim2(); o++) {
                             int oq = a->row_idx_[o][q];
                             sum += a->get(oq,rs) * b->get(p,o);
                        }
                        A2d_[pq][rs] = (alpha * sum) + (beta * A2d_[pq][rs]);
                   }
              }
         }
    }
    }

    // C(pq,rs) = \sum_{o} A(po,rs) B(o,q)
    else if (target_x == 2 && target_y == 1) {
    #pragma omp parallel for
    for (int p = 0; p < d1_; p++) {
         for (int q = 0; q < d2_; q++) {
              int pq = row_idx_[p][q];
              for (int r = 0; r < d3_; r++) {
                   for (int s = 0; s < d4_; s++) {
                        int rs = col_idx_[r][s];
                        float sum = 0.0;
                        for (int o = 0; o < b->dim1(); o++) {
                             int po = a->row_idx_[p][o];
                             sum += a->get(po,rs) * b->get(o,q);
                        }
                        A2d_[pq][rs] = (alpha * sum) + (beta * A2d_[pq][rs]);
                   }
              }
         }
    }
    }

    // C(pq,rs) = \sum_{o} A(po,rs) B(q,o)
    else if (target_x == 2 && target_y == 2) {
    #pragma omp parallel for
    for (int p = 0; p < d1_; p++) {
         for (int q = 0; q < d2_; q++) {
              int pq = row_idx_[p][q];
              for (int r = 0; r < d3_; r++) {
                   for (int s = 0; s < d4_; s++) {
                        int rs = col_idx_[r][s];
                        float sum = 0.0;
                        for (int o = 0; o < b->dim2(); o++) {
                             int po = a->row_idx_[p][o];
                             sum += a->get(po,rs) * b->get(q,o);
                        }
                        A2d_[pq][rs] = (alpha * sum) + (beta * A2d_[pq][rs]);
                   }
              }
         }
    }
    }

    // C(pq,rs) = \sum_{o} A(pq,os) B(o,r)
    else if (target_x == 3 && target_y == 1) {
    #pragma omp parallel for
    for (int p = 0; p < d1_; p++) {
         for (int q = 0; q < d2_; q++) {
              int pq = row_idx_[p][q];
              for (int r = 0; r < d3_; r++) {
                   for (int s = 0; s < d4_; s++) {
                        int rs = col_idx_[r][s];
                        float sum = 0.0;
                        for (int o = 0; o < b->dim1(); o++) {
                             int os = a->col_idx_[o][s];
                             sum += a->get(pq,os) * b->get(o,r);
                        }
                        A2d_[pq][rs] = (alpha * sum) + (beta * A2d_[pq][rs]);
                   }
              }
         }
    }
    }

    // C(pq,rs) = \sum_{o} A(pq,os) B(r,o)
    else if (target_x == 3 && target_y == 2) {
    #pragma omp parallel for
    for (int p = 0; p < d1_; p++) {
         for (int q = 0; q < d2_; q++) {
              int pq = row_idx_[p][q];
              for (int r = 0; r < d3_; r++) {
                   for (int s = 0; s < d4_; s++) {
                        int rs = col_idx_[r][s];
                        float sum = 0.0;
                        for (int o = 0; o < b->dim2(); o++) {
                             int os = a->col_idx_[o][s];
                             sum += a->get(pq,os) * b->get(r,o);
                        }
                        A2d_[pq][rs] = (alpha * sum) + (beta * A2d_[pq][rs]);
                   }
              }
         }
    }
    }
    */

    else {
        outfile->Printf("\tcontract424: Unrecognized targets! \n");
    }

}  //

void Tensor2f::contract442(int target_a, int target_b, const SharedTensor2f &a, const SharedTensor2f &b, float alpha,
                           float beta) {
    char ta;
    char tb;
    int lda, ldb, ldc;
    int m, n, k;

    // C(p,q) = \sum_{rst} A(pr,st) B(qr,st)
    if (target_a == 1 && target_b == 1) {
        ta = 'n';
        tb = 't';
        m = dim1_;
        n = dim2_;
        k = a->d2_ * a->d3_ * a->d4_;
        lda = k;
        ldb = k;
        ldc = n;

        if (m && n && k) {
            C_SGEMM(ta, tb, m, n, k, alpha, &(a->A2d_[0][0]), lda, &(b->A2d_[0][0]), ldb, beta, &(A2d_[0][0]), ldc);
        }
    }

    // C(p,q) = \sum_{rst} A(pr,st) B(rq,st)
    else if (target_a == 1 && target_b == 2) {
        ta = 'n';
        tb = 't';
        m = dim1_;
        n = dim2_;
        k = a->d2_ * a->d3_ * a->d4_;
        lda = k;
        ldb = k;
        ldc = n;

        SharedTensor2f temp = SharedTensor2f(new Tensor2f("temp", b->d2_, b->d1_, b->d3_, b->d4_));
        temp->sort(2134, b, 1.0, 0.0);
        if (m && n && k) {
            C_SGEMM(ta, tb, m, n, k, alpha, &(a->A2d_[0][0]), lda, &(temp->A2d_[0][0]), ldb, beta, &(A2d_[0][0]), ldc);
        }
        temp.reset();
    }

    // C(p,q) = \sum_{rst} A(pr,st) B(rs,qt)
    else if (target_a == 1 && target_b == 3) {
        ta = 'n';
        tb = 'n';
        m = dim1_;
        n = dim2_;
        k = a->d2_ * a->d3_ * a->d4_;
        lda = k;
        ldb = n;
        ldc = n;

        SharedTensor2f temp = SharedTensor2f(new Tensor2f("temp", b->d1_, b->d2_, b->d4_, b->d3_));
        temp->sort(1243, b, 1.0, 0.0);
        if (m && n && k) {
            C_SGEMM(ta, tb, m, n, k, alpha, &(a->A2d_[0][0]), lda, &(temp->A2d_[0][0]), ldb, beta, &(A2d_[0][0]), ldc);
        }
        temp.reset();
    }

    // C(p,q) = \sum_{rst} A(pr,st) B(rs,tq)
    else if (target_a == 1 && target_b == 4) {
        ta = 'n';
        tb = 'n';
        m = dim1_;
        n = dim2_;
        k = a->d2_ * a->d3_ * a->d4_;
        lda = k;
        ldb = n;
        ldc = n;
        if (m && n && k) {
            C_SGEMM(ta, tb, m, n, k, alpha, &(a->A2d_[0][0]), lda, &(b->A2d_[0][0]), ldb, beta, &(A2d_[0][0]), ldc);
        }
    }

    // C(p,q) = \sum_{rst} A(rs,tp) B(rs,tq) = \sum_{rst} X(p,rst) B(rst,q)
    else if (target_a == 4 && target_b == 4) {
        ta = 'n';
        tb = 'n';
        m = dim1_;
        n = dim2_;
        k = a->d1_ * a->d2_ * a->d3_;
        lda = k;
        ldb = n;
        ldc = n;

        SharedTensor2f temp = SharedTensor2f(new Tensor2f("temp", a->d4_, a->d1_, a->d2_, a->d3_));
        temp->sort(4123, a, 1.0, 0.0);

        if (m && n && k) {
            C_SGEMM(ta, tb, m, n, k, alpha, &(temp->A2d_[0][0]), lda, &(b->A2d_[0][0]), ldb, beta, &(A2d_[0][0]), ldc);
        }
        temp.reset();
    }

    // C(p,q) = \sum_{rst} A(rs,tp) B(rs,qt)
    else if (target_a == 4 && target_b == 3) {
        ta = 'n';
        tb = 'n';
        m = dim1_;
        n = dim2_;
        k = a->d2_ * a->d3_ * a->d1_;
        lda = k;
        ldb = n;
        ldc = n;

        SharedTensor2f temp1 = SharedTensor2f(new Tensor2f("temp", a->d4_, a->d1_, a->d2_, a->d3_));
        SharedTensor2f temp2 = SharedTensor2f(new Tensor2f("temp", b->d1_, b->d2_, b->d4_, b->d3_));
        temp1->sort(4123, a, 1.0, 0.0);
        temp2->sort(1243, b, 1.0, 0.0);
        if (m && n && k) {
            C_SGEMM(ta, tb, m, n, k, alpha, &(temp1->A2d_[0][0]), lda, &(temp2->A2d_[0][0]), ldb, beta, &(A2d_[0][0]),
                    ldc);
        }
        temp1.reset();
        temp2.reset();
    }

    // C(p,q) = \sum_{rst} A(rp,st) B(rq,st)
    else if (target_a == 2 && target_b == 2) {
        ta = 'n';
        tb = 't';
        m = dim1_;
        n = dim2_;
        k = a->d1_ * a->d3_ * a->d4_;
        lda = k;
        ldb = k;
        ldc = n;

        SharedTensor2f temp1 = SharedTensor2f(new Tensor2f("temp1", a->d2_, a->d1_, a->d3_, a->d4_));
        SharedTensor2f temp2 = SharedTensor2f(new Tensor2f("temp2", b->d2_, b->d1_, b->d3_, b->d4_));
        temp1->sort(2134, a, 1.0, 0.0);
        temp2->sort(2134, b, 1.0, 0.0);
        if (m && n && k) {
            C_SGEMM(ta, tb, m, n, k, alpha, &(temp1->A2d_[0][0]), lda, &(temp2->A2d_[0][0]), ldb, beta, &(A2d_[0][0]),
                    ldc);
        }
        temp1.reset();
        temp2.reset();
    }

    // C(p,q) = \sum_{rst} A(rs,pt) B(rs,qt)
    else if (target_a == 3 && target_b == 3) {
        ta = 't';
        tb = 'n';
        m = dim1_;
        n = dim2_;
        k = a->d1_ * a->d2_ * a->d4_;
        lda = m;
        ldb = n;
        ldc = n;

        SharedTensor2f temp1 = SharedTensor2f(new Tensor2f("temp1", a->d1_, a->d2_, a->d4_, a->d3_));
        SharedTensor2f temp2 = SharedTensor2f(new Tensor2f("temp2", b->d1_, b->d2_, b->d4_, b->d3_));
        temp1->sort(1243, a, 1.0, 0.0);
        temp2->sort(1243, b, 1.0, 0.0);
        if (m && n && k) {
            C_SGEMM(ta, tb, m, n, k, alpha, &(temp1->A2d_[0][0]), lda, &(temp2->A2d_[0][0]), ldb, beta, &(A2d_[0][0]),
                    ldc);
        }
        temp1.reset();
        temp2.reset();
    }

    else {
        outfile->Printf("contract442: Unrecognized targets!");
    }

}  //

void Tensor2f::gemv(bool transa, const SharedTensor2f &a, const SharedTensor1f &b, float alpha, float beta) {
    char ta = transa ? 't' : 'n';
    int m, n, k, incx, incy, lda;

    m = a->dim1();
    n = a->dim2();
    lda = n;
    incx = 1;  // increments in elements of b vector
    incy = 1;  // increments in elements of A2d_

    if (m && n) {
        C_SGEMV(ta, m, n, alpha, &(a->A2d_[0][0]), lda, b->A1d_, incx, beta, &(A2d_[0][0]), incy);
    }
}  //

void Tensor2f::gemv(bool transa, const SharedTensor2f &a, const SharedTensor2f &b, float alpha, float beta) {
    char ta = transa ? 't' : 'n';
    int m, n, k, incx, incy, lda;

    m = a->dim1();
    n = a->dim2();
    lda = n;
    incx = 1;  // increments in elements of b vector
    incy = 1;  // increments in elements of A2d_

    if (m && n) {
        C_SGEMV(ta, m, n, alpha, &(a->A2d_[0][0]), lda, &(b->A2d_[0][0]), incx, beta, &(A2d_[0][0]), incy);
    }
}  //

// void Tensor2f::davidson(int n_eigval, const SharedTensor2f &eigvectors, const SharedTensor1f &eigvalues, float cutoff,
//                         int print) {
//     david(A2d_, dim1_, n_eigval, eigvalues->A1d_, eigvectors->A2d_, cutoff, print);

// }  //

void Tensor2f::add(const SharedTensor2f &a) {
    size_t length = (size_t)dim1_ * (size_t)dim2_;
    C_SAXPY(length, 1.0, a->A2d_[0], 1, A2d_[0], 1);
}  //

void Tensor2f::add(float **a) {
    size_t length = (size_t)dim1_ * (size_t)dim2_;
    C_SAXPY(length, 1.0, a[0], 1, A2d_[0], 1);
}  //

void Tensor2f::add(float alpha, const SharedTensor2f &Adum) {
    SharedTensor2f temp = SharedTensor2f(new Tensor2f(Adum->dim1_, Adum->dim2_));
    temp->copy(Adum);
    temp->scale(alpha);
    add(temp);
}  //

void Tensor2f::add(int i, int j, float value) { A2d_[i][j] += value; }  //

void Tensor2f::subtract(const SharedTensor2f &a) {
    size_t length = (size_t)dim1_ * (size_t)dim2_;
    C_SAXPY(length, -1.0, a->A2d_[0], 1, A2d_[0], 1);
}  //

void Tensor2f::subtract(int i, int j, float value) { A2d_[i][j] -= value; }  //

void Tensor2f::axpy(float **a, float alpha) {
    size_t length = (size_t)dim1_ * (size_t)dim2_;
    C_SAXPY(length, alpha, a[0], 1, A2d_[0], 1);
}  //

void Tensor2f::axpy(const SharedTensor2f &a, float alpha) {
    size_t length = (size_t)dim1_ * (size_t)dim2_;
    C_SAXPY(length, alpha, a->A2d_[0], 1, A2d_[0], 1);
}  //

void Tensor2f::axpy(size_t length, int inc_a, const SharedTensor2f &a, int inc_2d, float alpha) {
    C_SAXPY(length, alpha, a->A2d_[0], inc_a, A2d_[0], inc_2d);
}  //

void Tensor2f::axpy(size_t length, int start_a, int inc_a, const SharedTensor2f &A, int start_2d, int inc_2d,
                    float alpha) {
    C_SAXPY(length, alpha, A->A2d_[0] + start_a, inc_a, A2d_[0] + start_2d, inc_2d);
}  //

float Tensor2f::norm() {
    float value = 0.0;
    size_t length = (size_t)dim1_ * (size_t)dim2_;
    value = C_SNRM2(length, A2d_[0], 1);
    return value;
}  //

float **Tensor2f::transpose2() {
    float **temp = block_matrix_float   (dim2_, dim1_);
    memset(temp[0], 0, sizeof(float) * dim1_ * dim2_);
#pragma omp parallel for
    for (int i = 0; i < dim2_; ++i) {
        for (int j = 0; j < dim1_; ++j) {
            temp[i][j] = A2d_[j][i];
        }
    }

    return temp;
}  //

SharedTensor2f Tensor2f::transpose() {
    SharedTensor2f temp = SharedTensor2f(new Tensor2f(dim2_, dim1_));
#pragma omp parallel for
    for (int i = 0; i < dim2_; ++i) {
        for (int j = 0; j < dim1_; ++j) {
            temp->A2d_[i][j] = A2d_[j][i];
        }
    }

    return temp;
}  //

void Tensor2f::trans(const SharedTensor2f &A) {
#pragma omp parallel for
    for (int i = 0; i < dim1_; ++i) {
        for (int j = 0; j < dim2_; ++j) {
            A2d_[i][j] = A->A2d_[j][i];
        }
    }

}  //

void Tensor2f::trans(float **A) {
#pragma omp parallel for
    for (int i = 0; i < dim1_; ++i) {
        for (int j = 0; j < dim2_; ++j) {
            A2d_[i][j] = A[j][i];
        }
    }

}  //

void Tensor2f::copy(float **a) {
    // size_t size = dim1_ * dim2_ * sizeof(float);
    // if (size) memcpy(&(A2d_[0][0]), &(a[0][0]), size);
    size_t length;
    length = (size_t)dim1_ * (size_t)dim2_;
    C_SCOPY(length, a[0], 1, A2d_[0], 1);
}

void Tensor2f::copy(const SharedTensor2f &Adum) {
    // Make sure that matrices are in the same size
    bool same = true;
    if (dim2_ != Adum->dim2_ || dim1_ != Adum->dim1_) same = false;

    if (same == false) {
        release();
        dim1_ = Adum->dim1_;
        dim2_ = Adum->dim2_;
        memalloc();
    }

    // If matrices are in the same size
    size_t length;
    length = (size_t)dim1_ * (size_t)dim2_;
    if (dim1_ != 0 && dim2_ != 0) {
        // memcpy(A2d_[0], Adum->A2d_[0], dim1_ * dim2_ * sizeof(float));
        C_SCOPY(length, Adum->A2d_[0], 1, A2d_[0], 1);
    }
}  //

void Tensor2f::copy(size_t length, const SharedTensor2f &A, int inc_a, int inc_2d) {
    C_SCOPY(length, A->A2d_[0], inc_a, A2d_[0], inc_2d);
}  //

void Tensor2f::copy(const SharedTensor2f &A, int start) {
    memcpy(A2d_[0], A->A2d_[0] + start, dim1_ * dim2_ * sizeof(float));
}  //

void Tensor2f::pcopy(const SharedTensor2f &A, int dim_copy, int dim_skip) {
    float *temp = new float[dim_copy];
    int syc = 0;
    // A[m] is getting the pointer to the m-th row of A.
    // A[0]+m is getting the pointer to the m-th element of A.
    for (int i = 0; i < dim1_ * dim2_; i += dim_copy) {
        memcpy(temp, A->A2d_[0] + syc, dim_copy * sizeof(float));
        memcpy(A2d_[0] + i, temp, dim_copy * sizeof(float));
        syc += dim_copy + dim_skip;
    }
    delete[] temp;

}  //

void Tensor2f::pcopy(const SharedTensor2f &A, int dim_copy, int dim_skip, int start) {
    float *temp = new float[dim_copy];
    int syc = 0;
    // A[m] is getting the pointer to the m-th row of A.
    // A[0]+m is getting the pointer to the m-th element of A.
    for (int i = 0; i < dim1_ * dim2_; i += dim_copy) {
        memcpy(temp, A->A2d_[0] + start + syc, dim_copy * sizeof(float));
        memcpy(A2d_[0] + i, temp, dim_copy * sizeof(float));
        syc += dim_copy + dim_skip;
    }
    delete[] temp;

}  //


void Tensor2f::double2float(const SharedTensor2d &D ) {
    #pragma omp parallel for
    for (int i = 0; i < dim1_; ++i) {
        for (int j = 0; j < dim2_; ++j) {
            A2d_[i][j] = static_cast<float>(D->get(i,j));
            // A2d_[i][j] = static_cast<float>(D->A2d_[j][i]);
        }
    }

    // size_t length;
    // length = (size_t)dim1_ * (size_t)dim2_;
    // for (int i = 0; i < length ; i += 1) {
    //      A2d_[i]=static_cast<float>(D->get(i,A2d_[i]);
        // A2d_[i]=(float)(D->A2d_[i]);
        // }
}  //



// void Tensor2f::diagonalize(const SharedTensor2f &eigvectors, const SharedTensor1f &eigvalues, float cutoff) {
//     sq_rsp(dim1_, dim2_, A2d_, eigvalues->A1d_, 1, eigvectors->A2d_, cutoff);

// }  //

// void Tensor2f::cdsyev(char jobz, char uplo, const SharedTensor2f &eigvectors, const SharedTensor1f &eigvalues) {
//     if (dim1_) {
//         int lwork = 3 * dim2_;
//         float **work = block_matrix_float   (dim1_, lwork);
//         memset(work[0], 0.0, sizeof(float) * dim1_ * lwork);
//         C_DSYEV(jobz, uplo, dim1_, &(A2d_[0][0]), dim2_, eigvalues->A1d_, &(work[0][0]), lwork);
//         free_block_float(work);
//     }
// }  //

// void Tensor2f::cdgesv(const SharedTensor1f &Xvec) {
//     if (dim1_) {
//         int errcod;
//         int *ipiv = new int[dim1_];
//         memset(ipiv, 0, sizeof(int) * dim1_);
//         errcod = 0;
//         errcod = C_SGESV(dim1_, 1, &(A2d_[0][0]), dim2_, &(ipiv[0]), Xvec->A1d_, dim2_);
//         delete[] ipiv;
//     }
// }  //

// void Tensor2f::cdgesv(const SharedTensor1f &Xvec, int errcod) {
//     if (dim1_) {
//         int *ipiv = new int[dim1_];
//         memset(ipiv, 0, sizeof(int) * dim1_);
//         errcod = 0;
//         errcod = C_SGESV(dim1_, 1, &(A2d_[0][0]), dim2_, &(ipiv[0]), Xvec->A1d_, dim2_);
//         delete[] ipiv;
//     }
// }  //

// void Tensor2f::cdgesv(float *Xvec) {
//     if (dim1_) {
//         int errcod;
//         int *ipiv = new int[dim1_];
//         memset(ipiv, 0, sizeof(int) * dim1_);
//         errcod = 0;
//         errcod = C_SGESV(dim1_, 1, &(A2d_[0][0]), dim2_, &(ipiv[0]), Xvec, dim2_);
//         delete[] ipiv;
//     }
// }  //

// void Tensor2f::cdgesv(float *Xvec, int errcod) {
//     if (dim1_) {
//         int *ipiv = new int[dim1_];
//         memset(ipiv, 0, sizeof(int) * dim1_);
//         errcod = 0;
//         errcod = C_SGESV(dim1_, 1, &(A2d_[0][0]), dim2_, &(ipiv[0]), Xvec, dim2_);
//         delete[] ipiv;
//     }
// }  //

// void Tensor2f::lineq_flin(const SharedTensor1f &Xvec, float *det) {
//     if (dim1_) {
//         flin(A2d_, Xvec->A1d_, dim1_, 1, det);
//     }
// }  //

// void Tensor2f::lineq_flin(float *Xvec, float *det) {
//     if (dim1_) {
//         flin(A2d_, Xvec, dim1_, 1, det);
//     }
// }  //

// void Tensor2f::lineq_pople(const SharedTensor1f &Xvec, int num_vecs, float cutoff) {
//     if (dim1_) {
//         pople(A2d_, Xvec->A1d_, dim1_, num_vecs, cutoff, "outfile", 0);
//     }
// }  //

// void Tensor2f::lineq_pople(float *Xvec, int num_vecs, float cutoff) {
//     if (dim1_) {
//         pople(A2d_, Xvec, dim1_, num_vecs, cutoff, "outfile", 0);
//     }
// }  //

void Tensor2f::level_shift(float value) {
#pragma omp parallel for
    for (int i = 0; i < dim1_; ++i) {
        subtract(i, i, value);
    }

}  //

void Tensor2f::outer_product(const SharedTensor1f &x, const SharedTensor1f &y) {
#pragma omp parallel for
    for (int i = 0; i < x->dim1_; i++) {
        for (int j = 0; j < y->dim1_; j++) {
            A2d_[i][j] = x->A1d_[i] * y->A1d_[j];
        }
    }

}  //

// TODO:
// DGER compute the rank-one update of a general matrix: A <-- A + alpha * x * yT
// dger(m, n, alpha, x, incx, y, incy, a, lda);

void Tensor2f::scale(float a) {
    // size_t size;
    size_t size;
    size = (size_t)dim1_ * (size_t)dim2_;
    if (size) C_SSCAL(size, a, &(A2d_[0][0]), 1);
}  //

void Tensor2f::scale_row(int m, float a) { C_SSCAL((size_t)dim1_, a, &(A2d_[m][0]), 1); }  //

void Tensor2f::scale_column(int n, float a) { C_SSCAL((size_t)dim2_, a, &(A2d_[0][n]), dim1_); }  //

void Tensor2f::identity() {
    zero();
    for (int i = 0; i < dim1_; ++i) A2d_[i][i] = 1.0;
}  //

float Tensor2f::trace() {
    float value = 0.0;
    for (int i = 0; i < dim1_; ++i) value += A2d_[i][i];
    return value;
}  //

void Tensor2f::transform(const SharedTensor2f &a, const SharedTensor2f &transformer) {
    SharedTensor2f temp = SharedTensor2f(new Tensor2f(a->dim1_, transformer->dim2_));
    temp->gemm(false, false, a, transformer, 1.0, 0.0);
    gemm(true, false, transformer, temp, 1.0, 0.0);
}  //

void Tensor2f::back_transform(const SharedTensor2f &a, const SharedTensor2f &transformer) {
    SharedTensor2f temp = SharedTensor2f(new Tensor2f(a->dim1_, transformer->dim2_));
    temp->gemm(false, true, a, transformer, 1.0, 0.0);
    gemm(false, false, transformer, temp, 1.0, 0.0);
}  //

void Tensor2f::back_transform(const SharedTensor2f &a, const SharedTensor2f &transformer, float alpha, float beta) {
    SharedTensor2f temp = SharedTensor2f(new Tensor2f(a->dim1_, transformer->dim2_));
    temp->gemm(false, true, a, transformer, 1.0, 0.0);
    gemm(false, false, transformer, temp, alpha, beta);
}  //

void Tensor2f::pseudo_transform(const SharedTensor2f &a, const SharedTensor2f &transformer) {
    SharedTensor2f temp = SharedTensor2f(new Tensor2f(a->dim1_, transformer->dim2_));
    temp->gemm(false, false, a, transformer, 1.0, 0.0);
    gemm(false, false, transformer, temp, 1.0, 0.0);
}  //

void Tensor2f::triple_gemm(const SharedTensor2f &a, const SharedTensor2f &b, const SharedTensor2f &c) {
    if (a->dim2_ == b->dim1_ && b->dim2_ == c->dim1_ && a->dim1_ == dim1_ && c->dim2_ == dim2_) {
        SharedTensor2f bc = SharedTensor2f(new Tensor2f(b->dim1_, c->dim2_));
        bc->gemm(false, false, b, c, 1.0, 0.0);
        gemm(false, false, a, bc, 1.0, 0.0);
    } else {
        outfile->Printf("\n Warning!!! Matrix dimensions do NOT match in triple_gemm().\n");
    }

}  //

float Tensor2f::vector_dot(float **rhs) {
    float value = 0.0;
    // size_t size = dim1_ * dim2_;
    size_t size;
    size = (size_t)dim1_ * (size_t)dim2_;
    if (size) value += C_SDOT(size, (&A2d_[0][0]), 1, &(rhs[0][0]), 1);
    return value;
}  //

float Tensor2f::vector_dot(const SharedTensor2f &rhs) {
    float value = 0.0;
    // size_t size = dim1_ * dim2_;
    size_t size;
    size = (size_t)dim1_ * (size_t)dim2_;
    if (size) value += C_SDOT(size, (&A2d_[0][0]), 1, &(rhs->A2d_[0][0]), 1);
    return value;
}  //

void Tensor2f::write(std::shared_ptr<psi::PSIO> psio, size_t fileno) {
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno))
        already_open = true;
    else
        psio->open(fileno, PSIO_OPEN_OLD);
    psio->write_entry(fileno, const_cast<char *>(name_.c_str()), (char *)A2d_[0], sizeof(float) * dim1_ * dim2_);
    if (!already_open) psio->close(fileno, 1);  // Close and keep
}  //

void Tensor2f::write(std::shared_ptr<psi::PSIO> psio, size_t fileno, psio_address start, psio_address *end) {
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno))
        already_open = true;
    else
        psio->open(fileno, PSIO_OPEN_OLD);
    size_t size_ = (size_t)dim1_ * dim2_ * sizeof(float);
    psio->write(fileno, const_cast<char *>(name_.c_str()), (char *)A2d_[0], size_, start, end);
    if (!already_open) psio->close(fileno, 1);  // Close and keep
}  //

void Tensor2f::write(psi::PSIO *const psio, size_t fileno) {
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno))
        already_open = true;
    else
        psio->open(fileno, PSIO_OPEN_OLD);
    psio->write_entry(fileno, const_cast<char *>(name_.c_str()), (char *)A2d_[0], sizeof(float) * dim1_ * dim2_);
    if (!already_open) psio->close(fileno, 1);  // Close and keep
}  //

void Tensor2f::write(psi::PSIO *const psio, size_t fileno, psio_address start, psio_address *end) {
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno))
        already_open = true;
    else
        psio->open(fileno, PSIO_OPEN_OLD);
    size_t size_ = (size_t)dim1_ * dim2_ * sizeof(float);
    psio->write(fileno, const_cast<char *>(name_.c_str()), (char *)A2d_[0], size_, start, end);
    if (!already_open) psio->close(fileno, 1);  // Close and keep
}  //

void Tensor2f::write(psi::PSIO &psio, size_t fileno) { write(&psio, fileno); }  //

void Tensor2f::write(psi::PSIO &psio, size_t fileno, psio_address start, psio_address *end) {
    write(&psio, fileno, start, end);
}  //

void Tensor2f::write(std::shared_ptr<psi::PSIO> psio, const std::string &filename, size_t fileno) {
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno))
        already_open = true;
    else
        psio->open(fileno, PSIO_OPEN_OLD);
    psio->write_entry(fileno, const_cast<char *>(filename.c_str()), (char *)A2d_[0], sizeof(float) * dim1_ * dim2_);
    if (!already_open) psio->close(fileno, 1);  // Close and keep
}  //

void Tensor2f::write(std::shared_ptr<psi::PSIO> psio, size_t fileno, bool three_index, bool symm) {
    // Form Lower triangular part
    if (three_index && symm) {
        int ntri_col = 0.5 * d2_ * (d2_ + 1);
        SharedTensor2f temp = SharedTensor2f(new Tensor2f("temp", d1_, ntri_col));
#pragma omp parallel for
        for (int R = 0; R < d1_; R++) {
            for (int p = 0; p < d2_; p++) {
                for (int q = 0; q < d3_; q++) {
                    int pq = col_idx_[p][q];
                    int pq_sym = index2(p, q);
                    temp->set(R, pq_sym, A2d_[R][pq]);
                }
            }
        }

        // Check to see if the file is open
        bool already_open = false;
        if (psio->open_check(fileno))
            already_open = true;
        else
            psio->open(fileno, PSIO_OPEN_OLD);
        psio->write_entry(fileno, const_cast<char *>(name_.c_str()), (char *)temp->A2d_[0],
                          sizeof(float) * dim1_ * ntri_col);
        if (!already_open) psio->close(fileno, 1);  // Close and keep
        temp.reset();
    }

    else {
        // Check to see if the file is open
        bool already_open = false;
        if (psio->open_check(fileno))
            already_open = true;
        else
            psio->open(fileno, PSIO_OPEN_OLD);
        psio->write_entry(fileno, const_cast<char *>(name_.c_str()), (char *)A2d_[0], sizeof(float) * dim1_ * dim2_);
        if (!already_open) psio->close(fileno, 1);  // Close and keep
    }

}  //

void Tensor2f::write(std::shared_ptr<psi::PSIO> psio, const std::string &filename, size_t fileno, bool three_index,
                     bool symm) {
    // Form Lower triangular part
    if (three_index && symm) {
        int ntri_col = 0.5 * d2_ * (d2_ + 1);
        SharedTensor2f temp = SharedTensor2f(new Tensor2f("temp", d1_, ntri_col));
#pragma omp parallel for
        for (int R = 0; R < d1_; R++) {
            for (int p = 0; p < d2_; p++) {
                for (int q = 0; q < d3_; q++) {
                    int pq = col_idx_[p][q];
                    int pq_sym = index2(p, q);
                    temp->set(R, pq_sym, A2d_[R][pq]);
                }
            }
        }

        // Check to see if the file is open
        bool already_open = false;
        if (psio->open_check(fileno))
            already_open = true;
        else
            psio->open(fileno, PSIO_OPEN_OLD);
        psio->write_entry(fileno, const_cast<char *>(filename.c_str()), (char *)temp->A2d_[0],
                          sizeof(float) * dim1_ * ntri_col);
        if (!already_open) psio->close(fileno, 1);  // Close and keep
        temp.reset();
    }

    else {
        // Check to see if the file is open
        bool already_open = false;
        if (psio->open_check(fileno))
            already_open = true;
        else
            psio->open(fileno, PSIO_OPEN_OLD);
        psio->write_entry(fileno, const_cast<char *>(filename.c_str()), (char *)A2d_[0],
                          sizeof(float) * dim1_ * dim2_);
        if (!already_open) psio->close(fileno, 1);  // Close and keep
    }

}  //

void Tensor2f::write_symm(std::shared_ptr<psi::PSIO> psio, size_t fileno) {
    // Form Lower triangular part
    int ntri_col = 0.5 * dim1_ * (dim1_ + 1);
    SharedTensor1f temp = SharedTensor1f(new Tensor1f("temp", ntri_col));
#pragma omp parallel for
    for (int p = 0; p < dim1_; p++) {
        for (int q = 0; q <= p; q++) {
            int pq = index2(p, q);
            temp->set(pq, A2d_[p][q]);
        }
    }

    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno))
        already_open = true;
    else
        psio->open(fileno, PSIO_OPEN_OLD);
    psio->write_entry(fileno, const_cast<char *>(name_.c_str()), (char *)&(temp->A1d_[0]), sizeof(float) * ntri_col);
    if (!already_open) psio->close(fileno, 1);  // Close and keep
    temp.reset();

}  //

void Tensor2f::write_anti_symm(std::shared_ptr<psi::PSIO> psio, size_t fileno) {
    // Form Lower triangular part
    int ntri_row, ntri_col;
    if (dim1_ > 1) {
        ntri_row = 0.5 * d1_ * (d1_ - 1);
    } else if (dim1_ == 1) {
        ntri_row = 1;
    }
    if (dim2_ > 1) {
        ntri_col = 0.5 * d3_ * (d3_ - 1);
    } else if (dim2_ == 1) {
        ntri_col = 1;
    }
    SharedTensor2f temp = SharedTensor2f(new Tensor2f("temp", ntri_row, ntri_col));
#pragma omp parallel for
    for (int p = 1; p < d1_; p++) {
        for (int q = 0; q < p; q++) {
            int pq = row_idx_[p][q];
            int pq2 = idx_asym(p, q);
            for (int r = 1; r < d3_; r++) {
                for (int s = 0; s < r; s++) {
                    int rs = col_idx_[r][s];
                    int rs2 = idx_asym(r, s);
                    temp->set(pq2, rs2, A2d_[pq][rs]);
                }
            }
        }
    }

    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno))
        already_open = true;
    else
        psio->open(fileno, PSIO_OPEN_OLD);
    psio->write_entry(fileno, const_cast<char *>(name_.c_str()), (char *)temp->A2d_[0],
                      sizeof(float) * ntri_row * ntri_col);
    if (!already_open) psio->close(fileno, 1);  // Close and keep
    temp.reset();

}  //

void Tensor2f::read(psi::PSIO *psio, size_t fileno) {
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno))
        already_open = true;
    else
        psio->open(fileno, PSIO_OPEN_OLD);
    psio->read_entry(fileno, const_cast<char *>(name_.c_str()), (char *)A2d_[0], sizeof(float) * dim1_ * dim2_);
    if (!already_open) psio->close(fileno, 1);  // Close and keep
}

void Tensor2f::read(psi::PSIO *psio, size_t fileno, psio_address start, psio_address *end) {
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno))
        already_open = true;
    else
        psio->open(fileno, PSIO_OPEN_OLD);
    size_t size_ = (size_t)dim1_ * dim2_ * sizeof(float);
    psio->read(fileno, const_cast<char *>(name_.c_str()), (char *)A2d_[0], size_, start, end);
    if (!already_open) psio->close(fileno, 1);  // Close and keep
}

void Tensor2f::read(std::shared_ptr<psi::PSIO> psio, size_t fileno) {
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno))
        already_open = true;
    else
        psio->open(fileno, PSIO_OPEN_OLD);
    psio->read_entry(fileno, const_cast<char *>(name_.c_str()), (char *)A2d_[0], sizeof(float) * dim1_ * dim2_);
    if (!already_open) psio->close(fileno, 1);  // Close and keep
}

void Tensor2f::read(std::shared_ptr<psi::PSIO> psio, size_t fileno, psio_address start, psio_address *end) {
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno))
        already_open = true;
    else
        psio->open(fileno, PSIO_OPEN_OLD);
    size_t size_ = (size_t)dim1_ * dim2_ * sizeof(float);
    psio->read(fileno, const_cast<char *>(name_.c_str()), (char *)A2d_[0], size_, start, end);
    if (!already_open) psio->close(fileno, 1);  // Close and keep
}

void Tensor2f::read(psi::PSIO &psio, size_t fileno) { read(&psio, fileno); }  //

void Tensor2f::read(psi::PSIO &psio, size_t fileno, psio_address start, psio_address *end) {
    read(&psio, fileno, start, end);
}  //

void Tensor2f::read(std::shared_ptr<psi::PSIO> psio, size_t fileno, bool three_index, bool symm) {
    // Form Lower triangular part
    if (three_index && symm) {
        int ntri_col = 0.5 * d2_ * (d2_ + 1);
        SharedTensor2f temp = SharedTensor2f(new Tensor2f("temp", d1_, ntri_col));

        // Check to see if the file is open
        bool already_open = false;
        if (psio->open_check(fileno))
            already_open = true;
        else
            psio->open(fileno, PSIO_OPEN_OLD);
        psio->read_entry(fileno, const_cast<char *>(name_.c_str()), (char *)temp->A2d_[0],
                         sizeof(float) * dim1_ * ntri_col);
        if (!already_open) psio->close(fileno, 1);  // Close and keep

#pragma omp parallel for
        for (int R = 0; R < d1_; R++) {
            for (int p = 0; p < d2_; p++) {
                for (int q = 0; q < d3_; q++) {
                    int pq = col_idx_[p][q];
                    int pq_sym = index2(p, q);
                    A2d_[R][pq] = temp->get(R, pq_sym);
                }
            }
        }
        temp.reset();
    }

    else {
        // Check to see if the file is open
        bool already_open = false;
        if (psio->open_check(fileno))
            already_open = true;
        else
            psio->open(fileno, PSIO_OPEN_OLD);
        psio->read_entry(fileno, const_cast<char *>(name_.c_str()), (char *)A2d_[0], sizeof(float) * dim1_ * dim2_);
        if (!already_open) psio->close(fileno, 1);  // Close and keep
    }

}  //

void Tensor2f::read_symm(std::shared_ptr<psi::PSIO> psio, size_t fileno) {
    // Form Lower triangular part
    int ntri_col = 0.5 * dim1_ * (dim1_ + 1);
    SharedTensor1f temp = SharedTensor1f(new Tensor1f("temp", ntri_col));

    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno))
        already_open = true;
    else
        psio->open(fileno, PSIO_OPEN_OLD);
    psio->read_entry(fileno, const_cast<char *>(name_.c_str()), (char *)&(temp->A1d_[0]), sizeof(float) * ntri_col);
    if (!already_open) psio->close(fileno, 1);  // Close and keep

#pragma omp parallel for
    for (int p = 0; p < dim1_; p++) {
        for (int q = 0; q <= p; q++) {
            int pq = index2(p, q);
            A2d_[p][q] = temp->get(pq);
            A2d_[q][p] = temp->get(pq);
        }
    }
    temp.reset();
}  //

void Tensor2f::read_anti_symm(std::shared_ptr<psi::PSIO> psio, size_t fileno) {
    // Form Lower triangular part
    int ntri_row, ntri_col;
    if (dim1_ > 1) {
        ntri_row = 0.5 * d1_ * (d1_ - 1);
    } else if (dim1_ == 1) {
        ntri_row = 1;
    }
    if (dim2_ > 1) {
        ntri_col = 0.5 * d3_ * (d3_ - 1);
    } else if (dim2_ == 1) {
        ntri_col = 1;
    }

    SharedTensor2f temp = SharedTensor2f(new Tensor2f("temp", ntri_row, ntri_col));

    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno))
        already_open = true;
    else
        psio->open(fileno, PSIO_OPEN_OLD);
    psio->read_entry(fileno, const_cast<char *>(name_.c_str()), (char *)temp->A2d_[0],
                     sizeof(float) * ntri_row * ntri_col);
    if (!already_open) psio->close(fileno, 1);  // Close and keep

#pragma omp parallel for
    for (int p = 1; p < d1_; p++) {
        for (int q = 0; q < p; q++) {
            int pq = row_idx_[p][q];
            int qp = row_idx_[q][p];
            int pq2 = idx_asym(p, q);
            for (int r = 1; r < d3_; r++) {
                for (int s = 0; s < r; s++) {
                    int rs = col_idx_[r][s];
                    int sr = col_idx_[s][r];
                    int rs2 = idx_asym(r, s);
                    float value = temp->get(pq2, rs2);
                    A2d_[pq][rs] = value;
                    A2d_[pq][sr] = -1.0 * value;
                    A2d_[qp][rs] = -1.0 * value;
                    A2d_[qp][sr] = value;
                }
            }
        }
    }
    temp.reset();

}  //

// bool Tensor2f::read(PSIO *psio, int itap, const char *label, int dim) {
//     int ntri = 0.5 * dim * (dim + 1);
//     float *mybuffer = init_array(ntri);
//     memset(mybuffer, 0, sizeof(float) * ntri);
//     IWL::read_one(psio, itap, label, mybuffer, ntri, 0, 0, "outfile");

//     float **Asq = block_matrix_float    (dim, dim);
//     memset(Asq[0], 0, sizeof(float) * dim * dim);
//     tri_to_sq(mybuffer, Asq, dim);
//     free(mybuffer);

//     set(Asq);
//     free_block_float(Asq);
//     return true;
// }  //

// bool Tensor2f::read(std::shared_ptr<psi::PSIO> psio, int itap, const char *label, int dim) {
//     int ntri = 0.5 * dim * (dim + 1);
//     float *mybuffer = init_array(ntri);
//     memset(mybuffer, 0, sizeof(float) * ntri);
//     IWL::read_one(psio.get(), itap, label, mybuffer, ntri, 0, 0, "outfile");

//     float **Asq = block_matrix_float    (dim, dim);
//     memset(Asq[0], 0, sizeof(float) * dim * dim);
//     tri_to_sq(mybuffer, Asq, dim);
//     free(mybuffer);

//     set(Asq);
//     free_block_float(Asq);
//     return true;
// }  //

void Tensor2f::save(std::shared_ptr<psi::PSIO> psio, size_t fileno) {
    write(psio, fileno);
    release();
}  //

void Tensor2f::save(psi::PSIO *const psio, size_t fileno) {
    write(psio, fileno);
    release();
}  //

void Tensor2f::save(psi::PSIO &psio, size_t fileno) {
    write(&psio, fileno);
    release();
}  //

void Tensor2f::load(std::shared_ptr<psi::PSIO> psio, size_t fileno, std::string name, int d1, int d2) {
    init(name, d1, d2);
    read(psio, fileno);
}  //

void Tensor2f::load(psi::PSIO *const psio, size_t fileno, std::string name, int d1, int d2) {
    init(name, d1, d2);
    read(psio, fileno);
}  //

void Tensor2f::load(psi::PSIO &psio, size_t fileno, std::string name, int d1, int d2) {
    init(name, d1, d2);
    read(&psio, fileno);
}  //

void Tensor2f::mywrite(const std::string &filename) {
    // write binary data
    std::ofstream OutFile;
    OutFile.open(const_cast<char *>(filename.c_str()), std::ios::out | std::ios::binary);
    OutFile.write((char *)A2d_[0], dim1_ * dim2_ * sizeof(float));
    OutFile.close();
}  //

void Tensor2f::mywrite(int fileno) {
    std::ostringstream convert;
    convert << fileno;
    std::string scr = PSIOManager::shared_object()->get_default_path();
    std::string pid_ = psio_getpid();
    // std::string fname = scr + "psi_dfocc." + convert.str();
    std::string fname = scr + "psi." + pid_ + "." + convert.str();

    // write binary data
    std::ofstream OutFile;
    OutFile.open(const_cast<char *>(fname.c_str()), std::ios::out | std::ios::binary);
    OutFile.write((char *)A2d_[0], dim1_ * dim2_ * sizeof(float));
    OutFile.close();
}  //

void Tensor2f::mywrite(int fileno, bool append) {
    std::ostringstream convert;
    convert << fileno;
    std::string scr = PSIOManager::shared_object()->get_default_path();
    std::string pid_ = psio_getpid();
    // std::string fname = scr + "psi_dfocc." + convert.str();
    std::string fname = scr + "psi." + pid_ + "." + convert.str();

    // write binary data
    std::ofstream OutFile;
    if (append)
        OutFile.open(const_cast<char *>(fname.c_str()), std::ios::out | std::ios::binary | std::ios::app);
    else
        OutFile.open(const_cast<char *>(fname.c_str()), std::ios::out | std::ios::binary);
    OutFile.write((char *)A2d_[0], dim1_ * dim2_ * sizeof(float));
    OutFile.close();
}  //

void Tensor2f::myread(const std::string &filename) {
    // read binary data
    std::ifstream InFile;
    InFile.open(const_cast<char *>(filename.c_str()), std::ios::in | std::ios::binary);
    InFile.read((char *)A2d_[0], dim1_ * dim2_ * sizeof(float));
    InFile.close();

}  //

void Tensor2f::myread(int fileno) {
    std::ostringstream convert;
    convert << fileno;
    std::string scr = PSIOManager::shared_object()->get_default_path();
    std::string pid_ = psio_getpid();
    // std::string fname = scr + "psi_dfocc." + convert.str();
    std::string fname = scr + "psi." + pid_ + "." + convert.str();

    // read binary data
    std::ifstream InFile;
    InFile.open(const_cast<char *>(fname.c_str()), std::ios::in | std::ios::binary);
    InFile.read((char *)A2d_[0], dim1_ * dim2_ * sizeof(float));
    InFile.close();

}  //

void Tensor2f::myread(int fileno, bool append) {
    std::ostringstream convert;
    convert << fileno;
    std::string scr = PSIOManager::shared_object()->get_default_path();
    std::string pid_ = psio_getpid();
    // std::string fname = scr + "psi_dfocc." + convert.str();
    std::string fname = scr + "psi." + pid_ + "." + convert.str();

    // read binary data
    std::ifstream InFile;
    if (append)
        InFile.open(const_cast<char *>(fname.c_str()), std::ios::in | std::ios::binary | std::ios::app);
    else
        InFile.open(const_cast<char *>(fname.c_str()), std::ios::in | std::ios::binary);
    InFile.read((char *)A2d_[0], dim1_ * dim2_ * sizeof(float));
    InFile.close();

}  //

void Tensor2f::myread(int fileno, size_t start) {
    std::ostringstream convert;
    convert << fileno;
    std::string scr = PSIOManager::shared_object()->get_default_path();
    std::string pid_ = psio_getpid();
    std::string fname = scr + "psi." + pid_ + "." + convert.str();
    // std::string fname = scr + "psi_dfocc." + convert.str();

    // read binary data
    std::ifstream InFile;
    InFile.open(const_cast<char *>(fname.c_str()), std::ios::in | std::ios::binary);
    InFile.seekg(start, std::ios::beg);
    InFile.read((char *)A2d_[0], dim1_ * dim2_ * sizeof(float));
    InFile.close();

}  //

float **Tensor2f::to_block_matrix () {
    float **temp = block_matrix_float   (dim1_, dim2_);
    // memcpy(&(temp[0][0]), &(A2d_[0][0]), dim1_ * dim2_ * sizeof(float));
    size_t length;
    length = (size_t)dim1_ * (size_t)dim2_;
    C_SCOPY(length, A2d_[0], 1, temp[0], 1);
    return temp;
}  //

// float *Tensor2f::to_lower_triangle() {
//     if (dim1_ != dim2_) return NULL;
//     int ntri = 0.5 * dim1_ * (dim1_ + 1);
//     float *tri = new float[ntri];
//     float **temp = to_block_matrix    ();
//     sq_to_tri(temp, tri, dim1_);
//     free_block_float(temp);
//     return tri;
// }  //

void Tensor2f::to_shared_matrix(SharedMatrix A) {
#pragma omp parallel for
    for (int i = 0; i < dim1_; ++i) {
        for (int j = 0; j < dim2_; ++j) {
            A->set(0, i, j, A2d_[i][j]);
        }
    }
}  //

void Tensor2f::to_matrix(SharedMatrix A) {
#pragma omp parallel for
    for (int i = 0; i < dim1_; ++i) {
        for (int j = 0; j < dim2_; ++j) {
            A->set(i, j, A2d_[i][j]);
        }
    }
}  //

void Tensor2f::to_pointer(float *A) {
#pragma omp parallel for
    for (int i = 0; i < dim1_; ++i) {
        for (int j = 0; j < dim2_; ++j) {
            int ij = j + (i * dim2_);
            A[ij] = A2d_[i][j];
        }
    }
}  //

void Tensor2f::mgs() {
    float rmgs1, rmgs2;
    //#pragma omp parallel for
    for (int k = 0; k < dim1_; k++) {  // loop-1
        rmgs1 = 0.0;

        for (int i = 0; i < dim1_; i++) {
            rmgs1 += A2d_[i][k] * A2d_[i][k];
        }

        rmgs1 = std::sqrt(rmgs1);

        for (int i = 0; i < dim1_; i++) {
            A2d_[i][k] /= rmgs1;
        }

        for (int j = (k + 1); j < dim1_; j++) {  // loop-2
            rmgs2 = 0.0;

            for (int i = 0; i < dim1_; i++) {
                rmgs2 += A2d_[i][k] * A2d_[i][j];
            }

            for (int i = 0; i < dim1_; i++) {
                A2d_[i][j] -= rmgs2 * A2d_[i][k];
            }
        }  // end 2
    }      // end 1

}  //

// void Tensor2f::gs() {
//     if (dim1_ != 0 && dim2_ != 0) {
//         schmidt(A2d_, dim1_, dim2_, "outfile");
//     }
// }  //

float *Tensor2f::row_vector(int n) {
    float *temp = new float[dim2_];
    memset(temp, 0, dim2_ * sizeof(float));
    for (int i = 0; i < dim2_; i++) temp[i] = A2d_[n][i];
    return temp;
}  //

float *Tensor2f::column_vector(int n) {
    float *temp = new float[dim1_];
    memset(temp, 0, dim1_ * sizeof(float));
    for (int i = 0; i < dim1_; i++) temp[i] = A2d_[i][n];
    return temp;
}  //

void Tensor2f::sort(int sort_type, const SharedTensor2f &A, float alpha, float beta) {
    int d1 = A->d1_;
    int d2 = A->d2_;
    int d3 = A->d3_;
    int d4 = A->d4_;

    if (sort_type == 1243) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = A->row_idx_[p][q];
                for (int r = 0; r < d3; r++) {
                    for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int sr = col_idx_[s][r];
                        A2d_[pq][sr] = (alpha * A->A2d_[pq][rs]) + (beta * A2d_[pq][sr]);
                    }
                }
            }
        }
    }

    else if (sort_type == 1324) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = A->row_idx_[p][q];
                for (int r = 0; r < d3; r++) {
                    for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int pr = row_idx_[p][r];
                        int qs = col_idx_[q][s];
                        A2d_[pr][qs] = (alpha * A->A2d_[pq][rs]) + (beta * A2d_[pr][qs]);
                    }
                }
            }
        }
    }

    else if (sort_type == 1342) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = A->row_idx_[p][q];
                for (int r = 0; r < d3; r++) {
                    for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int pr = row_idx_[p][r];
                        int sq = col_idx_[s][q];
                        A2d_[pr][sq] = (alpha * A->A2d_[pq][rs]) + (beta * A2d_[pr][sq]);
                    }
                }
            }
        }
    }

    else if (sort_type == 1423) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = A->row_idx_[p][q];
                for (int r = 0; r < d3; r++) {
                    for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int ps = row_idx_[p][s];
                        int qr = col_idx_[q][r];
                        A2d_[ps][qr] = (alpha * A->A2d_[pq][rs]) + (beta * A2d_[ps][qr]);
                    }
                }
            }
        }
    }

    else if (sort_type == 1432) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = A->row_idx_[p][q];
                for (int r = 0; r < d3; r++) {
                    int rq = col_idx_[r][q];
                    for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int ps = row_idx_[p][s];
                        A2d_[ps][rq] = (alpha * A->A2d_[pq][rs]) + (beta * A2d_[ps][rq]);
                    }
                }
            }
        }
    }

    else if (sort_type == 2134) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = A->row_idx_[p][q];
                for (int r = 0; r < d3; r++) {
                    for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int qp = row_idx_[q][p];
                        A2d_[qp][rs] = (alpha * A->A2d_[pq][rs]) + (beta * A2d_[qp][rs]);
                    }
                }
            }
        }
    }

    else if (sort_type == 2143) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = A->row_idx_[p][q];
                for (int r = 0; r < d3; r++) {
                    for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int qp = row_idx_[q][p];
                        int sr = col_idx_[s][r];
                        A2d_[qp][sr] = (alpha * A->A2d_[pq][rs]) + (beta * A2d_[qp][sr]);
                    }
                }
            }
        }
    }

    else if (sort_type == 2314) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = A->row_idx_[p][q];
                for (int r = 0; r < d3; r++) {
                    for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int qr = row_idx_[q][r];
                        int ps = col_idx_[p][s];
                        A2d_[qr][ps] = (alpha * A->A2d_[pq][rs]) + (beta * A2d_[qr][ps]);
                    }
                }
            }
        }
    }

    else if (sort_type == 2341) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = A->row_idx_[p][q];
                for (int r = 0; r < d3; r++) {
                    for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int qr = row_idx_[q][r];
                        int sp = col_idx_[s][p];
                        A2d_[qr][sp] = (alpha * A->A2d_[pq][rs]) + (beta * A2d_[qr][sp]);
                    }
                }
            }
        }
    }

    else if (sort_type == 2413) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = A->row_idx_[p][q];
                for (int r = 0; r < d3; r++) {
                    for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int qs = row_idx_[q][s];
                        int pr = col_idx_[p][r];
                        A2d_[qs][pr] = (alpha * A->A2d_[pq][rs]) + (beta * A2d_[qs][pr]);
                    }
                }
            }
        }
    }

    else if (sort_type == 2431) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = A->row_idx_[p][q];
                for (int r = 0; r < d3; r++) {
                    for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int qs = row_idx_[q][s];
                        int rp = col_idx_[r][p];
                        A2d_[qs][rp] = (alpha * A->A2d_[pq][rs]) + (beta * A2d_[qs][rp]);
                    }
                }
            }
        }
    }

    else if (sort_type == 3124) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = A->row_idx_[p][q];
                for (int r = 0; r < d3; r++) {
                    for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int rp = row_idx_[r][p];
                        int qs = col_idx_[q][s];
                        A2d_[rp][qs] = (alpha * A->A2d_[pq][rs]) + (beta * A2d_[rp][qs]);
                    }
                }
            }
        }
    }

    else if (sort_type == 3142) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = A->row_idx_[p][q];
                for (int r = 0; r < d3; r++) {
                    for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int rp = row_idx_[r][p];
                        int sq = col_idx_[s][q];
                        A2d_[rp][sq] = (alpha * A->A2d_[pq][rs]) + (beta * A2d_[rp][sq]);
                    }
                }
            }
        }
    }

    else if (sort_type == 3214) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = A->row_idx_[p][q];
                for (int r = 0; r < d3; r++) {
                    int rq = row_idx_[r][q];
                    for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int ps = col_idx_[p][s];
                        A2d_[rq][ps] = (alpha * A->A2d_[pq][rs]) + (beta * A2d_[rq][ps]);
                    }
                }
            }
        }
    }

    else if (sort_type == 3241) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = A->row_idx_[p][q];
                for (int r = 0; r < d3; r++) {
                    for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int rq = row_idx_[r][q];
                        int sp = col_idx_[s][p];
                        A2d_[rq][sp] = (alpha * A->A2d_[pq][rs]) + (beta * A2d_[rq][sp]);
                    }
                }
            }
        }
    }

    else if (sort_type == 3412) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = A->row_idx_[p][q];
                for (int r = 0; r < d3; r++) {
                    for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        A2d_[rs][pq] = (alpha * A->A2d_[pq][rs]) + (beta * A2d_[rs][pq]);
                    }
                }
            }
        }
    }

    else if (sort_type == 3421) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = A->row_idx_[p][q];
                for (int r = 0; r < d3; r++) {
                    for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int qp = col_idx_[q][p];
                        A2d_[rs][qp] = (alpha * A->A2d_[pq][rs]) + (beta * A2d_[rs][qp]);
                    }
                }
            }
        }
    }

    else if (sort_type == 4123) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = A->row_idx_[p][q];
                for (int r = 0; r < d3; r++) {
                    for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int sp = row_idx_[s][p];
                        int qr = col_idx_[q][r];
                        A2d_[sp][qr] = (alpha * A->A2d_[pq][rs]) + (beta * A2d_[sp][qr]);
                    }
                }
            }
        }
    }

    else if (sort_type == 4132) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = A->row_idx_[p][q];
                for (int r = 0; r < d3; r++) {
                    for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int sp = row_idx_[s][p];
                        int rq = col_idx_[r][q];
                        A2d_[sp][rq] = (alpha * A->A2d_[pq][rs]) + (beta * A2d_[sp][rq]);
                    }
                }
            }
        }
    }

    else if (sort_type == 4213) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = A->row_idx_[p][q];
                for (int r = 0; r < d3; r++) {
                    for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int sq = row_idx_[s][q];
                        int pr = col_idx_[p][r];
                        A2d_[sq][pr] = (alpha * A->A2d_[pq][rs]) + (beta * A2d_[sq][pr]);
                    }
                }
            }
        }
    }

    else if (sort_type == 4231) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = A->row_idx_[p][q];
                for (int r = 0; r < d3; r++) {
                    for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int sq = row_idx_[s][q];
                        int rp = col_idx_[r][p];
                        A2d_[sq][rp] = (alpha * A->A2d_[pq][rs]) + (beta * A2d_[sq][rp]);
                    }
                }
            }
        }
    }

    else if (sort_type == 4312) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = A->row_idx_[p][q];
                for (int r = 0; r < d3; r++) {
                    for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int sr = row_idx_[s][r];
                        A2d_[sr][pq] = (alpha * A->A2d_[pq][rs]) + (beta * A2d_[sr][pq]);
                    }
                }
            }
        }
    }

    else if (sort_type == 4321) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = A->row_idx_[p][q];
                for (int r = 0; r < d3; r++) {
                    for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int sr = row_idx_[s][r];
                        int qp = col_idx_[q][p];
                        A2d_[sr][qp] = (alpha * A->A2d_[pq][rs]) + (beta * A2d_[sr][qp]);
                    }
                }
            }
        }
    }

    else {
        outfile->Printf("\tUnrecognized sort type!\n");
        throw PSIEXCEPTION("Unrecognized sort type!");
    }

}  //

void Tensor2f::sort3a(int sort_type, int d1, int d2, int d3, const SharedTensor2f &A, float alpha, float beta) {
    if (sort_type == 132) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                for (int r = 0; r < d3; r++) {
                    int rq = q + (r * d2);
                    int qr = r + (q * d3);
                    A2d_[p][rq] = (alpha * A->A2d_[p][qr]) + (beta * A2d_[p][rq]);
                }
            }
        }
    }

    else {
        outfile->Printf("\tUnrecognized sort type!\n");
        throw PSIEXCEPTION("Unrecognized sort type!");
    }

}  //

void Tensor2f::sort3b(int sort_type, int d1, int d2, int d3, const SharedTensor2f &A, float alpha, float beta) {
    if (sort_type == 132) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = q + (p * d2);
                for (int r = 0; r < d3; r++) {
                    int pr = r + (p * d3);
                    A2d_[pr][q] = (alpha * A->A2d_[pq][r]) + (beta * A2d_[pr][q]);
                }
            }
        }
    }

    else if (sort_type == 213) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = q + (p * d2);
                int qp = p + (q * d1);
                for (int r = 0; r < d3; r++) {
                    A2d_[qp][r] = (alpha * A->A2d_[pq][r]) + (beta * A2d_[qp][r]);
                }
            }
        }
    }

    else if (sort_type == 312) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = q + (p * d2);
                for (int r = 0; r < d3; r++) {
                    int rp = p + (r * d1);
                    A2d_[rp][q] = (alpha * A->A2d_[pq][r]) + (beta * A2d_[rp][q]);
                }
            }
        }
    }

    else if (sort_type == 231) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = q + (p * d2);
                for (int r = 0; r < d3; r++) {
                    int qr = r + (q * d3);
                    A2d_[qr][p] = (alpha * A->A2d_[pq][r]) + (beta * A2d_[qr][p]);
                }
            }
        }
    }

    else if (sort_type == 321) {
#pragma omp parallel for
        for (int p = 0; p < d1; p++) {
            for (int q = 0; q < d2; q++) {
                int pq = q + (p * d2);
                for (int r = 0; r < d3; r++) {
                    int rq = q + (r * d2);
                    A2d_[rq][p] = (alpha * A->A2d_[pq][r]) + (beta * A2d_[rq][p]);
                }
            }
        }
    }

    else {
        outfile->Printf("\tUnrecognized sort type!\n");
        throw PSIEXCEPTION("Unrecognized sort type!");
    }

}  //

void Tensor2f::apply_denom(int frzc, int occ, const SharedTensor2f &fock) {
    int aocc = d1_;
    int avir = d3_;

#pragma omp parallel for
    for (int i = 0; i < aocc; i++) {
        float di = fock->A2d_[i + frzc][i + frzc];
        for (int j = 0; j < aocc; j++) {
            float dij = di + fock->A2d_[j + frzc][j + frzc];
            int ij = row_idx_[i][j];
            for (int a = 0; a < avir; a++) {
                float dija = dij - fock->A2d_[a + occ][a + occ];
                for (int b = 0; b < avir; b++) {
                    float dijab = dija - fock->A2d_[b + occ][b + occ];
                    int ab = col_idx_[a][b];
                    A2d_[ij][ab] /= dijab;
                }
            }
        }
    }
}  //

void Tensor2f::apply_denom_os(int frzc, int occA, int occB, const SharedTensor2f &fockA, const SharedTensor2f &fockB) {
    int aoccA = d1_;
    int aoccB = d2_;
    int avirA = d3_;
    int avirB = d4_;

#pragma omp parallel for
    for (int i = 0; i < aoccA; i++) {
        float di = fockA->A2d_[i + frzc][i + frzc];
        for (int j = 0; j < aoccB; j++) {
            float dij = di + fockB->A2d_[j + frzc][j + frzc];
            int ij = row_idx_[i][j];
            for (int a = 0; a < avirA; a++) {
                float dija = dij - fockA->A2d_[a + occA][a + occA];
                for (int b = 0; b < avirB; b++) {
                    float dijab = dija - fockB->A2d_[b + occB][b + occB];
                    int ab = col_idx_[a][b];
                    A2d_[ij][ab] /= dijab;
                }
            }
        }
    }
}  //

void Tensor2f::apply_denom_chem(int frzc, int occ, const SharedTensor2f &fock) {
    int aocc = d1_;
    int avir = d2_;

#pragma omp parallel for
    for (int i = 0; i < aocc; i++) {
        float di = fock->A2d_[i + frzc][i + frzc];
        for (int a = 0; a < avir; a++) {
            float dia = di - fock->A2d_[a + occ][a + occ];
            int ia = row_idx_[i][a];
            for (int j = 0; j < aocc; j++) {
                float diaj = dia + fock->A2d_[j + frzc][j + frzc];
                for (int b = 0; b < avir; b++) {
                    float diajb = diaj - fock->A2d_[b + occ][b + occ];
                    int jb = col_idx_[j][b];
                    A2d_[ia][jb] /= diajb;
                }
            }
        }
    }
}  //

void Tensor2f::reg_denom(int frzc, int occ, const SharedTensor2f &fock, float reg) {
    int aocc = d1_;
    int avir = d3_;

#pragma omp parallel for
    for (int i = 0; i < aocc; i++) {
        float di = fock->A2d_[i + frzc][i + frzc] - reg;
        for (int j = 0; j < aocc; j++) {
            float dij = di + fock->A2d_[j + frzc][j + frzc];
            int ij = row_idx_[i][j];
            for (int a = 0; a < avir; a++) {
                float dija = dij - fock->A2d_[a + occ][a + occ];
                for (int b = 0; b < avir; b++) {
                    float dijab = dija - fock->A2d_[b + occ][b + occ];
                    int ab = col_idx_[a][b];
                    A2d_[ij][ab] /= dijab;
                }
            }
        }
    }
}  //

void Tensor2f::reg_denom_os(int frzc, int occA, int occB, const SharedTensor2f &fockA, const SharedTensor2f &fockB,
                            float reg) {
    int aoccA = d1_;
    int aoccB = d2_;
    int avirA = d3_;
    int avirB = d4_;

#pragma omp parallel for
    for (int i = 0; i < aoccA; i++) {
        float di = fockA->A2d_[i + frzc][i + frzc] - reg;
        for (int j = 0; j < aoccB; j++) {
            float dij = di + fockB->A2d_[j + frzc][j + frzc];
            int ij = row_idx_[i][j];
            for (int a = 0; a < avirA; a++) {
                float dija = dij - fockA->A2d_[a + occA][a + occA];
                for (int b = 0; b < avirB; b++) {
                    float dijab = dija - fockB->A2d_[b + occB][b + occB];
                    int ab = col_idx_[a][b];
                    A2d_[ij][ab] /= dijab;
                }
            }
        }
    }
}  //

void Tensor2f::reg_denom_chem(int frzc, int occ, const SharedTensor2f &fock, float reg) {
    int aocc = d1_;
    int avir = d2_;

#pragma omp parallel for
    for (int i = 0; i < aocc; i++) {
        float di = fock->A2d_[i + frzc][i + frzc] - reg;
        for (int a = 0; a < avir; a++) {
            float dia = di - fock->A2d_[a + occ][a + occ];
            int ia = row_idx_[i][a];
            for (int j = 0; j < aocc; j++) {
                float diaj = dia + fock->A2d_[j + frzc][j + frzc];
                for (int b = 0; b < avir; b++) {
                    float diajb = diaj - fock->A2d_[b + occ][b + occ];
                    int jb = col_idx_[j][b];
                    A2d_[ia][jb] /= diajb;
                }
            }
        }
    }
}  //

void Tensor2f::dirprd(const SharedTensor2f &a, const SharedTensor2f &b) {
#pragma omp parallel for
    for (int i = 0; i < dim1_; i++) {
        for (int j = 0; j < dim2_; j++) {
            A2d_[i][j] = a->get(i, j) * b->get(i, j);
        }
    }
}  //

void Tensor2f::dirprd123(const SharedTensor1f &a, const SharedTensor2f &b, float alpha, float beta) {
    int d1 = dim1_;
    int d2 = b->dim1();
    int d3 = b->dim2();
#pragma omp parallel for
    for (int i = 0; i < d1; i++) {
        for (int j = 0; j < d2; j++) {
            for (int k = 0; k < d3; k++) {
                int jk = k + (j * d3);
                // A2d_[i][jk] = a->get(i) * b->get(j,k);
                A2d_[i][jk] = (alpha * a->get(i) * b->get(j, k)) + (beta * A2d_[i][jk]);
            }
        }
    }
}  //

void Tensor2f::dirprd123(bool transb, const SharedTensor1f &a, const SharedTensor2f &b, float alpha, float beta) {
    if (transb) {
        int d1 = dim1_;
        int d2 = b->dim2();
        int d3 = b->dim1();
#pragma omp parallel for
        for (int i = 0; i < d1; i++) {
            for (int j = 0; j < d2; j++) {
                for (int k = 0; k < d3; k++) {
                    int jk = k + (j * d3);
                    A2d_[i][jk] = (alpha * a->get(i) * b->get(k, j)) + (beta * A2d_[i][jk]);
                }
            }
        }
    }

    else {
        int d1 = dim1_;
        int d2 = b->dim1();
        int d3 = b->dim2();
#pragma omp parallel for
        for (int i = 0; i < d1; i++) {
            for (int j = 0; j < d2; j++) {
                for (int k = 0; k < d3; k++) {
                    int jk = k + (j * d3);
                    // A2d_[i][jk] = a->get(i) * b->get(j,k);
                    A2d_[i][jk] = (alpha * a->get(i) * b->get(j, k)) + (beta * A2d_[i][jk]);
                }
            }
        }
    }
}  //

void Tensor2f::dirprd112(const SharedTensor1f &a, const SharedTensor1f &b) {
#pragma omp parallel for
    for (int i = 0; i < dim1_; i++) {
        for (int j = 0; j < dim2_; j++) {
            A2d_[i][j] = a->get(i) * b->get(j);
        }
    }
}  //

void Tensor2f::dirprd112(const SharedTensor1f &a, const SharedTensor1f &b, float alpha, float beta) {
#pragma omp parallel for
    for (int i = 0; i < dim1_; i++) {
        for (int j = 0; j < dim2_; j++) {
            A2d_[i][j] = (alpha * a->get(i) * b->get(j)) + (beta * A2d_[i][j]);
        }
    }
}  //

void Tensor2f::dirprd224(const SharedTensor2f &a, const SharedTensor2f &b) {
#pragma omp parallel for
    for (int i = 0; i < d1_; i++) {
        for (int j = 0; j < d2_; j++) {
            int ij = row_idx_[i][j];
            for (int k = 0; k < d3_; k++) {
                for (int l = 0; l < d4_; l++) {
                    int kl = col_idx_[k][l];
                    A2d_[ij][kl] = a->get(i, j) * b->get(k, l);
                }
            }
        }
    }
}  //

void Tensor2f::dirprd224(const SharedTensor2f &a, const SharedTensor2f &b, float alpha, float beta) {
#pragma omp parallel for
    for (int i = 0; i < d1_; i++) {
        for (int j = 0; j < d2_; j++) {
            int ij = row_idx_[i][j];
            for (int k = 0; k < d3_; k++) {
                for (int l = 0; l < d4_; l++) {
                    int kl = col_idx_[k][l];
                    A2d_[ij][kl] = (alpha * a->get(i, j) * b->get(k, l)) + (beta * A2d_[ij][kl]);
                }
            }
        }
    }
}  //

float *Tensor2f::to_vector(const SharedTensor2i &pair_idx) {
    float *temp = new float[dim1_ * dim2_];
#pragma omp parallel for
    for (int i = 0; i < dim1_; i++) {
        for (int j = 0; j < dim2_; j++) {
            int ij = pair_idx->get(i, j);
            temp[ij] = A2d_[i][j];
        }
    }
    return temp;
}  //

float *Tensor2f::to_vector() {
    float *temp = new float[dim1_ * dim2_];
#pragma omp parallel for
    for (int i = 0; i < dim1_; i++) {
        for (int j = 0; j < dim2_; j++) {
            int ij = (i * dim2_) + j;
            temp[ij] = A2d_[i][j];
        }
    }
    return temp;
}  //

float Tensor2f::rms() {
    float summ = 0.0;
#pragma omp parallel for
    for (int i = 0; i < dim1_; ++i) {
        for (int j = 0; j < dim2_; ++j) {
            summ += A2d_[i][j] * A2d_[i][j];
        }
    }
    summ = std::sqrt(summ / (dim1_ * dim2_));

    return summ;
}  //

float Tensor2f::rms(const SharedTensor2f &a) {
    float summ = 0.0;
#pragma omp parallel for
    for (int i = 0; i < dim1_; ++i) {
        for (int j = 0; j < dim2_; ++j) {
            summ += (A2d_[i][j] - a->A2d_[i][j]) * (A2d_[i][j] - a->A2d_[i][j]);
        }
    }
    summ = std::sqrt(summ / (dim1_ * dim2_));

    return summ;
}  //

void Tensor2f::set_act_oo(int aocc, const SharedTensor2f &a) {
#pragma omp parallel for
    for (int i = 0; i < aocc; i++) {
        for (int j = 0; j < aocc; j++) {
            A2d_[i][j] = a->get(i, j);
        }
    }
}  //

void Tensor2f::set_act_oo(int frzc, int aocc, const SharedTensor2f &a) {
#pragma omp parallel for
    for (int i = 0; i < aocc; i++) {
        for (int j = 0; j < aocc; j++) {
            A2d_[i + frzc][j + frzc] = a->get(i, j);
        }
    }
}  //

void Tensor2f::set_oo(const SharedTensor2f &a) {
    int occ = a->dim1();
#pragma omp parallel for
    for (int i = 0; i < occ; i++) {
        for (int j = 0; j < occ; j++) {
            A2d_[i][j] = a->get(i, j);
        }
    }
}  //

void Tensor2f::set_act_ov(int frzc, const SharedTensor2f &A) {
    int aocc = A->dim1();
    int avir = A->dim2();
#pragma omp parallel for
    for (int i = 0; i < aocc; i++) {
        for (int a = 0; a < avir; a++) {
            A2d_[i + frzc][a] = A->get(i, a);
        }
    }
}  //

void Tensor2f::set_act_vo(int frzc, const SharedTensor2f &A) {
    int avir = A->dim1();
    int aocc = A->dim2();
#pragma omp parallel for
    for (int a = 0; a < avir; a++) {
        for (int i = 0; i < aocc; i++) {
            A2d_[a][i + frzc] = A->get(a, i);
        }
    }
}  //

void Tensor2f::set_act_vv(int occ, int avir, const SharedTensor2f &A) {
#pragma omp parallel for
    for (int a = 0; a < avir; a++) {
        for (int b = 0; b < avir; b++) {
            A2d_[a + occ][b + occ] = A->get(a, b);
        }
    }
}  //

void Tensor2f::set_act_vv(const SharedTensor2f &A) {
    int avir = A->dim1();
#pragma omp parallel for
    for (int a = 0; a < avir; a++) {
        for (int b = 0; b < avir; b++) {
            A2d_[a][b] = A->get(a, b);
        }
    }
}  //

void Tensor2f::set_vv(int occ, const SharedTensor2f &A) {
    int vir = A->dim1();
#pragma omp parallel for
    for (int a = 0; a < vir; a++) {
        for (int b = 0; b < vir; b++) {
            A2d_[a + occ][b + occ] = A->get(a, b);
        }
    }
}  //

void Tensor2f::set_ov(const SharedTensor2f &A) {
    int occ = A->dim1();
    int vir = A->dim2();
#pragma omp parallel for
    for (int i = 0; i < occ; i++) {
        for (int a = 0; a < vir; a++) {
            A2d_[i][a + occ] = A->get(i, a);
        }
    }
}  //

void Tensor2f::set_vo(const SharedTensor2f &A) {
    int vir = A->dim1();
    int occ = A->dim2();
#pragma omp parallel for
    for (int a = 0; a < vir; a++) {
        for (int i = 0; i < occ; i++) {
            A2d_[a + occ][i] = A->get(a, i);
        }
    }
}  //

void Tensor2f::add_oo(const SharedTensor2f &A, float alpha, float beta) {
    int occ = A->dim1();
#pragma omp parallel for
    for (int i = 0; i < occ; i++) {
        for (int j = 0; j < occ; j++) {
            A2d_[i][j] = (alpha * A->get(i, j)) + (beta * A2d_[i][j]);
        }
    }
}  //

void Tensor2f::add_vv(int occ, const SharedTensor2f &A, float alpha, float beta) {
    int vir = A->dim1();
#pragma omp parallel for
    for (int a = 0; a < vir; a++) {
        for (int b = 0; b < vir; b++) {
            A2d_[a + occ][b + occ] = (alpha * A->get(a, b)) + (beta * A2d_[a + occ][b + occ]);
        }
    }
}  //

void Tensor2f::add_ov(const SharedTensor2f &A, float alpha, float beta) {
    int occ = A->dim1();
    int vir = A->dim2();
#pragma omp parallel for
    for (int i = 0; i < occ; i++) {
        for (int a = 0; a < vir; a++) {
            A2d_[i][a + occ] = (alpha * A->get(i, a)) + (beta * A2d_[i][a + occ]);
        }
    }
}  //

void Tensor2f::add_vo(const SharedTensor2f &A, float alpha, float beta) {
    int vir = A->dim1();
    int occ = A->dim2();
#pragma omp parallel for
    for (int a = 0; a < vir; a++) {
        for (int i = 0; i < occ; i++) {
            A2d_[a + occ][i] = (alpha * A->get(a, i)) + (beta * A2d_[a + occ][i]);
        }
    }
}  //

void Tensor2f::add_aocc_fc(const SharedTensor2f &A, float alpha, float beta) {
    int aocc = A->dim1();
    int frzc = A->dim2();
#pragma omp parallel for
    for (int i = 0; i < aocc; i++) {
        for (int j = 0; j < frzc; j++) {
            A2d_[i + frzc][j] = (alpha * A->get(i, j)) + (beta * A2d_[i + frzc][j]);
        }
    }
}  //

void Tensor2f::add_fc_aocc(const SharedTensor2f &A, float alpha, float beta) {
    int frzc = A->dim1();
    int aocc = A->dim2();
#pragma omp parallel for
    for (int i = 0; i < frzc; i++) {
        for (int j = 0; j < aocc; j++) {
            A2d_[i][j + frzc] = (alpha * A->get(i, j)) + (beta * A2d_[i][j + frzc]);
        }
    }
}  //

void Tensor2f::set3_oo(const SharedTensor2f &A) {
    int naux = A->d1_;
    int occ = A->d2_;
#pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
        for (int i = 0; i < occ; i++) {
            for (int j = 0; j < occ; j++) {
                int ij = A->col_idx_[i][j];
                int oo = col_idx_[i][j];
                A2d_[Q][oo] = A->get(Q, ij);
            }
        }
    }
}  //

void Tensor2f::set3_act_oo(int frzc, const SharedTensor2f &A) {
    int naux = A->d1_;
    int aoccA = A->d2_;
    int aoccB = A->d3_;
    int occA = d2_;
    int occB = d3_;
#pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
        for (int i = 0; i < aoccA; i++) {
            for (int j = 0; j < aoccB; j++) {
                int ij = A->col_idx_[i][j];
                int oo = ((i + frzc) * occB) + j + frzc;
                A2d_[Q][oo] = A->get(Q, ij);
            }
        }
    }
}  //

void Tensor2f::add3_oo(const SharedTensor2f &A, float alpha, float beta) {
    int naux = A->d1_;
    int occ = A->d2_;
#pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
        for (int i = 0; i < occ; i++) {
            for (int j = 0; j < occ; j++) {
                int ij = A->col_idx_[i][j];
                int oo = col_idx_[i][j];
                A2d_[Q][oo] = (alpha * A->get(Q, ij)) + (beta * A2d_[Q][oo]);
            }
        }
    }
}  //

void Tensor2f::set3_act_ov(int frzc, int aocc, int avir, int vir, const SharedTensor2f &A) {
    int naux = dim1_;
#pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
        for (int i = 0; i < aocc; i++) {
            for (int a = 0; a < avir; a++) {
                int ia = (i * avir) + a;
                int ov = ((i + frzc) * vir) + a;
                A2d_[Q][ov] = A->get(Q, ia);
            }
        }
    }
}  //

void Tensor2f::set3_ov(const SharedTensor2f &A) {
    int naux = dim1_;
    int occ = A->d2_;
    int vir = A->d3_;
#pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
        for (int i = 0; i < occ; i++) {
            for (int a = 0; a < vir; a++) {
                int ia = A->col_idx_[i][a];
                int ov = col_idx_[i][a + occ];
                A2d_[Q][ov] = A->get(Q, ia);
            }
        }
    }
}  //

void Tensor2f::set3_vo(const SharedTensor2f &A) {
    int naux = dim1_;
    int vir = A->d2_;
    int occ = A->d3_;
#pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
        for (int a = 0; a < vir; a++) {
            for (int i = 0; i < occ; i++) {
                int ai = A->col_idx_[a][i];
                int vo = col_idx_[a + occ][i];
                A2d_[Q][vo] = A->get(Q, ai);
            }
        }
    }
}  //

void Tensor2f::set3_vv(const SharedTensor2f &A, int occ) {
    int naux = dim1_;
    int vir = A->d2_;
#pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
        for (int a = 0; a < vir; a++) {
            for (int b = 0; b < vir; b++) {
                int ab = A->col_idx_[a][b];
                int vv = col_idx_[a + occ][b + occ];
                A2d_[Q][vv] = A->get(Q, ab);
            }
        }
    }
}  //

void Tensor2f::set3_act_vv(const SharedTensor2f &A) {
    int naux = dim1_;
    int avir = A->d2_;
#pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
        for (int a = 0; a < avir; a++) {
            for (int b = 0; b < avir; b++) {
                int ab = A->col_idx_[a][b];
                int vv = col_idx_[a][b];
                A2d_[Q][vv] = A->get(Q, ab);
            }
        }
    }
}  //

void Tensor2f::swap_3index_col(const SharedTensor2f &A) {
    int d1 = A->d1_;
    int d2 = A->d2_;
    int d3 = A->d3_;

#pragma omp parallel for
    for (int Q = 0; Q < d1; Q++) {
        for (int p = 0; p < d2; p++) {
            for (int q = 0; q < d3; q++) {
                int pq = A->col_idx_[p][q];
                int qp = col_idx_[q][p];
                A2d_[Q][qp] = A->A2d_[Q][pq];
            }
        }
    }
}  //

void Tensor2f::form_oo(const SharedTensor2f &A) {
    int occ = dim1_;
#pragma omp parallel for
    for (int i = 0; i < occ; i++) {
        for (int j = 0; j < occ; j++) {
            A2d_[i][j] = A->get(i, j);
        }
    }
}  //

void Tensor2f::form_act_oo(int frzc, const SharedTensor2f &A) {
    int occ = dim1_;
#pragma omp parallel for
    for (int i = 0; i < occ; i++) {
        for (int j = 0; j < occ; j++) {
            A2d_[i][j] = A->get(i + frzc, j + frzc);
        }
    }
}  //

void Tensor2f::form_vv(int occ, const SharedTensor2f &A) {
    int vir = dim1_;
#pragma omp parallel for
    for (int a = 0; a < vir; a++) {
        for (int b = 0; b < vir; b++) {
            A2d_[a][b] = A->get(a + occ, b + occ);
        }
    }
}  //

void Tensor2f::form_act_vv(int occ, const SharedTensor2f &A) {
    int vir = dim1_;
#pragma omp parallel for
    for (int a = 0; a < vir; a++) {
        for (int b = 0; b < vir; b++) {
            A2d_[a][b] = A->get(a + occ, b + occ);
        }
    }
}  //

void Tensor2f::form_vo(const SharedTensor2f &A) {
    int vir = dim1_;
    int occ = dim2_;
#pragma omp parallel for
    for (int a = 0; a < vir; a++) {
        for (int i = 0; i < occ; i++) {
            A2d_[a][i] = A->get(a + occ, i);
        }
    }
}  //

void Tensor2f::form_vo(int occ, const SharedTensor2f &A) {
    int vir = dim1_;
    int occ2 = dim2_;
#pragma omp parallel for
    for (int a = 0; a < vir; a++) {
        for (int i = 0; i < occ2; i++) {
            A2d_[a][i] = A->get(a + occ, i);
        }
    }
}  //

void Tensor2f::form_act_vo(int frzc, const SharedTensor2f &A) {
    int vir = dim1_;
    int occ = dim2_;
#pragma omp parallel for
    for (int a = 0; a < vir; a++) {
        for (int i = 0; i < occ; i++) {
            A2d_[a][i] = A->get(a + occ, i + frzc);
        }
    }
}  //

void Tensor2f::form_act_vo(int frzc, int occ, const SharedTensor2f &A) {
    int vir = dim1_;
    int occ2 = dim2_;
#pragma omp parallel for
    for (int a = 0; a < vir; a++) {
        for (int i = 0; i < occ2; i++) {
            A2d_[a][i] = A->get(a + occ, i + frzc);
        }
    }
}  //

void Tensor2f::form_ov(const SharedTensor2f &A) {
    int vir = dim2_;
    int occ = dim1_;
#pragma omp parallel for
    for (int i = 0; i < occ; i++) {
        for (int a = 0; a < vir; a++) {
            A2d_[i][a] = A->get(i, a + occ);
        }
    }
}  //

void Tensor2f::form_ov(int occ, const SharedTensor2f &A) {
    int vir = dim2_;
    int occ2 = dim1_;
#pragma omp parallel for
    for (int i = 0; i < occ2; i++) {
        for (int a = 0; a < vir; a++) {
            A2d_[i][a] = A->get(i, a + occ);
        }
    }
}  //

void Tensor2f::form_act_ov(int frzc, const SharedTensor2f &A) {
    int vir = dim2_;
    int occ = dim1_;
#pragma omp parallel for
    for (int i = 0; i < occ; i++) {
        for (int a = 0; a < vir; a++) {
            A2d_[i][a] = A->get(i + frzc, a + occ);
        }
    }
}  //

void Tensor2f::form_act_ov(int frzc, int occ, const SharedTensor2f &A) {
    int vir = dim2_;
    int occ2 = dim1_;
#pragma omp parallel for
    for (int i = 0; i < occ2; i++) {
        for (int a = 0; a < vir; a++) {
            A2d_[i][a] = A->get(i + frzc, a + occ);
        }
    }
}  //

void Tensor2f::form_ooAB(const SharedTensor2f &A) {
    int occA = dim1_;
    int occB = dim2_;
#pragma omp parallel for
    for (int i = 0; i < occA; i++) {
        for (int j = 0; j < occB; j++) {
            A2d_[i][j] = A->get(i, j);
        }
    }
}  //

void Tensor2f::form_b_ij(int frzc, const SharedTensor2f &A) {
#pragma omp parallel for
    for (int Q = 0; Q < d1_; Q++) {
        for (int i = 0; i < d2_; i++) {
            for (int j = 0; j < d3_; j++) {
                int ij = col_idx_[i][j];
                int oo = A->col_idx_[i + frzc][j + frzc];
                A2d_[Q][ij] = A->get(Q, oo);
            }
        }
    }
}  //

void Tensor2f::form_b_ia(int frzc, const SharedTensor2f &A) {
#pragma omp parallel for
    for (int Q = 0; Q < d1_; Q++) {
        for (int i = 0; i < d2_; i++) {
            for (int a = 0; a < d3_; a++) {
                int ia = col_idx_[i][a];
                int ov = A->col_idx_[i + frzc][a];
                A2d_[Q][ia] = A->get(Q, ov);
            }
        }
    }
}  //

void Tensor2f::form_b_ab(const SharedTensor2f &A) {
#pragma omp parallel for
    for (int Q = 0; Q < d1_; Q++) {
        for (int a = 0; a < d2_; a++) {
            for (int b = 0; b < d3_; b++) {
                int ab = col_idx_[a][b];
                int vv = A->col_idx_[a][b];
                A2d_[Q][ab] = A->get(Q, vv);
            }
        }
    }
}  //

void Tensor2f::form_b_kl(const SharedTensor2f &A) {
    int naux = d1_;
    int aocc = d2_;
    int frzc = d3_;
#pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
        for (int i = 0; i < aocc; i++) {
            for (int j = 0; j < frzc; j++) {
                int ij = A->col_idx_[i + frzc][j];
                int oo = col_idx_[i][j];
                A2d_[Q][oo] = A->get(Q, ij);
            }
        }
    }
}  //

void Tensor2f::form_b_ki(const SharedTensor2f &A) {
    int naux = d1_;
    int aocc = d2_;
    int occ = d3_;
    int frzc = occ - aocc;
#pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
        for (int i = 0; i < aocc; i++) {
            for (int j = 0; j < occ; j++) {
                int ij = A->col_idx_[i + frzc][j];
                int oo = col_idx_[i][j];
                A2d_[Q][oo] = A->get(Q, ij);
            }
        }
    }
}  //

void Tensor2f::form_b_ka(const SharedTensor2f &A) {
    int naux = d1_;
    int aocc = d2_;
    int vir = d3_;
    int occ = A->d2_;
    int frzc = occ - aocc;
#pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
        for (int i = 0; i < aocc; i++) {
            for (int j = 0; j < vir; j++) {
                int ij = A->col_idx_[i + frzc][j];
                int oo = col_idx_[i][j];
                A2d_[Q][oo] = A->get(Q, ij);
            }
        }
    }
}  //

void Tensor2f::form_b_li(const SharedTensor2f &A) {
    int naux = d1_;
    int frzc = d2_;
    int occ = d3_;
#pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
        for (int i = 0; i < frzc; i++) {
            for (int j = 0; j < occ; j++) {
                int ij = A->col_idx_[i][j];
                int oo = col_idx_[i][j];
                A2d_[Q][oo] = A->get(Q, ij);
            }
        }
    }
}  //

void Tensor2f::form_b_il(const SharedTensor2f &A) {
    int naux = d1_;
    int occ = d2_;
    int frzc = d3_;
#pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
        for (int i = 0; i < occ; i++) {
            for (int j = 0; j < frzc; j++) {
                int ij = A->col_idx_[i][j];
                int oo = col_idx_[i][j];
                A2d_[Q][oo] = A->get(Q, ij);
            }
        }
    }
}  //

void Tensor2f::form_b_la(const SharedTensor2f &A) {
    int naux = d1_;
    int frzc = d2_;
    int vir = d3_;
#pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
        for (int i = 0; i < frzc; i++) {
            for (int j = 0; j < vir; j++) {
                int ij = A->col_idx_[i][j];
                int oo = col_idx_[i][j];
                A2d_[Q][oo] = A->get(Q, ij);
            }
        }
    }
}  //

void Tensor2f::symmetrize() {
    SharedTensor2f temp = SharedTensor2f(new Tensor2f(dim2_, dim1_));
    temp = transpose();
    add(temp);
    scale(0.5);
    temp.reset();
    /*
    #pragma omp parallel for
    for (int i=0; i<dim2_; ++i) {
        for (int j=0; j<dim1_; ++j) {
             A2d_[i][j] = 0.5 * (A2d_[i][j] + A2d_[j][i]);
        }
    }
    */

}  //

void Tensor2f::symmetrize(const SharedTensor2f &A) {
#pragma omp parallel for
    for (int i = 0; i < dim1_; ++i) {
        for (int j = 0; j < dim2_; ++j) {
            A2d_[i][j] = 0.5 * (A->A2d_[i][j] + A->A2d_[j][i]);
        }
    }

}  //

void Tensor2f::symmetrize3(const SharedTensor2f &A) {
    SharedTensor2f temp = SharedTensor2f(new Tensor2f("temp", d1_, d3_, d2_));
    temp->swap_3index_col(A);
    add(temp);
    scale(0.5);
    temp.reset();
}  //

void Tensor2f::symm_packed(const SharedTensor2f &A) {
// Form symetric packed 3-index
#pragma omp parallel for
    for (int R = 0; R < A->d1_; R++) {
        for (int p = 0; p < A->d2_; p++) {
            for (int q = 0; q <= p; q++) {
                int pq = A->col_idx_[p][q];
                int pq_sym = index2(p, q);
                float perm = (p == q ? 1.0 : 2.0);
                A2d_[R][pq_sym] = perm * A->get(R, pq);
            }
        }
    }

}  //

void Tensor2f::ltm(const SharedTensor2f &A) {
// Form Lower triangular part
#pragma omp parallel for
    for (int R = 0; R < A->d1_; R++) {
        for (int p = 0; p < A->d2_; p++) {
            for (int q = 0; q <= p; q++) {
                int pq = A->col_idx_[p][q];
                int pq_sym = index2(p, q);
                A2d_[R][pq_sym] = A->get(R, pq);
            }
        }
    }

}  //

void Tensor2f::expand23(int d1, int d2, int d3, const SharedTensor2f &A) {
// Convert Lower triangular to full tensor
#pragma omp parallel for
    for (int p = 0; p < d1; p++) {
        for (int q = 0; q < d2; q++) {
            for (int r = 0; r < d3; r++) {
                int pq = (p * d2) + q;
                int qr = index2(q, r);
                A2d_[pq][r] = A->get(p, qr);
            }
        }
    }

}  //

void Tensor2f::set_row(const SharedTensor2f &A, int n) {
#pragma omp parallel for
    for (int i = 0; i < d3_; i++) {
        for (int j = 0; j < d4_; j++) {
            int ij = col_idx_[i][j];
            A2d_[n][ij] = A->get(i, j);
        }
    }
}  //

void Tensor2f::set_column(const SharedTensor2f &A, int n) {
#pragma omp parallel for
    for (int i = 0; i < d1_; i++) {
        for (int j = 0; j < d2_; j++) {
            int ij = row_idx_[i][j];
            A2d_[ij][n] = A->get(i, j);
        }
    }

    /*
      size_t length;
      length = (size_t)dim1_ * (size_t)dim2_;
      C_SCOPY(length, A->A2d_[0], 1, &(A2d_[0][n]), 1);
    */
}  //

void Tensor2f::get_row(const SharedTensor2f &A, int n) {
#pragma omp parallel for
    for (int i = 0; i < A->d3_; i++) {
        for (int j = 0; j < A->d4_; j++) {
            int ij = A->col_idx_[i][j];
            A2d_[i][j] = A->get(n, ij);
        }
    }
}  //

void Tensor2f::get_column(const SharedTensor2f &A, int n) {
#pragma omp parallel for
    for (int i = 0; i < A->d1_; i++) {
        for (int j = 0; j < A->d2_; j++) {
            int ij = A->row_idx_[i][j];
            A2d_[i][j] = A->get(ij, n);
        }
    }

    /*
      size_t length;
      length = (size_t)dim1_ * (size_t)dim2_;
      C_SCOPY(length, &(A->A2d_[0][n]), 1, A2d_[0], 1);
    */
}  //

void Tensor2f::add2row(const SharedTensor2f &A, int n) {
#pragma omp parallel for
    for (int i = 0; i < d3_; i++) {
        for (int j = 0; j < d4_; j++) {
            int ij = col_idx_[i][j];
            A2d_[n][ij] += A->get(i, j);
        }
    }
}  //

void Tensor2f::add2col(const SharedTensor2f &A, int n) {
#pragma omp parallel for
    for (int i = 0; i < d1_; i++) {
        for (int j = 0; j < d2_; j++) {
            int ij = row_idx_[i][j];
            A2d_[ij][n] += A->get(i, j);
        }
    }
}  //

void Tensor2f::symm4(const SharedTensor2f &a) {
#pragma omp parallel for
    for (int i = 0; i < a->d1_; i++) {
        for (int j = 0; j <= i; j++) {
            int ij = a->row_idx_[i][j];
            int ji = a->row_idx_[j][i];
            int ij2 = index2(i, j);
            for (int k = 0; k < a->d3_; k++) {
                for (int l = 0; l <= k; l++) {
                    int kl = a->col_idx_[k][l];
                    int kl2 = index2(k, l);
                    A2d_[ij2][kl2] = 0.5 * (a->get(ij, kl) + a->get(ji, kl));
                }
            }
        }
    }
}  //

void Tensor2f::symm_col4(const SharedTensor2f &a) {
#pragma omp parallel for
    for (int i = 0; i < a->d1_; i++) {
        for (int j = 0; j <= i; j++) {
            int ij = a->row_idx_[i][j];
            int ij2 = index2(i, j);
            for (int k = 0; k < a->d3_; k++) {
                for (int l = 0; l <= k; l++) {
                    int kl = a->col_idx_[k][l];
                    int lk = a->col_idx_[l][k];
                    int kl2 = index2(k, l);
                    A2d_[ij2][kl2] = 0.5 * (a->get(ij, kl) + a->get(ij, lk));
                }
            }
        }
    }
}  //

void Tensor2f::antisymm4(const SharedTensor2f &a) {
#pragma omp parallel for
    for (int i = 0; i < a->d1_; i++) {
        for (int j = 0; j <= i; j++) {
            int ij = a->row_idx_[i][j];
            int ji = a->row_idx_[j][i];
            int ij2 = index2(i, j);
            for (int k = 0; k < a->d3_; k++) {
                for (int l = 0; l <= k; l++) {
                    int kl = a->col_idx_[k][l];
                    int kl2 = index2(k, l);
                    A2d_[ij2][kl2] = 0.5 * (a->get(ij, kl) - a->get(ji, kl));
                }
            }
        }
    }
}  //

void Tensor2f::antisymm_col4(const SharedTensor2f &a) {
#pragma omp parallel for
    for (int i = 0; i < a->d1_; i++) {
        for (int j = 0; j <= i; j++) {
            int ij = a->row_idx_[i][j];
            int ij2 = index2(i, j);
            for (int k = 0; k < a->d3_; k++) {
                for (int l = 0; l <= k; l++) {
                    int kl = a->col_idx_[k][l];
                    int lk = a->col_idx_[l][k];
                    int kl2 = index2(k, l);
                    A2d_[ij2][kl2] = 0.5 * (a->get(ij, kl) - a->get(ij, lk));
                }
            }
        }
    }
}  //

void Tensor2f::symm_row_packed4(const SharedTensor2f &a) {
#pragma omp parallel for
    for (int i = 0; i < a->d1_; i++) {
        for (int j = 0; j <= i; j++) {
            int ij = a->row_idx_[i][j];
            int ji = a->row_idx_[j][i];
            int ij2 = index2(i, j);
            float perm = (i == j ? 1.0 : 2.0);
            for (int k = 0; k < a->d3_; k++) {
                for (int l = 0; l <= k; l++) {
                    int kl = a->col_idx_[k][l];
                    int kl2 = index2(k, l);
                    A2d_[ij2][kl2] = 0.5 * perm * (a->get(ij, kl) + a->get(ji, kl));
                }
            }
        }
    }
}  //

void Tensor2f::symm_col_packed4(const SharedTensor2f &a) {
#pragma omp parallel for
    for (int i = 0; i < a->d1_; i++) {
        for (int j = 0; j <= i; j++) {
            int ij = a->row_idx_[i][j];
            int ji = a->row_idx_[j][i];
            int ij2 = index2(i, j);
            for (int k = 0; k < a->d3_; k++) {
                for (int l = 0; l <= k; l++) {
                    int kl = a->col_idx_[k][l];
                    int kl2 = index2(k, l);
                    float perm = (k == l ? 1.0 : 2.0);
                    A2d_[ij2][kl2] = 0.5 * perm * (a->get(ij, kl) + a->get(ji, kl));
                }
            }
        }
    }
}  //

void Tensor2f::antisymm_row_packed4(const SharedTensor2f &a) {
#pragma omp parallel for
    for (int i = 0; i < a->d1_; i++) {
        for (int j = 0; j <= i; j++) {
            int ij = a->row_idx_[i][j];
            int ji = a->row_idx_[j][i];
            int ij2 = index2(i, j);
            float perm = (i == j ? 1.0 : 2.0);
            for (int k = 0; k < a->d3_; k++) {
                for (int l = 0; l <= k; l++) {
                    int kl = a->col_idx_[k][l];
                    int kl2 = index2(k, l);
                    A2d_[ij2][kl2] = 0.5 * perm * (a->get(ij, kl) - a->get(ji, kl));
                }
            }
        }
    }
}  //

void Tensor2f::antisymm_col_packed4(const SharedTensor2f &a) {
#pragma omp parallel for
    for (int i = 0; i < a->d1_; i++) {
        for (int j = 0; j <= i; j++) {
            int ij = a->row_idx_[i][j];
            int ji = a->row_idx_[j][i];
            int ij2 = index2(i, j);
            for (int k = 0; k < a->d3_; k++) {
                for (int l = 0; l <= k; l++) {
                    int kl = a->col_idx_[k][l];
                    int kl2 = index2(k, l);
                    float perm = (k == l ? 1.0 : 2.0);
                    A2d_[ij2][kl2] = 0.5 * perm * (a->get(ij, kl) - a->get(ji, kl));
                }
            }
        }
    }
}  //

void Tensor2f::tei_cs1_anti_symm(const SharedTensor2f &J, const SharedTensor2f &K) {
    sort(1243, K, -1.0, 0.0);
    axpy(J, 2.0);
}  //

void Tensor2f::tei_cs2_anti_symm(const SharedTensor2f &J, const SharedTensor2f &K) {
    sort(2134, K, -1.0, 0.0);
    axpy(J, 2.0);
}  //

void Tensor2f::tei_cs3_anti_symm(const SharedTensor2f &J, const SharedTensor2f &K) {
    sort(1432, K, -1.0, 0.0);
    axpy(J, 2.0);
}  //

void Tensor2f::tei_cs4_anti_symm(const SharedTensor2f &J, const SharedTensor2f &K) {
    sort(3214, K, -1.0, 0.0);
    axpy(J, 2.0);
}  //

void Tensor2f::P_ijab(const SharedTensor2f &A) {
    // iajb --> ijab
    sort(1324, A, 1.0, 1.0);

    // jaib --> ijab
    sort(3124, A, -1.0, 1.0);

    // ibja --> ijab
    sort(1342, A, -1.0, 1.0);

    // jbia --> ijab
    sort(3142, A, 1.0, 1.0);

}  //

void Tensor2f::cont444(int t_a1, int t_a2, int f_a1, int f_a2, const SharedTensor2f &A, int t_b1, int t_b2, int f_b1,
                       int f_b2, const SharedTensor2f &B, float alpha, float beta) {
    char ta, tb;
    int nca, ncb, ncc;
    int m, n, k;
    int r1, r2, c1, c2;
    int rr1, rr2, cc1, cc2;
    int dim_t, dim_u;

    // C(pq,rs) = \sum_{tu} A(pq,tu) B(tu,rs)
    if (t_a1 == 3 && t_a2 == 4 && t_b1 == 1 && t_b2 == 2) {
        ta = 'n';
        tb = 'n';
        m = dim1_;
        n = dim2_;
        k = A->dim2_;
        nca = k;  // lda
        ncb = n;  // ldb
        ncc = n;  // ldc

        C_SGEMM(ta, tb, m, n, k, alpha, &(A->A2d_[0][0]), nca, &(B->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
    }

    else {
        // r1
        if (t_a1 == 1) {
            r1 = t_a1;
            dim_t = A->d1_;
        } else if (t_a2 == 1) {
            r1 = t_a2;
            dim_u = A->d1_;
        } else if (f_a1 == 1)
            r1 = f_a1;
        else if (f_a2 == 1)
            r1 = f_a2;

        // r2
        if (t_a1 == 2) {
            r2 = t_a1;
            dim_t = A->d2_;
        } else if (t_a2 == 2) {
            r2 = t_a2;
            dim_u = A->d2_;
        } else if (f_a1 == 2)
            r2 = f_a1;
        else if (f_a2 == 2)
            r2 = f_a2;

        // c1
        if (t_a1 == 3) {
            c1 = t_a1;
            dim_t = A->d3_;
        } else if (t_a2 == 3) {
            c1 = t_a2;
            dim_u = A->d3_;
        } else if (f_a1 == 3)
            c1 = f_a1;
        else if (f_a2 == 3)
            c1 = f_a2;

        // c2
        if (t_a1 == 4) {
            c2 = t_a1;
            dim_t = A->d4_;
        } else if (t_a2 == 4) {
            c2 = t_a2;
            dim_u = A->d4_;
        } else if (f_a1 == 4)
            c2 = f_a1;
        else if (f_a2 == 4)
            c2 = f_a2;

        // outfile->Printf("\tDimensions of A: %2d, %2d, %2d, %2d  \n", r1,r2,c1,c2);

        // Sort A(..,..) to A(pq,tu)
        SharedTensor2f temp1 = SharedTensor2f(new Tensor2f("temp1", d1_, d2_, dim_t, dim_u));
#pragma omp parallel for
        for (int p = 0; p < d1_; p++) {
            for (int q = 0; q < d2_; q++) {
                int pq = temp1->row_idx_[p][q];
                for (int t = 0; t < dim_t; t++) {
                    for (int u = 0; u < dim_u; u++) {
                        int tu = temp1->col_idx_[t][u];

                        if (r1 == f_a1)
                            rr1 = p;
                        else if (r1 == f_a2)
                            rr1 = q;
                        else if (r1 == t_a1)
                            rr1 = t;
                        else if (r1 == t_a2)
                            rr1 = u;

                        if (r2 == f_a1)
                            rr2 = p;
                        else if (r2 == f_a2)
                            rr2 = q;
                        else if (r2 == t_a1)
                            rr2 = t;
                        else if (r2 == t_a2)
                            rr2 = u;

                        if (c1 == f_a1)
                            cc1 = p;
                        else if (c1 == f_a2)
                            cc1 = q;
                        else if (c1 == t_a1)
                            cc1 = t;
                        else if (c1 == t_a2)
                            cc1 = u;

                        if (c2 == f_a1)
                            cc2 = p;
                        else if (c2 == f_a2)
                            cc2 = q;
                        else if (c2 == t_a1)
                            cc2 = t;
                        else if (c2 == t_a2)
                            cc2 = u;

                        int row = rr2 + (rr1 * A->d2_);
                        int col = cc2 + (cc1 * A->d4_);

                        temp1->A2d_[pq][tu] = A->A2d_[row][col];
                    }
                }
            }
        }
        // temp1->print();

        // r1
        if (t_b1 == 1)
            r1 = t_b1;
        else if (t_b2 == 1)
            r1 = t_b2;
        else if (f_b1 == 1)
            r1 = f_b1;
        else if (f_b2 == 1)
            r1 = f_b2;

        // r2
        if (t_b1 == 2)
            r2 = t_b1;
        else if (t_b2 == 2)
            r2 = t_b2;
        else if (f_b1 == 2)
            r2 = f_b1;
        else if (f_b2 == 2)
            r2 = f_b2;

        // c1
        if (t_b1 == 3)
            c1 = t_b1;
        else if (t_b2 == 3)
            c1 = t_b2;
        else if (f_b1 == 3)
            c1 = f_b1;
        else if (f_b2 == 3)
            c1 = f_b2;

        // c2
        if (t_b1 == 4)
            c2 = t_b1;
        else if (t_b2 == 4)
            c2 = t_b2;
        else if (f_b1 == 4)
            c2 = f_b1;
        else if (f_b2 == 4)
            c2 = f_b2;

        // outfile->Printf("\tDimensions of B: %2d, %2d, %2d, %2d  \n", r1,r2,c1,c2);

        // Sort B(..,..) to B(tu,rs)
        SharedTensor2f temp2 = SharedTensor2f(new Tensor2f("temp2", dim_t, dim_u, d3_, d4_));
#pragma omp parallel for
        for (int t = 0; t < dim_t; t++) {
            for (int u = 0; u < dim_u; u++) {
                int tu = temp2->row_idx_[t][u];
                for (int r = 0; r < d3_; r++) {
                    for (int s = 0; s < d4_; s++) {
                        int rs = temp2->col_idx_[r][s];

                        if (r1 == f_b1)
                            rr1 = r;
                        else if (r1 == f_b2)
                            rr1 = s;
                        else if (r1 == t_b1)
                            rr1 = t;
                        else if (r1 == t_b2)
                            rr1 = u;

                        if (r2 == f_b1)
                            rr2 = r;
                        else if (r2 == f_b2)
                            rr2 = s;
                        else if (r2 == t_b1)
                            rr2 = t;
                        else if (r2 == t_b2)
                            rr2 = u;

                        if (c1 == f_b1)
                            cc1 = r;
                        else if (c1 == f_b2)
                            cc1 = s;
                        else if (c1 == t_b1)
                            cc1 = t;
                        else if (c1 == t_b2)
                            cc1 = u;

                        if (c2 == f_b1)
                            cc2 = r;
                        else if (c2 == f_b2)
                            cc2 = s;
                        else if (c2 == t_b1)
                            cc2 = t;
                        else if (c2 == t_b2)
                            cc2 = u;

                        int row = rr2 + (rr1 * B->d2_);
                        int col = cc2 + (cc1 * B->d4_);

                        temp2->A2d_[tu][rs] = B->A2d_[row][col];
                    }
                }
            }
        }
        // temp2->print();

        ta = 'n';
        tb = 'n';
        m = dim1_;
        n = dim2_;
        k = temp1->dim2();
        nca = k;  // lda
        ncb = n;  // ldb
        ncc = n;  // ldc

        C_SGEMM(ta, tb, m, n, k, alpha, &(temp1->A2d_[0][0]), nca, &(temp2->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
        temp1.reset();
        temp2.reset();

    }  // end else

}  //

void Tensor2f::cont444(bool delete_a, int t_a1, int t_a2, int f_a1, int f_a2, SharedTensor2f &A, bool delete_b,
                       int t_b1, int t_b2, int f_b1, int f_b2, SharedTensor2f &B, float alpha, float beta) {
    char ta, tb;
    int nca, ncb, ncc;
    int m, n, k;
    int r1, r2, c1, c2;
    int rr1, rr2, cc1, cc2;
    int dim_t, dim_u;

    // C(pq,rs) = \sum_{tu} A(pq,tu) B(tu,rs)
    if (t_a1 == 3 && t_a2 == 4 && t_b1 == 1 && t_b2 == 2) {
        ta = 'n';
        tb = 'n';
        m = dim1_;
        n = dim2_;
        k = A->dim2_;
        nca = k;  // lda
        ncb = n;  // ldb
        ncc = n;  // ldc

        C_SGEMM(ta, tb, m, n, k, alpha, &(A->A2d_[0][0]), nca, &(B->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
    }

    else {
        // r1
        if (t_a1 == 1) {
            r1 = t_a1;
            dim_t = A->d1_;
        } else if (t_a2 == 1) {
            r1 = t_a2;
            dim_u = A->d1_;
        } else if (f_a1 == 1)
            r1 = f_a1;
        else if (f_a2 == 1)
            r1 = f_a2;

        // r2
        if (t_a1 == 2) {
            r2 = t_a1;
            dim_t = A->d2_;
        } else if (t_a2 == 2) {
            r2 = t_a2;
            dim_u = A->d2_;
        } else if (f_a1 == 2)
            r2 = f_a1;
        else if (f_a2 == 2)
            r2 = f_a2;

        // c1
        if (t_a1 == 3) {
            c1 = t_a1;
            dim_t = A->d3_;
        } else if (t_a2 == 3) {
            c1 = t_a2;
            dim_u = A->d3_;
        } else if (f_a1 == 3)
            c1 = f_a1;
        else if (f_a2 == 3)
            c1 = f_a2;

        // c2
        if (t_a1 == 4) {
            c2 = t_a1;
            dim_t = A->d4_;
        } else if (t_a2 == 4) {
            c2 = t_a2;
            dim_u = A->d4_;
        } else if (f_a1 == 4)
            c2 = f_a1;
        else if (f_a2 == 4)
            c2 = f_a2;

        // outfile->Printf("\tDimensions of A: %2d, %2d, %2d, %2d  \n", r1,r2,c1,c2);

        // Sort A(..,..) to A(pq,tu)
        SharedTensor2f temp1 = SharedTensor2f(new Tensor2f("temp1", d1_, d2_, dim_t, dim_u));
#pragma omp parallel for
        for (int p = 0; p < d1_; p++) {
            for (int q = 0; q < d2_; q++) {
                int pq = temp1->row_idx_[p][q];
                for (int t = 0; t < dim_t; t++) {
                    for (int u = 0; u < dim_u; u++) {
                        int tu = temp1->col_idx_[t][u];

                        if (r1 == f_a1)
                            rr1 = p;
                        else if (r1 == f_a2)
                            rr1 = q;
                        else if (r1 == t_a1)
                            rr1 = t;
                        else if (r1 == t_a2)
                            rr1 = u;

                        if (r2 == f_a1)
                            rr2 = p;
                        else if (r2 == f_a2)
                            rr2 = q;
                        else if (r2 == t_a1)
                            rr2 = t;
                        else if (r2 == t_a2)
                            rr2 = u;

                        if (c1 == f_a1)
                            cc1 = p;
                        else if (c1 == f_a2)
                            cc1 = q;
                        else if (c1 == t_a1)
                            cc1 = t;
                        else if (c1 == t_a2)
                            cc1 = u;

                        if (c2 == f_a1)
                            cc2 = p;
                        else if (c2 == f_a2)
                            cc2 = q;
                        else if (c2 == t_a1)
                            cc2 = t;
                        else if (c2 == t_a2)
                            cc2 = u;

                        int row = rr2 + (rr1 * A->d2_);
                        int col = cc2 + (cc1 * A->d4_);

                        temp1->A2d_[pq][tu] = A->A2d_[row][col];
                    }
                }
            }
        }
        // temp1->print();
        if (delete_a) A.reset();

        // r1
        if (t_b1 == 1)
            r1 = t_b1;
        else if (t_b2 == 1)
            r1 = t_b2;
        else if (f_b1 == 1)
            r1 = f_b1;
        else if (f_b2 == 1)
            r1 = f_b2;

        // r2
        if (t_b1 == 2)
            r2 = t_b1;
        else if (t_b2 == 2)
            r2 = t_b2;
        else if (f_b1 == 2)
            r2 = f_b1;
        else if (f_b2 == 2)
            r2 = f_b2;

        // c1
        if (t_b1 == 3)
            c1 = t_b1;
        else if (t_b2 == 3)
            c1 = t_b2;
        else if (f_b1 == 3)
            c1 = f_b1;
        else if (f_b2 == 3)
            c1 = f_b2;

        // c2
        if (t_b1 == 4)
            c2 = t_b1;
        else if (t_b2 == 4)
            c2 = t_b2;
        else if (f_b1 == 4)
            c2 = f_b1;
        else if (f_b2 == 4)
            c2 = f_b2;

        // outfile->Printf("\tDimensions of B: %2d, %2d, %2d, %2d  \n", r1,r2,c1,c2);

        // Sort B(..,..) to B(tu,rs)
        SharedTensor2f temp2 = SharedTensor2f(new Tensor2f("temp2", dim_t, dim_u, d3_, d4_));
#pragma omp parallel for
        for (int t = 0; t < dim_t; t++) {
            for (int u = 0; u < dim_u; u++) {
                int tu = temp2->row_idx_[t][u];
                for (int r = 0; r < d3_; r++) {
                    for (int s = 0; s < d4_; s++) {
                        int rs = temp2->col_idx_[r][s];

                        if (r1 == f_b1)
                            rr1 = r;
                        else if (r1 == f_b2)
                            rr1 = s;
                        else if (r1 == t_b1)
                            rr1 = t;
                        else if (r1 == t_b2)
                            rr1 = u;

                        if (r2 == f_b1)
                            rr2 = r;
                        else if (r2 == f_b2)
                            rr2 = s;
                        else if (r2 == t_b1)
                            rr2 = t;
                        else if (r2 == t_b2)
                            rr2 = u;

                        if (c1 == f_b1)
                            cc1 = r;
                        else if (c1 == f_b2)
                            cc1 = s;
                        else if (c1 == t_b1)
                            cc1 = t;
                        else if (c1 == t_b2)
                            cc1 = u;

                        if (c2 == f_b1)
                            cc2 = r;
                        else if (c2 == f_b2)
                            cc2 = s;
                        else if (c2 == t_b1)
                            cc2 = t;
                        else if (c2 == t_b2)
                            cc2 = u;

                        int row = rr2 + (rr1 * B->d2_);
                        int col = cc2 + (cc1 * B->d4_);

                        temp2->A2d_[tu][rs] = B->A2d_[row][col];
                    }
                }
            }
        }
        // temp2->print();
        if (delete_b) B.reset();

        ta = 'n';
        tb = 'n';
        m = dim1_;
        n = dim2_;
        k = temp1->dim2();
        nca = k;  // lda
        ncb = n;  // ldb
        ncc = n;  // ldc

        C_SGEMM(ta, tb, m, n, k, alpha, &(temp1->A2d_[0][0]), nca, &(temp2->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
        temp1.reset();
        temp2.reset();

    }  // end else

}  //

void Tensor2f::cont444(std::string idx_c, std::string idx_a, std::string idx_b, bool delete_a, bool delete_b,
                       SharedTensor2f &A, SharedTensor2f &B, float alpha, float beta) {
    char ta, tb;
    int nca, ncb, ncc;
    int m, n, k;
    int r1, r2, c1, c2;
    int rr1, rr2, cc1, cc2;
    int dim_t, dim_u;
    int t_a1, t_a2, f_a1, f_a2;
    int t_b1, t_b2, f_b1, f_b2;

    // Expected order: C(pq,rs) = \sum_{tu} A(pq,tu) B(tu,rs)
    /*
    for (int i = 0; i < idx_a.size(); i++) {
         outfile->Printf("\tIDX_A[%1d]: %c \n", i, idx_a[i]);
    }
    */

    // Find dummy & free indices for A
    // f_a1
    if (idx_a[0] == idx_c[0])
        f_a1 = 1;
    else if (idx_a[1] == idx_c[0])
        f_a1 = 2;
    else if (idx_a[2] == idx_c[0])
        f_a1 = 3;
    else if (idx_a[3] == idx_c[0])
        f_a1 = 4;

    // f_a2
    if (idx_a[0] == idx_c[1])
        f_a2 = 1;
    else if (idx_a[1] == idx_c[1])
        f_a2 = 2;
    else if (idx_a[2] == idx_c[1])
        f_a2 = 3;
    else if (idx_a[3] == idx_c[1])
        f_a2 = 4;

    // t_a1
    if (f_a1 == 1 && f_a2 == 2) {
        t_a1 = 3;
        t_a2 = 4;
    } else if (f_a1 == 1 && f_a2 == 3) {
        t_a1 = 2;
        t_a2 = 4;
    } else if (f_a1 == 1 && f_a2 == 4) {
        t_a1 = 2;
        t_a2 = 3;
    } else if (f_a1 == 2 && f_a2 == 3) {
        t_a1 = 1;
        t_a2 = 4;
    } else if (f_a1 == 2 && f_a2 == 4) {
        t_a1 = 1;
        t_a2 = 3;
    } else if (f_a1 == 3 && f_a2 == 4) {
        t_a1 = 1;
        t_a2 = 2;
    } else if (f_a1 == 2 && f_a2 == 1) {
        t_a1 = 3;
        t_a2 = 4;
    } else if (f_a1 == 3 && f_a2 == 1) {
        t_a1 = 2;
        t_a2 = 4;
    } else if (f_a1 == 4 && f_a2 == 1) {
        t_a1 = 2;
        t_a2 = 3;
    } else if (f_a1 == 3 && f_a2 == 2) {
        t_a1 = 1;
        t_a2 = 4;
    } else if (f_a1 == 4 && f_a2 == 2) {
        t_a1 = 1;
        t_a2 = 3;
    } else if (f_a1 == 4 && f_a2 == 3) {
        t_a1 = 1;
        t_a2 = 2;
    }
    // outfile->Printf("\tf_a1, f_a2, t_a1, t_a2: %1d, %1d, %1d, %1d  \n", f_a1,f_a2,t_a1,t_a2);

    // Find dummy & free indices for B
    // f_b1
    if (idx_b[0] == idx_c[2])
        f_b1 = 1;
    else if (idx_b[1] == idx_c[2])
        f_b1 = 2;
    else if (idx_b[2] == idx_c[2])
        f_b1 = 3;
    else if (idx_b[3] == idx_c[2])
        f_b1 = 4;

    // f_b2
    if (idx_b[0] == idx_c[3])
        f_b2 = 1;
    else if (idx_b[1] == idx_c[3])
        f_b2 = 2;
    else if (idx_b[2] == idx_c[3])
        f_b2 = 3;
    else if (idx_b[3] == idx_c[3])
        f_b2 = 4;

    // t_b1
    if (idx_b[0] == idx_a[t_a1 - 1])
        t_b1 = 1;
    else if (idx_b[1] == idx_a[t_a1 - 1])
        t_b1 = 2;
    else if (idx_b[2] == idx_a[t_a1 - 1])
        t_b1 = 3;
    else if (idx_b[3] == idx_a[t_a1 - 1])
        t_b1 = 4;

    // t_b2
    if (idx_b[0] == idx_a[t_a2 - 1])
        t_b2 = 1;
    else if (idx_b[1] == idx_a[t_a2 - 1])
        t_b2 = 2;
    else if (idx_b[2] == idx_a[t_a2 - 1])
        t_b2 = 3;
    else if (idx_b[3] == idx_a[t_a2 - 1])
        t_b2 = 4;
    // outfile->Printf("\tf_b1, f_b2, t_b1, t_b2: %1d, %1d, %1d, %1d  \n", f_b1,f_b2,t_b1,t_b2);

    // Figure out A
    // r1
    if (t_a1 == 1) {
        r1 = t_a1;
        dim_t = A->d1_;
    } else if (t_a2 == 1) {
        r1 = t_a2;
        dim_u = A->d1_;
    } else if (f_a1 == 1)
        r1 = f_a1;
    else if (f_a2 == 1)
        r1 = f_a2;

    // r2
    if (t_a1 == 2) {
        r2 = t_a1;
        dim_t = A->d2_;
    } else if (t_a2 == 2) {
        r2 = t_a2;
        dim_u = A->d2_;
    } else if (f_a1 == 2)
        r2 = f_a1;
    else if (f_a2 == 2)
        r2 = f_a2;

    // c1
    if (t_a1 == 3) {
        c1 = t_a1;
        dim_t = A->d3_;
    } else if (t_a2 == 3) {
        c1 = t_a2;
        dim_u = A->d3_;
    } else if (f_a1 == 3)
        c1 = f_a1;
    else if (f_a2 == 3)
        c1 = f_a2;

    // c2
    if (t_a1 == 4) {
        c2 = t_a1;
        dim_t = A->d4_;
    } else if (t_a2 == 4) {
        c2 = t_a2;
        dim_u = A->d4_;
    } else if (f_a1 == 4)
        c2 = f_a1;
    else if (f_a2 == 4)
        c2 = f_a2;

    // outfile->Printf("\tDimensions of A: %2d, %2d, %2d, %2d  \n", r1,r2,c1,c2);

    // Sort A(..,..) to A(pq,tu)
    SharedTensor2f temp1 = SharedTensor2f(new Tensor2f("temp1", d1_, d2_, dim_t, dim_u));
#pragma omp parallel for
    for (int p = 0; p < d1_; p++) {
        for (int q = 0; q < d2_; q++) {
            int pq = temp1->row_idx_[p][q];
            for (int t = 0; t < dim_t; t++) {
                for (int u = 0; u < dim_u; u++) {
                    int tu = temp1->col_idx_[t][u];

                    if (r1 == f_a1)
                        rr1 = p;
                    else if (r1 == f_a2)
                        rr1 = q;
                    else if (r1 == t_a1)
                        rr1 = t;
                    else if (r1 == t_a2)
                        rr1 = u;

                    if (r2 == f_a1)
                        rr2 = p;
                    else if (r2 == f_a2)
                        rr2 = q;
                    else if (r2 == t_a1)
                        rr2 = t;
                    else if (r2 == t_a2)
                        rr2 = u;

                    if (c1 == f_a1)
                        cc1 = p;
                    else if (c1 == f_a2)
                        cc1 = q;
                    else if (c1 == t_a1)
                        cc1 = t;
                    else if (c1 == t_a2)
                        cc1 = u;

                    if (c2 == f_a1)
                        cc2 = p;
                    else if (c2 == f_a2)
                        cc2 = q;
                    else if (c2 == t_a1)
                        cc2 = t;
                    else if (c2 == t_a2)
                        cc2 = u;

                    int row = rr2 + (rr1 * A->d2_);
                    int col = cc2 + (cc1 * A->d4_);

                    temp1->A2d_[pq][tu] = A->A2d_[row][col];
                }
            }
        }
    }
    // temp1->print();
    if (delete_a) A.reset();

    // Figure out B
    // r1
    if (t_b1 == 1)
        r1 = t_b1;
    else if (t_b2 == 1)
        r1 = t_b2;
    else if (f_b1 == 1)
        r1 = f_b1;
    else if (f_b2 == 1)
        r1 = f_b2;

    // r2
    if (t_b1 == 2)
        r2 = t_b1;
    else if (t_b2 == 2)
        r2 = t_b2;
    else if (f_b1 == 2)
        r2 = f_b1;
    else if (f_b2 == 2)
        r2 = f_b2;

    // c1
    if (t_b1 == 3)
        c1 = t_b1;
    else if (t_b2 == 3)
        c1 = t_b2;
    else if (f_b1 == 3)
        c1 = f_b1;
    else if (f_b2 == 3)
        c1 = f_b2;

    // c2
    if (t_b1 == 4)
        c2 = t_b1;
    else if (t_b2 == 4)
        c2 = t_b2;
    else if (f_b1 == 4)
        c2 = f_b1;
    else if (f_b2 == 4)
        c2 = f_b2;

    // outfile->Printf("\tDimensions of B: %2d, %2d, %2d, %2d  \n", r1,r2,c1,c2);

    // Sort B(..,..) to B(tu,rs)
    SharedTensor2f temp2 = SharedTensor2f(new Tensor2f("temp2", dim_t, dim_u, d3_, d4_));
#pragma omp parallel for
    for (int t = 0; t < dim_t; t++) {
        for (int u = 0; u < dim_u; u++) {
            int tu = temp2->row_idx_[t][u];
            for (int r = 0; r < d3_; r++) {
                for (int s = 0; s < d4_; s++) {
                    int rs = temp2->col_idx_[r][s];

                    if (r1 == f_b1)
                        rr1 = r;
                    else if (r1 == f_b2)
                        rr1 = s;
                    else if (r1 == t_b1)
                        rr1 = t;
                    else if (r1 == t_b2)
                        rr1 = u;

                    if (r2 == f_b1)
                        rr2 = r;
                    else if (r2 == f_b2)
                        rr2 = s;
                    else if (r2 == t_b1)
                        rr2 = t;
                    else if (r2 == t_b2)
                        rr2 = u;

                    if (c1 == f_b1)
                        cc1 = r;
                    else if (c1 == f_b2)
                        cc1 = s;
                    else if (c1 == t_b1)
                        cc1 = t;
                    else if (c1 == t_b2)
                        cc1 = u;

                    if (c2 == f_b1)
                        cc2 = r;
                    else if (c2 == f_b2)
                        cc2 = s;
                    else if (c2 == t_b1)
                        cc2 = t;
                    else if (c2 == t_b2)
                        cc2 = u;

                    int row = rr2 + (rr1 * B->d2_);
                    int col = cc2 + (cc1 * B->d4_);

                    temp2->A2d_[tu][rs] = B->A2d_[row][col];
                }
            }
        }
    }
    // temp2->print();
    if (delete_b) B.reset();

    ta = 'n';
    tb = 'n';
    m = dim1_;
    n = dim2_;
    k = temp1->dim2();
    nca = k;  // lda
    ncb = n;  // ldb
    ncc = n;  // ldc

    C_SGEMM(ta, tb, m, n, k, alpha, &(temp1->A2d_[0][0]), nca, &(temp2->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
    temp1.reset();
    temp2.reset();

}  //

void Tensor2f::cont444(std::string idx_c, std::string idx_a, std::string idx_b, SharedTensor2f &A, SharedTensor2f &B,
                       float alpha, float beta) {
    char ta, tb;
    int nca, ncb, ncc;
    int m, n, k;
    int r1, r2, c1, c2;
    int rr1, rr2, cc1, cc2;
    int dim_t, dim_u;
    int t_a1, t_a2, f_a1, f_a2;
    int t_b1, t_b2, f_b1, f_b2;
    int d1_a, d2_a, d3_a, d4_a;  // Dimensions of sorted A tensor
    int d1_b, d2_b, d3_b, d4_b;
    int sort_a, sort_b;

    // Find free indices (pq) for A
    // f_a1
    if (idx_a[0] == idx_c[0]) {
        f_a1 = 1;
        d1_a = A->d1_;
    } else if (idx_a[1] == idx_c[0]) {
        f_a1 = 2;
        d1_a = A->d2_;
    } else if (idx_a[2] == idx_c[0]) {
        f_a1 = 3;
        d1_a = A->d3_;
    } else if (idx_a[3] == idx_c[0]) {
        f_a1 = 4;
        d1_a = A->d4_;
    }

    // f_a2
    if (idx_a[0] == idx_c[1]) {
        f_a2 = 1;
        d2_a = A->d1_;
    } else if (idx_a[1] == idx_c[1]) {
        f_a2 = 2;
        d2_a = A->d2_;
    } else if (idx_a[2] == idx_c[1]) {
        f_a2 = 3;
        d2_a = A->d3_;
    } else if (idx_a[3] == idx_c[1]) {
        f_a2 = 4;
        d2_a = A->d4_;
    }

    // Find target indices (tu) for A
    // t_a1
    if (f_a1 == 1 && f_a2 == 2) {
        t_a1 = 3;
        t_a2 = 4;
        d3_a = A->d3_;
        d4_a = A->d4_;
    } else if (f_a1 == 1 && f_a2 == 3) {
        t_a1 = 2;
        t_a2 = 4;
        d3_a = A->d2_;
        d4_a = A->d4_;
    } else if (f_a1 == 1 && f_a2 == 4) {
        t_a1 = 2;
        t_a2 = 3;
        d3_a = A->d2_;
        d4_a = A->d3_;
    } else if (f_a1 == 2 && f_a2 == 3) {
        t_a1 = 1;
        t_a2 = 4;
        d3_a = A->d1_;
        d4_a = A->d4_;
    } else if (f_a1 == 2 && f_a2 == 4) {
        t_a1 = 1;
        t_a2 = 3;
        d3_a = A->d1_;
        d4_a = A->d3_;
    } else if (f_a1 == 3 && f_a2 == 4) {
        t_a1 = 1;
        t_a2 = 2;
        d3_a = A->d1_;
        d4_a = A->d2_;
    } else if (f_a1 == 2 && f_a2 == 1) {
        t_a1 = 3;
        t_a2 = 4;
        d3_a = A->d3_;
        d4_a = A->d4_;
    } else if (f_a1 == 3 && f_a2 == 1) {
        t_a1 = 2;
        t_a2 = 4;
        d3_a = A->d2_;
        d4_a = A->d4_;
    } else if (f_a1 == 4 && f_a2 == 1) {
        t_a1 = 2;
        t_a2 = 3;
        d3_a = A->d2_;
        d4_a = A->d3_;
    } else if (f_a1 == 3 && f_a2 == 2) {
        t_a1 = 1;
        t_a2 = 4;
        d3_a = A->d1_;
        d4_a = A->d4_;
    } else if (f_a1 == 4 && f_a2 == 2) {
        t_a1 = 1;
        t_a2 = 3;
        d3_a = A->d1_;
        d4_a = A->d3_;
    } else if (f_a1 == 4 && f_a2 == 3) {
        t_a1 = 1;
        t_a2 = 2;
        d3_a = A->d1_;
        d4_a = A->d2_;
    }
    // outfile->Printf("\tf_a1, f_a2, t_a1, t_a2: %1d, %1d, %1d, %1d  \n", f_a1,f_a2,t_a1,t_a2);

    // Sort A(..,..) to A(pq,tu)
    sort_a = (f_a1 * 1000) + (f_a2 * 100) + (t_a1 * 10) + t_a2;
    SharedTensor2f temp1 = SharedTensor2f(new Tensor2f("temp1", d1_a, d2_a, d3_a, d4_a));
    temp1->sort(sort_a, A, 1.0, 0.0);
    A.reset();
    // temp1->print();

    // Find target (tu) & free (rs) indices for B
    // f_b1
    if (idx_b[0] == idx_c[2]) {
        f_b1 = 1;
        d3_b = B->d1_;
    } else if (idx_b[1] == idx_c[2]) {
        f_b1 = 2;
        d3_b = B->d2_;
    } else if (idx_b[2] == idx_c[2]) {
        f_b1 = 3;
        d3_b = B->d3_;
    } else if (idx_b[3] == idx_c[2]) {
        f_b1 = 4;
        d3_b = B->d4_;
    }

    // f_b2
    if (idx_b[0] == idx_c[3]) {
        f_b2 = 1;
        d4_b = B->d1_;
    } else if (idx_b[1] == idx_c[3]) {
        f_b2 = 2;
        d4_b = B->d2_;
    } else if (idx_b[2] == idx_c[3]) {
        f_b2 = 3;
        d4_b = B->d3_;
    } else if (idx_b[3] == idx_c[3]) {
        f_b2 = 4;
        d4_b = B->d4_;
    }

    // t_b1
    if (idx_b[0] == idx_a[t_a1 - 1]) {
        t_b1 = 1;
        d1_b = B->d1_;
    } else if (idx_b[1] == idx_a[t_a1 - 1]) {
        t_b1 = 2;
        d1_b = B->d2_;
    } else if (idx_b[2] == idx_a[t_a1 - 1]) {
        t_b1 = 3;
        d1_b = B->d3_;
    } else if (idx_b[3] == idx_a[t_a1 - 1]) {
        t_b1 = 4;
        d1_b = B->d4_;
    }

    // t_b2
    if (idx_b[0] == idx_a[t_a2 - 1]) {
        t_b2 = 1;
        d2_b = B->d1_;
    } else if (idx_b[1] == idx_a[t_a2 - 1]) {
        t_b2 = 2;
        d2_b = B->d2_;
    } else if (idx_b[2] == idx_a[t_a2 - 1]) {
        t_b2 = 3;
        d2_b = B->d3_;
    } else if (idx_b[3] == idx_a[t_a2 - 1]) {
        t_b2 = 4;
        d2_b = B->d4_;
    }
    // outfile->Printf("\tf_b1, f_b2, t_b1, t_b2: %1d, %1d, %1d, %1d  \n", f_b1,f_b2,t_b1,t_b2);

    // Sort B(..,..) to B(tu,rs)
    sort_b = (t_b1 * 1000) + (t_b2 * 100) + (f_b1 * 10) + f_b2;
    SharedTensor2f temp2 = SharedTensor2f(new Tensor2f("temp2", d1_b, d2_b, d3_b, d4_b));
    temp2->sort(sort_b, B, 1.0, 0.0);
    B.reset();
    // temp2->print();

    ta = 'n';
    tb = 'n';
    m = dim1_;
    n = dim2_;
    k = temp1->dim2();
    nca = k;  // lda
    ncb = n;  // ldb
    ncc = n;  // ldc

    C_SGEMM(ta, tb, m, n, k, alpha, &(temp1->A2d_[0][0]), nca, &(temp2->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
    temp1.reset();
    temp2.reset();
}  //

void Tensor2f::cont343(std::string idx_c, std::string idx_a, std::string idx_b, bool delete_b, SharedTensor2f &A,
                       SharedTensor2f &B, float alpha, float beta) {
    char ta, tb;
    int nca, ncb, ncc;
    int m, n, k;
    int r1, r2, c1, c2;
    int rr1, rr2, cc1, cc2;
    int t_b1, t_b2, f_b1, f_b2;

    // Expected order: C(Q,pq) = \sum_{rs} A(Q,rs) B(rs,pq)

    // Find dummy & free indices for B
    // f_b1
    if (idx_b[0] == idx_c[0])
        f_b1 = 1;
    else if (idx_b[1] == idx_c[0])
        f_b1 = 2;
    else if (idx_b[2] == idx_c[0])
        f_b1 = 3;
    else if (idx_b[3] == idx_c[0])
        f_b1 = 4;

    // f_b2
    if (idx_b[0] == idx_c[1])
        f_b2 = 1;
    else if (idx_b[1] == idx_c[1])
        f_b2 = 2;
    else if (idx_b[2] == idx_c[1])
        f_b2 = 3;
    else if (idx_b[3] == idx_c[1])
        f_b2 = 4;

    // t_b1
    if (idx_b[0] == idx_a[0])
        t_b1 = 1;
    else if (idx_b[1] == idx_a[0])
        t_b1 = 2;
    else if (idx_b[2] == idx_a[0])
        t_b1 = 3;
    else if (idx_b[3] == idx_a[0])
        t_b1 = 4;

    // t_b2
    if (idx_b[0] == idx_a[1])
        t_b2 = 1;
    else if (idx_b[1] == idx_a[1])
        t_b2 = 2;
    else if (idx_b[2] == idx_a[1])
        t_b2 = 3;
    else if (idx_b[3] == idx_a[1])
        t_b2 = 4;

    // Figure out B
    // r1
    if (t_b1 == 1)
        r1 = t_b1;
    else if (t_b2 == 1)
        r1 = t_b2;
    else if (f_b1 == 1)
        r1 = f_b1;
    else if (f_b2 == 1)
        r1 = f_b2;

    // r2
    if (t_b1 == 2)
        r2 = t_b1;
    else if (t_b2 == 2)
        r2 = t_b2;
    else if (f_b1 == 2)
        r2 = f_b1;
    else if (f_b2 == 2)
        r2 = f_b2;

    // c1
    if (t_b1 == 3)
        c1 = t_b1;
    else if (t_b2 == 3)
        c1 = t_b2;
    else if (f_b1 == 3)
        c1 = f_b1;
    else if (f_b2 == 3)
        c1 = f_b2;

    // c2
    if (t_b1 == 4)
        c2 = t_b1;
    else if (t_b2 == 4)
        c2 = t_b2;
    else if (f_b1 == 4)
        c2 = f_b1;
    else if (f_b2 == 4)
        c2 = f_b2;

    // outfile->Printf("\tDimensions of B: %2d, %2d, %2d, %2d  \n", r1,r2,c1,c2);

    // Sort B(..,..) to B(rs,pq)
    SharedTensor2f temp = SharedTensor2f(new Tensor2f("temp", A->d2_, A->d3_, d2_, d3_));
#pragma omp parallel for
    for (int r = 0; r < A->d2_; r++) {
        for (int s = 0; s < A->d3_; s++) {
            int rs = temp->row_idx_[r][s];
            for (int p = 0; p < d2_; p++) {
                for (int q = 0; q < d3_; q++) {
                    int pq = temp->col_idx_[p][q];

                    if (r1 == f_b1)
                        rr1 = p;
                    else if (r1 == f_b2)
                        rr1 = q;
                    else if (r1 == t_b1)
                        rr1 = r;
                    else if (r1 == t_b2)
                        rr1 = s;

                    if (r2 == f_b1)
                        rr2 = p;
                    else if (r2 == f_b2)
                        rr2 = q;
                    else if (r2 == t_b1)
                        rr2 = r;
                    else if (r2 == t_b2)
                        rr2 = s;

                    if (c1 == f_b1)
                        cc1 = p;
                    else if (c1 == f_b2)
                        cc1 = q;
                    else if (c1 == t_b1)
                        cc1 = r;
                    else if (c1 == t_b2)
                        cc1 = s;

                    if (c2 == f_b1)
                        cc2 = p;
                    else if (c2 == f_b2)
                        cc2 = q;
                    else if (c2 == t_b1)
                        cc2 = r;
                    else if (c2 == t_b2)
                        cc2 = s;

                    int row = rr2 + (rr1 * B->d2_);
                    int col = cc2 + (cc1 * B->d4_);

                    temp->A2d_[rs][pq] = B->A2d_[row][col];
                }
            }
        }
    }
    // temp->print();
    if (delete_b) B.reset();

    ta = 'n';
    tb = 'n';
    m = d1_;
    n = d2_ * d3_;
    k = temp->dim1();
    nca = k;  // lda
    ncb = n;  // ldb
    ncc = n;  // ldc

    C_SGEMM(ta, tb, m, n, k, alpha, &(A->A2d_[0][0]), nca, &(temp->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
    temp.reset();

}  //

void Tensor2f::cont442(std::string idx_c, std::string idx_a, std::string idx_b, bool delete_a, bool delete_b,
                       SharedTensor2f &A, SharedTensor2f &B, float alpha, float beta) {
    char ta, tb;
    int nca, ncb, ncc;
    int m, n, k;
    int r1, r2, c1, c2;
    int rr1, rr2, cc1, cc2;
    int dim_r, dim_s, dim_t;
    int t_a1, t_a2, t_a3, f_a1;
    int t_b1, t_b2, t_b3, f_b1;

    // Expected order: C(pq) = \sum_{rst} A(pr,st) B(rs,tq)

    // Find dummy & free indices for A
    // f_a1
    if (idx_a[0] == idx_c[0])
        f_a1 = 1;
    else if (idx_a[1] == idx_c[0])
        f_a1 = 2;
    else if (idx_a[2] == idx_c[0])
        f_a1 = 3;
    else if (idx_a[3] == idx_c[0])
        f_a1 = 4;

    // t_a1, t_a2, t_a3
    if (f_a1 == 1) {
        t_a1 = 2;
        t_a2 = 3;
        t_a3 = 4;
    } else if (f_a1 == 2) {
        t_a1 = 1;
        t_a2 = 3;
        t_a3 = 4;
    } else if (f_a1 == 3) {
        t_a1 = 1;
        t_a2 = 2;
        t_a3 = 4;
    } else if (f_a1 == 4) {
        t_a1 = 1;
        t_a2 = 2;
        t_a3 = 3;
    }

    // Find dummy & free indices for B
    // f_b1
    if (idx_b[0] == idx_c[1])
        f_b1 = 1;
    else if (idx_b[1] == idx_c[1])
        f_b1 = 2;
    else if (idx_b[2] == idx_c[1])
        f_b1 = 3;
    else if (idx_b[3] == idx_c[1])
        f_b1 = 4;

    // t_b1
    if (idx_b[0] == idx_a[t_a1 - 1])
        t_b1 = 1;
    else if (idx_b[1] == idx_a[t_a1 - 1])
        t_b1 = 2;
    else if (idx_b[2] == idx_a[t_a1 - 1])
        t_b1 = 3;
    else if (idx_b[3] == idx_a[t_a1 - 1])
        t_b1 = 4;

    // t_b2
    if (idx_b[0] == idx_a[t_a2 - 1])
        t_b2 = 1;
    else if (idx_b[1] == idx_a[t_a2 - 1])
        t_b2 = 2;
    else if (idx_b[2] == idx_a[t_a2 - 1])
        t_b2 = 3;
    else if (idx_b[3] == idx_a[t_a2 - 1])
        t_b2 = 4;

    // t_b3
    if (idx_b[0] == idx_a[t_a3 - 1])
        t_b3 = 1;
    else if (idx_b[1] == idx_a[t_a3 - 1])
        t_b3 = 2;
    else if (idx_b[2] == idx_a[t_a3 - 1])
        t_b3 = 3;
    else if (idx_b[3] == idx_a[t_a3 - 1])
        t_b3 = 4;

    // Figure out A
    // r1
    if (t_a1 == 1) {
        r1 = t_a1;
        dim_r = A->d1_;
    } else if (t_a2 == 1) {
        r1 = t_a2;
        dim_s = A->d1_;
    } else if (t_a3 == 1) {
        r1 = t_a3;
        dim_t = A->d1_;
    } else if (f_a1 == 1)
        r1 = f_a1;

    // r2
    if (t_a1 == 2) {
        r2 = t_a1;
        dim_r = A->d2_;
    } else if (t_a2 == 2) {
        r2 = t_a2;
        dim_s = A->d2_;
    } else if (t_a3 == 2) {
        r2 = t_a3;
        dim_t = A->d2_;
    } else if (f_a1 == 2)
        r2 = f_a1;

    // c1
    if (t_a1 == 3) {
        c1 = t_a1;
        dim_r = A->d3_;
    } else if (t_a2 == 3) {
        c1 = t_a2;
        dim_s = A->d3_;
    } else if (t_a3 == 3) {
        c1 = t_a3;
        dim_t = A->d3_;
    } else if (f_a1 == 3)
        c1 = f_a1;

    // c2
    if (t_a1 == 4) {
        c2 = t_a1;
        dim_r = A->d4_;
    } else if (t_a2 == 4) {
        c2 = t_a2;
        dim_s = A->d4_;
    } else if (t_a3 == 4) {
        c2 = t_a3;
        dim_t = A->d4_;
    } else if (f_a1 == 3)
        c2 = f_a1;

    // Sort A(..,..) to A(pr,st)
    SharedTensor2f temp1 = SharedTensor2f(new Tensor2f("temp1", dim1_, dim_r, dim_s, dim_t));
#pragma omp parallel for
    for (int p = 0; p < dim1_; p++) {
        for (int r = 0; r < dim_r; r++) {
            int pr = temp1->row_idx_[p][r];
            for (int s = 0; s < dim_s; s++) {
                for (int t = 0; t < dim_t; t++) {
                    int st = temp1->col_idx_[s][t];

                    if (r1 == f_a1)
                        rr1 = p;
                    else if (r1 == t_a1)
                        rr1 = r;
                    else if (r1 == t_a2)
                        rr1 = s;
                    else if (r1 == t_a3)
                        rr1 = t;

                    if (r2 == f_a1)
                        rr2 = p;
                    else if (r2 == t_a1)
                        rr2 = r;
                    else if (r2 == t_a2)
                        rr2 = s;
                    else if (r2 == t_a3)
                        rr2 = t;

                    if (c1 == f_a1)
                        cc1 = p;
                    else if (c1 == t_a1)
                        cc1 = r;
                    else if (c1 == t_a2)
                        cc1 = s;
                    else if (c1 == t_a3)
                        cc1 = t;

                    if (c2 == f_a1)
                        cc2 = p;
                    else if (c2 == t_a1)
                        cc2 = r;
                    else if (c2 == t_a2)
                        cc2 = s;
                    else if (c2 == t_a3)
                        cc2 = t;

                    int row = rr2 + (rr1 * A->d2_);
                    int col = cc2 + (cc1 * A->d4_);

                    temp1->A2d_[pr][st] = A->A2d_[row][col];
                }
            }
        }
    }
    // temp1->print();
    if (delete_a) A.reset();

    // Figure out B
    // r1
    if (t_b1 == 1)
        r1 = t_b1;
    else if (t_b2 == 1)
        r1 = t_b2;
    else if (t_b3 == 1)
        r1 = t_b3;
    else if (f_b1 == 1)
        r1 = f_b1;

    // r2
    if (t_b1 == 2)
        r2 = t_b1;
    else if (t_b2 == 2)
        r2 = t_b2;
    else if (t_b3 == 2)
        r2 = t_b3;
    else if (f_b1 == 2)
        r2 = f_b1;

    // c1
    if (t_b1 == 3)
        c1 = t_b1;
    else if (t_b2 == 3)
        c1 = t_b2;
    else if (t_b3 == 3)
        c1 = t_b3;
    else if (f_b1 == 3)
        c1 = f_b1;

    // c2
    if (t_b1 == 4)
        c2 = t_b1;
    else if (t_b2 == 4)
        c2 = t_b2;
    else if (t_b3 == 4)
        c2 = t_b3;
    else if (f_b1 == 4)
        c2 = f_b1;

    // Sort B(..,..) to B(rs,tq)
    SharedTensor2f temp2 = SharedTensor2f(new Tensor2f("temp2", dim_r, dim_s, dim_t, dim2_));
#pragma omp parallel for
    for (int r = 0; r < dim_r; r++) {
        for (int s = 0; s < dim_s; s++) {
            int rs = temp2->row_idx_[r][s];
            for (int t = 0; t < dim_t; t++) {
                for (int q = 0; q < dim2_; q++) {
                    int tq = temp2->col_idx_[t][q];

                    if (r1 == t_b1)
                        rr1 = r;
                    else if (r1 == t_b2)
                        rr1 = s;
                    else if (r1 == t_b3)
                        rr1 = t;
                    else if (r1 == f_b1)
                        rr1 = q;

                    if (r2 == t_b1)
                        rr2 = r;
                    else if (r2 == t_b2)
                        rr2 = s;
                    else if (r2 == t_b3)
                        rr2 = t;
                    else if (r2 == f_b1)
                        rr2 = q;

                    if (c1 == t_b1)
                        cc1 = r;
                    else if (c1 == t_b2)
                        cc1 = s;
                    else if (c1 == t_b3)
                        cc1 = t;
                    else if (c1 == f_b1)
                        cc1 = q;

                    if (c2 == t_b1)
                        cc2 = r;
                    else if (c2 == t_b2)
                        cc2 = s;
                    else if (c2 == t_b3)
                        cc2 = t;
                    else if (c2 == f_b1)
                        cc2 = q;

                    int row = rr2 + (rr1 * B->d2_);
                    int col = cc2 + (cc1 * B->d4_);

                    temp2->A2d_[rs][tq] = B->A2d_[row][col];
                }
            }
        }
    }
    // temp2->print();
    if (delete_b) B.reset();

    ta = 'n';
    tb = 'n';
    m = dim1_;
    n = dim2_;
    k = dim_r * dim_s * dim_t;
    nca = k;  // lda
    ncb = n;  // ldb
    ncc = n;  // ldc

    C_SGEMM(ta, tb, m, n, k, alpha, &(temp1->A2d_[0][0]), nca, &(temp2->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
    temp1.reset();
    temp2.reset();

}  //

void Tensor2f::cont424(std::string idx_c, std::string idx_a, std::string idx_b, bool delete_a, SharedTensor2f &A,
                       SharedTensor2f &B, float alpha, float beta) {
    char ta, tb;
    int nca, ncb, ncc;
    int m, n, k;
    int r1, r2, c1, c2;
    int rr1, rr2, cc1, cc2;
    int dim_t;
    int t_a1, f_a1, f_a2, f_a3;
    int t_b1, f_b1;

    // Expected order: C(pq,rs) = \sum_{t} A(pq,rt) B(t,s)

    // Find dummy & free indices for A
    // f_a1
    if (idx_a[0] == idx_c[0])
        f_a1 = 1;
    else if (idx_a[1] == idx_c[0])
        f_a1 = 2;
    else if (idx_a[2] == idx_c[0])
        f_a1 = 3;
    else if (idx_a[3] == idx_c[0])
        f_a1 = 4;

    // f_a2
    if (idx_a[0] == idx_c[1])
        f_a2 = 1;
    else if (idx_a[1] == idx_c[1])
        f_a2 = 2;
    else if (idx_a[2] == idx_c[1])
        f_a2 = 3;
    else if (idx_a[3] == idx_c[1])
        f_a2 = 4;

    // f_a3
    if (idx_a[0] == idx_c[2])
        f_a3 = 1;
    else if (idx_a[1] == idx_c[2])
        f_a3 = 2;
    else if (idx_a[2] == idx_c[2])
        f_a3 = 3;
    else if (idx_a[3] == idx_c[2])
        f_a3 = 4;

    // t_a1
    if (idx_a[0] != idx_c[0] && idx_a[0] != idx_c[1] && idx_a[0] != idx_c[2] && idx_a[0] != idx_c[3])
        t_a1 = 1;
    else if (idx_a[1] != idx_c[0] && idx_a[1] != idx_c[1] && idx_a[1] != idx_c[2] && idx_a[1] != idx_c[3])
        t_a1 = 2;
    else if (idx_a[2] != idx_c[0] && idx_a[2] != idx_c[1] && idx_a[2] != idx_c[2] && idx_a[2] != idx_c[3])
        t_a1 = 3;
    else if (idx_a[3] != idx_c[0] && idx_a[3] != idx_c[1] && idx_a[3] != idx_c[2] && idx_a[3] != idx_c[3])
        t_a1 = 4;
    // outfile->Printf("\tf_a1, f_a2, f_a3, t_a1: %1d, %1d, %1d, %1d  \n", f_a1,f_a2,f_a3,t_a1);

    // Find dummy & free indices for B
    // f_b1 and t_b1
    if (idx_b[0] == idx_c[3]) {
        f_b1 = 1;
        t_b1 = 2;
        dim_t = B->dim2();
    } else if (idx_b[1] == idx_c[3]) {
        f_b1 = 2;
        t_b1 = 1;
        dim_t = B->dim1();
    }
    // outfile->Printf("\tf_b1, t_b1: %1d, %1d \n", f_b1,t_b1);

    // Figure out A
    // r1
    if (t_a1 == 1)
        r1 = t_a1;
    else if (f_a1 == 1)
        r1 = f_a1;
    else if (f_a2 == 1)
        r1 = f_a2;
    else if (f_a3 == 1)
        r1 = f_a3;

    // r2
    if (t_a1 == 2)
        r2 = t_a1;
    else if (f_a1 == 2)
        r2 = f_a1;
    else if (f_a2 == 2)
        r2 = f_a2;
    else if (f_a3 == 2)
        r2 = f_a3;

    // c1
    if (t_a1 == 3)
        c1 = t_a1;
    else if (f_a1 == 3)
        c1 = f_a1;
    else if (f_a2 == 3)
        c1 = f_a2;
    else if (f_a3 == 3)
        c1 = f_a3;

    // c2
    if (t_a1 == 4)
        c2 = t_a1;
    else if (f_a1 == 4)
        c2 = f_a1;
    else if (f_a2 == 4)
        c2 = f_a2;
    else if (f_a3 == 4)
        c2 = f_a3;

    // Sort A(..,..) to A(pq,rt)
    SharedTensor2f temp = SharedTensor2f(new Tensor2f("temp", d1_, d2_, d3_, dim_t));
#pragma omp parallel for
    for (int p = 0; p < d1_; p++) {
        for (int q = 0; q < d2_; q++) {
            int pq = temp->row_idx_[p][q];
            for (int r = 0; r < d3_; r++) {
                for (int t = 0; t < dim_t; t++) {
                    int rt = temp->col_idx_[r][t];

                    if (r1 == f_a1)
                        rr1 = p;
                    else if (r1 == f_a2)
                        rr1 = q;
                    else if (r1 == f_a3)
                        rr1 = r;
                    else if (r1 == t_a1)
                        rr1 = t;

                    if (r2 == f_a1)
                        rr2 = p;
                    else if (r2 == f_a2)
                        rr2 = q;
                    else if (r2 == f_a3)
                        rr2 = r;
                    else if (r2 == t_a1)
                        rr2 = t;

                    if (c1 == f_a1)
                        cc1 = p;
                    else if (c1 == f_a2)
                        cc1 = q;
                    else if (c1 == f_a3)
                        cc1 = r;
                    else if (c1 == t_a1)
                        cc1 = t;

                    if (c2 == f_a1)
                        cc2 = p;
                    else if (c2 == f_a2)
                        cc2 = q;
                    else if (c2 == f_a3)
                        cc2 = r;
                    else if (c2 == t_a1)
                        cc2 = t;

                    int row = rr2 + (rr1 * A->d2_);
                    int col = cc2 + (cc1 * A->d4_);

                    temp->A2d_[pq][rt] = A->A2d_[row][col];
                }
            }
        }
    }
    if (delete_a) A.reset();

    ta = 'n';
    if (t_b1 == 1)
        tb = 'n';
    else
        tb = 't';
    m = d1_ * d2_ * d3_;
    n = d4_;
    k = temp->d4_;
    nca = k;  // lda
    // ldb
    if (tb == 't')
        ncb = k;
    else
        ncb = n;
    ncc = n;  // ldc

    C_SGEMM(ta, tb, m, n, k, alpha, &(temp->A2d_[0][0]), nca, &(B->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
    temp.reset();
}  //

void Tensor2f::cont244(std::string idx_c, std::string idx_a, std::string idx_b, bool delete_b, SharedTensor2f &A,
                       SharedTensor2f &B, float alpha, float beta) {
    char ta, tb;
    int nca, ncb, ncc;
    int m, n, k;
    int r1, r2, c1, c2;
    int rr1, rr2, cc1, cc2;
    int dim_t;
    int t_a1, f_a1;
    int t_b1, f_b1, f_b2, f_b3;

    // Expected order: C(pq,rs) = \sum_{t} A(p,t) B(tq,rs)

    // Find dummy & free indices for B
    // f_b1
    if (idx_b[0] == idx_c[1])
        f_b1 = 1;
    else if (idx_b[1] == idx_c[1])
        f_b1 = 2;
    else if (idx_b[2] == idx_c[1])
        f_b1 = 3;
    else if (idx_b[3] == idx_c[1])
        f_b1 = 4;

    // f_b2
    if (idx_b[0] == idx_c[2])
        f_b2 = 1;
    else if (idx_b[1] == idx_c[2])
        f_b2 = 2;
    else if (idx_b[2] == idx_c[2])
        f_b2 = 3;
    else if (idx_b[3] == idx_c[2])
        f_b2 = 4;

    // f_b3
    if (idx_b[0] == idx_c[3])
        f_b3 = 1;
    else if (idx_b[1] == idx_c[3])
        f_b3 = 2;
    else if (idx_b[2] == idx_c[3])
        f_b3 = 3;
    else if (idx_b[3] == idx_c[3])
        f_b3 = 4;

    // t_b1
    if (idx_b[0] != idx_c[0] && idx_b[0] != idx_c[1] && idx_b[0] != idx_c[2] && idx_b[0] != idx_c[3])
        t_b1 = 1;
    else if (idx_b[1] != idx_c[0] && idx_b[1] != idx_c[1] && idx_b[1] != idx_c[2] && idx_b[1] != idx_c[3])
        t_b1 = 2;
    else if (idx_b[2] != idx_c[0] && idx_b[2] != idx_c[1] && idx_b[2] != idx_c[2] && idx_b[2] != idx_c[3])
        t_b1 = 3;
    else if (idx_b[3] != idx_c[0] && idx_b[3] != idx_c[1] && idx_b[3] != idx_c[2] && idx_b[3] != idx_c[3])
        t_b1 = 4;

    // Find dummy & free indices for A
    // f_a1 and t_a1
    if (idx_a[0] == idx_c[0]) {
        f_a1 = 1;
        t_a1 = 2;
        dim_t = A->dim2();
    } else if (idx_a[1] == idx_c[0]) {
        f_a1 = 2;
        t_a1 = 1;
        dim_t = A->dim1();
    }

    // Figure out B
    // r1
    if (t_b1 == 1)
        r1 = t_b1;
    else if (f_b1 == 1)
        r1 = f_b1;
    else if (f_b2 == 1)
        r1 = f_b2;
    else if (f_b3 == 1)
        r1 = f_b3;

    // r2
    if (t_b1 == 2)
        r2 = t_b1;
    else if (f_b1 == 2)
        r2 = f_b1;
    else if (f_b2 == 2)
        r2 = f_b2;
    else if (f_b3 == 2)
        r2 = f_b3;

    // c1
    if (t_b1 == 3)
        c1 = t_b1;
    else if (f_b1 == 3)
        c1 = f_b1;
    else if (f_b2 == 3)
        c1 = f_b2;
    else if (f_b3 == 3)
        c1 = f_b3;

    // c2
    if (t_b1 == 4)
        c2 = t_b1;
    else if (f_b1 == 4)
        c2 = f_b1;
    else if (f_b2 == 4)
        c2 = f_b2;
    else if (f_b3 == 4)
        c2 = f_b3;

    // Sort B(..,..) to B(tq,rs)
    SharedTensor2f temp = SharedTensor2f(new Tensor2f("temp", dim_t, d2_, d3_, d4_));
#pragma omp parallel for
    for (int t = 0; t < dim_t; t++) {
        for (int q = 0; q < d2_; q++) {
            int tq = temp->row_idx_[t][q];
            for (int r = 0; r < d3_; r++) {
                for (int s = 0; s < d4_; s++) {
                    int rs = temp->col_idx_[r][s];

                    if (r1 == t_b1)
                        rr1 = t;
                    else if (r1 == f_b1)
                        rr1 = q;
                    else if (r1 == f_b2)
                        rr1 = r;
                    else if (r1 == f_b3)
                        rr1 = s;

                    if (r2 == t_b1)
                        rr2 = t;
                    else if (r2 == f_b1)
                        rr2 = q;
                    else if (r2 == f_b2)
                        rr2 = r;
                    else if (r2 == f_b3)
                        rr2 = s;

                    if (c1 == t_b1)
                        cc1 = t;
                    else if (c1 == f_b1)
                        cc1 = q;
                    else if (c1 == f_b2)
                        cc1 = r;
                    else if (c1 == f_b3)
                        cc1 = s;

                    if (c2 == t_b1)
                        cc2 = t;
                    else if (c2 == f_b1)
                        cc2 = q;
                    else if (c2 == f_b2)
                        cc2 = r;
                    else if (c2 == f_b3)
                        cc2 = s;

                    int row = rr2 + (rr1 * B->d2_);
                    int col = cc2 + (cc1 * B->d4_);

                    temp->A2d_[tq][rs] = B->A2d_[row][col];
                }
            }
        }
    }
    if (delete_b) B.reset();

    if (t_a1 == 2)
        ta = 'n';
    else
        ta = 't';
    tb = 'n';
    m = d1_;
    n = d2_ * d3_ * d4_;
    k = temp->d1_;
    // nca = transa ? m : k; // lda
    // ncb = transb ? k : n; // ldb
    if (ta == 't')
        nca = m;
    else
        nca = k;
    ncb = n;
    ncc = n;  // ldc

    C_SGEMM(ta, tb, m, n, k, alpha, &(A->A2d_[0][0]), nca, &(temp->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
    temp.reset();

}  //

void Tensor2f::cont233(std::string idx_c, std::string idx_a, std::string idx_b, SharedTensor2f &A, SharedTensor2f &B,
                       float alpha, float beta) {
    char ta, tb;
    int lda, ldb, ldc;
    int m, n, k;
    int dim_r;
    int t_a1, f_a1;
    int t_b1, f_b1;

    // Expected order: C(Q,pq) = \sum_{r} A(p,r) B(Q,rq)

    // Find dummy & free indices for A
    // f_a1 and t_a1
    if (idx_a[0] == idx_c[0]) {
        f_a1 = 1;
        t_a1 = 2;
        dim_r = A->dim2();
    } else if (idx_a[1] == idx_c[0]) {
        f_a1 = 2;
        t_a1 = 1;
        dim_r = A->dim1();
    }

    // Find dummy & free indices for B
    // f_b1 and t_b1
    if (idx_b[0] == idx_c[1]) {
        f_b1 = 1;
        t_b1 = 2;
    } else if (idx_b[1] == idx_c[1]) {
        f_b1 = 2;
        t_b1 = 1;
    }

    m = d2_;
    n = d3_;

    if (t_a1 == 2)
        ta = 'n';
    else
        ta = 't';

    if (t_b1 == 1)
        tb = 'n';
    else
        tb = 't';

    if (ta == 'n')
        k = A->dim2();
    else
        k = A->dim1();

    // lda = transa ? m : k;
    if (ta == 't')
        lda = m;
    else
        lda = k;

    // ldb = transb ? k : n;
    if (tb == 't')
        ldb = k;
    else
        ldb = n;

    ldc = n;

#pragma omp parallel for
    for (int Q = 0; Q < dim1_; Q++) {
        C_SGEMM(ta, tb, m, n, k, alpha, A->A2d_[0], lda, B->A2d_[Q], ldb, beta, A2d_[Q], ldc);
    }

}  //

void Tensor2f::cont323(std::string idx_c, std::string idx_a, std::string idx_b, bool delete_a, SharedTensor2f &A,
                       SharedTensor2f &B, float alpha, float beta) {
    char ta, tb;
    int nca, ncb, ncc;
    int m, n, k;
    int dim_r;
    int t_a1, f_a1;
    int t_b1, f_b1;
    int r1, c1;
    int rr1, cc1;

    // Expected order: C(Q,pq) = \sum_{r} A(Q,pr) B(r,q)

    // Find dummy & free indices for A
    // f_a1 and t_a1
    if (idx_a[0] == idx_c[0]) {
        f_a1 = 1;
        t_a1 = 2;
    } else if (idx_a[1] == idx_c[0]) {
        f_a1 = 2;
        t_a1 = 1;
    }

    // Find dummy & free indices for B
    // f_b1 and t_b1
    if (idx_b[0] == idx_c[1]) {
        f_b1 = 1;
        t_b1 = 2;
        dim_r = B->dim2();
    } else if (idx_b[1] == idx_c[1]) {
        f_b1 = 2;
        t_b1 = 1;
        dim_r = B->dim1();
    }

    // Figure out A
    // r1
    if (t_a1 == 1)
        r1 = t_a1;
    else if (f_a1 == 1)
        r1 = f_a1;

    // c1
    if (t_a1 == 2)
        c1 = t_a1;
    else if (f_a1 == 2)
        c1 = f_a1;

    // Sort A(Q,..) to A(Q,pr)
    SharedTensor2f temp = SharedTensor2f(new Tensor2f("temp", d1_, d2_, dim_r));
#pragma omp parallel for
    for (int Q = 0; Q < dim1_; Q++) {
        for (int p = 0; p < d2_; p++) {
            for (int r = 0; r < dim_r; r++) {
                int pr = r + (p * dim_r);

                if (r1 == f_a1)
                    rr1 = p;
                else if (r1 == t_a1)
                    rr1 = r;

                if (c1 == f_a1)
                    cc1 = p;
                else if (c1 == t_a1)
                    cc1 = r;

                int col = cc1 + (rr1 * A->d3_);

                temp->A2d_[Q][pr] = A->A2d_[Q][col];
            }
        }
    }
    if (delete_a) A.reset();

    m = d1_ * d2_;
    n = d3_;
    k = dim_r;
    ta = 'n';
    if (t_b1 == 1)
        tb = 'n';
    else
        tb = 't';
    nca = k;  // lda
              // ldb = transb ? k : n;
    if (tb == 't')
        ncb = k;
    else
        ncb = n;
    ncc = n;

    C_SGEMM(ta, tb, m, n, k, alpha, &(temp->A2d_[0][0]), nca, &(B->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
    temp.reset();
}  //

void Tensor2f::cont332(std::string idx_c, std::string idx_a, std::string idx_b, bool delete_a, bool delete_b,
                       SharedTensor2f &A, SharedTensor2f &B, float alpha, float beta) {
    char ta, tb;
    int nca, ncb, ncc;  // number of columns
    int m, n, k;
    int dim_r;
    int t_a1, f_a1;
    int t_b1, f_b1;
    int ra, ca, rb, cb;
    int rra, cca, rrb, ccb;

    // Expected order: C(pq) = \sum_{Qr} A(Q,rp) B(Q,rq)

    // Find dummy & free indices for A
    // f_a1 and t_a1
    if (idx_a[0] == idx_c[0]) {
        f_a1 = 1;
        t_a1 = 2;
        dim_r = A->d3_;
    } else if (idx_a[1] == idx_c[0]) {
        f_a1 = 2;
        t_a1 = 1;
        dim_r = A->d2_;
    }

    // Find dummy & free indices for B
    // f_b1 and t_b1
    if (idx_b[0] == idx_c[1]) {
        f_b1 = 1;
        t_b1 = 2;
    } else if (idx_b[1] == idx_c[1]) {
        f_b1 = 2;
        t_b1 = 1;
    }

    // Figure out A
    // ra
    if (t_a1 == 1)
        ra = t_a1;
    else if (f_a1 == 1)
        ra = f_a1;

    // ca
    if (t_a1 == 2)
        ca = t_a1;
    else if (f_a1 == 2)
        ca = f_a1;

    // Sort A(Q,..) to A(Q,rp)
    SharedTensor2f temp1 = SharedTensor2f(new Tensor2f("temp1", A->d1_, dim_r, dim1_));
#pragma omp parallel for
    for (int Q = 0; Q < A->d1_; Q++) {
        for (int r = 0; r < dim_r; r++) {
            for (int p = 0; p < dim1_; p++) {
                int rp = p + (r * dim1_);

                if (ra == t_a1)
                    rra = r;
                else if (ra == f_a1)
                    rra = p;

                if (ca == t_a1)
                    cca = r;
                else if (ca == f_a1)
                    cca = p;

                int col = cca + (rra * A->d3_);

                temp1->A2d_[Q][rp] = A->A2d_[Q][col];
            }
        }
    }
    if (delete_a) A.reset();

    // Figure out B
    // rb
    if (t_b1 == 1)
        rb = t_b1;
    else if (f_b1 == 1)
        rb = f_b1;

    // cb
    if (t_b1 == 2)
        cb = t_b1;
    else if (f_b1 == 2)
        cb = f_b1;

    // Sort B(Q,..) to B(Q,rq)
    SharedTensor2f temp2 = SharedTensor2f(new Tensor2f("temp2", B->d1_, dim_r, dim2_));
#pragma omp parallel for
    for (int Q = 0; Q < B->d1_; Q++) {
        for (int r = 0; r < dim_r; r++) {
            for (int q = 0; q < dim2_; q++) {
                int rq = q + (r * dim2_);

                if (rb == t_b1)
                    rrb = r;
                else if (rb == f_b1)
                    rrb = q;

                if (cb == t_b1)
                    ccb = r;
                else if (cb == f_b1)
                    ccb = q;

                int col = ccb + (rrb * B->d3_);

                temp2->A2d_[Q][rq] = B->A2d_[Q][col];
            }
        }
    }
    if (delete_b) B.reset();

    m = dim1_;
    n = dim2_;
    k = A->d1_ * dim_r;
    ta = 't';
    tb = 'n';
    nca = m;
    ncb = n;
    ncc = n;

    C_SGEMM(ta, tb, m, n, k, alpha, &(temp1->A2d_[0][0]), nca, &(temp2->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
    temp2.reset();
    temp2.reset();
}  //

float Tensor2f::get_max_element() {
    float value = 0.0;
#pragma omp parallel for
    for (int i = 0; i < dim1_; i++) {
        for (int j = 0; j < dim2_; j++) {
            if (std::fabs(A2d_[i][j]) > value) value = std::fabs(A2d_[i][j]);
        }
    }
    return value;
}  //

/********************************************************************************************/
/************************** 3d array ********************************************************/
/********************************************************************************************/
Tensor3f::Tensor3f(int d1, int d2, int d3) {
    A3d_ = NULL;
    dim1_ = d1;
    dim2_ = d2;
    dim3_ = d3;
    memalloc();
}  //

Tensor3f::Tensor3f(std::string name, int d1, int d2, int d3) {
    A3d_ = NULL;
    dim1_ = d1;
    dim2_ = d2;
    dim3_ = d3;
    name_ = name;
    memalloc();
}  //

Tensor3f::Tensor3f() {
    A3d_ = NULL;
    dim1_ = 0;
    dim2_ = 0;
    dim3_ = 0;

}  //

Tensor3f::~Tensor3f() { release(); }  //

void Tensor3f::memalloc() {
    if (A3d_) release();
    A3d_ = init_3d_array_float(dim1_, dim2_, dim3_);
    zero();
}  //

void Tensor3f::init(int d1, int d2, int d3) {
    dim1_ = d1;
    dim2_ = d2;
    dim3_ = d3;
    memalloc();
}  //

void Tensor3f::init(std::string name, int d1, int d2, int d3) {
    dim1_ = d1;
    dim2_ = d2;
    dim3_ = d3;
    name_ = name;
    memalloc();
}  //

void Tensor3f::zero() { memset(&(A3d_[0][0][0]), 0, sizeof(float) * dim1_ * dim2_ * dim3_); }  //

void Tensor3f::print() {
    if (name_.length()) outfile->Printf("\n ## %s ##\n", name_.c_str());
    for (int i = 0; i < dim1_; i++) {
        outfile->Printf("\n Irrep: %d\n", i + 1);
        print_mat(A3d_[i], dim2_, dim3_, "outfile");
    }

}  //

void Tensor3f::release() {
    if (!A3d_) return;
    free_3d_array(A3d_, dim1_, dim2_);
    A3d_ = NULL;
}  //

void Tensor3f::set(int h, int i, int j, float value) { A3d_[h][i][j] = value; }  //

float Tensor3f::get(int h, int i, int j) { return A3d_[h][i][j]; }  //


/********************************************************************************************/
/********************************************************************************************/
}  // namespace dfoccwave
}  // namespace psi
