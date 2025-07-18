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

// Latest revision on April 28, 2013.
#include <stdio.h>
#include "psi4/libqt/qt.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.hpp"
#include "arrays.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libmints/matrix.h"

namespace psi {
namespace dfoccwave {

/********************************************************************************************/
/************************** 1d array ********************************************************/
/********************************************************************************************/
Array1d::Array1d(int d1) {
    A1d_ = NULL;
    dim1_ = d1;
    memalloc();
}  //

Array1d::Array1d(std::string name, int d1) {
    A1d_ = NULL;
    dim1_ = d1;
    name_ = name;
    memalloc();
}  //

Array1d::Array1d() {
    A1d_ = NULL;
    dim1_ = 0;

}  //

Array1d::~Array1d() { release(); }  //

Array1d *Array1d::generate(int d1) { return new Array1d(d1); }

Array1d *Array1d::generate(std::string name, int d1) { return new Array1d(name, d1); }

void Array1d::memalloc() {
    if (A1d_) release();
    A1d_ = new double[dim1_];
}  //

void Array1d::init(int d1) {
    dim1_ = d1;
    if (A1d_) release();
    A1d_ = new double[dim1_];
}  //

void Array1d::init(std::string name, int d1) {
    dim1_ = d1;
    name_ = name;
    if (A1d_) release();
    A1d_ = new double[dim1_];
}  //

void Array1d::zero() { memset(A1d_, 0, sizeof(double) * dim1_); }  //

void Array1d::print() {
    if (name_.length()) outfile->Printf("\n ## %s ##\n", name_.c_str());
    for (int p = 0; p < dim1_; p++) {
        outfile->Printf(" %3d %10.7f \n", p, A1d_[p]);
    }

}  //

void Array1d::print(std::string out) {
    std::shared_ptr<psi::PsiOutStream> printer =
        (out == "outfile" ? outfile : std::shared_ptr<PsiOutStream>(new PsiOutStream(out)));
    if (name_.length()) printer->Printf("\n ## %s ##\n", name_.c_str());
    for (int p = 0; p < dim1_; p++) {
        printer->Printf(" %3d %10.7f \n", p, A1d_[p]);
    }

}  //

void Array1d::release() {
    if (!A1d_) return;
    delete[] A1d_;
    A1d_ = NULL;
}  //

void Array1d::set(int i, double value) { A1d_[i] = value; }  //

void Array1d::set(double *vec) {
    for (int i = 0; i < dim1_; ++i) A1d_[i] = vec[i];
}  //

void Array1d::set(const Array1d *vec) {
    for (int i = 0; i < dim1_; ++i) A1d_[i] = vec->A1d_[i];
}  //

double Array1d::get(int i) { return A1d_[i]; }  //

void Array1d::add(const Array1d *Adum) {
    double *lhs, *rhs;
    size_t size = dim1_;
    if (size) {
        lhs = A1d_;
        rhs = Adum->A1d_;
        for (size_t ij = 0; ij < size; ++ij) {
            *lhs += *rhs;
            lhs++;
            rhs++;
        }
    }
}  //

void Array1d::add(int i, double value) { A1d_[i] += value; }  //

void Array1d::subtract(const Array1d *Adum) {
    double *lhs, *rhs;
    size_t size = dim1_;
    if (size) {
        lhs = A1d_;
        rhs = Adum->A1d_;
        for (size_t ij = 0; ij < size; ++ij) {
            *lhs -= *rhs;
            lhs++;
            rhs++;
        }
    }
}  //

void Array1d::subtract(int i, double value) { A1d_[i] -= value; }  //

double Array1d::rms() {
    double summ = 0.0;
    for (int i = 0; i < dim1_; ++i) summ += A1d_[i] * A1d_[i];
    summ = sqrt(summ) / dim1_;

    return summ;
}  //

double Array1d::rms(const Array1d *Atemp) {
    double summ = 0.0;
    for (int i = 0; i < dim1_; ++i) summ += (A1d_[i] - Atemp->A1d_[i]) * (A1d_[i] - Atemp->A1d_[i]);
    summ = sqrt(summ) / dim1_;

    return summ;
}  //

double Array1d::dot(const Array1d *y) {
    double value = 0.0;
    int incx = 1;
    int incy = 1;
    if (dim1_ == y->dim1_) value = C_DDOT(dim1_, A1d_, incx, y->A1d_, incy);
    return value;
}  //

void Array1d::gbmv(bool transa, double alpha, const Array2d *a, const Array1d *b, double beta) {
    char ta = transa ? 't' : 'n';
    int m, n, k, kl, ku, incx, incy, lda;

    m = a->dim1_;
    n = a->dim2_;
    k = b->dim1_;
    kl = m - 1;  // # of subdiagonal of matrix A, at most kl = m - 1
    ku = n - 1;  // # of superdiagonal of matrix A, at most ku = n - 1
    lda = kl + ku + 1;
    incx = 1;  // increments in elements of b vector
    incy = 1;  // increments in elements of A1d_

    // A1d_ = alpha * A * b + beta, where A is a band matrix
    if (m && n) {
        C_DGBMV(ta, m, n, kl, ku, alpha, &(a->A2d_[0][0]), lda, b->A1d_, incx, beta, A1d_, incy);
    }
}  //

void Array1d::gemv(bool transa, double alpha, const Array2d *a, const Array1d *b, double beta) {
    char ta = transa ? 't' : 'n';
    int m, n, k, kl, ku, incx, incy, lda;

    m = a->dim1_;
    n = a->dim2_;
    k = b->dim1_;
    lda = m;
    incx = 1;  // increments in elements of b vector
    incy = 1;  // increments in elements of A1d_

    // A1d_ = alpha * A * b + beta, where A is a general matrix
    if (m && n) {
        C_DGEMV(ta, m, n, alpha, &(a->A2d_[0][0]), lda, b->A1d_, incx, beta, A1d_, incy);
    }
}  //

double Array1d::xay(const Array2d *a, const Array1d *y) {
    double value = 0.0;
    Array1d *ay = new Array1d(a->dim1_);
    ay->zero();
    ay->gemv(false, 1.0, a, y, 0.0);
    value = dot(ay);
    delete ay;
    return value;
}  //

void Array1d::scale(double a) {
    size_t size = dim1_;
    if (size) C_DSCAL(size, a, A1d_, 1);
}  //

void Array1d::copy(double *x) {
    size_t size;
    size = dim1_ * sizeof(double);
    if (size) memcpy(&(A1d_[0]), &(x[0]), size);
}  //

void Array1d::copy(const Array1d *x) {
    size_t size;
    size = dim1_ * sizeof(double);
    if (size) memcpy(&(A1d_[0]), &(x->A1d_[0]), size);
}  //

void Array1d::row_vector(Array2d *A, int n) {
    int dim = A->dim2();
    for (int i = 0; i < dim; i++) A1d_[i] = A->get(n, i);
}  //

void Array1d::column_vector(Array2d *A, int n) {
    int dim = A->dim1();
    for (int i = 0; i < dim; i++) A1d_[i] = A->get(i, n);
}  //

void Array1d::dirprd(Array1d *a, Array1d *b) {
    int dima = a->dim1();
    int dimb = b->dim1();

    if (dima == dimb && dima == dim1_) {
        for (int i = 0; i < dim1_; i++) A1d_[i] = a->get(i) * b->get(i);
    } else
        throw SanityCheckError("Vector dimensions do NOT match!", __FILE__, __LINE__);
}  //

/********************************************************************************************/
/************************** 2d array ********************************************************/
/********************************************************************************************/
Array2d::Array2d(int d1, int d2) {
    A2d_ = NULL;
    dim1_ = d1;
    dim2_ = d2;
    memalloc();
}  //

Array2d::Array2d(std::string name, int d1, int d2) {
    A2d_ = NULL;
    dim1_ = d1;
    dim2_ = d2;
    name_ = name;
    memalloc();
}  //

Array2d::Array2d() {
    A2d_ = NULL;
    dim1_ = 0;
    dim2_ = 0;

}  //

Array2d::Array2d(psi::PSIO *psio, size_t fileno, std::string name, int d1, int d2) {
    A2d_ = NULL;
    dim1_ = d1;
    dim2_ = d2;
    name_ = name;
    memalloc();
    read(psio, fileno);
}

Array2d::Array2d(std::shared_ptr<psi::PSIO> psio, size_t fileno, std::string name, int d1, int d2) {
    A2d_ = NULL;
    dim1_ = d1;
    dim2_ = d2;
    name_ = name;
    memalloc();
    read(psio, fileno);
}

Array2d::Array2d(psi::PSIO &psio, size_t fileno, std::string name, int d1, int d2) {
    A2d_ = NULL;
    dim1_ = d1;
    dim2_ = d2;
    name_ = name;
    memalloc();
    read(&psio, fileno);
}  //

Array2d::~Array2d() { release(); }  //

Array2d *Array2d::generate(int d1, int d2) { return new Array2d(d1, d2); }

Array2d *Array2d::generate(std::string name, int d1, int d2) { return new Array2d(name, d1, d2); }

void Array2d::memalloc() {
    if (A2d_) release();
    A2d_ = block_matrix(dim1_, dim2_);
}  //

void Array2d::init(int d1, int d2) {
    dim1_ = d1;
    dim2_ = d2;
    if (A2d_) release();
    A2d_ = block_matrix(dim1_, dim2_);
}  //

void Array2d::init(std::string name, int d1, int d2) {
    dim1_ = d1;
    dim2_ = d2;
    name_ = name;
    if (A2d_) release();
    A2d_ = block_matrix(dim1_, dim2_);
}  //

void Array2d::zero() { memset(A2d_[0], 0, sizeof(double) * dim1_ * dim2_); }  //

void Array2d::zero_diagonal() {
    if (dim1_ == dim2_) {
        for (int i = 0; i < dim1_; i++) A2d_[i][i] = 0.0;
    }
}  //

void Array2d::print() {
    if (A2d_) {
        if (name_.length()) outfile->Printf("\n ## %s ##\n", name_.c_str());
        print_mat(A2d_, dim1_, dim2_, "outfile");
    }
}  //

void Array2d::print(std::string out) {
    std::shared_ptr<psi::PsiOutStream> printer =
        (out == "outfile" ? outfile : std::shared_ptr<PsiOutStream>(new PsiOutStream(out)));
    if (A2d_) {
        if (name_.length()) printer->Printf("\n ## %s ##\n", name_.c_str());
        print_mat(A2d_, dim1_, dim2_, out);
    }
}  //

void Array2d::release() {
    if (!A2d_) return;
    free_block(A2d_);
    A2d_ = NULL;
}  //

void Array2d::set(int i, int j, double value) { A2d_[i][j] = value; }  //

void Array2d::set(double **A) {
    if (A == NULL) return;
    for (int i = 0; i < dim1_; ++i) {
        for (int j = 0; j < dim2_; ++j) {
            A2d_[i][j] = A[i][j];
        }
    }
}  //

void Array2d::set(Array2d *A) {
    if (A == NULL) return;
    for (int i = 0; i < dim1_; ++i) {
        for (int j = 0; j < dim2_; ++j) {
            A2d_[i][j] = A->A2d_[i][j];
        }
    }
}  //

void Array2d::set(SharedMatrix A) {
    for (int i = 0; i < dim1_; ++i) {
        for (int j = 0; j < dim2_; ++j) {
            A2d_[i][j] = A->get(0, i, j);
        }
    }
}  //

double Array2d::get(int i, int j) { return A2d_[i][j]; }  //

void Array2d::gemm(bool transa, bool transb, const Array2d *a, const Array2d *b, double alpha, double beta) {
    char ta = transa ? 't' : 'n';
    char tb = transb ? 't' : 'n';
    int m, n, k, nca, ncb, ncc;

    m = dim1_;
    n = dim2_;
    k = transa ? a->dim1_ : a->dim2_;
    nca = transa ? m : k;
    ncb = transb ? k : n;
    ncc = n;

    if (m && n && k) {
        C_DGEMM(ta, tb, m, n, k, alpha, &(a->A2d_[0][0]), nca, &(b->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
    }
}  //

void Array2d::contract(bool transa, bool transb, int m, int n, int k, const Array2d *a, const Array2d *b, double alpha,
                       double beta) {
    char ta = transa ? 't' : 'n';
    char tb = transb ? 't' : 'n';
    int nca, ncb, ncc;

    nca = transa ? m : k;
    ncb = transb ? k : n;
    ncc = n;

    if (m && n && k) {
        C_DGEMM(ta, tb, m, n, k, alpha, &(a->A2d_[0][0]), nca, &(b->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
    }
}  //

void Array2d::contract323(bool transa, bool transb, int m, int n, const Array2d *a, const Array2d *b, double alpha,
                          double beta) {
    char ta = transa ? 't' : 'n';
    char tb = transb ? 't' : 'n';
    int k, nca, ncb, ncc;

    k = transb ? b->dim2_ : b->dim1_;
    nca = transa ? m : k;
    ncb = transb ? k : n;
    ncc = n;

    if (m && n && k) {
        for (int Q = 0; Q < dim1_; Q++) {
            C_DGEMM(ta, tb, m, n, k, alpha, a->A2d_[Q], nca, b->A2d_[0], ncb, beta, A2d_[Q], ncc);
        }
    }
}  //

void Array2d::contract233(bool transa, bool transb, int m, int n, const Array2d *a, const Array2d *b, double alpha,
                          double beta) {
    char ta = transa ? 't' : 'n';
    char tb = transb ? 't' : 'n';
    int k, lda, ldb, ldc;

    k = transa ? a->dim1_ : a->dim2_;
    lda = transa ? m : k;
    ldb = transb ? k : n;
    ldc = n;

    if (m && n && k) {
        for (int Q = 0; Q < dim1_; Q++) {
            C_DGEMM(ta, tb, m, n, k, alpha, a->A2d_[0], lda, b->A2d_[Q], ldb, beta, A2d_[Q], ldc);
        }
    }
}  //

void Array2d::add(const Array2d *Adum) {
    double *lhs, *rhs;
    size_t size = dim1_ * dim2_;
    if (size) {
        lhs = A2d_[0];
        rhs = Adum->A2d_[0];
        for (size_t ij = 0; ij < size; ++ij) {
            *lhs += *rhs;
            lhs++;
            rhs++;
        }
    }
}  //

void Array2d::add(double alpha, const Array2d *Adum) {
    Array2d *temp;
    temp = new Array2d(Adum->dim1_, Adum->dim2_);
    temp->copy(Adum);
    temp->scale(alpha);
    add(temp);
    delete temp;
}  //

void Array2d::add(int i, int j, double value) { A2d_[i][j] += value; }  //

void Array2d::subtract(const Array2d *Adum) {
    double *lhs, *rhs;
    size_t size = dim1_ * dim2_;
    if (size) {
        lhs = A2d_[0];
        rhs = Adum->A2d_[0];
        for (size_t ij = 0; ij < size; ++ij) {
            *lhs -= *rhs;
            lhs++;
            rhs++;
        }
    }
}  //

void Array2d::subtract(int i, int j, double value) { A2d_[i][j] -= value; }  //

Array2d *Array2d::transpose() {
    Array2d *temp;
    temp = new Array2d(dim2_, dim1_);
    temp->zero();

    for (int i = 0; i < dim2_; ++i) {
        for (int j = 0; j < dim1_; ++j) {
            temp->set(i, j, A2d_[j][i]);
        }
    }

    return temp;
}  //

void Array2d::copy(const Array2d *Adum) {
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
    if (dim1_ != 0 && dim2_ != 0) {
        memcpy(A2d_[0], Adum->A2d_[0], dim1_ * dim2_ * sizeof(double));
    }
}  //

void Array2d::copy(double **a) {
    size_t size = dim1_ * dim2_ * sizeof(double);
    if (size) memcpy(&(A2d_[0][0]), &(a[0][0]), size);
}

void Array2d::level_shift(double value) {
    for (int i = 0; i < dim1_; ++i) {
        subtract(i, i, value);
    }

}  //

void Array2d::outer_product(const Array1d *x, const Array1d *y) {
    for (int i = 0; i < x->dim1_; i++) {
        for (int j = 0; j < y->dim1_; j++) {
            A2d_[i][j] = x->A1d_[i] * y->A1d_[j];
        }
    }

}  //

// TODO:
// DGER compute the rank-one update of a general matrix: A <-- A + alpha * x * yT
// dger(m, n, alpha, x, incx, y, incy, a, lda);

void Array2d::scale(double a) {
    size_t size;
    size = dim1_ * dim2_;
    if (size) C_DSCAL(size, a, &(A2d_[0][0]), 1);
}  //

void Array2d::scale_row(int m, double a) { C_DSCAL(dim1_, a, &(A2d_[m][0]), 1); }  //

void Array2d::scale_column(int n, double a) { C_DSCAL(dim2_, a, &(A2d_[0][n]), dim1_); }  //

void Array2d::identity() {
    zero();
    for (int i = 0; i < dim1_; ++i) A2d_[i][i] = 1.0;
}  //

double Array2d::trace() {
    double value = 0.0;
    for (int i = 0; i < dim1_; ++i) value += A2d_[i][i];
    return value;
}  //

void Array2d::transform(const Array2d *a, const Array2d *transformer) {
    Array2d *temp = new Array2d(a->dim1_, transformer->dim2_);
    temp->zero();
    temp->gemm(false, false, a, transformer, 1.0, 0.0);
    gemm(true, false, transformer, temp, 1.0, 0.0);
    delete temp;
}  //

void Array2d::back_transform(const Array2d *a, const Array2d *transformer) {
    Array2d *temp = new Array2d(a->dim1_, transformer->dim2_);
    temp->zero();
    temp->gemm(false, true, a, transformer, 1.0, 0.0);
    gemm(false, false, transformer, temp, 1.0, 0.0);
    delete temp;
}  //

void Array2d::pseudo_transform(const Array2d *a, const Array2d *transformer) {
    Array2d *temp = new Array2d(a->dim1_, transformer->dim2_);
    temp->zero();
    temp->gemm(false, false, a, transformer, 1.0, 0.0);
    gemm(false, false, transformer, temp, 1.0, 0.0);
    delete temp;
}  //

void Array2d::triple_gemm(const Array2d *a, const Array2d *b, const Array2d *c) {
    if (a->dim2_ == b->dim1_ && b->dim2_ == c->dim1_ && a->dim1_ == dim1_ && c->dim2_ == dim2_) {
        Array2d *bc = new Array2d(b->dim1_, c->dim2_);
        bc->zero();
        bc->gemm(false, false, b, c, 1.0, 0.0);
        gemm(false, false, a, bc, 1.0, 0.0);
        delete bc;
    } else {
        outfile->Printf("\n Warning!!! Matrix dimensions do NOT match in triple_gemm().\n");
    }

}  //

double Array2d::vector_dot(Array2d *rhs) {
    double value = 0.0;
    size_t size = dim1_ * dim2_;
    if (size) value += C_DDOT(size, (&A2d_[0][0]), 1, &(rhs->A2d_[0][0]), 1);
    return value;
}  //

double Array2d::vector_dot(double **rhs) {
    double value = 0.0;
    size_t size = dim1_ * dim2_;
    if (size) value += C_DDOT(size, (&A2d_[0][0]), 1, &(rhs[0][0]), 1);
    return value;
}  //

void Array2d::write(std::shared_ptr<psi::PSIO> psio, size_t fileno) {
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno))
        already_open = true;
    else
        psio->open(fileno, PSIO_OPEN_OLD);
    psio->write_entry(fileno, const_cast<char *>(name_.c_str()), (char *)A2d_[0], sizeof(double) * dim1_ * dim2_);
    if (!already_open) psio->close(fileno, 1);  // Close and keep
}  //

void Array2d::write(psi::PSIO *const psio, size_t fileno) {
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno))
        already_open = true;
    else
        psio->open(fileno, PSIO_OPEN_OLD);
    psio->write_entry(fileno, const_cast<char *>(name_.c_str()), (char *)A2d_[0], sizeof(double) * dim1_ * dim2_);
    if (!already_open) psio->close(fileno, 1);  // Close and keep
}  //

void Array2d::write(psi::PSIO &psio, size_t fileno) { write(&psio, fileno); }  //

void Array2d::read(psi::PSIO *psio, size_t fileno) {
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno))
        already_open = true;
    else
        psio->open(fileno, PSIO_OPEN_OLD);
    psio->read_entry(fileno, const_cast<char *>(name_.c_str()), (char *)A2d_[0], sizeof(double) * dim1_ * dim2_);
    if (!already_open) psio->close(fileno, 1);  // Close and keep
}

void Array2d::read(std::shared_ptr<psi::PSIO> psio, size_t fileno) {
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno))
        already_open = true;
    else
        psio->open(fileno, PSIO_OPEN_OLD);
    psio->read_entry(fileno, const_cast<char *>(name_.c_str()), (char *)A2d_[0], sizeof(double) * dim1_ * dim2_);
    if (!already_open) psio->close(fileno, 1);  // Close and keep
}

void Array2d::read(psi::PSIO &psio, size_t fileno) { read(&psio, fileno); }  //

bool Array2d::read(PSIO *psio, int itap, const char *label, int dim) {
    int ntri = 0.5 * dim * (dim + 1);
    double *mybuffer = init_array(ntri);
    memset(mybuffer, 0, sizeof(double) * ntri);
    IWL::read_one(psio, itap, label, mybuffer, ntri, 0, 0, "outfile");

    double **Asq = block_matrix(dim, dim);
    memset(Asq[0], 0, sizeof(double) * dim * dim);
    tri_to_sq(mybuffer, Asq, dim);
    free(mybuffer);

    set(Asq);
    free_block(Asq);
    return true;
}  //

bool Array2d::read(std::shared_ptr<psi::PSIO> psio, int itap, const char *label, int dim) {
    int ntri = 0.5 * dim * (dim + 1);
    double *mybuffer = init_array(ntri);
    memset(mybuffer, 0, sizeof(double) * ntri);
    IWL::read_one(psio.get(), itap, label, mybuffer, ntri, 0, 0, "outfile");

    double **Asq = block_matrix(dim, dim);
    memset(Asq[0], 0, sizeof(double) * dim * dim);
    tri_to_sq(mybuffer, Asq, dim);
    free(mybuffer);

    set(Asq);
    free_block(Asq);
    return true;
}  //

void Array2d::save(std::shared_ptr<psi::PSIO> psio, size_t fileno) {
    write(psio, fileno);
    release();
}  //

void Array2d::save(psi::PSIO *const psio, size_t fileno) {
    write(psio, fileno);
    release();
}  //

void Array2d::save(psi::PSIO &psio, size_t fileno) {
    write(&psio, fileno);
    release();
}  //

void Array2d::load(std::shared_ptr<psi::PSIO> psio, size_t fileno, std::string name, int d1, int d2) {
    init(name, d1, d2);
    read(psio, fileno);
}  //

void Array2d::load(psi::PSIO *const psio, size_t fileno, std::string name, int d1, int d2) {
    init(name, d1, d2);
    read(psio, fileno);
}  //

void Array2d::load(psi::PSIO &psio, size_t fileno, std::string name, int d1, int d2) {
    init(name, d1, d2);
    read(&psio, fileno);
}  //

double **Array2d::to_block_matrix() {
    double **temp = block_matrix(dim1_, dim2_);
    memcpy(&(temp[0][0]), &(A2d_[0][0]), dim1_ * dim2_ * sizeof(double));
    return temp;
}  //

double *Array2d::to_lower_triangle() {
    if (dim1_ != dim2_) return NULL;
    int ntri = 0.5 * dim1_ * (dim1_ + 1);
    double *tri = new double[ntri];
    double **temp = to_block_matrix();
    sq_to_tri(temp, tri, dim1_);
    free_block(temp);
    return tri;
}  //

void Array2d::to_shared_matrix(SharedMatrix A) {
    for (int i = 0; i < dim1_; ++i) {
        for (int j = 0; j < dim2_; ++j) {
            A->set(0, i, j, A2d_[i][j]);
        }
    }
}  //

void Array2d::mgs() {
    double rmgs1, rmgs2;

    for (int k = 0; k < dim1_; k++) {  // loop-1
        rmgs1 = 0.0;

        for (int i = 0; i < dim1_; i++) {
            rmgs1 += A2d_[i][k] * A2d_[i][k];
        }

        rmgs1 = sqrt(rmgs1);

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

void Array2d::gs() {
    if (dim1_ != 0 && dim2_ != 0) {
        schmidt(A2d_, dim1_, dim2_, "outfile");
    }
}  //

double *Array2d::row_vector(int n) {
    double *temp = new double[dim2_];
    memset(temp, 0, dim2_ * sizeof(double));
    for (int i = 0; i < dim2_; i++) temp[i] = A2d_[n][i];
    return temp;
}  //

double *Array2d::column_vector(int n) {
    double *temp = new double[dim1_];
    memset(temp, 0, dim1_ * sizeof(double));
    for (int i = 0; i < dim1_; i++) temp[i] = A2d_[i][n];
    return temp;
}  //

void Array2d::sort1432(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx,
                       Array2i *row_pair_idx2, Array2i *col_pair_idx2) {
    for (int p = 0; p < d1; p++) {
        for (int q = 0; q < d2; q++) {
            int pq = row_pair_idx->A2i_[p][q];
            for (int r = 0; r < d3; r++) {
                for (int s = 0; s < d4; s++) {
                    int rs = col_pair_idx->A2i_[r][s];
                    int ps = row_pair_idx2->A2i_[p][s];
                    int rq = col_pair_idx2->A2i_[r][q];
                    A2d_[ps][rq] = A->A2d_[pq][rs];
                }
            }
        }
    }
}  //

void Array2d::sort2134(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx,
                       Array2i *row_pair_idx2) {
    for (int p = 0; p < d1; p++) {
        for (int q = 0; q < d2; q++) {
            int pq = row_pair_idx->A2i_[p][q];
            for (int r = 0; r < d3; r++) {
                for (int s = 0; s < d4; s++) {
                    int rs = col_pair_idx->A2i_[r][s];
                    int qp = row_pair_idx2->A2i_[p][s];
                    A2d_[qp][rs] = A->A2d_[pq][rs];
                }
            }
        }
    }
}  //

void Array2d::sort1243(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx,
                       Array2i *col_pair_idx2) {
    for (int p = 0; p < d1; p++) {
        for (int q = 0; q < d2; q++) {
            int pq = row_pair_idx->A2i_[p][q];
            for (int r = 0; r < d3; r++) {
                for (int s = 0; s < d4; s++) {
                    int rs = col_pair_idx->A2i_[r][s];
                    int sr = col_pair_idx2->A2i_[s][r];
                    A2d_[pq][sr] = A->A2d_[pq][rs];
                }
            }
        }
    }
}  //

void Array2d::sort2413(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx,
                       Array2i *row_pair_idx2, Array2i *col_pair_idx2) {
    for (int p = 0; p < d1; p++) {
        for (int q = 0; q < d2; q++) {
            int pq = row_pair_idx->A2i_[p][q];
            for (int r = 0; r < d3; r++) {
                for (int s = 0; s < d4; s++) {
                    int rs = col_pair_idx->A2i_[r][s];
                    int qs = row_pair_idx2->A2i_[q][s];
                    int pr = col_pair_idx2->A2i_[p][r];
                    A2d_[qs][pr] = A->A2d_[pq][rs];
                }
            }
        }
    }
}  //

void Array2d::sort2143(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx,
                       Array2i *row_pair_idx2, Array2i *col_pair_idx2) {
    for (int p = 0; p < d1; p++) {
        for (int q = 0; q < d2; q++) {
            int pq = row_pair_idx->A2i_[p][q];
            for (int r = 0; r < d3; r++) {
                for (int s = 0; s < d4; s++) {
                    int rs = col_pair_idx->A2i_[r][s];
                    int qp = row_pair_idx2->A2i_[q][p];
                    int sr = col_pair_idx2->A2i_[s][r];
                    A2d_[qp][sr] = A->A2d_[pq][rs];
                }
            }
        }
    }
}  //

void Array2d::sort4231(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx,
                       Array2i *row_pair_idx2, Array2i *col_pair_idx2) {
    for (int p = 0; p < d1; p++) {
        for (int q = 0; q < d2; q++) {
            int pq = row_pair_idx->A2i_[p][q];
            for (int r = 0; r < d3; r++) {
                for (int s = 0; s < d4; s++) {
                    int rs = col_pair_idx->A2i_[r][s];
                    int sq = row_pair_idx2->A2i_[s][q];
                    int rp = col_pair_idx2->A2i_[r][p];
                    A2d_[sq][rp] = A->A2d_[pq][rs];
                }
            }
        }
    }
}  //

void Array2d::sort3142(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx,
                       Array2i *row_pair_idx2, Array2i *col_pair_idx2) {
    for (int p = 0; p < d1; p++) {
        for (int q = 0; q < d2; q++) {
            int pq = row_pair_idx->A2i_[p][q];
            for (int r = 0; r < d3; r++) {
                for (int s = 0; s < d4; s++) {
                    int rs = col_pair_idx->A2i_[r][s];
                    int rp = row_pair_idx2->A2i_[r][p];
                    int sq = col_pair_idx2->A2i_[s][q];
                    A2d_[rp][sq] = A->A2d_[pq][rs];
                }
            }
        }
    }
}  //

void Array2d::sortp1432(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx,
                        Array2i *row_pair_idx2, Array2i *col_pair_idx2) {
    for (int p = 0; p < d1; p++) {
        for (int q = 0; q < d2; q++) {
            int pq = row_pair_idx->A2i_[p][q];
            for (int r = 0; r < d3; r++) {
                for (int s = 0; s < d4; s++) {
                    int rs = col_pair_idx->A2i_[r][s];
                    int ps = row_pair_idx2->A2i_[p][s];
                    int rq = col_pair_idx2->A2i_[r][q];
                    A2d_[ps][rq] += A->A2d_[pq][rs];
                }
            }
        }
    }
}  //

void Array2d::sortp2134(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx,
                        Array2i *row_pair_idx2) {
    for (int p = 0; p < d1; p++) {
        for (int q = 0; q < d2; q++) {
            int pq = row_pair_idx->A2i_[p][q];
            for (int r = 0; r < d3; r++) {
                for (int s = 0; s < d4; s++) {
                    int rs = col_pair_idx->A2i_[r][s];
                    int qp = row_pair_idx2->A2i_[p][s];
                    A2d_[qp][rs] += A->A2d_[pq][rs];
                }
            }
        }
    }
}  //

void Array2d::sortp1243(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx,
                        Array2i *col_pair_idx2) {
    for (int p = 0; p < d1; p++) {
        for (int q = 0; q < d2; q++) {
            int pq = row_pair_idx->A2i_[p][q];
            for (int r = 0; r < d3; r++) {
                for (int s = 0; s < d4; s++) {
                    int rs = col_pair_idx->A2i_[r][s];
                    int sr = col_pair_idx2->A2i_[s][r];
                    A2d_[pq][sr] += A->A2d_[pq][rs];
                }
            }
        }
    }
}  //

void Array2d::sortp2413(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx,
                        Array2i *row_pair_idx2, Array2i *col_pair_idx2) {
    for (int p = 0; p < d1; p++) {
        for (int q = 0; q < d2; q++) {
            int pq = row_pair_idx->A2i_[p][q];
            for (int r = 0; r < d3; r++) {
                for (int s = 0; s < d4; s++) {
                    int rs = col_pair_idx->A2i_[r][s];
                    int qs = row_pair_idx2->A2i_[q][s];
                    int pr = col_pair_idx2->A2i_[p][r];
                    A2d_[qs][pr] += A->A2d_[pq][rs];
                }
            }
        }
    }
}  //

void Array2d::sortp2143(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx,
                        Array2i *row_pair_idx2, Array2i *col_pair_idx2) {
    for (int p = 0; p < d1; p++) {
        for (int q = 0; q < d2; q++) {
            int pq = row_pair_idx->A2i_[p][q];
            for (int r = 0; r < d3; r++) {
                for (int s = 0; s < d4; s++) {
                    int rs = col_pair_idx->A2i_[r][s];
                    int qp = row_pair_idx2->A2i_[q][p];
                    int sr = col_pair_idx2->A2i_[s][r];
                    A2d_[qp][sr] += A->A2d_[pq][rs];
                }
            }
        }
    }
}  //

void Array2d::sortp4231(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx,
                        Array2i *row_pair_idx2, Array2i *col_pair_idx2) {
    for (int p = 0; p < d1; p++) {
        for (int q = 0; q < d2; q++) {
            int pq = row_pair_idx->A2i_[p][q];
            for (int r = 0; r < d3; r++) {
                for (int s = 0; s < d4; s++) {
                    int rs = col_pair_idx->A2i_[r][s];
                    int sq = row_pair_idx2->A2i_[s][q];
                    int rp = col_pair_idx2->A2i_[r][p];
                    A2d_[sq][rp] += A->A2d_[pq][rs];
                }
            }
        }
    }
}  //

void Array2d::sortp3142(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx,
                        Array2i *row_pair_idx2, Array2i *col_pair_idx2) {
    for (int p = 0; p < d1; p++) {
        for (int q = 0; q < d2; q++) {
            int pq = row_pair_idx->A2i_[p][q];
            for (int r = 0; r < d3; r++) {
                for (int s = 0; s < d4; s++) {
                    int rs = col_pair_idx->A2i_[r][s];
                    int rp = row_pair_idx2->A2i_[r][p];
                    int sq = col_pair_idx2->A2i_[s][q];
                    A2d_[rp][sq] += A->A2d_[pq][rs];
                }
            }
        }
    }
}  //

/********************************************************************************************/
/************************** 3d array ********************************************************/
/********************************************************************************************/
Array3d::Array3d(int d1, int d2, int d3) {
    A3d_ = NULL;
    dim1_ = d1;
    dim2_ = d2;
    dim3_ = d3;
    memalloc();
}  //

Array3d::Array3d(std::string name, int d1, int d2, int d3) {
    A3d_ = NULL;
    dim1_ = d1;
    dim2_ = d2;
    dim3_ = d3;
    name_ = name;
    memalloc();
}  //

Array3d::Array3d() {
    A3d_ = NULL;
    dim1_ = 0;
    dim2_ = 0;
    dim3_ = 0;

}  //

Array3d::~Array3d() { release(); }  //

Array3d *Array3d::generate(int d1, int d2, int d3) { return new Array3d(d1, d2, d3); }

Array3d *Array3d::generate(std::string name, int d1, int d2, int d3) { return new Array3d(name, d1, d2, d3); }

void Array3d::memalloc() {
    if (A3d_) release();
    A3d_ = (double ***)malloc(sizeof(double ***) * dim1_);
    for (int i = 0; i < dim1_; i++) {
        A3d_[i] = block_matrix(dim2_, dim3_);
    }
}  //

void Array3d::init(int d1, int d2, int d3) {
    dim1_ = d1;
    dim2_ = d2;
    dim3_ = d3;
    if (A3d_) release();
    A3d_ = (double ***)malloc(sizeof(double ***) * dim1_);
    for (int i = 0; i < dim1_; i++) {
        A3d_[i] = block_matrix(dim2_, dim3_);
    }
}  //

void Array3d::init(std::string name, int d1, int d2, int d3) {
    dim1_ = d1;
    dim2_ = d2;
    dim3_ = d3;
    name_ = name;
    if (A3d_) release();
    A3d_ = (double ***)malloc(sizeof(double ***) * dim1_);
    for (int i = 0; i < dim1_; i++) {
        A3d_[i] = block_matrix(dim2_, dim3_);
    }
}  //

void Array3d::zero() {
    for (int i = 0; i < dim1_; i++) {
        memset(&(A3d_[i][0][0]), 0, sizeof(double) * dim2_ * dim3_);
    }
}  //

void Array3d::print() {
    if (name_.length()) outfile->Printf("\n ## %s ##\n", name_.c_str());
    for (int i = 0; i < dim1_; i++) {
        outfile->Printf("\n Irrep: %d\n", i + 1);
        print_mat(A3d_[i], dim2_, dim3_, "outfile");
    }

}  //

void Array3d::release() {
    if (!A3d_) return;
    for (int i = 0; i < dim1_; i++) {
        free_block(A3d_[i]);
    }
    A3d_ = NULL;
}  //

void Array3d::set(int h, int i, int j, double value) { A3d_[h][i][j] = value; }  //

double Array3d::get(int h, int i, int j) { return A3d_[h][i][j]; }  //

/********************************************************************************************/
/************************** 1i array ********************************************************/
/********************************************************************************************/
Array1i::Array1i(int d1) {
    A1i_ = NULL;
    dim1_ = d1;
    memalloc();
}  //

Array1i::Array1i(std::string name, int d1) {
    A1i_ = NULL;
    dim1_ = d1;
    name_ = name;
    memalloc();
}  //

Array1i::Array1i() {
    A1i_ = NULL;
    dim1_ = 0;

}  //

Array1i::~Array1i() { release(); }  //

Array1i *Array1i::generate(int d1) { return new Array1i(d1); }

Array1i *Array1i::generate(std::string name, int d1) { return new Array1i(name, d1); }

void Array1i::memalloc() {
    if (A1i_) release();
    A1i_ = new int[dim1_];
}  //

void Array1i::init(int d1) {
    dim1_ = d1;
    if (A1i_) release();
    A1i_ = new int[dim1_];
}  //

void Array1i::init(std::string name, int d1) {
    dim1_ = d1;
    name_ = name;
    if (A1i_) release();
    A1i_ = new int[dim1_];
}  //

void Array1i::zero() { memset(A1i_, 0, sizeof(int) * dim1_); }  //

void Array1i::print() {
    if (name_.length()) outfile->Printf("\n ## %s ##\n", name_.c_str());
    for (int p = 0; p < dim1_; p++) {
        outfile->Printf(" %3d %3d \n", p, A1i_[p]);
    }

}  //

void Array1i::release() {
    if (!A1i_) return;
    delete[] A1i_;
    A1i_ = NULL;
}  //

void Array1i::set(int i, int value) { A1i_[i] = value; }  //

int Array1i::get(int i) { return A1i_[i]; }  //

void Array1i::add(const Array1i *Adum) {
    int *lhs, *rhs;
    size_t size = dim1_;
    if (size) {
        lhs = A1i_;
        rhs = Adum->A1i_;
        for (size_t ij = 0; ij < size; ++ij) {
            *lhs += *rhs;
            lhs++;
            rhs++;
        }
    }
}  //

void Array1i::add(int i, int value) { A1i_[i] += value; }  //

void Array1i::subtract(const Array1i *Adum) {
    int *lhs, *rhs;
    size_t size = dim1_;
    if (size) {
        lhs = A1i_;
        rhs = Adum->A1i_;
        for (size_t ij = 0; ij < size; ++ij) {
            *lhs -= *rhs;
            lhs++;
            rhs++;
        }
    }
}  //

void Array1i::subtract(int i, int value) { A1i_[i] -= value; }  //

/********************************************************************************************/
/************************** 2i array ********************************************************/
/********************************************************************************************/
Array2i::Array2i(int d1, int d2) {
    A2i_ = NULL;
    dim1_ = d1;
    dim2_ = d2;
    memalloc();
}  //

Array2i::Array2i(std::string name, int d1, int d2) {
    A2i_ = NULL;
    dim1_ = d1;
    dim2_ = d2;
    name_ = name;
    memalloc();
}  //

Array2i::Array2i() {
    A2i_ = NULL;
    dim1_ = 0;
    dim2_ = 0;

}  //

Array2i::~Array2i() { release(); }  //

Array2i *Array2i::generate(int d1, int d2) { return new Array2i(d1, d2); }

Array2i *Array2i::generate(std::string name, int d1, int d2) { return new Array2i(name, d1, d2); }

void Array2i::memalloc() {
    if (A2i_) release();
    A2i_ = init_int_matrix(dim1_, dim2_);
}  //

void Array2i::init(int d1, int d2) {
    dim1_ = d1;
    dim2_ = d2;
    if (A2i_) release();
    A2i_ = init_int_matrix(dim1_, dim2_);
}  //

void Array2i::init(std::string name, int d1, int d2) {
    dim1_ = d1;
    dim2_ = d2;
    name_ = name;
    if (A2i_) release();
    A2i_ = init_int_matrix(dim1_, dim2_);
}  //

void Array2i::zero() { memset(A2i_[0], 0, sizeof(int) * dim1_ * dim2_); }  //

void Array2i::zero_diagonal() {
    if (dim1_ == dim2_) {
        for (int i = 0; i < dim1_; i++) A2i_[i][i] = 0.0;
    }
}  //

void Array2i::print() {
    if (name_.length()) outfile->Printf("\n ## %s ##\n", name_.c_str());
    print_int_mat(A2i_, dim1_, dim2_, "outfile");

}  //

void Array2i::print(std::string out) {
    std::shared_ptr<psi::PsiOutStream> printer =
        (out == "outfile" ? outfile : std::shared_ptr<PsiOutStream>(new PsiOutStream(out)));
    if (name_.length()) printer->Printf("\n ## %s ##\n", name_.c_str());
    print_int_mat(A2i_, dim1_, dim2_, out);
}  //

void Array2i::release() {
    if (!A2i_) return;
    free_int_matrix(A2i_);
    A2i_ = NULL;
}  //

void Array2i::set(int i, int j, int value) { A2i_[i][j] = value; }  //

void Array2i::set(int **A) {
    if (A == NULL) return;
    for (int i = 0; i < dim1_; ++i) {
        for (int j = 0; j < dim2_; ++j) {
            A2i_[i][j] = A[i][j];
        }
    }
}  //

double Array2i::get(int i, int j) { return A2i_[i][j]; }  //

void Array2i::add(const Array2i *Adum) {
    int *lhs, *rhs;
    size_t size = dim1_ * dim2_;
    if (size) {
        lhs = A2i_[0];
        rhs = Adum->A2i_[0];
        for (size_t ij = 0; ij < size; ++ij) {
            *lhs += *rhs;
            lhs++;
            rhs++;
        }
    }
}  //

void Array2i::add(int i, int j, int value) { A2i_[i][j] += value; }  //

void Array2i::subtract(const Array2i *Adum) {
    int *lhs, *rhs;
    size_t size = dim1_ * dim2_;
    if (size) {
        lhs = A2i_[0];
        rhs = Adum->A2i_[0];
        for (size_t ij = 0; ij < size; ++ij) {
            *lhs -= *rhs;
            lhs++;
            rhs++;
        }
    }
}  //

void Array2i::subtract(int i, int j, int value) { A2i_[i][j] -= value; }  //

Array2i *Array2i::transpose() {
    Array2i *temp;
    temp = new Array2i(dim2_, dim1_);
    temp->zero();

    for (int i = 0; i < dim2_; ++i) {
        for (int j = 0; j < dim1_; ++j) {
            temp->set(i, j, A2i_[j][i]);
        }
    }

    return temp;
}  //

void Array2i::copy(const Array2i *Adum) {
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
    if (dim1_ != 0 && dim2_ != 0) {
        memcpy(A2i_[0], Adum->A2i_[0], dim1_ * dim2_ * sizeof(int));
    }
}  //

void Array2i::copy(int **a) {
    size_t size = dim1_ * dim2_ * sizeof(int);
    if (size) memcpy(&(A2i_[0][0]), &(a[0][0]), size);
}

void Array2i::identity() {
    zero();
    for (int i = 0; i < dim1_; ++i) A2i_[i][i] = 1.0;
}  //

int Array2i::trace() {
    int value = 0;
    for (int i = 0; i < dim1_; ++i) value += A2i_[i][i];
    return value;
}  //

int **Array2i::to_int_matrix() {
    int **temp = init_int_matrix(dim1_, dim2_);
    memcpy(&(temp[0][0]), &(A2i_[0][0]), dim1_ * dim2_ * sizeof(int));
    return temp;
}  //

/********************************************************************************************/
/************************** 3i array ********************************************************/
/********************************************************************************************/
Array3i::Array3i(int d1, int d2, int d3) {
    A3i_ = NULL;
    dim1_ = d1;
    dim2_ = d2;
    dim3_ = d3;
    memalloc();
}  //

Array3i::Array3i(std::string name, int d1, int d2, int d3) {
    A3i_ = NULL;
    dim1_ = d1;
    dim2_ = d2;
    dim3_ = d3;
    name_ = name;
    memalloc();
}  //

Array3i::Array3i() {
    A3i_ = NULL;
    dim1_ = 0;
    dim2_ = 0;
    dim3_ = 0;

}  //

Array3i::~Array3i() { release(); }  //

Array3i *Array3i::generate(int d1, int d2, int d3) { return new Array3i(d1, d2, d3); }  //

Array3i *Array3i::generate(std::string name, int d1, int d2, int d3) { return new Array3i(name, d1, d2, d3); }  //

void Array3i::memalloc() {
    if (A3i_) release();
    A3i_ = (int ***)malloc(sizeof(int ***) * dim1_);
    for (int i = 0; i < dim1_; i++) {
        A3i_[i] = init_int_matrix(dim2_, dim3_);
    }
}  //

void Array3i::init(int d1, int d2, int d3) {
    dim1_ = d1;
    dim2_ = d2;
    dim3_ = d3;
    if (A3i_) release();
    A3i_ = (int ***)malloc(sizeof(int ***) * dim1_);
    for (int i = 0; i < dim1_; i++) {
        A3i_[i] = init_int_matrix(dim2_, dim3_);
    }
}  //

void Array3i::init(std::string name, int d1, int d2, int d3) {
    dim1_ = d1;
    dim2_ = d2;
    dim3_ = d3;
    name_ = name;
    if (A3i_) release();
    A3i_ = (int ***)malloc(sizeof(int ***) * dim1_);
    for (int i = 0; i < dim1_; i++) {
        A3i_[i] = init_int_matrix(dim2_, dim3_);
    }
}  //

void Array3i::zero() {
    for (int i = 0; i < dim1_; i++) {
        memset(&(A3i_[i][0][0]), 0, sizeof(int) * dim2_ * dim3_);
    }
}  //

void Array3i::print() {
    if (name_.length()) outfile->Printf("\n ## %s ##\n", name_.c_str());
    for (int i = 0; i < dim1_; i++) {
        outfile->Printf("\n Irrep: %d\n", i + 1);
        print_int_mat(A3i_[i], dim2_, dim3_, "outfile");
    }

}  //

void Array3i::release() {
    if (!A3i_) return;
    for (int i = 0; i < dim1_; i++) {
        free_int_matrix(A3i_[i]);
    }
    A3i_ = NULL;
}  //

void Array3i::set(int h, int i, int j, int value) { A3i_[h][i][j] = value; }  //

int Array3i::get(int h, int i, int j) { return A3i_[h][i][j]; }  //

/********************************************************************************************/
/********************************************************************************************/
}  // namespace dfoccwave
}  // namespace psi
