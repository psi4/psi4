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

// Latest revision on April 25, 2013.

#include "arrays.h"

#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libqt/qt.h"

#include <cstdio>
#include <cmath>

using namespace psi;
namespace psi {
namespace occwave {

/********************************************************************************************/
/************************** 1d array ********************************************************/
/********************************************************************************************/
Array1d::Array1d(int d1) {
    A1d_ = nullptr;
    dim1_ = d1;
    memalloc();
}  //

Array1d::Array1d(std::string name, int d1) {
    A1d_ = nullptr;
    dim1_ = d1;
    name_ = name;
    memalloc();
}  //

Array1d::Array1d() {
    A1d_ = nullptr;
    dim1_ = 0;

}  //

Array1d::~Array1d() { release(); }  //

Array1d* Array1d::generate(int d1) { return new Array1d(d1); }

Array1d* Array1d::generate(std::string name, int d1) { return new Array1d(name, d1); }

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

void Array1d::print(std::string out_fname) {
    auto mode = std::ostream::app;
    auto printer = out_fname == "outfile" ? outfile : std::make_shared<PsiOutStream>(out_fname, mode);
    if (name_.length()) printer->Printf("\n ## %s ##\n", name_.c_str());
    for (int p = 0; p < dim1_; p++) {
        printer->Printf(" %3d %10.7f \n", p, A1d_[p]);
    }
}  //

void Array1d::release() {
    if (!A1d_) return;
    delete[] A1d_;
    A1d_ = nullptr;
}  //

void Array1d::set(int i, double value) { A1d_[i] = value; }  //

void Array1d::set(double* vec) {
    for (int i = 0; i < dim1_; ++i) A1d_[i] = vec[i];
}  //

void Array1d::set(const Array1d* vec) {
    for (int i = 0; i < dim1_; ++i) A1d_[i] = vec->A1d_[i];
}  //

double Array1d::get(int i) { return A1d_[i]; }  //

void Array1d::add(const Array1d* Adum) {
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

void Array1d::subtract(const Array1d* Adum) {
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
    summ = std::sqrt(summ) / dim1_;

    return summ;
}  //

double Array1d::rms(const Array1d* Atemp) {
    double summ = 0.0;
    for (int i = 0; i < dim1_; ++i) summ += (A1d_[i] - Atemp->A1d_[i]) * (A1d_[i] - Atemp->A1d_[i]);
    summ = std::sqrt(summ) / dim1_;

    return summ;
}  //

double Array1d::dot(const Array1d* y) {
    double value = 0.0;
    int incx = 1;
    int incy = 1;
    if (dim1_ == y->dim1_) value = C_DDOT(dim1_, A1d_, incx, y->A1d_, incy);
    return value;
}  //

void Array1d::gemv(bool transa, double alpha, const Array2d* a, const Array1d* b, double beta) {
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

double Array1d::xay(const Array2d* a, const Array1d* y) {
    double value = 0.0;
    auto* ay = new Array1d(a->dim1_);
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

void Array1d::copy(double* x) {
    size_t size;
    size = dim1_ * sizeof(double);
    if (size) memcpy(&(A1d_[0]), &(x[0]), size);
}  //

void Array1d::copy(const Array1d* x) {
    size_t size;
    size = dim1_ * sizeof(double);
    if (size) memcpy(&(A1d_[0]), &(x->A1d_[0]), size);
}  //

void Array1d::dirprd(Array1d* a, Array1d* b) {
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
    A2d_ = nullptr;
    dim1_ = d1;
    dim2_ = d2;
    memalloc();
}  //

Array2d::Array2d(std::string name, int d1, int d2) {
    A2d_ = nullptr;
    dim1_ = d1;
    dim2_ = d2;
    name_ = name;
    memalloc();
}  //

Array2d::Array2d() {
    A2d_ = nullptr;
    dim1_ = 0;
    dim2_ = 0;

}  //

Array2d::~Array2d() { release(); }  //

Array2d* Array2d::generate(int d1, int d2) { return new Array2d(d1, d2); }

Array2d* Array2d::generate(std::string name, int d1, int d2) { return new Array2d(name, d1, d2); }

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
    if (name_.length()) outfile->Printf("\n ## %s ##\n", name_.c_str());
    print_mat(A2d_, dim1_, dim2_, "outfile");

}  //

void Array2d::print(std::string out_fname) {
    auto mode = std::ostream::app;
    auto printer = out_fname == "outfile" ? outfile : std::make_shared<PsiOutStream>(out_fname, mode);
    if (name_.length()) printer->Printf("\n ## %s ##\n", name_.c_str());
    print_mat(A2d_, dim1_, dim2_, out_fname);
}  //

void Array2d::release() {
    if (!A2d_) return;
    free_block(A2d_);
    A2d_ = nullptr;
}  //

void Array2d::set(int i, int j, double value) { A2d_[i][j] = value; }  //

void Array2d::set(double** A) {
    if (A == nullptr) return;
    for (int i = 0; i < dim1_; ++i) {
        for (int j = 0; j < dim2_; ++j) {
            A2d_[i][j] = A[i][j];
        }
    }
}  //

double Array2d::get(int i, int j) { return A2d_[i][j]; }  //

void Array2d::gemm(bool transa, bool transb, double alpha, const Array2d* a, const Array2d* b, double beta) {
    char ta = transa ? 't' : 'n';
    char tb = transb ? 't' : 'n';
    int m, n, k, nca, ncb, ncc;

    m = dim1_;
    n = dim2_;
    k = a->dim2_;
    nca = transa ? m : k;
    ncb = transb ? k : n;
    ncc = n;

    if (m && n && k) {
        C_DGEMM(ta, tb, m, n, k, alpha, &(a->A2d_[0][0]), nca, &(b->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
    }
}  //

void Array2d::add(const Array2d* Adum) {
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

void Array2d::add(int i, int j, double value) { A2d_[i][j] += value; }  //

void Array2d::subtract(const Array2d* Adum) {
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

Array2d* Array2d::transpose() {
    Array2d* temp;
    temp = new Array2d(dim2_, dim1_);
    temp->zero();

    for (int i = 0; i < dim2_; ++i) {
        for (int j = 0; j < dim1_; ++j) {
            temp->set(i, j, A2d_[j][i]);
        }
    }

    return temp;
}  //

void Array2d::copy(const Array2d* Adum) {
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
        memcpy(&(A2d_[0]), &(Adum->A2d_[0]), dim1_ * dim2_ * sizeof(double));
    }
}  //

void Array2d::copy(double** a) {
    size_t size = dim1_ * dim2_ * sizeof(double);
    if (size) memcpy(&(A2d_[0][0]), &(a[0][0]), size);
}

void Array2d::cdgesv(Array1d* Xvec, int errcod) {
    if (dim1_) {
        int* ipiv = init_int_array(dim1_);
        memset(ipiv, 0, sizeof(int) * dim1_);
        errcod = 0;
        errcod = C_DGESV(dim1_, 1, &(A2d_[0][0]), dim2_, &(ipiv[0]), Xvec->A1d_, dim2_);
        free(ipiv);
    }
}  //

void Array2d::lineq_flin(Array1d* Xvec, double* det) {
    if (dim1_) {
        flin(A2d_, Xvec->A1d_, dim1_, 1, det);
    }
}  //

void Array2d::lineq_pople(Array1d* Xvec, int num_vecs, double cutoff) {
    if (dim1_) {
        pople(A2d_, Xvec->A1d_, dim1_, num_vecs, cutoff, "outfile", 0);
    }
}  //

// TODO:
// DGER compute the rank-one update of a general matrix: A <-- A + alpha * x * yT
// dger(m, n, alpha, x, incx, y, incy, a, lda);

/*
void Array2d::write(psi::PSIO* psio, size_t fileno)
{
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno)) already_open = true;
    else psio->open(fileno, PSIO_OPEN_OLD);}
    psio->write_entry(fileno, const_cast<char*>(name_.c_str()), (char*)A2d_[0], sizeof(double) * dim1_ * dim2_);
    if (!already_open) psio->close(fileno, 1);     // Close and keep
}//

void Array2d::write(shared_ptr<psi::PSIO> psio, size_t fileno)
{
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno)) already_open = true;
    else psio->open(fileno, PSIO_OPEN_NEW);}
    psio->write_entry(fileno, const_cast<char*>(name_.c_str()), (char*)A2d_[0], sizeof(double) * dim1_ * dim2_);
    if (!already_open) psio->close(fileno, 1);     // Close and keep
}//

void Array2d::write(psi::PSIO& psio, size_t fileno)
{
    write(&psio, fileno);
}//

void Array2d::read(psi::PSIO* psio, size_t fileno)
{
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno)) already_open = true;
    else psio->open(fileno, PSIO_OPEN_NEW);
    psio->read_entry(fileno, const_cast<char*>(name_.c_str()), (char*)A2d_[0], sizeof(double) * dim1_ * dim2_);
    if (!already_open) psio->close(fileno, 1);     // Close and keep
}

void Array2d::read(shared_ptr<psi::PSIO> psio, size_t fileno)
{
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno)) already_open = true;
    else psio->open(fileno, PSIO_OPEN_NEW);
    psio->read_entry(fileno, const_cast<char*>(name_.c_str()), (char*)A2d_[0], sizeof(double) * dim1_ * dim2_);
    if (!already_open) psio->close(fileno, 1);     // Close and keep
}

void Array2d::read(psi::PSIO& psio, size_t fileno)
{
    read(&psio, fileno);
}//

bool Array2d::read(PSIO* psio, int itap, const char *label, int dim)
{
    double *mybuffer = init_array(ntri);
    memset(mybuffer, 0, sizeof(double)*ntri);
    IWL::read_one(psio, itap, label, mybuffer, ntri, 0, 0, outfile);

    double **Asq = block_matrix(dim, dim);
    memset(Asq[0], 0, sizeof(double)*dim*dim);
    tri_to_sq(mybuffer,Asq,dim);
    free(mybuffer);

    set(Asq);
    free_block(Asq);
    return true;
}//

bool Array2d::read(shared_ptr<psi::PSIO> psio, int itap, const char *label, int dim)
{
    double *mybuffer = init_array(ntri);
    memset(mybuffer, 0, sizeof(double)*ntri);
    IWL::read_one(psio.get(), itap, label, mybuffer, ntri, 0, 0, outfile);

    double **Asq = block_matrix(dim, dim);
    memset(Asq[0], 0, sizeof(double)*dim*dim);
    tri_to_sq(mybuffer,Asq,dim);
    free(mybuffer);

    set(Asq);
    free_block(Asq);
    return true;
}//
*/

/********************************************************************************************/
/************************** 3d array ********************************************************/
/********************************************************************************************/
Array3d::Array3d(int d1, int d2, int d3) {
    A3d_ = nullptr;
    dim1_ = d1;
    dim2_ = d2;
    dim3_ = d3;
    memalloc();
}  //

Array3d::Array3d(std::string name, int d1, int d2, int d3) {
    A3d_ = nullptr;
    dim1_ = d1;
    dim2_ = d2;
    dim3_ = d3;
    name_ = name;
    memalloc();
}  //

Array3d::Array3d() {
    A3d_ = nullptr;
    dim1_ = 0;
    dim2_ = 0;
    dim3_ = 0;

}  //

Array3d::~Array3d() { release(); }  //

Array3d* Array3d::generate(int d1, int d2, int d3) { return new Array3d(d1, d2, d3); }

Array3d* Array3d::generate(std::string name, int d1, int d2, int d3) { return new Array3d(name, d1, d2, d3); }

void Array3d::memalloc() {
    if (A3d_) release();
    A3d_ = (double***)malloc(sizeof(double***) * dim1_);
    for (int i = 0; i < dim1_; i++) {
        A3d_[i] = block_matrix(dim2_, dim3_);
    }
}  //

void Array3d::init(int d1, int d2, int d3) {
    dim1_ = d1;
    dim2_ = d2;
    dim3_ = d3;
    if (A3d_) release();
    A3d_ = (double***)malloc(sizeof(double***) * dim1_);
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
    A3d_ = (double***)malloc(sizeof(double***) * dim1_);
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
    A3d_ = nullptr;
}  //

void Array3d::set(int h, int i, int j, double value) { A3d_[h][i][j] = value; }  //

double Array3d::get(int h, int i, int j) { return A3d_[h][i][j]; }  //

/********************************************************************************************/
/************************** 1i array ********************************************************/
/********************************************************************************************/
Array1i::Array1i(int d1) {
    A1i_ = nullptr;
    dim1_ = d1;
    memalloc();
}  //

Array1i::Array1i(std::string name, int d1) {
    A1i_ = nullptr;
    dim1_ = d1;
    name_ = name;
    memalloc();
}  //

Array1i::Array1i() {
    A1i_ = nullptr;
    dim1_ = 0;

}  //

Array1i::~Array1i() { release(); }  //

Array1i* Array1i::generate(int d1) { return new Array1i(d1); }

Array1i* Array1i::generate(std::string name, int d1) { return new Array1i(name, d1); }

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
    A1i_ = nullptr;
}  //

void Array1i::set(int i, int value) { A1i_[i] = value; }  //

int Array1i::get(int i) { return A1i_[i]; }  //

void Array1i::add(const Array1i* Adum) {
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

void Array1i::subtract(const Array1i* Adum) {
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
/************************** 3i array ********************************************************/
/********************************************************************************************/
Array3i::Array3i(int d1, int d2, int d3) {
    A3i_ = nullptr;
    dim1_ = d1;
    dim2_ = d2;
    dim3_ = d3;
    memalloc();
}  //

Array3i::Array3i(std::string name, int d1, int d2, int d3) {
    A3i_ = nullptr;
    dim1_ = d1;
    dim2_ = d2;
    dim3_ = d3;
    name_ = name;
    memalloc();
}  //

Array3i::Array3i() {
    A3i_ = nullptr;
    dim1_ = 0;
    dim2_ = 0;
    dim3_ = 0;

}  //

Array3i::~Array3i() { release(); }  //

Array3i* Array3i::generate(int d1, int d2, int d3) { return new Array3i(d1, d2, d3); }  //

Array3i* Array3i::generate(std::string name, int d1, int d2, int d3) { return new Array3i(name, d1, d2, d3); }  //

void Array3i::memalloc() {
    if (A3i_) release();
    A3i_ = (int***)malloc(sizeof(int***) * dim1_);
    for (int i = 0; i < dim1_; i++) {
        A3i_[i] = init_int_matrix(dim2_, dim3_);
    }
}  //

void Array3i::init(int d1, int d2, int d3) {
    dim1_ = d1;
    dim2_ = d2;
    dim3_ = d3;
    if (A3i_) release();
    A3i_ = (int***)malloc(sizeof(int***) * dim1_);
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
    A3i_ = (int***)malloc(sizeof(int***) * dim1_);
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
    A3i_ = nullptr;
}  //

void Array3i::set(int h, int i, int j, int value) { A3i_[h][i][j] = value; }  //

int Array3i::get(int h, int i, int j) { return A3i_[h][i][j]; }  //

/********************************************************************************************/
/********************************************************************************************/
}
Vector::Vector(const Dimension& dimpi, const occwave::Array1d& array) {
    nirrep_ = dimpi.n();
    dimpi_ = dimpi;
    name_ = array.name();
    v_ = std::vector<double>(array.array(), array.array() + array.dim1());
    assign_pointer_offsets();
}
}  // End Namespaces
