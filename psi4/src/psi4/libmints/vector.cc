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

#include "psi4/libqt/qt.h"
#include "matrix.h"
#include "vector.h"
#include "dimension.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <algorithm>
#include <numeric>
#include <cstdlib>
#include <string.h>
#include <math.h>

namespace psi {

Vector::Vector() {
    dimpi_ = nullptr;
    nirrep_ = 0;
    name_ = "";
}

Vector::Vector(const Vector &c) {
    nirrep_ = c.nirrep_;
    dimpi_ = c.dimpi_;
    alloc();
    copy_from(c);
    name_ = c.name_;
}

Vector::Vector(int nirreps, int *dimpi) : dimpi_(nirreps) {
    nirrep_ = nirreps;
    dimpi_ = dimpi;
    alloc();
}

Vector::Vector(int dim) : dimpi_(1) {
    nirrep_ = 1;
    dimpi_[0] = dim;
    alloc();
}

Vector::Vector(const std::string &name, int nirreps, int *dimpi) : dimpi_(nirreps) {
    nirrep_ = nirreps;
    dimpi_ = new int[nirrep_];
    for (int h = 0; h < nirrep_; ++h) dimpi_[h] = dimpi[h];
    alloc();
    name_ = name;
}

Vector::Vector(const std::string &name, int dim) : dimpi_(1) {
    nirrep_ = 1;
    dimpi_[0] = dim;
    alloc();
    name_ = name;
}

Vector::Vector(const Dimension &v) {
    nirrep_ = v.n();
    dimpi_ = v;
    alloc();
    name_ = v.name();
}

Vector::Vector(const std::string &name, const Dimension &v) {
    nirrep_ = v.n();
    dimpi_ = v;
    alloc();
    name_ = name;
}

Vector::~Vector() { release(); }

SharedVector Vector::create(const std::string &name, const Dimension &dim) {
    return std::make_shared<Vector>(name, dim);
}

void Vector::init(int nirreps, int *dimpi) {
    dimpi_.init(nirreps);
    nirrep_ = nirreps;
    dimpi_ = dimpi;
    alloc();
}

void Vector::init(int nirreps, const int *dimpi, const std::string &name) {
    name_ = name;
    dimpi_.init(nirreps);
    dimpi_ = dimpi;
    alloc();
}

void Vector::init(const Dimension &v) {
    name_ = v.name();
    nirrep_ = v.n();
    dimpi_ = v;
    alloc();
}

Vector *Vector::clone() {
    Vector *temp = new Vector(dimpi_);
    temp->copy(this);
    return temp;
}

void Vector::alloc() {
    if (vector_.size()) release();

    int total = dimpi_.sum();
    v_.resize(total);

    std::fill(vector_.begin(), vector_.end(), (double *)0);
    std::fill(v_.begin(), v_.end(), 0.0);

    assign_pointer_offsets();
}

void Vector::assign_pointer_offsets() {
    // Resize just to be sure it's the correct size
    vector_.resize(dimpi_.n(), 0);

    size_t offset = 0;
    for (int h = 0; h < nirrep_; ++h) {
        if (dimpi_[h])
            vector_[h] = &(v_[0]) + offset;
        else
            vector_[h] = nullptr;
        offset += dimpi_[h];
    }
}

void Vector::release() {
    std::fill(vector_.begin(), vector_.end(), (double *)0);
    std::fill(v_.begin(), v_.end(), 0.0);
}

void Vector::copy_from(const Vector &other) {
    nirrep_ = other.dimpi_.n();
    dimpi_ = other.dimpi_;
    v_ = other.v_;
    assign_pointer_offsets();
}

void Vector::copy(const Vector *rhs) { copy_from(*rhs); }

void Vector::copy(const Vector &rhs) { copy_from(rhs); }

void Vector::set(double *vec) { std::copy(vec, vec + dimpi_.sum(), v_.begin()); }

SharedVector Vector::get_block(const Slice &slice) {
    // check if slice is within bounds
    for (int h = 0; h < nirrep_; h++) {
        if (slice.end()[h] > dimpi_[h]) {
            std::string msg =
                "Invalid call to Vector::get_block(): Slice is out of bounds. Irrep = " + std::to_string(h);
            throw PSIEXCEPTION(msg);
        }
    }
    const Dimension &slice_begin = slice.begin();
    Dimension slice_dim = slice.end() - slice.begin();
    auto block = std::make_shared<Vector>("Block", slice_dim);
    for (int h = 0; h < nirrep_; h++) {
        int max_p = slice_dim[h];
        for (int p = 0; p < max_p; p++) {
            double value = get(h, p + slice_begin[h]);
            block->set(h, p, value);
        }
    }
    return block;
}

void Vector::set_block(const Slice &slice, SharedVector block) {
    // check if slice is within bounds
    for (int h = 0; h < nirrep_; h++) {
        if (slice.end()[h] > dimpi_[h]) {
            std::string msg =
                "Invalid call to Vector::set_block(): Slice is out of bounds. Irrep = " + std::to_string(h);
            throw PSIEXCEPTION(msg);
        }
    }
    const Dimension &slice_begin = slice.begin();
    Dimension slice_dim = slice.end() - slice.begin();
    for (int h = 0; h < nirrep_; h++) {
        int max_p = slice_dim[h];
        for (int p = 0; p < max_p; p++) {
            double value = block->get(h, p);
            set(h, p + slice_begin[h], value);
        }
    }
}

void Vector::zero() { std::fill(v_.begin(), v_.end(), 0.0); }

void Vector::print(std::string out, const char *extra) const {
    std::shared_ptr<psi::PsiOutStream> printer =
        (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));
    int h;
    if (extra == nullptr) {
        printer->Printf("\n # %s #\n", name_.c_str());
    } else {
        printer->Printf("\n # %s %s #\n", name_.c_str(), extra);
    }
    for (h = 0; h < nirrep_; ++h) {
        printer->Printf(" Irrep: %d\n", h + 1);
        for (int i = 0; i < dimpi_[h]; ++i) printer->Printf("   %4d: %10.7f\n", i + 1, vector_[h][i]);
        printer->Printf("\n");
    }
}

double *Vector::to_block_vector() {
    size_t size = 0;
    for (int h = 0; h < nirrep_; ++h) size += dimpi_[h];

    double *temp = new double[size];
    size_t offset = 0;
    for (int h = 0; h < nirrep_; ++h) {
        for (int i = 0; i < dimpi_[h]; ++i) {
            temp[i + offset] = vector_[h][i];
        }
        offset += dimpi_[h];
    }

    return temp;
}

void Vector::gemv(bool transa, double alpha, Matrix *A, Vector *X, double beta) {
    char trans = transa ? 't' : 'n';

    for (int h = 0; h < nirrep_; ++h) {
        C_DGEMV(trans, A->rowspi_[h], A->colspi_[h], alpha, &(A->matrix_[h][0][0]), A->rowspi_[h], &(X->vector_[h][0]),
                1, beta, &(vector_[h][0]), 1);
    }
}

double Vector::vector_dot(const SharedVector &other) { return vector_dot(*other.get()); }
double Vector::vector_dot(const Vector &other) {
    if (v_.size() != other.v_.size()) {
        throw PSIEXCEPTION("Vector::vector_dot: Vector sizes do not match!");
    }

    return C_DDOT(v_.size(), v_.data(), 1, const_cast<double *>(other.v_.data()), 1);
}
double Vector::dot(Vector *X) {
    if (v_.size() != X->v_.size()) {
        throw PSIEXCEPTION("Vector::vector_dot: Vector sizes do not match!");
    }

    return C_DDOT(v_.size(), v_.data(), 1, X->v_.data(), 1);
}

double Vector::sum_of_squares() { return C_DDOT(v_.size(), v_.data(), 1, v_.data(), 1); }
double Vector::norm() { return sqrt(sum_of_squares()); }
double Vector::rms() { return sqrt(sum_of_squares() / v_.size()); }

void Vector::scale(const double &sc) { C_DSCAL(v_.size(), sc, v_.data(), 1); }

void Vector::add(const SharedVector &other) { axpy(1.0, *other.get()); }
void Vector::add(const Vector &other) { axpy(1.0, other); }

void Vector::subtract(const SharedVector &other) { axpy(-1.0, *other.get()); }
void Vector::subtract(const Vector &other) { axpy(-1.0, other); }

void Vector::axpy(double scale, const SharedVector &other) { axpy(scale, *other.get()); }
void Vector::axpy(double scale, const Vector &other) {
    if (v_.size() != other.v_.size()) {
        throw PSIEXCEPTION("Vector::axpy: Vector sizes do not match!");
    }

    C_DAXPY(v_.size(), scale, const_cast<double *>(other.v_.data()), 1, v_.data(), 1);
}

}  // namespace psi
