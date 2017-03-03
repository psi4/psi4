/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/libqt/qt.h"
#include "matrix.h"
#include "vector.h"
#include "dimension.h"
#include "psi4/libparallel/ParallelPrinter.h"
#include "psi4/libparallel/mpi_wrapper.h"
#include "psi4/libparallel/local.h"

#include "psi4/pybind11.h"

#include <algorithm>
#include <numeric>
#include <cstdlib>
#include <string.h>

namespace psi {

Vector::Vector()
{
    dimpi_ = NULL;
    nirrep_ = 0;
    name_ = "";
}

Vector::Vector(const Vector &c)
{
    nirrep_ = c.nirrep_;
    dimpi_ = c.dimpi_;
    alloc();
    copy_from(c);
    name_ = c.name_;
}

Vector::Vector(int nirreps, int *dimpi)
        : dimpi_(nirreps)
{
    nirrep_ = nirreps;
    dimpi_ = dimpi;
    alloc();
}

Vector::Vector(int dim)
        : dimpi_(1)
{
    nirrep_ = 1;
    dimpi_[0] = dim;
    alloc();
}

Vector::Vector(const std::string &name, int nirreps, int *dimpi)
        : dimpi_(nirreps)
{
    nirrep_ = nirreps;
    dimpi_ = new int[nirrep_];
    for (int h = 0; h < nirrep_; ++h)
        dimpi_[h] = dimpi[h];
    alloc();
    name_ = name;
}

Vector::Vector(const std::string &name, int dim)
        : dimpi_(1)
{
    nirrep_ = 1;
    dimpi_[0] = dim;
    alloc();
    name_ = name;
}

Vector::Vector(const Dimension &v)
{
    nirrep_ = v.n();
    dimpi_ = v;
    alloc();
    name_ = v.name();
}

Vector::Vector(const std::string &name, const Dimension &v)
{
    nirrep_ = v.n();
    dimpi_ = v;
    alloc();
    name_ = name;
}

Vector::~Vector()
{
    release();
}

SharedVector Vector::create(const std::string &name, const Dimension &dim)
{
    return SharedVector(new Vector(name, dim));
}

void Vector::init(int nirreps, int *dimpi)
{
    dimpi_.init(nirreps);
    nirrep_ = nirreps;
    dimpi_ = dimpi;
    alloc();
}

void Vector::init(int nirreps, const int *dimpi, const std::string &name)
{
    name_ = name;
    dimpi_.init(nirreps);
    dimpi_ = dimpi;
    alloc();
}

void Vector::init(const Dimension &v)
{
    name_ = v.name();
    nirrep_ = v.n();
    dimpi_ = v;
    alloc();
}

Vector *Vector::clone()
{
    Vector *temp = new Vector(dimpi_);
    temp->copy(this);
    return temp;
}

void Vector::alloc()
{
    if (vector_.size())
        release();

    int total = dimpi_.sum();
    v_.resize(total);

    std::fill(vector_.begin(), vector_.end(), (double *) 0);
    std::fill(v_.begin(), v_.end(), 0.0);

    assign_pointer_offsets();
}

void Vector::assign_pointer_offsets()
{
    // Resize just to be sure it's the correct size
    vector_.resize(dimpi_.n(), 0);

    size_t offset = 0;
    for (int h = 0; h < nirrep_; ++h) {
        if (dimpi_[h])
            vector_[h] = &(v_[0]) + offset;
        else
            vector_[h] = NULL;
        offset += dimpi_[h];
    }
}

void Vector::release()
{
    std::fill(vector_.begin(), vector_.end(), (double *) 0);
    std::fill(v_.begin(), v_.end(), 0.0);
}

void Vector::copy_from(const Vector &other)
{
    nirrep_ = other.dimpi_.n();
    dimpi_ = other.dimpi_;
    v_ = other.v_;
    assign_pointer_offsets();
}

void Vector::copy(const Vector *rhs)
{
    copy_from(*rhs);
}

void Vector::copy(const Vector &rhs)
{
    copy_from(rhs);
}

void Vector::set(double *vec)
{
    std::copy(vec, vec + dimpi_.sum(), v_.begin());
}

void Vector::zero()
{
    std::fill(v_.begin(), v_.end(), 0.0);
}

double Vector::pyget(const py::tuple &key)
{
    int h = 0, elem = 0;
    h = key[0].cast<int>();
    elem = key[1].cast<int>();

    return get(h, elem);
}

void Vector::pyset(const py::tuple &key, double value)
{
    int h = 0, elem = 0;
    h = key[0].cast<int>();
    elem = key[1].cast<int>();

    set(h, elem, value);
}

double Vector::pyget(int key)
{
    int h = 0, elem = key;
    return get(h, elem);
}

void Vector::pyset(int key, double value)
{
    int h = 0, elem = key;
    set(h, elem, value);
}

void Vector::print(std::string out, const char *extra) const
{
    std::shared_ptr <psi::PsiOutStream> printer = (out == "outfile" ? outfile :
                                                   std::shared_ptr<OutFile>(new OutFile(out)));
    int h;
    if (extra == NULL) {
        printer->Printf("\n # %s #\n", name_.c_str());
    } else {
        printer->Printf("\n # %s %s #\n", name_.c_str(), extra);
    }
    for (h = 0; h < nirrep_; ++h) {
        printer->Printf(" Irrep: %d\n", h + 1);
        for (int i = 0; i < dimpi_[h]; ++i)
            printer->Printf("   %4d: %10.7f\n", i + 1, vector_[h][i]);
        printer->Printf("\n");
    }
}

double *Vector::to_block_vector()
{
    size_t size = 0;
    for (int h = 0; h < nirrep_; ++h)
        size += dimpi_[h];

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
        C_DGEMV(trans, A->rowspi_[h], A->colspi_[h], alpha, &(A->matrix_[h][0][0]), A->rowspi_[h],
                &(X->vector_[h][0]), 1, beta, &(vector_[h][0]), 1);
    }
}

double Vector::vector_dot(const std::shared_ptr<Vector> &other) { return vector_dot(*other.get()); }
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

void Vector::add(const std::shared_ptr<Vector> &other) { axpy(1.0, *other.get()); }
void Vector::add(const Vector &other) { axpy(1.0, other); }

void Vector::subtract(const std::shared_ptr<Vector> &other) { axpy(-1.0, *other.get()); }
void Vector::subtract(const Vector &other) { axpy(-1.0, other); }

void Vector::axpy(double scale, const SharedVector &other) { axpy(scale, *other.get()); }
void Vector::axpy(double scale, const Vector &other) {
    if (v_.size() != other.v_.size()) {
        throw PSIEXCEPTION("Vector::axpy: Vector sizes do not match!");
    }

    C_DAXPY(v_.size(), scale, const_cast<double *>(other.v_.data()), 1, v_.data(), 1);
}

void Vector::send() {}

void Vector::recv() {}

void Vector::bcast(int)
{
    // Assume the user allocated the matrix to the correct size first.
    /*std::cout<<"Someone is calling the vector bcast"<<std::endl;
    for (int h=0; h<nirrep_; ++h) {
        if (dimpi_[h] > 0)
            WorldComm->bcast(vector_[h], dimpi_[h], broadcaster);
    }*/
}

void Vector::sum()
{
    ///RMR-See note in Matrix::sum()
}

std::vector<py::buffer_info> Vector::array_interface(){
    std::vector<py::buffer_info> ret;

    if (numpy_shape_.size()) {
        if (nirrep_ > 1){
            throw PSIEXCEPTION("Vector::array_interface numpy shape with more than one irrep is not valid.");
        }

        std::vector<size_t> shape(numpy_shape_.size());
        std::vector<size_t> strides(numpy_shape_.size());
        size_t current_stride = sizeof(double);

        for (size_t i = numpy_shape_.size(); i-- > 0;) {
            shape[i] = numpy_shape_[i];
            strides[i] = current_stride;
            current_stride *= numpy_shape_[i];
        }
        ret.push_back(py::buffer_info(pointer(0), sizeof(double),
                                      py::format_descriptor<double>::format(),
                                      numpy_shape_.size(),
                                      shape, strides));

    } else {
        for (size_t h = 0; h < nirrep_; h++) {
            ret.push_back(py::buffer_info(
                pointer(h), sizeof(double),
                py::format_descriptor<double>::format(), 1,
                {static_cast<size_t>(dim(h))}, {sizeof(double)}));
        }
    }
    return ret;
}

} // namespace psi
