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

#include "vector.h"

#include <algorithm>
#include <numeric>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsio/psio.hpp"

#include "dimension.h"
#include "matrix.h"
#include "psi4/occ/arrays.h"

namespace psi {

std::unique_ptr<Vector> Vector::clone() const {
    auto temp = std::make_unique<Vector>(dimpi_);
    temp->copy(*this);
    return temp;
}

void Vector::sort(const IntVector& idxs) {
    auto orig = clone();
    if (dimpi_ != idxs.dimpi()) {
        throw PSIEXCEPTION("Vector::sort Indexing vector and vector to sort must have the same dimension.");
    }
    // WARNING! Function also requires each irrep to be a permutation of 0, 1, 2...
    for (int h = 0; h < nirrep(); h++) {
        for (int i = 0; i < dimpi_[h]; i++) {
            set(h, i, orig->get(h, idxs.get(h, i)));
        }
    }
}

SharedVector Vector::get_block(const Slice &slice) {
    // check if slice is within bounds
    for (int h = 0; h < nirrep(); h++) {
        if (slice.end()[h] > dimpi_[h]) {
            std::string msg =
                "Invalid call to Vector::get_block(): Slice is out of bounds. Irrep = " + std::to_string(h);
            throw PSIEXCEPTION(msg);
        }
    }
    const Dimension &slice_begin = slice.begin();
    Dimension slice_dim = slice.end() - slice.begin();
    auto block = std::make_shared<Vector>("Block", slice_dim);
    for (int h = 0; h < nirrep(); h++) {
        int max_p = slice_dim[h];
        for (int p = 0; p < max_p; p++) {
            double value = get(h, p + slice_begin[h]);
            block->set(h, p, value);
        }
    }
    return block;
}

void Vector::set_block(const Slice &slice, const Vector& block) {
    // check if slice is within bounds
    for (int h = 0; h < nirrep(); h++) {
        if (slice.end()[h] > dimpi_[h]) {
            std::string msg =
                "Invalid call to Vector::set_block(): Slice is out of bounds. Irrep = " + std::to_string(h);
            throw PSIEXCEPTION(msg);
        }
    }
    const Dimension &slice_begin = slice.begin();
    Dimension slice_dim = slice.end() - slice.begin();
    for (int h = 0; h < nirrep(); h++) {
        int max_p = slice_dim[h];
        for (int p = 0; p < max_p; p++) {
            double value = block.get(h, p);
            set(h, p + slice_begin[h], value);
        }
    }
}

void Vector::zero() { std::fill(v_.begin(), v_.end(), 0.0); }

void Vector::print(std::string out) const {
    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));
    printer->Printf("\n # %s #\n", name_.c_str());
    for (int h = 0; h < nirrep(); ++h) {
        printer->Printf(" Irrep: %d\n", h + 1);
        for (int i = 0; i < dimpi_[h]; ++i) printer->Printf("   %4d: %20.15f\n", i + 1, vector_[h][i]);
        printer->Printf("\n");
    }
}

void Vector::gemv(bool transa, double alpha, const Matrix& A, const Vector& X, double beta) {
    char trans = transa ? 't' : 'n';

    for (int h = 0; h < nirrep(); ++h) {
        C_DGEMV(trans, A.rowspi(h), A.colspi(h), alpha, A.get_pointer(h), A.colspi(h), &(X.vector_[h][0]), 1, beta,
                &(vector_[h][0]), 1);
    }
}

double Vector::vector_dot(const Vector &other) {
    if (v_.size() != other.v_.size()) {
        throw PSIEXCEPTION("Vector::vector_dot: Vector sizes do not match!");
    }

    return C_DDOT(v_.size(), v_.data(), 1, const_cast<double *>(other.v_.data()), 1);
}

double Vector::sum_of_squares() { return C_DDOT(v_.size(), v_.data(), 1, v_.data(), 1); }

double Vector::norm() { return sqrt(sum_of_squares()); }

double Vector::rms() { return sqrt(sum_of_squares() / v_.size()); }

void Vector::scale(double sc) { C_DSCAL(v_.size(), sc, v_.data(), 1); }

void Vector::add(const Vector &other) { axpy(1.0, other); }

void Vector::subtract(const Vector &other) { axpy(-1.0, other); }

void Vector::axpy(double scale, const Vector &other) {
    if (v_.size() != other.v_.size()) {
        throw PSIEXCEPTION("Vector::axpy: Vector sizes do not match!");
    }

    C_DAXPY(v_.size(), scale, const_cast<double *>(other.v_.data()), 1, v_.data(), 1);
}

void Vector::save(psi::PSIO *const psio, size_t fileno) {
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno)) {
        already_open = true;
    } else {
        psio->open(fileno, PSIO_OPEN_OLD);
    }

    for (int h = 0; h < nirrep(); ++h) {
        std::string str(name_);
        str += " Irrep " + std::to_string(h);

        // Write the sub-blocks
        if (dimpi_[h] > 0)
            psio->write_entry(fileno, const_cast<char *>(str.c_str()), (char *)vector_[h],
                              sizeof(double) * dimpi_[h]);
    }
    if (!already_open) psio->close(fileno, 1);  // Close and keep
}

void Vector::load(psi::PSIO *const psio, size_t fileno) {
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno)) {
        already_open = true;
    } else {
        psio->open(fileno, PSIO_OPEN_OLD);
    }

    for (int h = 0; h < nirrep(); ++h) {
        std::string str(name_);
        str += " Irrep " + std::to_string(h);

        // Read the sub-blocks
        if (dimpi_[h] > 0)
            psio->read_entry(fileno, str.c_str(), (char *)vector_[h],
                             sizeof(double) * dimpi_[h]);
    }

    if (!already_open) psio->close(fileno, 1);  // Close and keep
}
}  // namespace psi
