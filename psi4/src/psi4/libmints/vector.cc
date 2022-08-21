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

#include "psi4/libpsio/psio.hpp"

#include "dimension.h"
#include "matrix.h"
#include "psi4/occ/arrays.h"

namespace psi {

void Vector::gemv(bool transa, double alpha, const Matrix& A, const Vector& X, double beta) {
    char trans = transa ? 't' : 'n';

    for (int h = 0; h < nirrep(); ++h) {
        C_DGEMV(trans, A.rowspi(h), A.colspi(h), alpha, A.get_pointer(h), A.colspi(h), &(X.vector_[h][0]), 1, beta,
                &(vector_[h][0]), 1);
    }
}

double Vector::vector_dot(const Vector &other) const {
    if (v_.size() != other.v_.size()) {
        throw PSIEXCEPTION("Vector::vector_dot: Vector sizes do not match!");
    }

    return C_DDOT(v_.size(), v_.data(), 1, const_cast<double *>(other.v_.data()), 1);
}

double Vector::sum_of_squares() const { return C_DDOT(v_.size(), v_.data(), 1, v_.data(), 1); }

double Vector::norm() const { return sqrt(sum_of_squares()); }

double Vector::rms() const { return sqrt(sum_of_squares() / v_.size()); }

void Vector::scale(double sc) { C_DSCAL(v_.size(), sc, v_.data(), 1); }

void Vector::add(const Vector &other) { axpy(1.0, other); }

void Vector::subtract(const Vector &other) { axpy(-1.0, other); }

void Vector::axpy(double scale, const Vector &other) {
    if (v_.size() != other.v_.size()) {
        throw PSIEXCEPTION("Vector::axpy: Vector sizes do not match!");
    }

    C_DAXPY(v_.size(), scale, const_cast<double *>(other.v_.data()), 1, v_.data(), 1);
}

void Vector::save(psi::PSIO *const psio, size_t fileno) const {
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
