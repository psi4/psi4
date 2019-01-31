/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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
#include "psi4/pragma.h"
PRAGMA_WARNING_PUSH
PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
#include <memory>
PRAGMA_WARNING_POP
#include "psi4/libpsio/psio.hpp"
#include "diisentry.h"
#include <cmath>
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include <sstream>

namespace psi {

DIISEntry::DIISEntry(std::string label, int ID, int orderAdded, int errorVectorSize, double *errorVector,
                     int vectorSize, double *vector, std::shared_ptr<PSIO> psio)
    : _vectorSize(vectorSize),
      _errorVectorSize(errorVectorSize),
      _vector(vector),
      _errorVector(errorVector),
      _ID(ID),
      _orderAdded(orderAdded),
      _label(label),
      _psio(psio) {
    double sumSQ = C_DDOT(_errorVectorSize, _errorVector, 1, _errorVector, 1);
    _rmsError = sqrt(sumSQ / _errorVectorSize);
    _dotProducts[_ID] = sumSQ;
    _knownDotProducts[_ID] = true;
    std::stringstream s;
    s << _label << ":entry " << ID;
    _label = s.str();
}

void DIISEntry::open_psi_file() {
    if (_psio->open_check(PSIF_LIBDIIS) == 0) {
        _psio->open(PSIF_LIBDIIS, PSIO_OPEN_OLD);
    }
}

void DIISEntry::close_psi_file() {
    if (_psio->open_check(PSIF_LIBDIIS) == 1) {
        _psio->close(PSIF_LIBDIIS, 1);
    }
}

void DIISEntry::dump_vector_to_disk() {
    std::string label = _label + " vector";
    open_psi_file();
    _psio->write_entry(PSIF_LIBDIIS, label.c_str(), (char *)_vector, _vectorSize * sizeof(double));
    free_vector_memory();
}

void DIISEntry::read_vector_from_disk() {
    if (_vector == nullptr) {
        _vector = new double[_vectorSize];
        std::string label = _label + " vector";
        open_psi_file();
        _psio->read_entry(PSIF_LIBDIIS, label.c_str(), (char *)_vector, _vectorSize * sizeof(double));
    }
}

void DIISEntry::dump_error_vector_to_disk() {
    std::string label = _label + " error";
    open_psi_file();
    _psio->write_entry(PSIF_LIBDIIS, label.c_str(), (char *)_errorVector, _errorVectorSize * sizeof(double));
    free_error_vector_memory();
}

void DIISEntry::read_error_vector_from_disk() {
    if (_errorVector == nullptr) {
        _errorVector = new double[_errorVectorSize];
        std::string label = _label + " error";
        open_psi_file();
        _psio->read_entry(PSIF_LIBDIIS, label.c_str(), (char *)_errorVector, _errorVectorSize * sizeof(double));
    }
}

void DIISEntry::free_vector_memory() {
    if (_vector) delete[] _vector;
    _vector = nullptr;
}

void DIISEntry::free_error_vector_memory() {
    if (_errorVector) delete[] _errorVector;
    _errorVector = nullptr;
}

DIISEntry::~DIISEntry() {
    if (_vector != nullptr) delete[] _vector;
    if (_errorVector != nullptr) delete[] _errorVector;
}

}  // namespace psi
