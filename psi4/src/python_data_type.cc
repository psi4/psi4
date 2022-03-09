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

#include "python_data_type.h"

namespace psi {

PythonDataType::PythonDataType() {}

PythonDataType::PythonDataType(const py::object &p) : python_object_(p) {}

PythonDataType::~PythonDataType() {}

std::string PythonDataType::type() const { return std::string("python"); }

void PythonDataType::assign(const py::object &p) {
    python_object_ = p;
    changed();
}

const py::object &PythonDataType::to_python() const { return python_object_; }

}  // namespace psi
