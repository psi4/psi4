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

#include "psi4/pybind11.h"

#include "psi4/libdiis/diisentry.h"
#include "psi4/libdiis/diismanager.h"

using namespace psi;
namespace py = pybind11;

void export_diis(py::module &m) {
    py::class_<DIISManager, std::shared_ptr<DIISManager> >(m, "DIISManager", "docstring")
        .def(py::init<>())
        .def("reset_subspace", &DIISManager::reset_subspace, "docstring")
        .def("delete_diis_file", &DIISManager::delete_diis_file, "docstring");
}
