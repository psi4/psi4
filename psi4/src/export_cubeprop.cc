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

#include "psi4/pybind11.h"

#include "psi4/libcubeprop/cubeprop.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"

using namespace psi;
namespace py = pybind11;
using namespace pybind11::literals;

void export_cubeprop(py::module& m) {
    py::class_<CubeProperties, std::shared_ptr<CubeProperties>>(m, "CubeProperties", "docstring")
        .def(py::init<std::shared_ptr<Wavefunction>>())
        .def("compute_density", &CubeProperties::compute_density, "Compute and dump a cube file for a density matrix",
             "D"_a, "key"_a)
        .def("compute_orbitals", &CubeProperties::compute_orbitals,
             "Compute and dump a cube file for a set of orbitals", "C"_a, "indices"_a, "labels"_a, "key"_a)
        .def("basisset", &CubeProperties::basisset, "Returns orbital/primary basis set associated with cubeprop.")
        .def("raw_compute_properties", &CubeProperties::raw_compute_properties,
             "Compute all relevant properties from options object specifications");
}
