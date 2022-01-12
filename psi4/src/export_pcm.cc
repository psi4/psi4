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

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsipcm/psipcm.h"

using namespace psi;
namespace py = pybind11;
using namespace pybind11::literals;

#ifdef USING_PCMSolver

void export_pcm(py::module& m) {
    py::class_<PCM, std::shared_ptr<PCM>> pcm(m, "PCM", "Class interfacing with PCMSolver");

    py::enum_<PCM::CalcType>(pcm, "CalcType")
        .value("Total", PCM::CalcType::Total)
        .value("NucAndEle", PCM::CalcType::NucAndEle)
        .value("EleOnly", PCM::CalcType::EleOnly);

    pcm.def(py::init<std::string, int, std::shared_ptr<BasisSet>>())
        .def("compute_PCM_terms", &PCM::compute_PCM_terms, "Compute PCM contributions to energy and Fock matrix", "D"_a,
             "type"_a)
        .def("compute_V", &PCM::compute_V, "Computes electronic PCM contributions due to first-order perturbed densities", "D"_a);
}
#endif
