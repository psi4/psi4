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

#include "psi4/pybind11.h"

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpe/psipe.h"

#ifdef USING_cppe
using namespace psi;

void export_cppe(py::module& m) {
    py::class_<libcppe::BorderOptions, std::shared_ptr<libcppe::BorderOptions>> pe_border_options(
        m, "PeBorderOptions", "Border Options for PE library");
    py::enum_<libcppe::BorderType>(pe_border_options, "BorderType")
        .value("rem", libcppe::BorderType::rem)
        .value("redist", libcppe::BorderType::redist);
    pe_border_options.def(py::init<>())
        .def_readwrite("border_type", &libcppe::BorderOptions::border_type)
        .def_readwrite("rmin", &libcppe::BorderOptions::rmin)
        .def_readwrite("nredist", &libcppe::BorderOptions::nredist)
        .def_readwrite("redist_order", &libcppe::BorderOptions::redist_order)
        .def_readwrite("redist_pol", &libcppe::BorderOptions::redist_pol);

    py::class_<libcppe::PeOptions, std::shared_ptr<libcppe::PeOptions>> pe_options(m, "PeOptions",
                                                                                   "Options for PE library");
    pe_options.def(py::init<>())
        .def_readwrite("potfile", &libcppe::PeOptions::potfile)
        .def_readwrite("iso_pol", &libcppe::PeOptions::iso_pol)

        .def_readwrite("induced_thresh", &libcppe::PeOptions::induced_thresh)
        .def_readwrite("do_diis", &libcppe::PeOptions::do_diis)
        .def_readwrite("maxiter", &libcppe::PeOptions::maxiter)

        .def_readwrite("pe_border", &libcppe::PeOptions::pe_border)
        .def_readwrite("border_options", &libcppe::PeOptions::border_options);

    py::class_<PeState, std::shared_ptr<PeState>> pe(m, "PE", "Class interfacing with CPPE");
    py::enum_<PeState::CalcType>(pe, "CalcType")
        .value("total", PeState::CalcType::total)
        .value("electronic_only", PeState::CalcType::electronic_only);

    pe.def(py::init<libcppe::PeOptions, std::shared_ptr<BasisSet>>())
        .def("compute_pe_contribution", &PeState::compute_pe_contribution,
             "Compute PE contributions to energy and Fock matrix", py::arg("D"), py::arg("type"))
        .def("print_energy_summary", &PeState::print_energy_summary);
}
#endif
