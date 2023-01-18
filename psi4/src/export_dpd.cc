/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

#include "psi4/libdpd/dpd.h"
#include "psi4/libmints/matrix.h"

using namespace psi;
namespace py = pybind11;
using namespace pybind11::literals;

void export_dpd(py::module &m) {
    py::class_<dpdbuf4, std::shared_ptr<dpdbuf4>>(m, "dpdbuf4", "docstring")
        .def("axpy_matrix", &dpdbuf4::axpy_matrix, "Add 'a' times a Matrix to this.")
        .def("zero", &dpdbuf4::zero, "Fill all with entries.")
        .def("rowdim", [](dpdbuf4& buf) {
                std::vector<int> dim;
                for (int h = 0; h < buf.params->nirreps; ++h) {
                    dim.push_back(buf.params->rowtot[h]);
                }
                return Dimension(dim);
            }, "Return the dimensions of the row index.")
        .def("coldim", [](dpdbuf4& buf) {
                std::vector<int> dim;
                for (int h = 0; h < buf.params->nirreps; ++h) {
                    dim.push_back(buf.params->coltot[h]);
                }
                return Dimension(dim);
            }, "Return the dimensions of the column index.");

    py::class_<dpdfile2, std::shared_ptr<dpdfile2>>(m, "dpdfile2", "docstring")
        .def("axpy_matrix", &dpdfile2::axpy_matrix, "Add 'a' times a Matrix to this.")
        .def("zero", &dpdfile2::zero, "Fill all entries with zeroes.")
        .def("rowdim", [](dpdfile2& file) {
                std::vector<int> dim;
                for (int h = 0; h < file.params->nirreps; ++h) {
                    dim.push_back(file.params->rowtot[h]);
                }
                return Dimension(dim);
            }, "Return the dimensions of the row index.")
        .def("coldim", [](dpdfile2& file) {
                std::vector<int> dim;
                for (int h = 0; h < file.params->nirreps; ++h) {
                    dim.push_back(file.params->coltot[h]);
                }
                return Dimension(dim);
            }, "Return the dimensions of the column index.");
}
