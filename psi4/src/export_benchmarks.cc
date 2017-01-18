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

#include "psi4/libmints/benchmark.h"
#include "psi4/pybind11.h"

void export_benchmarks(py::module& m)
{
    m.def("benchmark_blas1",     &psi::benchmark_blas1, "docstring");
    m.def("benchmark_blas2",     &psi::benchmark_blas2, "docstring");
    m.def("benchmark_blas3",     &psi::benchmark_blas3, "docstring");
    m.def("benchmark_disk",      &psi::benchmark_disk, "docstring");
    m.def("benchmark_math",      &psi::benchmark_math, "docstring");
    m.def("benchmark_integrals", &psi::benchmark_integrals, "docstring");
}
