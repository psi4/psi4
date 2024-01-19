/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#include "psi4/libmints/benchmark.h"
#include "psi4/pybind11.h"

namespace py = pybind11;
using namespace pybind11::literals;

 void export_benchmarks(py::module& m) {
    m.def("benchmark_blas1", &psi::benchmark_blas1, "max_dim"_a, "min_time"_a,
          "Perform benchmark traverse of BLAS 1 routines. Use up to *max_dim* with each routine run at least *min_time* [s].");
    m.def("benchmark_blas2", &psi::benchmark_blas2, "max_dim"_a, "min_time"_a,
          "Perform benchmark traverse of BLAS 2 routines. Use up to *max_dim* with each routine run at least *min_time* [s].");
    m.def("benchmark_blas3", &psi::benchmark_blas3, "max_dim"_a, "min_time"_a, "nthread"_a = 1,
          "Perform benchmark traverse of BLAS 3 routines. Use up to *max_dim* with each routine run at least *min_time* [s] on *nthread*.");
    m.def("benchmark_disk", &psi::benchmark_disk, "max_dim"_a, "min_time"_a,
          "Perform benchmark of PSIO disk performance. Use up to *max_dim* with each routine run at least *min_time* [s].");
    m.def("benchmark_math", &psi::benchmark_math, "min_time"_a,
          "Perform benchmark of common double floating point operations including most of cmath. For each routine run at least *min_time* [s].");
    m.def("benchmark_integrals", &psi::benchmark_integrals, "max_am"_a, "min_time"_a,
          "Perform benchmark of psi integrals (of libmints type). Benchmark integrals called from different centers. For up to *max_am* with each shell combination run at least *min_time* [s].");
}
