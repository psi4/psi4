/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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

#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"

using namespace psi;
namespace py = pybind11;
using namespace pybind11::literals;

void export_misc(py::module &m) {
    m.def("timer_on", timer_on, "label"_a, "Start timer with *label*. Needs to be paired with :func:`psi4.core.timer_off`.");
    m.def("timer_off", timer_off, "label"_a, "Stop timer with *label*.");
    m.def("tstart", tstart, "Start module-level timer. Only one active at once.");
    m.def("tstop", tstop, "Stop module-level timer. Prints user, system, and total times to outfile.");
    m.def("collect_timers", collect_timers, "");
    m.def("clean_timers", clean_timers, "Reinitialize timers for independent ``timer.dat`` entries. Vital when earlier independent calc finished improperly.");
}
