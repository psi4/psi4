/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#ifdef USING_ecpint_RUNTIME
#include "psi4/libmints/ecpint_loader.h"
#endif

using namespace psi;
namespace py = pybind11;
using namespace pybind11::literals;

void export_misc(py::module &m) {
    m.def("timer_on", timer_on, "label"_a, "Start timer with *label*. Needs to be paired with :func:`psi4.core.timer_off`.");
    m.def("timer_off", timer_off, "label"_a, "Stop timer with *label*.");
    m.def("tstart", tstart, "Start module-level timer. Only one active at once.");
    m.def("tstop", tstop, "Stop module-level timer. Prints user, system, and total times to outfile.");
    m.def("clean_timers", clean_timers, "Reinitialize timers for independent ``timer.dat`` entries. Vital when earlier independent calc finished improperly.");
    m.def("get_timer_records", [](bool compact) {
        auto records = get_timer_records();

        py::dict out;

        for (const auto& rec : records) {
            py::dict d;

            d["wall_time"] = rec.wall_time;
            d["user_time"] = rec.user_time;
            d["system_time"] = rec.system_time;
            d["n_calls"] = rec.n_calls;

            if (!compact) {
                d["timer_name"] = rec.timer_name;

                if (rec.parent_id.empty()) {
                    d["parent_id"] = py::none();
                } else {
                    d["parent_id"] = rec.parent_id;
                }

                d["timer_path"] = rec.timer_path;
                d["level"] = rec.level;
            }

            out[py::str(rec.timer_id)] = d;
        }
        return out;
    },
    py::arg("compact") = false,
    "Get timing information as structured timer records. Returns a dictionary keyed by timer_id. "
    "Each value contains timing fields (wall_time, user_time, system_time, n_calls) and, when "
    "compact=False, hierarchy fields (parent_id, timer_name, timer_path, level).");

    m.def("ecpint_available", []() {
#ifdef USING_ecpint_RUNTIME
        // Runtime-optional build: libecpint is not linked into core.so, so the only way to
        // know whether ECP integrals can actually be computed is to attempt the dlopen and
        // see if it succeeds. This is the same check used internally by
        // IntegralFactory::ao_ecp() before touching any libecpint symbol.
        return ecpint_runtime::is_available();
#elif defined(USING_ecpint_STATIC)
        // Statically linked: always available if this build has ECP support at all.
        return true;
#else
        // Not compiled with ECP support at all.
        return false;
#endif
    }, "Returns whether ECP integrals (libecpint) are actually usable right now: true if compiled "
       "in statically, or (for runtime-optional builds) if libecpint could be dlopen'd. False if "
       "ECP support wasn't compiled in, or (runtime-optional builds) libecpint isn't installed.");
}
