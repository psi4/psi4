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

#include "psi4/pybind11.h"

#include "psi4/libpsio/psio.hpp"

using namespace psi;

void export_psio(py::module &m)
{
    py::class_<PSIO, std::shared_ptr<PSIO> >( m, "IO", "docstring" ).
        def( "state", &PSIO::state, "docstring" ).
        def( "open", &PSIO::open, "docstring" ).
        def( "close", &PSIO::close, "docstring" ).
        def( "rehash", &PSIO::rehash, "docstring" ).
        def( "open_check", &PSIO::open_check, "docstring" ).
        def( "tocclean", &PSIO::tocclean, "docstring" ).
        def( "tocprint", &PSIO::tocprint, "docstring" ).
        def( "tocwrite", &PSIO::tocwrite, "docstring" ).
        def( "set_pid", &PSIO::set_pid, "docstring" ).
        def_static("shared_object", &PSIO::shared_object, "docstring").
        def_static("get_default_namespace", &PSIO::get_default_namespace, "docstring").
        def_static("set_default_namespace", &PSIO::set_default_namespace,
            py::arg("ns"), "docstring").
        def_static("change_file_namespace", &PSIO::change_file_namespace,
            py::arg("fileno"), py::arg("ns1"), py::arg("ns2"), "docstring");

    py::class_<PSIOManager, std::shared_ptr<PSIOManager> >( m, "IOManager", "docstring" ).
        def_static("shared_object", &PSIOManager::shared_object, "docstring").
        def( "print_out", &PSIOManager::print_out, "docstring" ).
        def( "psiclean", &PSIOManager::psiclean, "docstring" ).
        def( "crashclean", &PSIOManager::crashclean, "docstring" ).
        def( "mark_file_for_retention", &PSIOManager::mark_file_for_retention, "docstring" ).
        def( "write_scratch_file", &PSIOManager::write_scratch_file, "docstring").
        def( "set_default_path", &PSIOManager::set_default_path, "docstring" ).
        def( "set_specific_path", &PSIOManager::set_specific_path, "docstring" ).
        def( "get_file_path", &PSIOManager::get_file_path, "docstring" ).
        def( "set_specific_retention", &PSIOManager::set_specific_retention, "docstring" ).
        def( "get_default_path", &PSIOManager::get_default_path, "docstring" );
}
