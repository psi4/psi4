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

#include "psi4/libpsio/psio.hpp"

using namespace psi;

void export_psio(py::module &m) {
    py::class_<PSIO, std::shared_ptr<PSIO> >(m, "IO", "docstring")
        .def("state", &PSIO::state, "Return 1 if PSIO library is activated")
        .def("open", &PSIO::open, "Open unit. Status can be PSIO_OPEN_OLD (if existing file is to be opened) or PSIO_OPEN_NEW if new file should be open", py::arg("unit"), py::arg("status"))
        .def("close", &PSIO::close, "Close unit. If keep == 0, will remove the file, else keep it", py::arg("unit"), py::arg("keep"))
        .def("rehash", &PSIO::rehash, "Sync up the object to the file on disk by closing and opening the file, if necessary", py::arg("unit"))
        .def("open_check", &PSIO::open_check, "Return 1 if unit is open", py::arg("unit"))
        .def("tocclean", &PSIO::tocclean, "Delete all TOC entries after the given key. If a blank key is given, the entire TOC will be wiped", py::arg("unit"), py::arg("key"))
        .def("tocprint", &PSIO::tocprint, "Print the table of contents for the given unit")
        .def("tocentry_exists", &PSIO::tocentry_exists, "Checks the TOC to see if a particular keyword exists there or not")
        .def("tocwrite", &PSIO::tocwrite, "Write the table of contents for passed file number")
        .def("getpid", &PSIO::getpid, "Lookup process id")
        .def("set_pid", &PSIO::set_pid, "Set process id", py::arg("pid"))
        .def_static("shared_object", &PSIO::shared_object, "Return the global shared object")
        .def_static("get_default_namespace", &PSIO::get_default_namespace, "Get the default namespace (for PREFIX.NAMESPACE.UNIT file numbering)")
        .def_static("set_default_namespace", &PSIO::set_default_namespace, "Set the current namespace (for PREFIX.NAMESPACE.UNIT file numbering)", py::arg("ns"))
        .def_static("change_file_namespace", &PSIO::change_file_namespace, "Change file number from ns1 to ns2",  py::arg("fileno"), py::arg("ns1"), py::arg("ns2"));

    py::class_<PSIOManager, std::shared_ptr<PSIOManager> >(m, "IOManager", "PSIOManager is a class designed to be used as a static object to track all PSIO operations in a given PSI4 computation")
        .def_static("shared_object", &PSIOManager::shared_object, "The one and (should be) only instance of PSIOManager for a PSI4 instance")
        .def("print_out", &PSIOManager::print_out, "Print the current status of PSI4 files")
        .def("psiclean", &PSIOManager::psiclean, "Execute the psiclean protocol, deleting all recorded files, except those currently marked for retention")
        .def("crashclean", &PSIOManager::crashclean, "Clean from disk-mirrored image after crash. NOT to be called during regular computation.")
        .def("mark_file_for_retention", &PSIOManager::mark_file_for_retention, "Mark a file to be retained after a psiclean operation, ie for use in a later computation", py::arg("full_path"), py::arg("retain"))
        .def("write_scratch_file", &PSIOManager::write_scratch_file, "Write a string to a temporary file.  The scratch file is opened and closed by this function.", py::arg("full_path"), py::arg("text"))
        .def("set_default_path", &PSIOManager::set_default_path, "Set the default path for files to be stored", py::arg("path"))
        .def("set_specific_path", &PSIOManager::set_specific_path, "Set the path for specific file numbers", py::arg("fileno"), py::arg("path"))
        .def("get_file_path", &PSIOManager::get_file_path, "Get the path for a specific file number", py::arg("fileno"))
        .def("set_specific_retention", &PSIOManager::set_specific_retention, "Set the specific file number to be retained", py::arg("fileno"), py::arg("retain"))
        .def("get_default_path", &PSIOManager::get_default_path, "Return the default path");
}
