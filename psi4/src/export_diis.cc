/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
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

#include "psi4/libdiis/diisentry.h"
#include "psi4/libdiis/diismanager.h"

using namespace psi;
namespace py = pybind11;
using namespace pybind11::literals;

void export_diis(py::module &m) {
    py::class_<DIISManager, std::shared_ptr<DIISManager> > diis(m, "DIISManager", "docstring");

    py::enum_<DIISManager::StoragePolicy>(diis, "StoragePolicy")
        .value("InCore", DIISManager::StoragePolicy::InCore)
        .value("OnDisk", DIISManager::StoragePolicy::OnDisk);

    py::enum_<DIISManager::RemovalPolicy>(diis, "RemovalPolicy")
        .value("LargestError", DIISManager::RemovalPolicy::LargestError)
        .value("OldestAdded", DIISManager::RemovalPolicy::OldestAdded);

    diis.def(py::init<>())
        .def(py::init<int, const std::string&, DIISManager::RemovalPolicy, DIISManager::StoragePolicy>(),
                "max_vectors"_a, "name"_a, "removal_policy"_a = DIISManager::RemovalPolicy::LargestError, "storage_policy"_a = DIISManager::StoragePolicy::InCore)
        .def("set_error_vector_size", [](DIISManager& diis, const SharedMatrix mat) {
                diis.set_error_vector_size(1, DIISEntry::InputType::Matrix, mat.get());
            })
        .def("set_vector_size", [](DIISManager& diis, const SharedMatrix mat) {
                diis.set_vector_size(1, DIISEntry::InputType::Matrix, mat.get());
            })
        .def("add_entry", [](DIISManager& diis, const SharedMatrix m1, const SharedMatrix m2) {
                diis.add_entry(2, m1.get(), m2.get());
            })
        .def("set_error_vector_size", [](DIISManager& diis, const SharedMatrix m1, const SharedMatrix m2) {
                diis.set_error_vector_size(2, DIISEntry::InputType::Matrix, m1.get(),  DIISEntry::InputType::Matrix, m2.get());
            })
        .def("set_vector_size", [](DIISManager& diis, const SharedMatrix m1, const SharedMatrix m2) {
                diis.set_vector_size(2, DIISEntry::InputType::Matrix, m1.get(),  DIISEntry::InputType::Matrix, m2.get());
             })
        .def("add_entry", [](DIISManager& diis, const SharedMatrix m1, const SharedMatrix m2, const SharedMatrix m3, const SharedMatrix m4) {
                diis.add_entry(4, m1.get(), m2.get(), m3.get(), m4.get());
            })
        .def("reset_subspace", &DIISManager::reset_subspace, "docstring")
        .def("delete_diis_file", &DIISManager::delete_diis_file, "docstring");

    py::class_<DIISEntry, std::shared_ptr<DIISEntry> > diis_entry(m, "DIISEntry", "docstring");
}
