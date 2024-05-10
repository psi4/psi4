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

#include "diismanager.h"

using namespace psi;

namespace psi {

/**
 *
 * @param maxSubspaceSize Maximum number of vectors allowed in the subspace
 * @param label: the base part of the label used to store the vectors to disk
 * @param removalPolicy: How to decide which vectors to remove when the subspace is full
 * @param storagePolicy: How to store the DIIS vectors
 */
DIISManager::DIISManager(int maxSubspaceSize, const std::string &label, RemovalPolicy removalPolicy,
                         StoragePolicy storagePolicy) {
    auto diis_file = py::module_::import("psi4").attr("driver").attr("diis");
    py::object pyRemovalPolicy, pyStoragePolicy;
    if (removalPolicy == RemovalPolicy::LargestError) {
        pyRemovalPolicy = diis_file.attr("RemovalPolicy").attr("LargestError");
    } else {
        pyRemovalPolicy = diis_file.attr("RemovalPolicy").attr("OldestAdded");
    }
    if (storagePolicy == StoragePolicy::InCore) {
        pyStoragePolicy = diis_file.attr("StoragePolicy").attr("InCore");
    } else {
        pyStoragePolicy = diis_file.attr("StoragePolicy").attr("OnDisk");
    }
    pydiis = diis_file.attr("DIIS")(maxSubspaceSize, label, pyRemovalPolicy, pyStoragePolicy);
}

int DIISManager::subspace_size() { return py::len(pydiis.attr("stored_vectors")); }

void DIISManager::delete_diis_file() {
    pydiis.attr("delete_diis_file")();
}

void DIISManager::reset_subspace() {
    pydiis.attr("reset_subspace")();
}

DIISManager::~DIISManager() {
}

}  // namespace psi
