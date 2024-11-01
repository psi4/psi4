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

#ifndef _PSI_SRC_LIB_LIBDIIS_DIISMANAGER_H_
#define _PSI_SRC_LIB_LIBDIIS_DIISMANAGER_H_

#include <vector>
#include <map>

#include "psi4/pragma.h"
#include "psi4/libmints/matrix.h"

#include "psi4/pybind11.h"

namespace psi {

class PSIO;

/**
   @brief The DIISManager class handles DIIS extrapolations.
 */

class PSI_API DIISManager {
   public:

    /**
     * @brief How the quantities are to be stored;
     *
     * OnDisk - Stored on disk, and retrieved when required
     * InCore - Stored in memory throughout
     */
    enum class StoragePolicy { InCore, OnDisk };
    /**
     * @brief How vectors are removed from the subspace, when required
     *
     * LargestError - The vector corresponding to the largest error is removed
     * OldestFirst - A first-in-first-out policy is used
     */
    enum class RemovalPolicy { LargestError, OldestAdded };

    DIISManager(int maxSubspaceSize, const std::string& label, RemovalPolicy = RemovalPolicy::LargestError, StoragePolicy = StoragePolicy::OnDisk);
    DIISManager() {}
    ~DIISManager();

    // Variadic templates to interface with Python.
    // If you're new to variadics, these allow multiple arguments.
    // MUST be implemented in header.
    template <typename ... types>
    void set_error_vector_size(types ... arrays) {
        pydiis.attr("set_error_vector_size")(arrays...);
    }
    template <typename ... types>
    void set_vector_size(types... arrays) {
        pydiis.attr("set_vector_size")(arrays...);
    }
    template <typename... types>
    bool extrapolate(types... arrays) {
        return py::len(pydiis.attr("extrapolate")(arrays...));
    }
    template <typename ... types>
    bool add_entry(types... arrays) {
        auto success = pydiis.attr("add_entry")(arrays...);
        return success.template cast<bool>();
    }

    void delete_diis_file();

    void reset_subspace();
    /// The number of vectors currently in the subspace
    int subspace_size();

  protected:

    py::object pydiis;

};

}  // namespace psi

#endif  // Header guard
