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

#ifndef _PSI_SRC_LIB_LIBDIIS_DIISMANAGER_H_
#define _PSI_SRC_LIB_LIBDIIS_DIISMANAGER_H_

#include <vector>
#include <map>

#include "psi4/pragma.h"
#include "psi4/libdiis/diisentry.h"
#include "psi4/libmints/matrix.h"

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
    DIISManager() { _maxSubspaceSize = 0; }
    ~DIISManager();

    // C-style variadic allows you to DIIS "direct sums" of quantities.
    // For instance, orbital-optimized theories can do a combined DIIS on orbital amplitudes and T2
    void set_error_vector_size(int numQuantities, ...);
    void set_vector_size(int numQuantities, ...);
    bool extrapolate(int numQuatities, ...);
    bool add_entry(int numQuatities, ...);

    // Wrappers for those who dislike variadic
    void set_error_vector_size(SharedMatrix error) { DIISManager::set_error_vector_size(1, error.get()); }

    void set_vector_size(SharedMatrix state) { DIISManager::set_vector_size(1, state.get()); }

    bool add_entry(SharedMatrix state, SharedMatrix error) {
        return DIISManager::add_entry(2, state.get(), error.get());
    }

    bool extrapolate(SharedMatrix extrapolated) { return DIISManager::extrapolate(1, extrapolated.get()); }

    int remove_entry();
    void reset_subspace();
    void delete_diis_file();
    /// The number of vectors currently in the subspace
    int subspace_size();

   protected:
    int get_next_entry_id();

    /// How the vectors are handled in memory
    StoragePolicy _storagePolicy;
    /// How vectors are removed from the subspace
    RemovalPolicy _removalPolicy;
    /// The maximum number of vectors allowed in the subspace
    int _maxSubspaceSize;
    /// The size of the error vector
    int _errorVectorSize;
    /// The size of the vector
    int _vectorSize;
    /// The number of components in the error vector
    int _numErrorVectorComponents;
    /// The number of components in the vector
    int _numVectorComponents;
    /// The counter that keeps track of how many entries have been added
    int _entryCount;
    /// The DIIS entries
    std::vector<DIISEntry> _subspace;
    /// The types used in building the vector and the error vector
    std::vector<DIISEntry::InputType> _componentTypes;
    /// The types used in the vector
    std::vector<size_t> _componentSizes;
    /// The label used in disk storage of the DIISEntry objects
    std::string _label;
    /// The PSIO object to use for I/O
    std::shared_ptr<PSIO> _psio;
};

}  // namespace psi

#endif  // Header guard
