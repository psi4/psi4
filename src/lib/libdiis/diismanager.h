/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#ifndef _PSI_SRC_LIB_LIBDIIS_DIISMANAGER_H_
#define _PSI_SRC_LIB_LIBDIIS_DIISMANAGER_H_

#include "diisentry.h"
#include <vector>
#include <map>

namespace boost {
template<class T>
class shared_ptr;
}

namespace psi{

class DIISEntry;
class PSIO;

  /**
     @brief The DIISManager class handles DIIS extrapolations.
   */

class DIISManager{
    public:
        /**
         * @brief How the quantities are to be stored;
         *
         * OnDisk - Stored on disk, and retrieved when required
         * InCore - Stored in memory throughout
         */
        enum StoragePolicy {InCore, OnDisk};
        /**
         * @brief How vectors are removed from the subspace, when required
         *
         * LargestError - The vector corresponding to the largest error is removed
         * OldestFirst - A first-in-first-out policy is used
         */
        enum RemovalPolicy {LargestError, OldestAdded};

        DIISManager(int maxSubspaceSize, const std::string& label,
                    RemovalPolicy = LargestError,
                    StoragePolicy = OnDisk);
        DIISManager() {_maxSubspaceSize = 0;}
        ~DIISManager();

        void set_error_vector_size(int numQuantities, ...);
        void set_vector_size(int numQuantities, ...);
        bool extrapolate(int numQuatities, ...);
        bool add_entry(int numQuatities, ...);
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
        std::vector<DIISEntry*> _subspace;
        /// The types used in building the vector and the error vector
        std::vector<DIISEntry::InputType> _componentTypes;
        /// The types used in the vector
        std::vector<size_t> _componentSizes;
        /// The label used in disk storage of the DIISEntry objects
        std::string _label;
        /// The PSIO object to use for I/O
        boost::shared_ptr<PSIO> _psio;
};

} // End namespace

#endif // Header guard
