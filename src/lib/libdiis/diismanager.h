#ifndef _PSI_SRC_LIB_LIBDIIS_DIISMANAGER_H_
#define _PSI_SRC_LIB_LIBDIIS_DIISMANAGER_H_

#include <libdpd/dpd.h>
#include <libpsio/psio.hpp>
#include <vector>
#include <stdarg.h>
#include <map>
#include "psifiles.h"
#include "diisentry.h"

using namespace psi;

namespace psi{ namespace libdiis{

  /**
     @Brief The DIISManager class handles DIIS extrapolations.
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
        
        DIISManager(int maxSubspaceSize, RemovalPolicy = LargestError,
                    StoragePolicy = InCore);
        DIISManager() {_maxSubspaceSize = 0;}
        ~DIISManager();

        void set_error_vector_size(int numQuantities, ...);
        void set_vector_size(int numQuantities, ...);
        void add_entry(int numQuatities, ...);
        void extrapolate(int numQuatities, ...);
        int remove_entry();
        /// The number of vectors currently in the subspace
        int subspace_size() {return _subspace.size();}
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
};

}} // End namespaces

// This is here so that files including this have clean(er) syntax
using namespace psi::libdiis;

#endif // Header guard
