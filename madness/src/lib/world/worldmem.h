#ifndef MADNESS_WORLD_WORLDMEM_H__INCLUDED
#define MADNESS_WORLD_WORLDMEM_H__INCLUDED

/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680


  $Id: world.cc 454 2008-01-26 03:08:15Z rjharrison $
*/

#include <madness_config.h>
#ifdef WORLD_GATHER_MEM_STATS
#include <new>
#endif // WORLD_GATHER_MEM_STATS
#include <cstddef>

namespace madness {

    /// Used to output memory statistics and control tracing, etc.

    /// There is currently only one instance of this class in worldmem.cc
    /// that is used by special versions of global new+delete enabled
    /// by compliling with WORLD_GATHER_MEM_STATS enabled.
    ///
    /// Default for max_mem_limit is unlimited.
    class WorldMemInfo {
#ifdef WORLD_GATHER_MEM_STATS
        friend void* ::operator new(size_t size) throw (std::bad_alloc);
        friend void ::operator delete(void *p) throw();
#endif
    private:
        /// Invoked when user pointer p is allocated with size bytes
        void do_new(void *p, std::size_t size);

        /// Invoked when user pointer p is deleted with size bytes
        void do_del(void *p, std::size_t size);

    public:
        static const unsigned long overhead = 4*sizeof(long) + 16;
        /// If you add new stats be sure that the initialization in worldmem.cc is OK
        unsigned long num_new_calls;   ///< Counts calls to new
        unsigned long num_del_calls;   ///< Counts calls to delete
        unsigned long cur_num_frags;   ///< Current number of allocated fragments
        unsigned long max_num_frags;   ///< Lifetime maximum number of allocated fragments
        unsigned long cur_num_bytes;   ///< Current amount of allocated memory in bytes
        unsigned long max_num_bytes;   ///< Lifetime maximum number of allocated bytes
        unsigned long max_mem_limit;   ///< if size+cur_num_bytes>max_mem_limit new will throw MadnessException
        bool trace;

        /// Prints memory use statistics to std::cout
        void print() const;

        /// Resets all counters to zero
        void reset();

        /// If trace is set true a message is printed for every new and delete
        void set_trace(bool trace) {
            this->trace = trace;
        }

        /// Set the maximum memory limit (trying to allocate more will throw MadnessException)
        void set_max_mem_limit(unsigned long max_mem_limit) {
            this->max_mem_limit = max_mem_limit;
        }
    };

    /// Returns pointer to internal structure
    WorldMemInfo* world_mem_info();
}
#endif // MADNESS_WORLD_WORLDMEM_H__INCLUDED
