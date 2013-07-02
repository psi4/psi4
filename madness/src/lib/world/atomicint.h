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

  $Id: atomicint.h 2384 2011-06-20 18:58:36Z xindongneu $
*/
#ifndef MADNESS_WORLD_ATOMICINT_H__INCLUDED
#define MADNESS_WORLD_ATOMICINT_H__INCLUDED

/// \file atomicint.h

/// \brief Implements AtomicInteger

#include <madness_config.h>

#ifdef HAVE_IBMBGP
#ifndef MADATOMIC_USE_GCC
#define MADATOMIC_USE_BGP
#endif
#elif defined(USE_X86_32_ASM) || defined(USE_X86_64_ASM)
#define MADATOMIC_USE_X86_ASM
#else
#define MADATOMIC_USE_GCC
#endif

#ifdef MADATOMIC_USE_BGP
#  include <bpcore/bgp_atomic_ops.h>
#elif defined(MADATOMIC_USE_AIX)
#  include <sys/atomic_op.h>
#elif defined(MADATOMIC_USE_GCC)
#  ifdef GCC_ATOMICS_IN_BITS
#    include <bits/atomicity.h>
#  else
#    include <ext/atomicity.h>
#  endif
#endif

namespace madness {

    /// An integer with atomic set, get, read+inc, read+dec, dec+test operations

    /// Only the default constructor is available and IT DOES NOT INITIALIZE THE VARIABLE.
    ///
    /// Conciously modeled after the TBB API to prepare for switching to it.
    class AtomicInt {
    private:

#ifdef MADATOMIC_USE_BGP
        typedef _BGP_Atomic atomic_int;
#else
        typedef volatile int atomic_int;
#endif
        atomic_int value;

        inline int exchange_and_add(int i) {
#ifdef MADATOMIC_USE_GCC
            return __gnu_cxx::__exchange_and_add(&value,i);
#elif defined(MADATOMIC_USE_X86_ASM)
            __asm__ __volatile__("lock; xaddl %0,%1" :"=r"(i) : "m"(value), "0"(i));
            return i;
#elif defined(MADATOMIC_USE_AIX)
            return fetch_and_add(&value,i);
#elif defined(MADATOMIC_USE_BGP)
            return _bgp_fetch_and_add(&value,i);
#else
#error ... atomic exchange_and_add operator must be implemented for this platform;
#endif
        }

    public:
        /// Returns the value of the counter with fence ensuring subsequent operations are not moved before the load
        operator int() const volatile {
#if defined(MADATOMIC_USE_BGP)
            int result = value.atom;
#else
	    int result = value;
#endif
	    // BARRIER to stop instructions migrating up
            __asm__ __volatile__ ("" : : : "memory");
            return result;
        }

        /// Sets the value of the counter with fence ensuring preceding operations are not moved after the store
        int operator=(int other) {
	    // BARRIER to stop instructions migrating down
            __asm__ __volatile__ ("" : : : "memory");
#if defined(MADATOMIC_USE_BGP)
            value.atom = other;
#else
            value = other;
#endif
            return other;
        }

        /// Sets the value of the counter with fences ensuring operations are not moved either side of the load+store
        AtomicInt& operator=(const AtomicInt& other) {
            *this = int(other);
            return *this;
        }

        /// Decrements the counter and returns the original value
        int operator--(int) {
            return exchange_and_add(-1);
        }

        /// Decrements the counter and returns the incremented value
        int operator--() {
            return exchange_and_add(-1) - 1;
        }

        /// Increments the counter and returns the original value
        int operator++(int) {
            return exchange_and_add(1);
        }

        /// Increments the counter and returns the incremented value
        int operator++() {
            return exchange_and_add(1) + 1;
        }

        /// Decrements the counter and returns true if the new value is zero
        bool dec_and_test() {
            return ((*this)-- == 1);
        }

#ifdef ATOMICINT_CAS
        /// Compare and swap

        /// Always returns original value; if (value == compare) value = newval.
        inline int compare_and_swap(int compare, int newval) {
#ifdef MADATOMIC_USE_GCC
            return __sync_val_compare_and_swap (&value, compare, newval);
#else
#error ... atomic exchange_and_add operator must be implemented for this platform;
#endif
        }
#endif

    }; // class AtomicInt

}
#endif // MADNESS_WORLD_ATOMICINT_H__INCLUDED
