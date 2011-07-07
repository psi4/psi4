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


  $Id: worldtime.h 2363 2011-06-12 02:58:06Z rjharrison $
*/


#ifndef MADNESS_WORLD_WORLDTIME_H__INCLUDED
#define MADNESS_WORLD_WORLDTIME_H__INCLUDED

#include <stdint.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <madness_config.h>

#ifdef _CRAY
#include <catamount/dclock.h>
#endif

// #ifdef UINT64_T
// typedef UINT64_T uint64_t;
// #endif


/// \file worldtime.h
/// \brief Wrappers around platform dependent timers and performance info


namespace madness {

    /// Returns the wall time in seconds relative to arbitrary origin

    /// As accurate and lightweight as we can get it, but may not
    /// be any better than the gettime of day system call.
    double wall_time();

    /// On some machines we have access to a cycle count

    /// For small intervals this is probably the most lightweight and accurate timer
    /// but may not be meaningful over long intervals due to O/S scheduling,
    /// migration to different cores, frequency shifts, etc.  On x86 uses rtdsc.
    /// Otherwise uses wall_time() in nanoseconds.
    static inline uint64_t cycle_count() {
        uint64_t x;
#ifdef HAVE_IBMBGP
     unsigned int rx, ry, rz;
     do
     {
         asm volatile ( "mftbu %0" : "=r"(rx) );
         asm volatile ( "mftb %0" : "=r"(ry) );
         asm volatile ( "mftbu %0" : "=r"(rz) );
     }
     while ( rx != rz );
     x = rx;
     x = (x << 32) | ry;
#elif defined(X86_32)
__asm__ volatile(".byte 0x0f, 0x31" : "=A"(x));
#elif defined(X86_64)
        unsigned int a,d;
__asm__ volatile("rdtsc" : "=a"(a), "=d"(d));
        x = ((uint64_t)a) | (((uint64_t)d)<<32);
#else
        x = wall_time()*1e9;
#endif
        return x;
    }

    /// Estimates frequency of the processor in Hz

    /// First call may take about 0.1s to execute.  Subsequent
    /// calls return value cached from the first call so does
    /// not respond to changing processor frequency.
    ///
    /// If cycle_count() is returning wall_time() in nanoseconds
    /// this will return 1GHz.
    ///
    /// If not available returns 0.
    double cpu_frequency();


    /// Returns the cpu time in seconds relative to arbitrary origin

    /// As accurate and lightweight as we can get it, but may not
    /// be any better than the clock system call.
    static inline double cpu_time() {
#if defined(X86_32) || defined(X86_64) || defined(HAVE_IBMBGP)
        static const double rfreq = 1.0/cpu_frequency();
        return cycle_count()*rfreq;
#elif defined(_CRAY)
        return dclock();
#else
        return double(clock())/CLOCKS_PER_SEC;
#endif
    }


    /// Do nothing and especially do not touch memory
    inline void cpu_relax() {
#if defined(X86_32) || defined(X86_64)
        asm volatile("rep;nop" : : : "memory");
#elif defined(HAVE_IBMBGP)
	for (int i=0; i<25; i++) {
	    asm volatile ("nop\n");
	}
#else
#endif
    }


    /// Sleep or spin for specified no. of microseconds

    /// Wrapper to ensure desired behavior (and what is that one might ask??)
    static inline void myusleep(int us) {
#if defined(HAVE_CRAYXT) || defined(HAVE_IBMBGP)
        double secs = us*1e-6;
        double start = cpu_time();
        while (cpu_time()-start < secs) {
            for (int i=0; i<100; ++i) cpu_relax();
        }
#else
        usleep(us);
#endif
    }
}


#endif // MADNESS_WORLD_WORLDTIME_H__INCLUDED
