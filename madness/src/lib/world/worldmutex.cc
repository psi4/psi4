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

  $Id: $
*/

#include <world/worldmutex.h>
#include <world/worldtime.h>

/// \file worldmutex.h
/// \brief Implements Mutex, MutexFair, Spinlock, ConditionVariable


namespace madness {

    void MutexWaiter::wait() {
        //#ifdef HAVE_CRAYXT
#ifdef USE_SPINLOCKS
        // The value of 300 below is purely empirical but apparently a substantial
        // backoff (or something that is a by-product of waiting) is necessary
        // to avoid clobbering the memory subsystem while spinning on the taskq.
        // The time is  "Time to  run 100000 chain of tasks" from running world.
        // 1000--> 2+us
        // 400 --> 1.7us
        // 300 --> 1.7us
        // 250 --> 2us
        // 200 --> 3.6us
        // 100 --> 40+us (ouch!)

        for (int i=0; i<300; ++i)  cpu_relax();
#else
        const unsigned int nspin  = 1000;    // Spin for 1,000 calls
        const unsigned int nsleep = 100000;  // Sleep 10us for 100,000 calls = 1s
        if (count++ < nspin) return;
        else if (count < nsleep) yield(10);
        else yield(10000);
#endif
    }

    RecursiveMutex::RecursiveMutex() {
        // Create recursive mutex attribute
        pthread_mutexattr_t attr;
        int result = pthread_mutexattr_init(&attr);
        if (result) MADNESS_EXCEPTION("RecursiveMutex attribute initialization failed.", result);
        result = pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
        if (result) MADNESS_EXCEPTION("RecursiveMutex attribute set type failed.", result);

        // Initialize the mutex
        result = pthread_mutex_init(&mutex, &attr);
        if (result) MADNESS_EXCEPTION("RecursiveMutex initialization failed.", result);

        // Destroy the mutex attribute
        result = pthread_mutexattr_destroy(&attr);
        if (result) MADNESS_EXCEPTION("RecursiveMutex initialization failed.", result);
    }




} // namespace madness

