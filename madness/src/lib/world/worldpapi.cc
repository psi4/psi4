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

  $Id$
*/
#include <madness_config.h>
#ifdef HAVE_PAPI
#include <papi.h>
#include <world/worldthread.h>
#include <world/worldpapi.h>
namespace madness {


    static int events[NUMEVENTS] = {PAPI_FP_OPS};
    static Mutex papi_mutex;
    static volatile long long total_values[NUMEVENTS] = {0};
    static volatile int th=0;

    void initialize_papi() {
        if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT)
            MADNESS_EXCEPTION("Could not init PAPI", 1);
        if (PAPI_thread_init(pthread_self) != PAPI_OK)
            MADNESS_EXCEPTION("Could not init PAPI thread API", 1);

        reset_papi_measurement();
    }

    void begin_papi_measurement() {
        PAPI_start_counters(events, NUMEVENTS);
    }

    void end_papi_measurement() {
        // Ignore PAPI errors since don't want to kill the calculation
        // at its very end
        long long values[NUMEVENTS];
        PAPI_stop_counters(values, NUMEVENTS);
        papi_mutex.lock();
        ++th;
        //std::cout << "PAPITHREAD " << values[0] << std::endl;
        for (int i=0; i<NUMEVENTS; ++i) total_values[i] += values[i];
        papi_mutex.unlock();
    }

    void reset_papi_measurement() {
        for (int i=0; i<NUMEVENTS; ++i) total_values[i] = 0;
    }

    const long long* get_papi_measurement() {
        return const_cast<long long*>(total_values);
    }
}

#endif
