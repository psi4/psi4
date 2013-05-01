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

#ifndef smartptr_timer_h_
#define smartptr_timer_h_

#include <string>
#include <map>
#include <iostream>
#include <utility>

#include "printstream.h"

namespace timer {


#define tstart(x) timer::Timer::start(x)
#define tstop(x) timer::Timer::stop(x)

class Timer {

    private:
        typedef std::pair<double, unsigned int> time_count_t;

        typedef std::map<std::string, time_count_t> active_map;

        typedef std::map<std::string, double> total_map;

        static active_map active_;

        static total_map totals_;

        static bool is_on_;
    
    public:
        
        /**
            Start the timer for a given region.  If this region is already
            active, it increases a "refcount" for the timer region.
            @param region  The name of the region to time
        */
        static void start(const std::string& region);

        /**
            Stop the timer for a given region. This decreases a "refcount"
            for the timer region.  If the refcount is 0, it stops the timer.
            @param region  The name of the region to time
        */
        static void stop(const std::string& region);

        /**
            @return unix time in seconds.
        */
        static double getTime();

        /**
            From a total time in seconds, figures out the number
            of hours, minutes, and seconds.
            @param total The total amount of time in seconds
            @param hours Reference returned
            @param minutes Reference returned
            @param seconds Referenced returned
        */
        static void getDisplay(
            double total,
            unsigned int& hours,
            unsigned int& minutes,
            double& seconds
        );

        /**
            Print time info for all of the regions
            @param os The stream to print to
        */
        static void print(std::ostream& os = std::cout);

        /**
            Activates the timer. Calls to stop and start will be registered.
        */
        static void turnOn();

        /**
            Deactivates the timer.  All calls to stop and start will
            be ignored.
        */
        static void turnOff();
        
};

}

#endif

