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

