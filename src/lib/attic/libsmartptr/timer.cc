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

#include <time.h>
#include <sys/time.h>
#include <cmath>

#include "timer.h"
#include "xmlarchive.h"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

using namespace timer;
using namespace std;

Timer::active_map Timer::active_;
Timer::total_map Timer::totals_;
bool Timer::is_on_ = true;

void
Timer::start(const std::string& region)
{
    if (!is_on_)
        return;

    active_map::iterator it = active_.find(region);
    if (it != active_.end())
    {
        active_[region].second += 1;
        return;
    }

    time_count_t time(getTime(), 1);
    active_[region] = time;
}

double
Timer::getTime()
{
    struct timeval tod;
    gettimeofday(&tod, 0);
    double val = tod.tv_sec + 1e-6 * tod.tv_usec;
    return val;
}

void
Timer::stop(const std::string& region)
{
    if (!is_on_)
        return;

    active_map::iterator it = active_.find(region);
    if (it == active_.end())
    {
        cerr << "Cannot find time region " << region << endl;
        abort();
    }

    //decrement the ref count
    it->second.second -= 1;

    //refcount not zero
    if (it->second.second > 0) return;

    double current = getTime();
    double start = it->second.first;
    double increment = current - start;

    totals_[region] += increment;
    active_.erase(it);
}

void
Timer::getDisplay(double total, unsigned int& hours, unsigned int& minutes, double& seconds)
{
    seconds = total;

    unsigned int sec_p_hour = 3600;
    unsigned int sec_p_min = 60;
    hours = (unsigned int) floor(total / sec_p_hour);
    seconds -= sec_p_hour * hours;

    minutes =  (unsigned int) floor(total / sec_p_min);
    seconds -= sec_p_min * minutes;
}

void
Timer::print(ostream& os)
{
    os << endl << "Timer Statistics" << endl; 
    for (total_map::iterator it(totals_.begin()); it != totals_.end(); ++it)
    {
        unsigned int hours; unsigned int mins; double secs;
        getDisplay(it->second, hours, mins, secs);
        os << stream_printf("%20s: %5d Hours, %2d Minutes, %9.4f Seconds", it->first.c_str(), hours, mins, secs) << endl;
    }
    os << endl;
}

void
Timer::turnOn()
{
    is_on_ = true;
}

void
Timer::turnOff()
{
    is_on_ = false;
}
