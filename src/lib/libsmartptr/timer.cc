#include <time.h>
#include <sys/time.h>
#include <cmath>

#include "timer.h"
#include "xmlarchive.h"

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
