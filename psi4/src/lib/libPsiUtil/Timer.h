/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */
#ifndef SRC_LIB_LIBPSIUTIL_TIMER_H_
#define SRC_LIB_LIBPSIUTIL_TIMER_H_
#include <sstream>
#include <vector>
#include <boost/chrono/chrono.hpp>
#include <boost/timer/timer.hpp>
#include <boost/shared_ptr.hpp>
#include "../libparallel2/LibParallel2.h"
namespace psi{
///I abuse the integers for MPI stuff, so I'm hardcoding them
enum TimeTypes{WALL=0,CPU=1,SYSTEM=2};

typedef boost::timer::nanosecond_type Time_t;

/** \brief A little wrapper class to hold the returned time values
 *
 *  Each time value contains the measurement, as well as the accuracy of
 *  the clock that made it.  The member Time_ is your measurement,
 *  Resolution_ is it's accuracy.  Both of these are actually std::pairs,
 *  where the first member is the actual value, in seconds, and the second
 *  is the the maximum number of meaningful decimal places.
 *
 *  Basic example:
 *  if Time_={9.012345678,3} and Resolution_={1e-3,3} it means your time
 *  is 9.012 s +/- 0.001 s.
 */
class TimeValue{
   public:
      TimeValue(TimeTypes Type);
      std::pair<double,int> Time_;
      std::pair<double,int> Resolution_;
      std::string PrintOut(const int i=0)const;
};

/** \brief A timer object (only called SmartTimer b/c Timer was taken)
 *
 *   libqt contains a timer already.  I don't like it because it's not
 *   super accurate and it doesn't amend itself to parallel timings. Hence
 *   this timer's creation.  Times are returned to you either as a pretty
 *   string suited to printing or as TimeValue objects.  If you are going
 *   to use the TimeValue objects take a peek at their documentation.
 *
 *   In general, the wall times generated with this object are accurate to
 *   around a nanosecond, whereas CPU and system times are accurate to
 *   around a millisecond.  Given that a double can hold 16 sig figs, this
 *   means that if your wall time exceeds about 7.6 years, or your
 *   CPU/System times exceed about 7.6 million years, you will loose a
 *   decimal place of accuracy; you've been warned.
 *
 *   Things to keep in mind:
 *     - Construction starts the timer
 *     - Calling start resets the timer
 *     - My default this timer times all MPI processes individually
       - Percent load imbalance, L (T's are times, assume N processes) is:
 *       \f[
 *        L=\left(\frac{T_{Max}}{\frac{1}{N}\sum_{i}^{N}T_i}-1\right)\times\ 100\%
 *       \f]
 *     - Although common usage is likely to do something like:
 *       \code
 *       Timer.stop();
 *       std::cout<<Timer.PrintOut()<<std::endl;
 *       \endcode
 *       and it would hence make sense to combine them, the print out
 *       necessarily requires some MPI functionality involving the
 *       synching of timers
 *
 */
class SmartTimer{
   private:
      ///The name of the timer
      std::stringstream Name_;
      ///Whether or not this time is unique for each MPI process
      bool IsParallel_;
      /** \brief Returns the percent load balance
       *
       *  It's assumed this is being called inside PrintOut and hence
       *  that the Times will be a 3*3*NProcs long vector.  Thinking
       *  of it as a rank 3 tensor, the first index is the process
       *  number (0...NProc-1), the second is the time type being indexed
       *  from 0 to 2 mapping to WALL, CPU, and SYSTEM respectively, the
       *  final index is also from 0 to 2 and is respectively: the time in
       *  seconds, the number of decimal places (cast to a double),
       *  and the resolution in seconds.
       */
      double ComputeLoadImbalance(TimeTypes Type,std::vector<double>& Times)const;
      ///The actual timer
      boost::timer::cpu_timer Timer_;
      ///The starting times
      boost::timer::cpu_times StartTimes_;
      ///The stopped times
      boost::timer::cpu_times StopTimes_;
      boost::shared_ptr<const LibParallel::Communicator> Comm_;
   public:
      ///Constructs a timer, with parallel statistics if desired
      SmartTimer(const std::string& Name,bool IsParallel=true);
      ///Starts (or restarts) the timer
      void Start();
      ///Stops the timer
      void Stop();
      ///Resumes timing (if currently timing does nothing)
      void Resume();
      ///Returns a desired time type in seconds
      TimeValue GetTime(TimeTypes Type)const;
      ///Returns true if the timer is running
      bool IsStopped()const;
      ///Returns the statistics in a pretty string
      std::string PrintOut()const;
};

}



#endif /* SRC_LIB_LIBPSIUTIL_TIMER_H_ */