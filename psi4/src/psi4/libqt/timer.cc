/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*!
** \file
** \brief Obtain user and system timings for blocks of code
** \ingroup QT
**
** TIMER.CC: These functions allow one to obtain user and system
** timings for arbitrary blocks of code.  If a code block is called
** repeatedly during the course of program execution, the timer
** functions will report the block's cumulative execution time and
** the number of calls. In addition, one may time multiple code blocks
** simultaneously, and even ``overlap'' timers.  Timing data is
** written to the file "timer.dat" at the end of timer execution,
** i.e., when timer_done() is called.
**
** To use the timer functions defined here:
**
** (1) Initialize the linked list of timers at the beginning of your
** program: timer_init();
**
** (2) Start a timer at the start of the block of code:
** timer_on("My Timer");
**
** (3) Stop the timer at the end of the block: timer_off("My Timer");
**
** (4) When all timer calls are complete, dump the linked list of
** timing data to the output file, "timer.dat": timer_done();
**
** NB this code uses system functions ctime(), time(), and times(),
** which may not quite be standard on all machines.
**
** T. Daniel Crawford, August 1999.
**
** Modified to use timeval structures for the module wall times, getting
** nanosecond precision on cumulated wall time instead of second.
** Useful to time integral computations where there can be millions
** to billion calls to functions.
**
** J. F. Gonthier, February 2016
*/

#include "psi4/libciomr/libciomr.h"
#include "psi4/libparallel/ParallelPrinter.h"
#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <sys/param.h>
#include <sys/times.h>
#include <unistd.h>

#include <chrono>
#include <list>
#include <map>
#include <string>
#include <utility>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
typedef int omp_lock_t;
#define omp_init_lock(lock_timer_p)                                            \
  do {                                                                         \
  } while (0)
#define omp_set_lock(lock_timer_p)                                             \
  do {                                                                         \
  } while (0)
#define omp_unset_lock(lock_timer_p)                                           \
  do {                                                                         \
  } while (0)
#define omp_destroy_lock(lock_timer_p)                                         \
  do {                                                                         \
  } while (0)
#endif

/* guess for HZ, if missing */
#ifndef HZ
#define HZ 60
#endif

namespace psi {

using clock = std::chrono::high_resolution_clock;

enum Timer_Status { OFF, ON, PARALLEL };

class Timer_thread {
private:
  Timer_Status status_;
  size_t n_calls_;
  clock::time_point wall_start_;
  clock::duration wtime_;

public:
  Timer_thread() {
    status_ = OFF;
    n_calls_ = 0;
    wtime_ = clock::duration::zero();
  }
  Timer_thread(Timer_Status status, size_t n_calls,
               clock::time_point wall_start, clock::duration wtime) {
    status_ = status;
    n_calls_ = n_calls;
    wall_start_ = wall_start;
    wtime_ = wtime;
  }
  bool turn_on() {
    if (status_ == ON)
      return true;
    status_ = ON;
    ++n_calls_;
    wall_start_ = clock::now();
    return false;
  }
  bool turn_off() {
    if (status_ == OFF)
      return true;
    status_ = OFF;
    wtime_ += clock::now() - wall_start_;
    return false;
  }
  Timer_Status get_status() const { return status_; }
  size_t get_n_calls() const { return n_calls_; }
  clock::time_point get_wall_start() const { return wall_start_; }
  clock::duration get_wtime() const { return wtime_; }
};

class Timer_Structure {
private:
  std::string key_;
  Timer_Status status_;
  size_t n_calls_;
  clock::time_point wall_start_;
  struct tms ontime_;
  double utime_;
  double stime_;
  clock::duration wtime_;
  std::vector<Timer_thread> thread_timers_;
  std::list<Timer_Structure> children_;
  Timer_Structure *parent_ptr_;

public:
  Timer_Structure(Timer_Structure *parent, std::string key)
      : parent_ptr_(parent), key_(key) {
    status_ = OFF;
    n_calls_ = 0;
    utime_ = 0;
    stime_ = 0;
    wtime_ = clock::duration::zero();
  }
  void turn_on(int thread_rank = 0) {
    size_t thread_size = 0;
    switch (status_) {
    case OFF:
      if (0 == thread_rank) {
        status_ = ON;
        ++n_calls_;
        times(&ontime_);
        wall_start_ = clock::now();
      } else {
        status_ = PARALLEL;
        Timer_thread thread_0_timer(OFF, n_calls_, wall_start_, wtime_);
        thread_timers_.push_back(thread_0_timer);
        thread_timers_.resize(thread_rank + 1);
        thread_timers_[thread_rank].turn_on();
      }
      break;
    case ON:
      if (0 == thread_rank) {
        std::string str = "Timer ";
        str += key_;
        str += " is already on.";
        throw PsiException(str, __FILE__, __LINE__);
      } else {
        status_ = PARALLEL;
        Timer_thread thread_0_timer(ON, n_calls_, wall_start_, wtime_);
        thread_timers_.push_back(thread_0_timer);
        thread_timers_.resize(thread_rank + 1);
        thread_timers_[thread_rank].turn_on();
      }
      break;
    case PARALLEL:
      thread_size = thread_timers_.size();
      if (thread_rank >= thread_size) {
        thread_timers_.resize(thread_rank + 1);
      }
      if (thread_timers_[thread_rank].turn_on()) {
        std::string str = "Timer ";
        str += key_;
        str += " on thread ";
        str += std::to_string(thread_rank);
        str += " is already on.";
        throw PsiException(str, __FILE__, __LINE__);
      }
      break;
    default:
      break;
    }
  }
  void turn_off(int thread_rank = 0) {
    size_t thread_size = 0;
    switch (status_) {
    case ON:
      if (0 == thread_rank) {
        tms offtime;
        status_ = OFF;
        times(&offtime);
        utime_ += ((double)(offtime.tms_utime - ontime_.tms_utime)) / HZ;
        stime_ += ((double)(offtime.tms_stime - ontime_.tms_stime)) / HZ;
        wtime_ += clock::now() - wall_start_;
      } else {
        std::string str = "Timer ";
        str += key_;
        str += " on thread ";
        str += std::to_string(thread_rank);
        str += " has never been turned on.";
        throw PsiException(str, __FILE__, __LINE__);
      }
      break;
    case OFF:
      if (0 == thread_rank) {
        std::string str = "Timer ";
        str += key_;
        str += " is already off.";
        throw PsiException(str, __FILE__, __LINE__);
      } else {
        std::string str = "Timer ";
        str += key_;
        str += " on thread ";
        str += std::to_string(thread_rank);
        str += " has never been turned on.";
        throw PsiException(str, __FILE__, __LINE__);
      }
      break;
    case PARALLEL:
      thread_size = thread_timers_.size();
      if (thread_rank >= thread_size) {
        std::string str = "Timer ";
        str += key_;
        str += " on thread ";
        str += std::to_string(thread_rank);
        str += " has never been turned on.";
        throw PsiException(str, __FILE__, __LINE__);
      }
      if (thread_timers_[thread_rank].turn_off()) {
        std::string str = "Timer ";
        str += key_;
        str += " on thread ";
        str += std::to_string(thread_rank);
        str += " is already off.";
        throw PsiException(str, __FILE__, __LINE__);
      }
      break;
    default:
      break;
    }
  }
  const std::string &get_key() const { return key_; }
  Timer_Status get_status() const { return status_; }
  size_t get_n_calls() const {
    if (status_ == PARALLEL) {
      size_t n_calls = 0;
      size_t thread_size = thread_timers_.size();
      for (size_t i = 0; i < thread_size; ++i) {
        n_calls += thread_timers_[i].get_n_calls();
      }
      return n_calls;
    }
    return n_calls_;
  }
  double get_utime() const { return utime_; }
  double get_stime() const { return stime_; }
  clock::time_point get_last_start_time() const { return wall_start_; }
  clock::duration get_wtime() const {
    if (status_ == PARALLEL) {
      clock::duration wtime = clock::duration::zero();
      size_t thread_size = thread_timers_.size();
      for (size_t i = 0; i < thread_size; ++i) {
        wtime += thread_timers_[i].get_wtime();
      }
      return wtime;
    }
    return wtime_;
  }
  const std::list<Timer_Structure> &get_children() const { return children_; }
  Timer_Structure &get_child(const std::string &key) {
    auto end_iter = children_.end();
    for (auto child = children_.begin(); child != end_iter; ++child) {
      if (child->get_key() == key) {
        return *child;
      }
    }
    children_.push_back(Timer_Structure(this, key));
    return children_.back();
  }
};

Timer_Structure root_timer(nullptr, "");
time_t timer_start, timer_end;
static omp_lock_t lock_timer;

/*!
** timer_init(): Initialize the linked list of timers
**
** \ingroup QT
*/
void timer_init(void) {
  extern time_t timer_start;
  omp_init_lock(&lock_timer);
  omp_set_lock(&lock_timer);
  extern Timer_Structure root_timer;
  timer_start = time(NULL);
  root_timer.turn_on();
  omp_unset_lock(&lock_timer);
}

/*!
** timer_done(): Close down all timers and write results to timer.dat
**
** \ingroup QT
*/
void timer_done(void) {
  extern time_t timer_start, timer_end;
  omp_set_lock(&lock_timer);
  extern Timer_Structure root_timer;
  root_timer.turn_off();
  // to be complete
  char *host;

  host = (char *)malloc(40 * sizeof(char));
  gethostname(host, 40);

  /* Dump the timing data to timer.dat and free the timers */
  std::shared_ptr<OutFile> printer(new OutFile("timer.dat", APPEND));
  printer->Printf("\n");
  printer->Printf("Host: %s\n", host);
  free(host);
  printer->Printf("\n");
  printer->Printf("Timers On : %s", ctime(&timer_start));
  timer_end = time(NULL);
  printer->Printf("Timers Off: %s", ctime(&timer_end));
  printer->Printf("\nWall Time:  %10.2f seconds\n\n",
                  std::chrono::duration_cast<std::chrono::duration<double>>(
                      root_timer.get_wtime())
                      .count());

  const std::list<Timer_Structure> &timer_list = root_timer.get_children();
  for (auto timer_iter = timer_list.begin(), end_iter = timer_list.end();
       timer_iter != end_iter; ++timer_iter) {
    double wtime = std::chrono::duration_cast<std::chrono::duration<double>>(
                       timer_iter->get_wtime())
                       .count();
    switch (timer_iter->get_status()) {
    case ON:
    //          printer->Printf("Warning: Timer hasn't turned off: ");
    case OFF:
      if (timer_iter->get_n_calls() > 1) {
        if (wtime < 10.0) {
          printer->Printf("%-12s: %10.2fu %10.2fs %10.6fw %6d calls\n",
                          timer_iter->get_key().c_str(),
                          timer_iter->get_utime(), timer_iter->get_stime(),
                          wtime, timer_iter->get_n_calls());
        } else {
          printer->Printf("%-12s: %10.2fu %10.2fs %10.2fw %6d calls\n",
                          timer_iter->get_key().c_str(),
                          timer_iter->get_utime(), timer_iter->get_stime(),
                          wtime, timer_iter->get_n_calls());
        }
      } else {
        if (wtime < 10.0) {
          printer->Printf("%-12s: %10.2fu %10.2fs %10.8fw %6d call\n",
                          timer_iter->get_key().c_str(),
                          timer_iter->get_utime(), timer_iter->get_stime(),
                          wtime, timer_iter->get_n_calls());
        } else {
          printer->Printf("%-12s: %10.2fu %10.2fs %10.2fw %6d call\n",
                          timer_iter->get_key().c_str(),
                          timer_iter->get_utime(), timer_iter->get_stime(),
                          wtime, timer_iter->get_n_calls());
        }
      }
      break;
    case PARALLEL:
      if (wtime < 10.0) {
        printer->Printf(
            "%-12s:       -.--u       -.--s       -.--w %10.6fp %6d calls\n",
            timer_iter->get_key().c_str(), wtime, timer_iter->get_n_calls());
      } else {
        printer->Printf(
            "%-12s:       -.--u       -.--s       -.--w %10.2fp %6d calls\n",
            timer_iter->get_key().c_str(), wtime, timer_iter->get_n_calls());
      }
    default:
      break;
    }
  }

  printer->Printf(
      "\n***********************************************************\n");

  omp_unset_lock(&lock_timer);
  omp_destroy_lock(&lock_timer);
}

/*!
** timer_on(): Turn on the timer with the name given as an argument.  Can
** be turned on and off, time will accumulate while on.
**
** \param key = Name of timer
**
** \ingroup QT
*/
void timer_on(const char *key, int thread_rank) {
  omp_set_lock(&lock_timer);
  extern Timer_Structure root_timer;
  std::string k(key);
  Timer_Structure &timer = root_timer.get_child(k);
  timer.turn_on(thread_rank);
  omp_unset_lock(&lock_timer);
}

/*!
** timer_off(): Turn off the timer with the name given as an argument.  Can
** be turned on and off, time will accumulate while on.
**
** \param key = Name of timer
**
** \ingroup QT
*/
void timer_off(const char *key, int thread_rank) {
  omp_set_lock(&lock_timer);
  extern Timer_Structure root_timer;
  std::string k(key);
  Timer_Structure &timer = root_timer.get_child(k);
  timer.turn_off(thread_rank);
  omp_unset_lock(&lock_timer);
}
}
