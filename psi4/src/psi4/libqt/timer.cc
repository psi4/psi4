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
**
** Modified to use std::chrono::high_resolution_clock for the module wall times.
** Modified implementation of timer to a Class of tree structured timers
** and a stack of timer to enable nested timer prints.
** Implemented timer for OpenMP parallism.
**
** Tianyuan Zhang, June 2017
*/

#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include <sys/time.h>
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
typedef int omp_lock_t;
#define omp_init_lock(lock_timer_p) \
    do {                            \
    } while (0)
#define omp_set_lock(lock_timer_p) \
    do {                           \
    } while (0)
#define omp_unset_lock(lock_timer_p) \
    do {                             \
    } while (0)
#define omp_destroy_lock(lock_timer_p) \
    do {                               \
    } while (0)
#endif

/* guess for HZ, if missing */
#ifndef HZ
#define HZ 60
#endif

namespace psi {

using clock = std::chrono::high_resolution_clock;


enum Timer_Status { OFF, ON, PARALLEL };

class Timer_Structure;

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
    Timer_thread(Timer_Status status, size_t n_calls, clock::time_point wall_start, clock::duration wtime) {
        status_ = status;
        n_calls_ = n_calls;
        wall_start_ = wall_start;
        wtime_ = wtime;
    }
    bool turn_on() {
        if (status_ == ON) return true;
        status_ = ON;
        ++n_calls_;
        wall_start_ = clock::now();
        return false;
    }
    bool turn_off() {
        if (status_ == OFF) return true;
        status_ = OFF;
        wtime_ += clock::now() - wall_start_;
        return false;
    }
    Timer_Status get_status() const { return status_; }
    void set_status(Timer_Status status) { status_ = status; }
    bool isOn() { return status_ == ON; }
    bool is_empty() {
        switch (status_) {
            case ON:
                return false;
            case OFF:
                if (n_calls_ != 0) return false;
                if (wtime_ != clock::duration::zero()) return false;
                return true;
            default:
                break;
        }
        return false;
    }
    size_t get_n_calls() const { return n_calls_; }
    void set_n_calls(size_t n_calls) { n_calls_ = n_calls; }
    clock::time_point get_wall_start() const { return wall_start_; }
    clock::duration get_wtime() const { return wtime_; }
    void set_wtime(clock::duration wtime) { wtime_ = wtime; }
    bool merge_move(Timer_thread *another) {
        if (another == this) return false;
        switch (status_) {
            case ON:
                if (another->status_ != OFF) {
                    return true;
                }
                break;
            case OFF:
                if (another->status_ == ON) {
                    status_ = ON;
                    another->status_ = OFF;
                    wall_start_ = another->wall_start_;
                }
                break;
            default:
                break;
        }
        n_calls_ += another->n_calls_;
        another->n_calls_ = 0;
        wtime_ += another->wtime_;
        another->wtime_ = clock::duration::zero();
        return false;
    }
    bool merge_move(Timer_Structure *another);
    bool clear() {
        if (status_ != OFF) {
            return true;
        }
        n_calls_ = 0;
        wtime_ = clock::duration::zero();
        return false;
    }
    Timer_thread &operator+=(const Timer_thread &rhs) {
        switch (status_) {
            case ON:
                --n_calls_;
                status_ = OFF;
            case OFF:
                switch (rhs.status_) {
                    case ON:
                        --n_calls_;
                    case OFF:
                        n_calls_ += rhs.n_calls_;
                        wtime_ += rhs.wtime_;
                        break;
                    default:
                        break;
                }
                break;
            default:
                break;
        }
        return *this;
    }
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
    Timer_Structure(Timer_Structure *parent, std::string key) : parent_ptr_(parent), key_(key) {
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
    void set_status(Timer_Status status) { status_ = status; }
    bool isOn(int thread_rank = 0) {
        size_t thread_size = 0;
        switch (status_) {
            case ON:
                if (thread_rank == 0) {
                    return true;
                }
                return false;
            case OFF:
                return false;
            case PARALLEL:
                thread_size = thread_timers_.size();
                if (thread_rank >= thread_size) {
                    return false;
                }
                return thread_timers_[thread_rank].isOn();
            default:
                break;
        }
        return false;
    }
    bool is_empty() {
        size_t thread_size;
        switch (status_) {
            case ON:
                return false;
            case OFF:
                if (n_calls_ != 0) return false;
                if (utime_ != 0) return false;
                if (stime_ != 0) return false;
                if (wtime_ != clock::duration::zero()) return false;
                break;
            case PARALLEL:
                thread_size = thread_timers_.size();
                for (size_t i = 0; i < thread_size; ++i) {
                    if (!thread_timers_[i].is_empty()) {
                        return false;
                    }
                }
                break;
            default:
                return false;
        }
        for (auto child_iter = children_.begin(), end_iter = children_.end(); child_iter != end_iter; ++child_iter) {
            if (!child_iter->is_empty()) {
                return false;
            }
        }
        return true;
    }
    bool onThreadOtherThan(int thread_rank) {
        size_t thread_size = 0;
        switch (status_) {
            case ON:
                if (thread_rank == 0) {
                    return false;
                }
                return true;
            case OFF:
                return false;
            case PARALLEL:
                thread_size = thread_timers_.size();
                for (int i = 0; i < thread_size; ++i) {
                    if (i != thread_rank and thread_timers_[thread_rank].isOn()) {
                        return true;
                    }
                }
                return false;
            default:
                break;
        }
        return false;
    }
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
    void set_n_calls(size_t n_calls) { n_calls_ = n_calls; }
    double get_utime() const { return utime_; }
    void set_utime(double utime) { utime_ = utime; }
    double get_stime() const { return stime_; }
    void set_stime(double stime) { stime_ = stime; }
    clock::time_point get_wall_start() const { return wall_start_; }
    void set_wall_start(clock::time_point wall_start) { wall_start_ = wall_start; }
    clock::duration get_wtime() const { return wtime_; }
    void set_wtime(clock::duration wtime) { wtime_ = wtime; }
    clock::duration get_total_wtime() const {
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
    const std::vector<Timer_thread> &get_threads() const { return thread_timers_; }
    const std::list<Timer_Structure> &get_children() const { return children_; }
    Timer_Structure *find_child(const std::string &key) {
        auto end_iter = children_.end();
        for (auto child = children_.begin(); child != end_iter; ++child) {
            if (child->get_key() == key) {
                return &(*child);
            }
        }
        return nullptr;
    }
    Timer_Structure *get_child(const std::string &key) {
        Timer_Structure *child = find_child(key);
        if (child != nullptr) {
            return child;
        }
        children_.push_back(Timer_Structure(this, key));
        return &(children_.back());
    }
    bool remove_child(const std::string &key) {
        auto end_iter = children_.end();
        for (auto child = children_.begin(); child != end_iter; ++child) {
            if (child->get_key() == key) {
                children_.erase(child);
                return true;
            }
        }
        return false;
    }
    bool remove_child(Timer_Structure *timer_ptr) {
        auto end_iter = children_.end();
        for (auto child = children_.begin(); child != end_iter; ++child) {
            if (&(*child) == timer_ptr) {
                children_.erase(child);
                return true;
            }
        }
        return false;
    }
    Timer_Structure *get_parent() const { return parent_ptr_; }
    bool merge_move(Timer_Structure *another, int thread_rank = 0) {
        if (another == this) return false;
        size_t thread_size = 0, another_thread_size = 0;
        switch (status_) {
            case ON:
            case OFF:
                if (thread_rank == 0) {
                    std::string str;
                    Timer_thread thread_0_timer(status_, n_calls_, wall_start_, wtime_);
                    switch (another->status_) {
                        case ON:
                            if (status_ == ON) {
                                str = "Timer ";
                                str += another->key_;
                                str += " is still on and cannot be merged.";
                                throw PsiException(str, __FILE__, __LINE__);
                            } else {
                                status_ = ON;
                                another->status_ = OFF;
                                wall_start_ = another->wall_start_;
                                ontime_ = another->ontime_;
                            }
                        case OFF:
                            n_calls_ += another->n_calls_;
                            another->n_calls_ = 0;
                            utime_ += another->utime_;
                            another->utime_ = 0;
                            stime_ += another->stime_;
                            another->stime_ = 0;
                            wtime_ += another->wtime_;
                            another->wtime_ = clock::duration::zero();
                            break;
                        case PARALLEL:
                            thread_timers_.push_back(thread_0_timer);
                            if (thread_timers_[0].merge_move(&(another->thread_timers_[0]))) {
                                std::string str = "Both timer with key ";
                                str += key_;
                                str += " on thread ";
                                str += std::to_string(thread_rank);
                                str += " are on and cannot be merged.";
                                throw PsiException(str, __FILE__, __LINE__);
                            }
                            status_ = PARALLEL;
                            break;
                        default:
                            break;
                    }

                } else {
                    std::string str;
                    Timer_thread thread_0_timer(status_, n_calls_, wall_start_, wtime_);
                    switch (another->status_) {
                        case ON:
                        case OFF:
                            str = "Timer ";
                            str += another->key_;
                            str += " on thread ";
                            str += std::to_string(thread_rank);
                            str += " has never been turned on.";
                            throw PsiException(str, __FILE__, __LINE__);
                            break;
                        case PARALLEL:
                            another_thread_size = another->thread_timers_.size();
                            if (another_thread_size <= thread_rank) {
                                str = "Timer ";
                                str += another->key_;
                                str += " on thread ";
                                str += std::to_string(thread_rank);
                                str += " has never been turned on.";
                                throw PsiException(str, __FILE__, __LINE__);
                            }
                            thread_timers_.push_back(thread_0_timer);
                            thread_timers_.resize(thread_rank + 1);
                            if (thread_timers_[thread_rank].merge_move(&(another->thread_timers_[thread_rank]))) {
                                std::string str = "Both timer with key ";
                                str += key_;
                                str += " on thread ";
                                str += std::to_string(thread_rank);
                                str += " are on and cannot be merged.";
                                throw PsiException(str, __FILE__, __LINE__);
                            }
                            status_ = PARALLEL;
                            break;
                        default:
                            break;
                    }
                }
                break;
            case PARALLEL:
                if (thread_rank == 0) {
                    std::string str;
                    switch (another->status_) {
                        case ON:
                        case OFF:
                            if (thread_timers_[0].merge_move(another)) {
                                std::string str = "Both timer with key ";
                                str += key_;
                                str += " on thread ";
                                str += std::to_string(thread_rank);
                                str += " are on and cannot be merged.";
                                throw PsiException(str, __FILE__, __LINE__);
                            }
                            break;
                        case PARALLEL:
                            if (thread_timers_[0].merge_move(&(another->thread_timers_[0]))) {
                                std::string str = "Both timer with key ";
                                str += key_;
                                str += " on thread ";
                                str += std::to_string(thread_rank);
                                str += " are on and cannot be merged.";
                                throw PsiException(str, __FILE__, __LINE__);
                            }
                            break;
                        default:
                            break;
                    }
                } else {
                    std::string str;
                    switch (another->status_) {
                        case ON:
                        case OFF:
                            str = "Timer ";
                            str += another->key_;
                            str += " on thread ";
                            str += std::to_string(thread_rank);
                            str += " has never been turned on.";
                            throw PsiException(str, __FILE__, __LINE__);
                            break;
                        case PARALLEL:
                            another_thread_size = another->thread_timers_.size();
                            if (another_thread_size <= thread_rank) {
                                str = "Timer ";
                                str += another->key_;
                                str += " on thread ";
                                str += std::to_string(thread_rank);
                                str += " has never been turned on.";
                                throw PsiException(str, __FILE__, __LINE__);
                            }
                            thread_size = thread_timers_.size();
                            if (thread_size <= thread_rank) {
                                thread_timers_.resize(thread_rank + 1);
                            }
                            if (thread_timers_[thread_rank].merge_move(&(another->thread_timers_[thread_rank]))) {
                                std::string str = "Both timer with key ";
                                str += key_;
                                str += " on thread ";
                                str += std::to_string(thread_rank);
                                str += " are on and cannot be merged.";
                                throw PsiException(str, __FILE__, __LINE__);
                            }
                            break;
                        default:
                            break;
                    }
                }
                break;
            default:
                break;
        }
        std::list<Timer_Structure *> to_be_removed_children;
        for (auto child_iter = another->children_.begin(), end_iter = another->children_.end(); child_iter != end_iter;
             ++child_iter) {
            if (get_child(child_iter->key_)->merge_move(&(*child_iter), thread_rank)) {
                to_be_removed_children.push_back(&(*child_iter));
            }
        }
        for (auto child_iter = to_be_removed_children.begin(), end_iter = to_be_removed_children.end();
             child_iter != end_iter; ++child_iter) {
            another->remove_child(*child_iter);
        }
        if (another->is_empty()) {
            return true;
        }
        return false;
    }
    void clear(int thread_rank = 0) {
        size_t thread_size = 0;
        switch (status_) {
            case ON:
                if (thread_rank == 0) {
                    std::string str = "Timer ";
                    str += key_;
                    str += " is still on and cannot be cleared.";
                    throw PsiException(str, __FILE__, __LINE__);
                } else {
                    std::string str = "Timer ";
                    str += key_;
                    str += " on thread ";
                    str += std::to_string(thread_rank);
                    str += " has never been turned on.";
                    throw PsiException(str, __FILE__, __LINE__);
                }
            case OFF:
                if (thread_rank == 0) {
                    n_calls_ = 0;
                    utime_ = 0;
                    stime_ = 0;
                    wtime_ = clock::duration::zero();
                } else {
                    std::string str = "Timer ";
                    str += key_;
                    str += " on thread ";
                    str += std::to_string(thread_rank);
                    str += " has never been turned on.";
                    throw PsiException(str, __FILE__, __LINE__);
                }
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
                if (thread_timers_[thread_rank].clear()) {
                    std::string str = "Timer ";
                    str += key_;
                    str += " on thread ";
                    str += std::to_string(thread_rank);
                    str += " is still on and cannot be cleared.";
                    throw PsiException(str, __FILE__, __LINE__);
                }
            default:
                break;
        }
    }
    Timer_Structure &operator+=(const Timer_Structure &rhs) {
        size_t thread_size, rhs_thread_size;
        switch (status_) {
            case ON:
                --n_calls_;
                status_ = OFF;
            case OFF:
                switch (rhs.status_) {
                    case ON:
                        --n_calls_;
                    case OFF:
                        n_calls_ += rhs.n_calls_;
                        utime_ += rhs.utime_;
                        stime_ += rhs.stime_;
                        wtime_ += rhs.wtime_;
                        break;
                    case PARALLEL:
                        status_ = PARALLEL;
                        rhs_thread_size = rhs.thread_timers_.size();
                        thread_timers_.push_back(Timer_thread(OFF, n_calls_, wall_start_, wtime_));
                        thread_timers_.resize(rhs_thread_size);
                        for (size_t i = 0; i < rhs_thread_size; ++i) {
                            thread_timers_[i] += rhs.thread_timers_[i];
                        }
                        break;
                    default:
                        break;
                }
                break;
            case PARALLEL:
                switch (rhs.status_) {
                    case ON:
                        thread_timers_[0].set_n_calls(thread_timers_[0].get_n_calls() - 1);
                    case OFF:
                        thread_timers_[0].set_n_calls(thread_timers_[0].get_n_calls() + rhs.n_calls_);
                        thread_timers_[0].set_wtime(thread_timers_[0].get_wtime() + rhs.wtime_);
                        break;
                    case PARALLEL:
                        thread_size = thread_timers_.size();
                        rhs_thread_size = rhs.thread_timers_.size();
                        if (thread_size < rhs_thread_size) {
                            thread_timers_.resize(rhs_thread_size);
                        }
                        for (size_t i = 0; i < rhs_thread_size; ++i) {
                            thread_timers_[i] += rhs.thread_timers_[i];
                        }
                        break;
                    default:
                        break;
                }
                break;
            default:
                break;
        }
        return *this;
    }
    std::list<Timer_Structure> summarize() {
        std::list<Timer_Structure> timers;
        for (auto child_iter = children_.begin(), end_child_iter = children_.end(); child_iter != end_child_iter;
             ++child_iter) {
            auto sum_iter = timers.begin(), end_sum_iter = timers.end();
            for (; sum_iter != end_sum_iter; ++sum_iter) {
                if (child_iter->key_ == sum_iter->key_) {
                    (*sum_iter) += (*child_iter);
                    break;
                }
            }
            if (sum_iter == end_sum_iter) {
                Timer_Structure temp(nullptr, child_iter->key_);
                temp += (*child_iter);
                timers.push_back(temp);
            }
            std::list<Timer_Structure> nested_timers = child_iter->summarize();
            for (auto nested_iter = nested_timers.begin(), end_nested_iter = nested_timers.end();
                 nested_iter != end_nested_iter; ++nested_iter) {
                auto sum_iter = timers.begin(), end_sum_iter = timers.end();
                for (; sum_iter != end_sum_iter; ++sum_iter) {
                    if (nested_iter->key_ == sum_iter->key_) {
                        (*sum_iter) += (*nested_iter);
                        break;
                    }
                }
                if (sum_iter == end_sum_iter) {
                    Timer_Structure temp(nullptr, nested_iter->key_);
                    temp += (*nested_iter);
                    timers.push_back(temp);
                }
            }
        }
        return timers;
    }
};

bool Timer_thread::merge_move(Timer_Structure *another) {
    switch (status_) {
        case ON:
            if (another->get_status() != OFF) {
                return true;
            }
            break;
        case OFF:
            if (another->get_status() == ON) {
                status_ = ON;
                another->set_status(OFF);
                wall_start_ = another->get_wall_start();
            }
            break;
        default:
            break;
    }
    n_calls_ += another->get_n_calls();
    another->set_n_calls(0);
    wtime_ += another->get_wtime();
    another->set_utime(0);
    another->set_stime(0);
    another->set_wtime(clock::duration::zero());
    return false;
}

Timer_Structure root_timer(nullptr, "");
std::vector<std::list<Timer_Structure *>> on_timers;
//Timer_Structure* parallel_start_ptr;
time_t timer_start, timer_end;
static omp_lock_t lock_timer;

std::string formatTimeNumberPrint(double time, size_t sig_fig = 4, size_t width = 10) {
    std::string s = std::to_string(time);
    size_t l;
    std::string append_front = "";
    size_t dot_index = s.find('.');
    if (dot_index != std::string::npos){
        if (dot_index >= sig_fig) {
            s = s.substr(0, dot_index);
        } else if (time >= 1.0) {
            s = s.substr(0, sig_fig + 1);
        } else {
            size_t sig_start = s.find_first_not_of('0',dot_index + 1);
            l = s.length();
            if (l - sig_start > sig_fig) {
                s = s.substr(0, sig_start + sig_fig);
            }
        }
    }
    l = s.length();
    if (l < width) {
        append_front.resize(width - l, ' ');
    }
    return append_front + s;
}

void print_timer(const Timer_Structure &timer, std::shared_ptr<PsiOutStream> printer, int align_key_width) {
    std::string key = timer.get_key();
    if (key.length() < align_key_width) {
        key.resize(align_key_width, ' ');
    }
    double wtime = std::chrono::duration_cast<std::chrono::duration<double>>(timer.get_total_wtime()).count();
    switch (timer.get_status()) {
        case ON:
        case OFF:
            printer->Printf("%s: %su %ss %sw %6d calls\n", key.c_str(),
                            formatTimeNumberPrint(timer.get_utime()).c_str(),
                            formatTimeNumberPrint(timer.get_stime()).c_str(),
                            formatTimeNumberPrint(wtime).c_str(),
                            timer.get_n_calls());
            break;
        case PARALLEL:
            printer->Printf("%s: %sp                         %6d calls\n", key.c_str(),
                            formatTimeNumberPrint(wtime).c_str(),
                            timer.get_n_calls());
        default:
            break;
    }
}

void print_nested_timer(const Timer_Structure &timer, std::shared_ptr<PsiOutStream> printer, std::string indent) {
    const std::list<Timer_Structure> &children = timer.get_children();
    for (auto child_iter = children.begin(), end_child_iter = children.end(); child_iter != end_child_iter;
         ++child_iter) {
        printer->Printf("%s", indent.c_str());
        print_timer(*child_iter, printer, 28 - indent.length());
        print_nested_timer(*child_iter, printer, indent + "| ");
    }
}

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
    extern std::vector<std::list<Timer_Structure *>> on_timers;
//    extern Timer_Structure* parallel_start_ptr;
//    parallel_start_ptr = nullptr;
    timer_start = time(NULL);
    root_timer.turn_on();
    std::list<Timer_Structure *> thread_0_list;
    thread_0_list.push_back(&root_timer);
    on_timers.push_back(thread_0_list);
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
    char *host;

    host = (char *)malloc(40 * sizeof(char));
    gethostname(host, 40);

    /* Dump the timing data to timer.dat and free the timers */
    std::shared_ptr<PsiOutStream> printer(new PsiOutStream("timer.dat",std::ostream::app));
    printer->Printf("\n");
    printer->Printf("Host: %s\n", host);
    free(host);
    printer->Printf("\n");
    printer->Printf("Timers On : %s", ctime(&timer_start));
    timer_end = time(NULL);
    printer->Printf("Timers Off: %s", ctime(&timer_end));
    printer->Printf("\nWall Time:  %10.2f seconds\n\n",
                    std::chrono::duration_cast<std::chrono::duration<double>>(root_timer.get_total_wtime()).count());

    const std::list<Timer_Structure> timer_list = root_timer.summarize();
    for (auto timer_iter = timer_list.begin(), end_iter = timer_list.end(); timer_iter != end_iter; ++timer_iter) {
        print_timer(*timer_iter, printer, 28);
    }

    printer->Printf("\n-----------------------------------------------------------\n");

    print_nested_timer(root_timer, printer, "");

    printer->Printf("\n***********************************************************\n");

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
void timer_on(std::string key, int thread_rank) {
    omp_set_lock(&lock_timer);
    extern std::vector<std::list<Timer_Structure *>> on_timers;
//    extern Timer_Structure* parallel_start_ptr;
    Timer_Structure *top_timer_ptr = nullptr;
    Timer_Structure *thread_0_top = on_timers[0].back();
    if (thread_rank == 0) {
        if (key == thread_0_top->get_key()) {
            thread_0_top->turn_on(thread_rank);
        } else {
            top_timer_ptr = thread_0_top->get_child(key);
            top_timer_ptr->turn_on(thread_rank);
            on_timers[0].push_back(top_timer_ptr);
        }
    } else {
        size_t thread_size = on_timers.size();
        if (thread_rank >= thread_size) {
            on_timers.resize(thread_rank + 1);
        }
        if (on_timers[thread_rank].empty()) {
            while (thread_0_top != nullptr) {
                if (thread_0_top->get_key() == key) {
                    top_timer_ptr = thread_0_top;
                    break;
                }
                else {
                    Timer_Structure *temp = thread_0_top->find_child(key);
                    if (temp != nullptr) {
                        top_timer_ptr = temp;
                        break;
                    }
                }
                thread_0_top = thread_0_top->get_parent();
            }
            if (top_timer_ptr == nullptr) {
                top_timer_ptr = on_timers[0].back()->get_child(key);
            }
            top_timer_ptr->turn_on(thread_rank);
            on_timers[thread_rank].push_back(top_timer_ptr);
        } else {
            top_timer_ptr = on_timers[thread_rank].back();
            if (key == top_timer_ptr->get_key()) {
                top_timer_ptr->turn_on(thread_rank);
            } else {
                top_timer_ptr = top_timer_ptr->get_child(key);
                top_timer_ptr->turn_on(thread_rank);
                on_timers[thread_rank].push_back(top_timer_ptr);
            }
        }
    }
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
void timer_off(std::string key, int thread_rank) {
    omp_set_lock(&lock_timer);
    extern std::vector<std::list<Timer_Structure *>> on_timers;
//    extern Timer_Structure* parallel_start_ptr;
    Timer_Structure *timer_ptr = nullptr;
    if (on_timers[thread_rank].empty()) {
        std::string str = "Timer ";
        str += key;
        str += " on thread ";
        str += std::to_string(thread_rank);
        str += " has never been turned on.";
        throw PsiException(str, __FILE__, __LINE__);
    }
    timer_ptr = on_timers[thread_rank].back();
    if (key == timer_ptr->get_key()) {
        timer_ptr->turn_off(thread_rank);
        on_timers[thread_rank].pop_back();
    } else {
        timer_ptr = nullptr;
        auto timer_iter = on_timers[thread_rank].end();
        --timer_iter;
        std::list<std::string> stack_keys;
        stack_keys.push_front((*timer_iter)->get_key());
        auto iter_begin = on_timers[thread_rank].begin();
        for (; timer_iter != iter_begin;) {
            --timer_iter;
            if ((*timer_iter)->get_key() == key) {
                timer_ptr = *timer_iter;
                break;
            } else {
                stack_keys.push_front((*timer_iter)->get_key());
            }
        }
        if (timer_ptr == nullptr) {
            std::string str = "Timer ";
            str += key;
            str += " on thread ";
            str += std::to_string(thread_rank);
            str += " is not on.";
            throw PsiException(str, __FILE__, __LINE__);
        }
        timer_ptr->turn_off(thread_rank);
        auto on_child_iter = timer_iter;
        ++on_child_iter;
        Timer_Structure *on_child_ptr = *(on_child_iter);
        Timer_Structure *parent_ptr = timer_ptr->get_parent();
        Timer_Structure *parent_on_child_ptr = parent_ptr->get_child(on_child_ptr->get_key());
        if (parent_on_child_ptr->merge_move(on_child_ptr, thread_rank)) {
            timer_ptr->remove_child(on_child_ptr);
        }
        on_timers[thread_rank].erase(timer_iter, on_timers[thread_rank].end());
        for (auto stack_iter = stack_keys.begin(), stack_end = stack_keys.end(); stack_iter != stack_end;
             ++stack_iter) {
            on_child_ptr = parent_ptr->get_child(*stack_iter);
            on_timers[thread_rank].push_back(on_child_ptr);
            parent_ptr = on_child_ptr;
        }
    }
    omp_unset_lock(&lock_timer);
}
}
