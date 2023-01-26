/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#ifdef _MSC_VER
#include <Winsock2.h>
#include <winsock.h>
#else
#include <unistd.h>
#endif

#include <algorithm>
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

#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
#include "psi4/times.h"

#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/exception.h"

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
    Timer_Structure(Timer_Structure *parent, const std::string &key) : parent_ptr_(parent), key_(key) {
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
    bool isOff() {
        size_t thread_size = 0;
        switch (status_) {
            case ON:
                return false;
            case OFF:
                return true;
            case PARALLEL:
                thread_size = thread_timers_.size();
                for (size_t thread_rank = 0; thread_rank < thread_size; ++thread_rank) {
                    if (thread_timers_[thread_rank].get_status() != OFF) {
                        return false;
                    }
                }
                return true;
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
    bool is_empty(int thread_rank) {
        if (thread_rank == 0) {
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
                    if (!thread_timers_[0].is_empty()) {
                        return false;
                    }
                    break;
                default:
                    return false;
            }
        } else {
            if (status_ == PARALLEL) {
                if (!thread_timers_[thread_rank].is_empty()) {
                    return false;
                }
            }
        }
        for (auto child_iter = children_.begin(), end_iter = children_.end(); child_iter != end_iter; ++child_iter) {
            if (!child_iter->is_empty()) {
                return false;
            }
        }
        return true;
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
    const std::list<Timer_Structure> &get_children() const { return children_; }
    bool all_children_off() {
        for (auto child_iter = children_.begin(), end_iter = children_.end(); child_iter != end_iter; ++child_iter) {
            if (!child_iter->isOff()) {
                return false;
            }
        }
        return true;
    }
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
    void set_parent(Timer_Structure *parent_ptr) { parent_ptr_ = parent_ptr; }
    bool merge_move(Timer_Structure *another, int thread_rank = 0) {
        if (another == this) return false;
        if (another->is_empty(thread_rank)) {
            return another->is_empty();
        }
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
                                str = "Both timer with key ";
                                str += another->key_;
                                str += " are on and cannot be merged.";
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
    void merge_move_all(Timer_Structure *another) {
        if (another == this) return;
        size_t thread_size = 0, another_thread_size = 0;
        std::string str;
        Timer_thread thread_0_timer(status_, n_calls_, wall_start_, wtime_);
        switch (status_) {
            case ON:
            case OFF:
                switch (another->status_) {
                    case ON:
                        if (status_ == ON) {
                            str = "Both timer with key ";
                            str += another->key_;
                            str += " are on and cannot be merged.";
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
                        another_thread_size = another->thread_timers_.size();
                        thread_timers_.push_back(thread_0_timer);
                        if (thread_timers_[0].merge_move(&(another->thread_timers_[0]))) {
                            str = "Both timer with key ";
                            str += key_;
                            str += " on thread 0";
                            str += " are on and cannot be merged.";
                            throw PsiException(str, __FILE__, __LINE__);
                        }
                        for (size_t i = 1; i < another_thread_size; ++i) {
                            thread_timers_.push_back(another->thread_timers_[i]);
                        }
                        status_ = PARALLEL;
                        break;
                    default:
                        break;
                }
                break;
            case PARALLEL:
                switch (another->status_) {
                    case ON:
                    case OFF:
                        if (thread_timers_[0].merge_move(another)) {
                            std::string str = "Both timer with key ";
                            str += key_;
                            str += " on thread ";
                            str += std::to_string(0);
                            str += " are on and cannot be merged.";
                            throw PsiException(str, __FILE__, __LINE__);
                        }
                        break;
                    case PARALLEL:
                        another_thread_size = another->thread_timers_.size();
                        thread_size = thread_timers_.size();
                        if (thread_size <= another_thread_size) {
                            for (size_t i = 0; i < thread_size; ++i) {
                                if (thread_timers_[i].merge_move(&(another->thread_timers_[i]))) {
                                    std::string str = "Both timer with key ";
                                    str += key_;
                                    str += " on thread ";
                                    str += std::to_string(i);
                                    str += " are on and cannot be merged.";
                                    throw PsiException(str, __FILE__, __LINE__);
                                }
                            }
                            for (size_t i = thread_size; i < another_thread_size; ++i) {
                                thread_timers_.push_back(another->thread_timers_[i]);
                            }
                        } else {
                            for (size_t i = 0; i < another_thread_size; ++i) {
                                if (thread_timers_[i].merge_move(&(another->thread_timers_[i]))) {
                                    std::string str = "Both timer with key ";
                                    str += key_;
                                    str += " on thread ";
                                    str += std::to_string(i);
                                    str += " are on and cannot be merged.";
                                    throw PsiException(str, __FILE__, __LINE__);
                                }
                            }
                        }
                        break;
                    default:
                        break;
                }
                break;
            default:
                break;
        }
        for (auto child_iter = another->children_.begin(), end_iter = another->children_.end(); child_iter != end_iter;
             ++child_iter) {
            get_child(child_iter->key_)->merge_move_all(&(*child_iter));
        }
        another->children_.clear();
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

Timer_Structure root_timer(nullptr, ""), parallel_timer(nullptr, "");
std::list<Timer_Structure *> ser_on_timers;
std::vector<std::list<Timer_Structure *>> par_on_timers;
std::time_t timer_start, timer_end;
bool skip_timers;
static omp_lock_t lock_timer;

void print_timer(const Timer_Structure &timer, std::shared_ptr<PsiOutStream> printer, int align_key_width) {
    std::string key = timer.get_key();
    if (key.length() < align_key_width) {
        key.resize(align_key_width, ' ');
    }
    double wtime = std::chrono::duration_cast<std::chrono::duration<double>>(timer.get_total_wtime()).count();
    switch (timer.get_status()) {
        case ON:
        case OFF:
            printer->Printf("%s: %10.3fu %10.3fs %10.3fw %6zu calls\n", key.c_str(), timer.get_utime(),
                            timer.get_stime(), wtime, timer.get_n_calls());
            break;
        case PARALLEL:
            printer->Printf("%s: %10.3fp                         %6zu calls\n", key.c_str(), wtime, timer.get_n_calls());
        default:
            break;
    }
}

void print_nested_timer(const Timer_Structure &timer, std::shared_ptr<PsiOutStream> printer,
                        const std::string &indent) {
    const std::list<Timer_Structure> &children = timer.get_children();
    for (auto child_iter = children.begin(), end_child_iter = children.end(); child_iter != end_child_iter;
         ++child_iter) {
        printer->Printf("%s", indent.c_str());
        print_timer(*child_iter, printer, 36 - indent.length());
        print_nested_timer(*child_iter, printer, indent + "| ");
    }
}

bool empty_parallel() {
    extern std::vector<std::list<Timer_Structure *>> par_on_timers;
    size_t on_timer_size = par_on_timers.size();
    for (size_t i = 0; i < on_timer_size; ++i) {
        if (par_on_timers[i].size() != 0) {
            return false;
        }
    }
    return true;
}

/*!
** timer_init(): Initialize the linked list of timers
**
** \ingroup QT
*/
void timer_init() {
    omp_init_lock(&lock_timer);
    omp_set_lock(&lock_timer);
    extern std::time_t timer_start;
    extern Timer_Structure root_timer;
    timer_start = std::time(nullptr);
    root_timer.turn_on();
    extern std::list<Timer_Structure *> ser_on_timers;
    ser_on_timers.push_back(&root_timer);
    extern bool skip_timers;
    skip_timers = false;
    omp_unset_lock(&lock_timer);
}

/*!
** timer_done(): Close down all timers and write results to timer.dat
**
** \ingroup QT
*/
void timer_done() {
    extern std::time_t timer_start, timer_end;
    omp_set_lock(&lock_timer);
    extern Timer_Structure root_timer;
    root_timer.turn_off();
    char *host;

    host = (char *)malloc(40 * sizeof(char));
    gethostname(host, 40);

    /* Dump the timing data to timer.dat and free the timers */
    auto mode = std::ostream::app;
    auto printer = std::make_shared<PsiOutStream>("timer.dat", mode);
    printer->Printf("\n");
    printer->Printf("Host: %s\n", host);
    free(host);
    printer->Printf("\n");
    // NOTE: ctime is not thread-safe and could potentially cause overwriting of the return strings. 
    printer->Printf("Timers On : %s", ctime(&timer_start)); // lgtm[cpp/potentially-dangerous-function] 
    timer_end = std::time(nullptr);
    printer->Printf("Timers Off: %s", ctime(&timer_end)); // lgtm[cpp/potentially-dangerous-function] 
    printer->Printf("\nWall Time:  %10.2f seconds\n\n",
                    std::chrono::duration_cast<std::chrono::duration<double>>(root_timer.get_total_wtime()).count());
    printer->Printf("                                                       Time (seconds)\n");
    printer->Printf("Module                               %12s%12s%12s%13s\n", "User", "System", "Wall", "Calls");

    const std::list<Timer_Structure> timer_list = root_timer.summarize();
    for (auto timer_iter = timer_list.begin(), end_iter = timer_list.end(); timer_iter != end_iter; ++timer_iter) {
        print_timer(*timer_iter, printer, 36);
    }

    printer->Printf("\n--------------------------------------------------------------------------------------\n");

    print_nested_timer(root_timer, printer, "");

    printer->Printf("\n**************************************************************************************\n");

    omp_unset_lock(&lock_timer);
    omp_destroy_lock(&lock_timer);
}

void clean_timers() {
    timer_done();
    Timer_Structure new_root_timer(nullptr, ""), new_parallel_timer(nullptr, "");
    extern Timer_Structure root_timer, parallel_timer;
    root_timer = new_root_timer;
    parallel_timer = new_parallel_timer;
    timer_init();
}

void start_skip_timers() {
    omp_set_lock(&lock_timer);
    extern bool skip_timers;
    skip_timers = true;
    omp_unset_lock(&lock_timer);
}

void stop_skip_timers() {
    omp_set_lock(&lock_timer);
    extern bool skip_timers;
    skip_timers = false;
    omp_unset_lock(&lock_timer);
}

/*!
** timer_on(): Turn on the timer with the name given as an argument.  Can
** be turned on and off, time will accumulate while on.
** Should only be called out of OpenMP parallel section.
**
** \param key = Name of timer
**
** \ingroup QT
*/
PSI_API void timer_on(const std::string &key) {
    omp_set_lock(&lock_timer);
    extern bool skip_timers;
    if (skip_timers) {
        omp_unset_lock(&lock_timer);
        return;
    }
    extern Timer_Structure parallel_timer;
    if (parallel_timer.get_parent() != nullptr) {
        std::string str = "Unable to turn on serial Timer ";
        str += key;
        str += " when parallel timers are not all off.";
        throw PsiException(str, __FILE__, __LINE__);
    }
    extern std::list<Timer_Structure *> ser_on_timers;
    Timer_Structure *top_timer_ptr = nullptr;
    Timer_Structure *top_timer = ser_on_timers.back();
    if (key == top_timer->get_key()) {
        top_timer->turn_on();
    } else {
        top_timer_ptr = top_timer->get_child(key);
        ser_on_timers.push_back(top_timer_ptr);
        top_timer_ptr->turn_on();
    }
    omp_unset_lock(&lock_timer);
}

/*!
** timer_off(): Turn off the timer with the name given as an argument.  Can
** be turned on and off, time will accumulate while on.
** Should only be called out of OpenMP parallel section.
**
** \param key = Name of timer
**
** \ingroup QT
*/
PSI_API void timer_off(const std::string &key) {
    omp_set_lock(&lock_timer);
    extern bool skip_timers;
    if (skip_timers) {
        omp_unset_lock(&lock_timer);
        return;
    }
    extern Timer_Structure parallel_timer;
    extern std::list<Timer_Structure *> ser_on_timers;
    if (parallel_timer.get_parent() != nullptr) {
        std::string str = "Unable to turn on serial Timer ";
        str += key;
        str += " when parallel timers are not all off.";
        throw PsiException(str, __FILE__, __LINE__);
    }
    Timer_Structure *timer_ptr = nullptr;
    timer_ptr = ser_on_timers.back();
    if (key == timer_ptr->get_key()) {
        timer_ptr->turn_off();
        ser_on_timers.pop_back();
    } else {
        timer_ptr = nullptr;
        auto timer_iter = ser_on_timers.end();
        --timer_iter;
        std::list<std::string> stack_keys;
        stack_keys.push_front((*timer_iter)->get_key());
        auto iter_begin = ser_on_timers.begin();
        while (timer_iter != iter_begin) {
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
            str += " is not on.";
            throw PsiException(str, __FILE__, __LINE__);
        }
        timer_ptr->turn_off();
        auto on_child_iter = timer_iter;
        ++on_child_iter;
        Timer_Structure *on_child_ptr = *(on_child_iter);
        Timer_Structure *parent_ptr = timer_ptr->get_parent();
        Timer_Structure *parent_on_child_ptr = parent_ptr->get_child(on_child_ptr->get_key());
        if (parent_on_child_ptr->merge_move(on_child_ptr)) {
            timer_ptr->remove_child(on_child_ptr);
        }
        ser_on_timers.erase(timer_iter, ser_on_timers.end());
        for (auto stack_iter = stack_keys.begin(), stack_end = stack_keys.end(); stack_iter != stack_end;
             ++stack_iter) {
            on_child_ptr = parent_ptr->get_child(*stack_iter);
            ser_on_timers.push_back(on_child_ptr);
            parent_ptr = on_child_ptr;
        }
    }
    omp_unset_lock(&lock_timer);
}

/*!
** parallel_timer_on(): Turn on the timer with the name given as an argument.  Can
** be turned on and off, time will accumulate while on.
** Should only be called in OpenMP parallel sections.
**
** \param key = Name of timer
**
** \ingroup QT
*/
void parallel_timer_on(const std::string &key, int thread_rank) {
    omp_set_lock(&lock_timer);
    extern bool skip_timers;
    if (skip_timers) {
        omp_unset_lock(&lock_timer);
        return;
    }
    extern std::list<Timer_Structure *> ser_on_timers;
    extern std::vector<std::list<Timer_Structure *>> par_on_timers;
    extern Timer_Structure parallel_timer;
    if (par_on_timers.size() <= thread_rank) {
        par_on_timers.resize(thread_rank + 1);
    }
    if (parallel_timer.get_parent() == nullptr) {
        parallel_timer.set_parent(ser_on_timers.back());
    }
    Timer_Structure *top_timer_ptr = nullptr;
    if (par_on_timers[thread_rank].empty()) {
        top_timer_ptr = parallel_timer.get_child(key);
        par_on_timers[thread_rank].push_back(top_timer_ptr);
        top_timer_ptr->turn_on(thread_rank);
    } else {
        top_timer_ptr = par_on_timers[thread_rank].back();
        if (key == top_timer_ptr->get_key()) {
            top_timer_ptr->turn_on(thread_rank);
        } else {
            top_timer_ptr = top_timer_ptr->get_child(key);
            par_on_timers[thread_rank].push_back(top_timer_ptr);
            top_timer_ptr->turn_on(thread_rank);
        }
    }
    omp_unset_lock(&lock_timer);
}

/*!
** parallel_timer_off(): Turn off the timer with the name given as an argument.  Can
** be turned on and off, time will accumulate while on.
** Should only be called in OpenMP parallel sections.
**
** \param key = Name of timer
**
** \ingroup QT
*/
void parallel_timer_off(const std::string &key, int thread_rank) {
    omp_set_lock(&lock_timer);
    extern bool skip_timers;
    if (skip_timers) {
        omp_unset_lock(&lock_timer);
        return;
    }
    extern std::vector<std::list<Timer_Structure *>> par_on_timers;
    extern Timer_Structure parallel_timer;
    if (par_on_timers[thread_rank].empty()) {
        std::string str = "Timer ";
        str += key;
        str += " on thread ";
        str += std::to_string(thread_rank);
        str += " has never been turned on.";
        throw PsiException(str, __FILE__, __LINE__);
    }
    Timer_Structure *timer_ptr = nullptr;
    timer_ptr = par_on_timers[thread_rank].back();
    if (key == timer_ptr->get_key()) {
        timer_ptr->turn_off(thread_rank);
        par_on_timers[thread_rank].pop_back();
    } else {
        timer_ptr = nullptr;
        auto timer_iter = par_on_timers[thread_rank].end();
        --timer_iter;
        std::list<std::string> stack_keys;
        stack_keys.push_front((*timer_iter)->get_key());
        auto iter_begin = par_on_timers[thread_rank].begin();
        while (timer_iter != iter_begin) {
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
        par_on_timers[thread_rank].erase(timer_iter, par_on_timers[thread_rank].end());
        for (auto stack_iter = stack_keys.begin(), stack_end = stack_keys.end(); stack_iter != stack_end;
             ++stack_iter) {
            on_child_ptr = parent_ptr->get_child(*stack_iter);
            par_on_timers[thread_rank].push_back(on_child_ptr);
            parent_ptr = on_child_ptr;
        }
    }
    if (parallel_timer.get_parent() != nullptr && empty_parallel()) {
        parallel_timer.get_parent()->merge_move_all(&parallel_timer);
        parallel_timer.set_parent(nullptr);
    }
    omp_unset_lock(&lock_timer);
}
}  // namespace psi
