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


  $Id: dqueue.h 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/


#ifndef MADNESS_WORLD_DQUEUE_H__INCLUDED
#define MADNESS_WORLD_DQUEUE_H__INCLUDED

#include <world/worldmutex.h>
#include <cstddef>
#include <utility>
#include <algorithm>
#include <iostream>
#include <stdint.h>

/// \file dqueue.h
/// \brief Implements DQueue

namespace madness {

    struct DQStats { // Dilly bar, blizzard, ...
        uint64_t npush_back;    //< #calls to push_back
        uint64_t npush_front;   //< #calls to push_front
        uint64_t npop_front;    //< #calls to pop_front
        uint64_t ngrow;         //< #calls to grow
        uint64_t nmax;          //< Lifetime max. entries in the queue

        DQStats()
                : npush_back(0), npush_front(0), npop_front(0), ngrow(0), nmax(0) {}
    };

    /// A thread safe, fast but simple doubled-ended queue.

    /// Since the point is speed, the implementation is a circular
    /// buffer rather than a linked list so as to avoid the new/del
    /// overhead.  It will grow as needed, but presently will not
    /// shrink.  Had to modify STL API to make things thread safe.
    ///
    /// It is now rather heavily specialized to its only use.
    template <typename T>
    class DQueue : private CONDITION_VARIABLE_TYPE {
        char pad[64]; // To put the lock and the data in separate cache lines
        volatile size_t n __attribute__((aligned(64)));        ///< Number of elements in the buffer
        volatile size_t sz;              ///< Current capacity
        volatile T* volatile buf;        ///< Actual buffer
        volatile int _front;  ///< Index of element at front of buffer
        volatile int _back;    ///< Index of element at back of buffer
        DQStats stats;

        void grow() {
            // ASSUME WE ALREADY HAVE THE MUTEX WHEN IN HERE
            ++(stats.ngrow);
            if (sz != n) MADNESS_EXCEPTION("assertion failure in dqueue::grow", static_cast<int>(sz));
            size_t oldsz = sz;
            if (sz < 32768)
                sz = 65536;
            else if (sz <= 1048576)
                sz *= 2;
            else
                sz += 1048576;
            volatile T* volatile nbuf = new T[sz];
            int lo = sz/2 - oldsz/2;
            for (int i=_front; i<int(oldsz); ++i,++lo) {
                nbuf[lo] = buf[i];
            }
            if (_front > 0) {
                for (int i=0; i<=_back; ++i,++lo) {
                    nbuf[lo] = buf[i];
                }
            }
            _front = sz/2 - oldsz/2;
            _back = _front + n - 1;
            buf = nbuf;
            //sanity_check();
        }

        void sanity_check() const {
            // ASSUME WE ALREADY HAVE THE MUTEX WHEN IN HERE
            int num = _back - _front + 1;
            if (num < 0) num += sz;
            if (num==int(sz) && n==0) num=0;
            if (num==0 && n==sz) num=sz;
            //if (long(n) != num) print("size",sz,"front",_front,"back",_back,"n",n,"num",num);
            MADNESS_ASSERT(long(n) == num);
        }

        void push_back_with_lock(const T& value) {
            size_t nn = n;
            size_t ss = sz;
            if (nn == ss) {
                grow();
                ss = sz;
            }
            ++nn;
            if (nn > stats.nmax) stats.nmax = nn;
            n = nn;

            int b = _back + 1;
            if (b >= int(ss)) b = 0;
            buf[b] = value;
            _back = b;
            ++(stats.npush_back);

            signal();
        }


    public:
        DQueue(size_t hint=200000) // was 32768
                : n(0)
                , sz(hint>2 ? hint : 2)
                , buf(new T[sz])
                , _front(sz/2)
                , _back(_front-1) {}

        virtual ~DQueue() {
            delete buf;
        }

        /// Insert value at front of queue
        void push_front(const T& value) {
            madness::ScopedMutex<CONDITION_VARIABLE_TYPE> obolus(this);
            //sanity_check();

            size_t nn = n;
            size_t ss = sz;
            if (nn == ss) {
                grow();
                ss = sz;
            }
            ++nn;
            if (nn > stats.nmax) stats.nmax = nn;
            n = nn;

            int f = _front - 1;
            if (f < 0) f = ss - 1;
            buf[f] = value;
            _front = f;
            ++(stats.npush_front);

            //sanity_check();
            signal();
            //broadcast();
        }

        /// Insert element at back of queue (default is just one copy)
        void push_back(const T& value, int ncopy=1) {
            madness::ScopedMutex<CONDITION_VARIABLE_TYPE> obolus(this);
            //sanity_check();
            while (ncopy--)
                push_back_with_lock(value);
            //sanity_check();
            //broadcast();
        }

        template <typename opT>
        void scan(opT& op) {
            madness::ScopedMutex<CONDITION_VARIABLE_TYPE> obolus(this);

            int f = _front;
            size_t nn = n;
            int size = int(sz);
            std::cout << "IN Q " << nn << std::endl;

            while (nn--) {
                T* p = const_cast<T*>(buf + f);
                if (!op(p)) break;
                ++f;
                if (f >= size) f = 0;
            }
        }

        /// Pop multiple values off the front of queue ... returns number popped ... might be zero

        /// r must refer to an array of dimension at least nmax ... you are presently
        /// given no more than max(size()/64,1) values ... arbitrary choice.
        ///
        /// multi-threaded tasks might cause fewer tasks to be taken
        int pop_front(int nmax, T* r, bool wait) {
            madness::ScopedMutex<CONDITION_VARIABLE_TYPE> obolus(this);

            size_t nn = n;

            if (nn==0 && wait) {
                while (n == 0) // !!! Must be n (memory) not nn (local copy)
                    CONDITION_VARIABLE_TYPE::wait();

                nn = n;
            }

            ++(stats.npop_front);
            if (nn) {
                size_t thesize = sz;
                //sanity_check();

                nmax = std::min(nmax,std::max(int(nn>>6),1));
                int retval; // Will return the number of items taken


                int f = _front;

                // Original loop was this
                //retval = nmax;
                //while (nmax--) {
                //    *r++ = buf[f++];
                //    if (f >= int(sz)) f = 0;
                //}

                // New loop includes checking for replicated multi-threaded task
                // ... take one task and then check that subsequent tasks differ
                nmax--;
                *r++ = buf[f++];
                if (f >= int(thesize)) f = 0;
                retval=1;
                while (nmax--) {
                    T ptr = buf[f];
                    if (ptr == *r) {
                        break;
                    }
                    else if (ptr) { // Null pointer indicates stolen task
                        *r++ = ptr;
                        ++f;
                        if (f >= int(thesize)) f = 0;
                        ++retval;
                    }
                }

                n = nn - retval;
                _front = f;

                //sanity_check();
                return retval;
            }
            else {
                return 0;
            }
        }

        /// Pop value off the front of queue
        std::pair<T,bool> pop_front(bool wait) {
            T r;
            int ngot = pop_front(1, &r, wait);
            return std::pair<T,bool>(r,ngot==1);
        }

        size_t size() const {
            return n;
        }

        bool empty() const {
            return n==0;
        }

        const DQStats& get_stats() const {
            return stats;
        }
    };

}

#endif // MADNESS_WORLD_DQUEUE_H__INCLUDED
