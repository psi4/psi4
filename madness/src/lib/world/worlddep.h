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


  $Id: worlddep.h 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/


#ifndef MADNESS_WORLD_WORLDDEP_H__INCLUDED
#define MADNESS_WORLD_WORLDDEP_H__INCLUDED

/// \file worlddep.h
/// \brief Defines DependencyInterface and CallbackInterface

#include <world/array.h>
#include <world/worldmutex.h>
#include <world/atomicint.h>

#include <typeinfo>

namespace madness {

    /// This class used for callbacks (e.g., for dependency tracking)
    class CallbackInterface {
    public:
        virtual void notify() = 0;

        virtual ~CallbackInterface() {}
    };


    /// Provides interface for tracking dependencies
    class DependencyInterface : public CallbackInterface, private Spinlock {
    private:
        AtomicInt ndepend;   ///< Counts dependencies

        static const int MAXCALLBACKS = 8;
        typedef Stack<CallbackInterface*,MAXCALLBACKS> callbackT;
        mutable volatile callbackT callbacks; ///< Called ONCE by dec() when ndepend==0

        void do_callbacks() const;

    public:
        DependencyInterface(int ndep = 0);


        /// Returns the number of unsatisfied dependencies
        int ndep() const;

        /// Returns true if ndepend == 0
        bool probe() const;


        /// Invoked by callbacks to notifiy of dependencies being satisfied
        void notify();


        /// Registers a callback for when ndepend=0

        /// If ndepend == 0, the callback is immediately invoked.
        void register_callback(CallbackInterface* callback);


        /// Increment the number of dependencies
        void inc();


        /// Decrement the number of dependencies and invoke callback if ndepend=0
        void dec();


        virtual ~DependencyInterface();
    };
}
#endif // MADNESS_WORLD_WORLDDEP_H__INCLUDED
