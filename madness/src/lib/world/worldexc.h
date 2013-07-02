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


  $Id: worldexc.h 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/


#ifndef MADNESS_WORLD_WORLDEXC_H__INCLUDED
#define MADNESS_WORLD_WORLDEXC_H__INCLUDED

/// \file worldexc.h
/// \brief Implements MadnessException

#include <iosfwd>
#include <exception>
#include <madness_config.h>

#ifndef MADNESS_DISPLAY_EXCEPTION_BREAK_MESSAGE
#define MADNESS_DISPLAY_EXCEPTION_BREAK_MESSAGE 1
#endif

namespace madness {

    /// Most exceptions thrown in MADNESS should be derived from these
    class MadnessException : public std::exception {
    public:
        const char* msg;
        const char* assertion;
        const int value;
        const int line;
        const char *function;
        const char *filename;

        // Capturing the line/function/filename info is best done with the macros below
        MadnessException(const char* msg, const char *assertion, int value,
                         int line, const char *function, const char *file)
                : msg(msg)
                , assertion(assertion)
                , value(value)
                , line(line)
                , function(function)
                , filename(file) {};

        virtual const char* what() const throw() {
            return msg;
        }

    };

    /// Print a MadnessException to the stream (for human consumption)

    /// Implemented in world.cc
    std::ostream& operator <<(std::ostream& out, const MadnessException& e);

    /// This function is executed just before a madness exception is thrown

    /// Implemented in world.cc
    void exception_break(bool);

#define MADNESS_EXCEPTION(msg,value)  { \
    madness::exception_break(MADNESS_DISPLAY_EXCEPTION_BREAK_MESSAGE); \
    throw madness::MadnessException(msg,0,value,__LINE__,__FUNCTION__,__FILE__); \
}


    /*
     * Default behaviour is MADNESS_ASSERTIONS throw a MADNESS exception
     *
     * Configure options are MADNESS_ASSERSIONS = THROW, ASSERT, DISABLE, ABORT
     *
     */

#ifdef MADNESS_ASSERTIONS_ABORT
#  define MADNESS_ASSERT(condition) \
     do {if (!(condition)) ((void (*)())0)();} while(0)
#endif

#ifdef MADNESS_ASSERTIONS_DISABLE
#  define MADNESS_ASSERT(condition)
#endif

#ifdef MADNESS_ASSERTIONS_ASSERT
#  include <cassert>
#  define MADNESS_ASSERT(condition) assert(condition)
#endif

#ifdef MADNESS_ASSERTIONS_THROW
#  define MADNESS_STRINGIZE(X) #X
#  define MADNESS_EXCEPTION_AT(F, L) MADNESS_STRINGIZE(F) "(" MADNESS_STRINGIZE(L) ")"
#  define MADNESS_ASSERT(condition) \
    do { \
        if (!(condition)) { \
            madness::exception_break(MADNESS_DISPLAY_EXCEPTION_BREAK_MESSAGE); \
            throw madness::MadnessException("MADNESS ASSERTION FAILED: " MADNESS_EXCEPTION_AT( __FILE__, __LINE__ ), \
                #condition,0,__LINE__,__FUNCTION__,__FILE__); \
        } \
    } while (0)
#endif

}

#endif // MADNESS_WORLD_WORLDEXC_H__INCLUDED
