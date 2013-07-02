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


  $Id: enable_if.h 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/

#ifndef MADNESS_WORLD_ENABLE_IF_H__INCLUDED
#define MADNESS_WORLD_ENABLE_IF_H__INCLUDED

//#include <world/typestuff.h>

namespace madness {

    /// enable_if_c from Boost for conditionally instantiating templates based on type

    /// Evaluates to \c returnT if \c B is true, otherwise to an invalid type expression
    /// which causes the template expression in which it is used to not be considered for
    /// overload resolution.
    template <bool B, class returnT = void>
    struct enable_if_c {
        typedef returnT type;
    };

    /// enable_if_c from Boost for conditionally instantiating templates based on type
    template <class returnT>
    struct enable_if_c<false, returnT> {};

    /// enable_if from Boost for conditionally instantiating templates based on type

    /// Evaluates to \c returnT if \c Cond::value is true, otherwise to an invalid type expression
    /// which causes the template expression in which it is used to not be considered for
    /// overload resolution.
    template <class Cond, class returnT = void>
    struct enable_if : public enable_if_c<Cond::value, returnT> {};

    /// disable_if from Boost for conditionally instantiating templates based on type

    /// Evaluates to \c returnT if \c Cond::value is false, otherwise to an invalid type expression
    /// which causes the template expression in which it is used to not be considered for
    /// overload resolution.
    template <bool B, class returnT = void>
    struct disable_if_c {
        typedef returnT type;
    };

    /// disable_if from Boost for conditionally instantiating templates based on type
    template <class returnT>
    struct disable_if_c<true, returnT> {};

    /// disable_if from Boost for conditionally instantiating templates based on type
    template <class Cond, class returnT = void>
    struct disable_if : public disable_if_c<Cond::value, returnT> {};

    /// lazy_enable_if_c from Boost for conditionally instantiating templates based on type

    /// Evaluates to \c returnT if \c B is true, otherwise to an invalid type expression
    /// which causes the template expression in which it is used to not be considered for
    /// overload resolution. This "lazy" version is used if \c T is only valid when
    /// B is true. Note: typename T::type is the return type and must be well formed.
    template <bool B, class returnT>
    struct lazy_enable_if_c {
      typedef typename returnT::type type;
    };

    /// lazy_enable_if_c from Boost for conditionally instantiating templates based on type
    template <class returnT>
    struct lazy_enable_if_c<false, returnT> { };

    /// lazy_enable_if from Boost for conditionally instantiating templates based on type

    /// Evaluates to \c returnT if \c Cond::value is true, otherwise to an invalid type expression
    /// which causes the template expression in which it is used to not be considered for
    /// overload resolution. This "lazy" version is used if \c returnT is only valid when
    /// Cond::value is true. Note: typename T::type is the return type and must be well formed.
    template <class Cond, class returnT>
    struct lazy_enable_if : public lazy_enable_if_c<Cond::value, returnT> { };

    /// lazy_disable_if_c from Boost for conditionally instantiating templates based on type

    /// Evaluates to \c returnT if \c B is false, otherwise to an invalid type expression
    /// which causes the template expression in which it is used to not be considered for
    /// overload resolution. This "lazy" version is used if \c returnT is only valid
    /// when B is false. Note: typename T::type is the return type and must be well formed.
    template <bool B, class returnT>
    struct lazy_disable_if_c {
      typedef typename returnT::type type;
    };

    /// lazy_disable_if_c from Boost for conditionally instantiating templates based on type
    template <class returnT>
    struct lazy_disable_if_c<true, returnT> {};

    /// lazy_disable_if from Boost for conditionally instantiating templates based on type

    /// Evaluates to \c returnT if \c Cond::value is false, otherwise to an invalid type expression
    /// which causes the template expression in which it is used to not be considered for
    /// overload resolution. This "lazy" version is used if \c returnT is only valid when
    /// Cond::value is false. Note: typename T::type is the return type and must be well formed.
    template <class Cond, class returnT>
    struct lazy_disable_if : public lazy_disable_if_c<Cond::value, returnT> {};

    /// enable_if_same (from Boost?) for conditionally instantiating templates if two types are equal

    /// Use example
    /// \code
    ///     template <class T> A(T& other, typename enable_if_same<A const,T>::type = 0) {
    /// \endcode
//    template <class T, class U, class returnT = void>
//    struct enable_if_same : public enable_if<std::is_same<T,U>, returnT> {};

} // namespace madness


/* Macros to make some of this stuff more readable */

/**

   \def ENABLE_IF(CONDITION,TYPEIFTRUE)
   \brief Macro to make enable_if<> template easier to use

   \def DISABLE_IF(CONDITION,TYPEIFTRUE)
   \brief Macro to make disable_if<> template easier to use


   \def DISABLE_IF(CONDITION,TYPEIFTRUE)
   \brief Macro to make enable_if<std::is_same< A , B > > template easier to use

*/

#define ENABLE_IF(CONDITION,TYPEIFTRUE)  typename madness::enable_if< CONDITION, TYPEIFTRUE >::type
#define DISABLE_IF(CONDITION,TYPEIFTRUE) typename madness::disable_if< CONDITION, TYPEIFTRUE >::type
#define ENABLE_IF_SAME(A,B,TYPEIFTRUE) typename madness::enable_if<std::is_same< A , B >, TYPEIFTRUE >::type
#define DISABLE_IF_SAME(A,B,TYPEIFTRUE) typename madness::disable_if<std::is_same< A , B >, TYPEIFTRUE >::type

#endif // MADNESS_WORLD_ENABLE_IF_H__INCLUDED
