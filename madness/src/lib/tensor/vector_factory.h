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


  $Id: vector_factory.h 1602 2009-12-27 19:53:06Z rjharrison $
*/


#ifndef MADNESS_TENSOR_VECTOR_FACTORY_H__INCLUDED
#define MADNESS_TENSOR_VECTOR_FACTORY_H__INCLUDED

#include <vector>

/// \file vector_factory.h
/// \brief Declares and implements factories for short vectors

/// \ingroup tensor

namespace madness {

    /// Returns a std::vector<T> initialized from the arguments
    template <typename T>
    inline std::vector<T> vector_factory(const T& v0) {
        std::vector<T> v(1);
        v[0] = v0;
        return v;
    }

    /// Returns a std::vector<T> initialized from the arguments
    template <typename T>
    inline std::vector<T> vector_factory(const T& v0, const T& v1) {
        std::vector<T> v(2);
        v[0] = v0;
        v[1] = v1;
        return v;
    }

    /// Returns a std::vector<T> initialized from the arguments
    template <typename T>
    inline std::vector<T> vector_factory(const T& v0, const T& v1,
                                         const T& v2) {
        std::vector<T> v(3);
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        return v;
    }

    /// Returns a std::vector<T> initialized from the arguments
    template <typename T>
    inline std::vector<T> vector_factory(const T& v0, const T& v1,
                                         const T& v2, const T& v3) {
        std::vector<T> v(4);
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        v[3] = v3;
        return v;
    }

    /// Returns a std::vector<T> initialized from the arguments
    template <typename T>
    inline std::vector<T> vector_factory(const T& v0, const T& v1,
                                         const T& v2, const T& v3,
                                         const T& v4) {
        std::vector<T> v(5);
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        v[3] = v3;
        v[4] = v4;
        return v;
    }

    /// Returns a std::vector<T> initialized from the arguments
    template <typename T>
    inline std::vector<T> vector_factory(const T& v0, const T& v1,
                                         const T& v2, const T& v3,
                                         const T& v4, const T& v5) {
        std::vector<T> v(6);
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        v[3] = v3;
        v[4] = v4;
        v[5] = v5;
        return v;
    }
}

#endif // MADNESS_TENSOR_VECTOR_FACTORY_H__INCLUDED
