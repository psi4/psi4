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


  $Id: vecar.h 1602 2009-12-27 19:53:06Z rjharrison $
*/


#ifndef MADNESS_WORLD_VECAR_H__INCLUDED
#define MADNESS_WORLD_VECAR_H__INCLUDED

/// \file vecar.h
/// \brief Implements archive wrapping an STL vector

#include <vector>
#include <world/archive.h>


namespace madness {
    namespace archive {

// With a bit of thought this could be generalized to several STL containers


        /// Wraps an archive around an STL vector for output
        class VectorOutputArchive : public BaseOutputArchive {
            mutable std::vector<unsigned char>* v;
        public:
            VectorOutputArchive(std::vector<unsigned char>& v, std::size_t hint=262144) : v(&v) {
                open(hint);
            };

            template <class T>
            inline
            typename madness::enable_if< madness::is_serializable<T>, void >::type
            store(const T* t, long n) const {
                const unsigned char* ptr = (unsigned char*) t;
                v->insert(v->end(),ptr,ptr+n*sizeof(T));
            }

            void open(std::size_t hint=262144) {
                v->clear();
                v->reserve(hint);
            };

            void close() {};

            void flush() {};
        };


        /// Wraps an archive around an STL vector for input
        class VectorInputArchive : public BaseInputArchive {
            mutable std::vector<unsigned char>* v;
            mutable std::size_t i;
        public:
            VectorInputArchive(std::vector<unsigned char>& v) : v(&v) , i(0) {}

            template <class T>
            inline
            typename madness::enable_if< madness::is_serializable<T>, void >::type
            load(T* t, long n) const {
                std::size_t m = n*sizeof(T);
                if (m+i >  v->size()) MADNESS_EXCEPTION("VectorInputArchive: reading past end", m+1);
                memcpy((unsigned char*) t, &((*v)[i]), m);
                i += m;
            }

            void open() {};

            void rewind() const {
                i=0;
            };

            std::size_t nbyte_avail() const {
                return v->size()-i;
            };

            void close() {}
        };
    }
}
#endif // MADNESS_WORLD_VECAR_H__INCLUDED
