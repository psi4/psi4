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


  $Id: bufar.h 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/


#ifndef MADNESS_WORLD_BUFAR_H__INCLUDED
#define MADNESS_WORLD_BUFAR_H__INCLUDED

/// \file bufar.h
/// \brief Implements an archive wrapping a memory buffer


#include <world/archive.h>
#include <world/print.h>
#include <cstring>


namespace madness {
    namespace archive {

        /// Wraps an archive around a memory buffer for output

        /// Type checking is disabled for efficiency.
        ///
        /// Throws MadnessException in case of buffer overflow
        ///
        /// Default constructor can be used to count stuff.
        class BufferOutputArchive : public BaseOutputArchive {
        private:
            unsigned char * const ptr;    // Buffer
            const std::size_t nbyte;      // Buffer size
            mutable std::size_t i;        // Current output location
            bool countonly;               // If true just count, don't copy
        public:
            BufferOutputArchive()
                    : ptr(0), nbyte(0), i(0), countonly(true) {}

            BufferOutputArchive(void* ptr, std::size_t nbyte)
                    : ptr((unsigned char *) ptr), nbyte(nbyte), i(0), countonly(false) {}

            template <class T>
            inline
            typename madness::enable_if< madness::is_serializable<T>, void >::type
            store(const T* t, long n) const {
                std::size_t m = n*sizeof(T);
                if (countonly) {
                    i += m;
                }
                else if (i+m > nbyte) {
                    madness::print("BufferOutputArchive:ptr,nbyte,i,n,m,i+m:",(void *)ptr,nbyte,i,n,m,i+m);
                    MADNESS_ASSERT(i+m<=nbyte);
                }
                else {
                    memcpy(ptr+i, t, m);
                    i += m;
                }
            }

            void open(std::size_t /*hint*/) {}

            void close() {}

            void flush() {}

            bool count_only() const { return countonly; }

            inline std::size_t size() const {
                return i;
            };
        };


        /// Convenience template computing the size of a buffer archive containing the arguments
        template <typename A, typename B, typename C, typename D, typename E, typename F, typename G, typename H, typename I, typename J>
        inline size_t bufar_size(const A& a, const B& b, const C& c, const D& d,
                                 const E& e, const F& f, const G& g, const H& h,
                                 const I& i, const J& j) {
            BufferOutputArchive count;
            count & a & b & c & d & e & f & g & h & i & j;
            return count.size();
        }


        /// Convenience template computing the size of a buffer archive containing the arguments
        template <typename A, typename B, typename C, typename D, typename E, typename F, typename G, typename H, typename I>
        inline size_t bufar_size(const A& a, const B& b, const C& c, const D& d,
                                 const E& e, const F& f, const G& g, const H& h,
                                 const I& i) {
            BufferOutputArchive count;
            count & a & b & c & d & e & f & g & h & i;
            return count.size();
        }


        /// Convenience template computing the size of a buffer archive containing the arguments
        template <typename A, typename B, typename C, typename D, typename E, typename F, typename G, typename H>
        inline size_t bufar_size(const A& a, const B& b, const C& c, const D& d,
                                 const E& e, const F& f, const G& g, const H& h) {
            BufferOutputArchive count;
            count & a & b & c & d & e & f & g & h;
            return count.size();
        }


        /// Convenience template computing the size of a buffer archive containing the arguments
        template <typename A, typename B, typename C, typename D, typename E, typename F, typename G>
        inline size_t bufar_size(const A& a, const B& b, const C& c, const D& d,
                                 const E& e, const F& f, const G& g) {
            BufferOutputArchive count;
            count & a & b & c & d & e & f & g;
            return count.size();
        }


        /// Convenience template computing the size of a buffer archive containing the arguments
        template <typename A, typename B, typename C, typename D, typename E, typename F>
        inline size_t bufar_size(const A& a, const B& b, const C& c, const D& d,
                                 const E& e, const F& f) {
            BufferOutputArchive count;
            count & a & b & c & d & e & f;
            return count.size();
        }


        /// Convenience template computing the size of a buffer archive containing the arguments
        template <typename A, typename B, typename C, typename D, typename E>
        inline size_t bufar_size(const A& a, const B& b, const C& c, const D& d,
                                 const E& e) {
            BufferOutputArchive count;
            count & a & b & c & d & e;
            return count.size();
        }


        /// Convenience template computing the size of a buffer archive containing the arguments
        template <typename A, typename B, typename C, typename D>
        inline size_t bufar_size(const A& a, const B& b, const C& c, const D& d) {
            BufferOutputArchive count;
            count & a & b & c & d;
            return count.size();
        }


        /// Convenience template computing the size of a buffer archive containing the arguments
        template <typename A, typename B, typename C>
        inline size_t bufar_size(const A& a, const B& b, const C& c) {
            BufferOutputArchive count;
            count & a & b & c;
            return count.size();
        }


        /// Convenience template computing the size of a buffer archive containing the arguments
        template <typename A, typename B>
        inline size_t bufar_size(const A& a, const B& b) {
            BufferOutputArchive count;
            count & a & b;
            return count.size();
        }


        /// Convenience template computing the size of a buffer archive containing the arguments
        template <typename A>
        inline size_t bufar_size(const A& a) {
            BufferOutputArchive count;
            count & a;
            return count.size();
        }


        /// Wraps an archive around a memory buffer for input

        /// Type checking is disabled for efficiency.
        ///
        /// Throws MadnessException in case of buffer overrun
        class BufferInputArchive : public BaseInputArchive {
        private:
            const unsigned char* const ptr;
            const std::size_t nbyte;
            mutable std::size_t i;

        public:
            BufferInputArchive(const void* ptr, std::size_t nbyte)
                    : ptr((const unsigned char *) ptr), nbyte(nbyte), i(0) {};

            template <class T>
            inline
            typename madness::enable_if< madness::is_serializable<T>, void >::type
            load(T* t, long n) const {
                std::size_t m = n*sizeof(T);
                MADNESS_ASSERT(m+i <=  nbyte);
                memcpy((unsigned char*) t, ptr+i, m);
                i += m;
            }

            void open() {};

            void rewind() const {
                i=0;
            };

            std::size_t nbyte_avail() const {
                return nbyte-i;
            };

            void close() {}
        };


        // No type checking over Buffer stream for efficiency
        template <class T>
        struct ArchivePrePostImpl<BufferOutputArchive,T> {
            static inline void preamble_store(const BufferOutputArchive& /*ar*/) {}
            static inline void postamble_store(const BufferOutputArchive& /*ar*/) {}
        };

        // No type checking over Buffer stream for efficiency
        template <class T>
        struct ArchivePrePostImpl<BufferInputArchive,T> {
            static inline void preamble_load(const BufferInputArchive& /*ar*/) {}
            static inline void postamble_load(const BufferInputArchive& /*ar*/) {}
        };
    }
}
#endif // MADNESS_WORLD_BUFAR_H__INCLUDED
