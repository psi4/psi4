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


  $Id$
*/


#ifndef MADNESS_WORLD_MPIAR_H__INCLUDED
#define MADNESS_WORLD_MPIAR_H__INCLUDED

#include <world/world.h>
#include <world/archive.h>
#include <world/vecar.h>

namespace madness {
    namespace archive {

        class MPIRawOutputArchive : public BaseOutputArchive {
            mutable World* world;
            ProcessID dest;
            int tag;
        public:
            MPIRawOutputArchive(World& world, const ProcessID& dest, int tag=SafeMPI::MPIAR_TAG)
                    : world(&world), dest(dest), tag(tag) {};

            template <class T>
            inline
            typename madness::enable_if< std::is_fundamental<T>, void >::type
            store(const T* t, long n) const {
                world->mpi.Send(t, n, dest, tag);
            }
        };

        class MPIRawInputArchive : public BaseInputArchive {
            mutable World* world;
            ProcessID src;
            int tag;
        public:
            MPIRawInputArchive(World& world, const ProcessID& src, int tag=SafeMPI::MPIAR_TAG)
                    : world(&world), src(src), tag(tag) {};

            template <class T>
            inline
            typename madness::enable_if< std::is_fundamental<T>, void >::type
            load(T* t, long n) const {
                world->mpi.Recv(t, n, src, tag);
            }
        };


        class MPIOutputArchive : public BaseOutputArchive {
            mutable World* world;
            ProcessID dest;
            int tag;
            const std::size_t bufsize;
            mutable std::vector<unsigned char> v;
            madness::archive::VectorOutputArchive var;
        public:
            MPIOutputArchive(World& world, const ProcessID& dest, int tag=SafeMPI::MPIAR_TAG)
                    : world(&world), dest(dest), tag(tag), bufsize(1024*1024), v(), var(v) {
                v.reserve(2*bufsize);
            };

            template <class T>
            inline
            typename madness::enable_if< std::is_fundamental<T>, void >::type
            store(const T* t, long n) const {
                if (v.size() > bufsize) flush();
                var.store(t,n);
                if (v.size() > bufsize) flush();
            }

            void flush() const {
                if (v.size()) {
                    world->mpi.Send(v.size(), dest, tag);
                    world->mpi.Send(&v[0], v.size(), dest, tag);
                    v.clear();
                    if (v.capacity() < 2*bufsize) v.reserve(2*bufsize); // ?? why ??
                }
            };

            void close() {
                flush();
            };

            ~MPIOutputArchive() {
                close();
            };
        };

        class MPIInputArchive : public BaseInputArchive {
            mutable World* world;
            ProcessID src;
            int tag;
            mutable std::vector<unsigned char> v;
            madness::archive::VectorInputArchive var;
        public:
            MPIInputArchive(World& world, const ProcessID& src, int tag=SafeMPI::MPIAR_TAG)
                    : world(&world), src(src), tag(tag), v(), var(v) {};

            template <class T>
            inline
            typename madness::enable_if< std::is_fundamental<T>, void >::type
            load(T* t, long n) const {
                if (!var.nbyte_avail()) {
                    var.rewind();
                    std::size_t m;
                    world->mpi.Recv(m, src, tag);
                    v.resize(m);
                    world->mpi.Recv(&v[0], m, src, tag);
                }
                var.load(t,n);
            }
        };

        // No type checking over MPI stream for efficiency
        template <class T>
        struct ArchivePrePostImpl<MPIRawOutputArchive,T> {
            static void preamble_store(const MPIRawOutputArchive& ar) {};
            static inline void postamble_store(const MPIRawOutputArchive& ar) {};
        };

        // No type checking over MPI stream for efficiency
        template <class T>
        struct ArchivePrePostImpl<MPIRawInputArchive,T> {
            static inline void preamble_load(const MPIRawInputArchive& ar) {};
            static inline void postamble_load(const MPIRawInputArchive& ar) {};
        };

        // No type checking over MPI stream for efficiency
        template <class T>
        struct ArchivePrePostImpl<MPIOutputArchive,T> {
            static void preamble_store(const MPIOutputArchive& ar) {};
            static inline void postamble_store(const MPIOutputArchive& ar) {};
        };

        // No type checking over MPI stream for efficiency
        template <class T>
        struct ArchivePrePostImpl<MPIInputArchive,T> {
            static inline void preamble_load(const MPIInputArchive& ar) {};
            static inline void postamble_load(const MPIInputArchive& ar) {};
        };


    }
}
#endif // MADNESS_WORLD_MPIAR_H__INCLUDED
