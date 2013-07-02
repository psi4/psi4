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
#ifndef MADNESS_WORLD_PARAR_H__INCLUDED
#define MADNESS_WORLD_PARAR_H__INCLUDED

/// \file parar.h
/// \brief Implements ParallelInputArchive and ParallelOutputArchive

#include <world/archive.h>
#include <world/binfsar.h>

#include <unistd.h>
#include <cstring>
#include <cstdio>

namespace madness {
    namespace archive {

        /// Objects that implement their own parallel archive interface should derive from this
        class ParallelSerializableObject {};

        /// Base class for input and output parallel archives

        /// Templated by the local archive (only tested for BinaryFstream(In/Out)putArchive).
        template <typename Archive>
        class BaseParallelArchive {
            World* world;       //< Yep, the world.
            mutable Archive ar; //< The local archive
            int nio;            //< Number of IO nodes (always includes node zero)
            bool do_fence;      //< If true (default) read/write of parallel objects fence before and after IO
            char fname[256];    //< Name of the archive
            int nclient;        //< Number of clients of this node including self.  Zero if not IO node.

        public:
            static const bool is_parallel_archive = true;

            BaseParallelArchive() : world(0), ar(), nio(0), do_fence(true) {}

            /// Returns the process doing IO for given node

            /// Currently assigned in round-robin-fashion to the first nio processes except
	    /// on IBM BG/P where use every 64'th
            ProcessID io_node(ProcessID rank) const {
#ifdef HAVE_IBMBGP
		return ((rank/64)%nio) * 64;
#else
                return rank%nio;
#endif
            }

            /// Returns the process doing IO for this node
            ProcessID my_io_node() const {
                MADNESS_ASSERT(world);
                return io_node(world->rank());
            }

            /// Returns the number of IO clients for this node including self (zero if not an IO node)
            int num_io_clients() const {
                MADNESS_ASSERT(world);
                return nclient;
            }

            /// Returns true if this node is doing physical IO
            bool is_io_node() const {
                MADNESS_ASSERT(world);
                return world->rank() == my_io_node();
            }

            /// Returns pointer to the world
            World* get_world() const {
                MADNESS_ASSERT(world);
                return world;
            }

            /// Opens the parallel archive

            /// When writing a new archive, the numer of writers
            /// specified is used.  When reading an existing archive,
            /// the number of ionodes is adjusted to to be the same as
            /// the number that wrote the original archive.  Presently
            /// we don't have logic to handle reading an archive using
            /// fewer processes originally used to write it.  If you
            /// want to fix this have a look in worlddc.h for the only
            /// spot that currently needs changing to make that work.
            ///
            /// The default number of IO nodes is one and there is an
            /// arbitrary maximum of 50 set. On IBM BG/P the maximum
	    /// is nproc/64.
            void open(World& world, const char* filename, int nwriter=1) {
                this->world = &world;
                nio = nwriter;
#ifdef HAVE_IBMBGP
		int maxio = (world.size()-1)/64 + 1;
#else
		int maxio = 50;
#endif
                if (nio > maxio) nio = maxio; // Sanity?
                if (nio > world.size()) nio = world.size();

                MADNESS_ASSERT(filename);
                MADNESS_ASSERT(strlen(filename)-1<sizeof(fname));
                strcpy(fname,filename); // Save the filename for later
                char buf[256];
                MADNESS_ASSERT(strlen(filename)+7 <= sizeof(buf));
                sprintf(buf, "%s.%5.5d", filename, world.rank());

                if (world.rank() == 0) {
                    ar.open(buf);
                    ar & nio; // read/write nio from/to the archive
                    MADNESS_ASSERT(nio <= world.size());
                }

                // Ensure all agree on value of nio that may also have changed if reading
                world.gop.broadcast(nio, 0);

                // Other reader/writers can now open the local archive
                if (is_io_node() && world.rank()) {
                    ar.open(buf);
                }

                // Count #client
                ProcessID me = world.rank();
                nclient=0;
                for (ProcessID p=0; p<world.size(); ++p) if (io_node(p) == me) ++nclient;

//                 if (is_io_node()) {
//                     madness::print("I am an IO node with",nclient,"clients and file",buf);
//                 }
//                 else {
//                     madness::print("I am a client served by",my_io_node(),fname);
//                 }
            }

            /// Returns true if the named, unopened archive exists on disk with read access ... collective
            static bool exists(World& world, const char* filename) {
                char buf[256];
                MADNESS_ASSERT(strlen(filename)+7 <= sizeof(buf));
                sprintf(buf, "%s.%5.5d", filename, world.rank());
                bool status;
                if (world.rank() == 0)
                    status = (access(buf, F_OK|R_OK) == 0);

                world.gop.broadcast(status);

                return status;
            }

            /// Closes the parallel archive
            void close() {
                MADNESS_ASSERT(world);
                if (is_io_node()) ar.close();
            }

            /// Returns a reference to local archive ... throws if not an IO node
            Archive& local_archive() const {
                MADNESS_ASSERT(world);
                MADNESS_ASSERT(is_io_node());
                return ar;
            }

            /// Same as world.gop.broadcast_serializable(obj, root)
            template <typename objT>
            void broadcast(objT& obj, ProcessID root) const {
                get_world()->gop.broadcast_serializable(obj, root);
            }

            /// Deletes the files associated with the archive of the given name

            /// Presently assumes a shared file system since process zero does the
            /// deleting
            static void remove(World& world, const char* filename) {
                if (world.rank() == 0) {
                    char buf[256];
                    MADNESS_ASSERT(strlen(filename)+7 <= sizeof(buf));
                    for (ProcessID p=0; p<world.size(); ++p) {
                        sprintf(buf, "%s.%5.5d", filename, p);
                        if (::remove(buf)) break;
                    }
                }
            }

            /// Removes the files associated with the current archive
            void remove() {
                MADNESS_ASSERT(world);
                remove(*world, fname);
            }


            bool dofence() const {
                return this->do_fence;
            }

            void set_dofence(bool dofence) {
                this->dofence = do_fence;
            }
        };


        /// An archive for storing local or parallel data wrapping BinaryFstreamOutputArchive

        /// Writes of process local objects only stores the data from process zero.
        ///
        /// Writes of parallel containers (presently only WorldContainer) store all data.
        ///
        /// Each of the server or IO nodes creates a
        /// BinaryFstreamOutputArchive with the name filename.rank.  Client
        /// processes send their data to servers in round-robin fashion.
        ///
        /// Process zero records the number of writers so that when the archive is opened
        /// for reading the number of readers is forced to match.
        class ParallelOutputArchive : public BaseParallelArchive<BinaryFstreamOutputArchive>, public BaseOutputArchive {
        public:
            ParallelOutputArchive() {}

            /// Creates a parallel archive for output with given base filename and number of IO nodes
            ParallelOutputArchive(World& world, const char* filename, int nio=1)  {
                open(world, filename, nio);
            }

            void flush() {
                if (is_io_node()) local_archive().flush();
            }
        };

        /// An archive for storing local or parallel data wrapping BinaryFstreamInputArchive

        /// Reads of process local objects loads the value originally stored by process zero
        /// which is then broadcast to all processes.
        ///
        /// Reads of parallel containers (presently only WorldContainer) load all data.
        ///
        /// The number of IO nodes or readers is presently ignored.  It is
        /// forced to be the same as the original number of writers and
        /// therefore you cannot presently read an archive from a parallel job
        /// with fewer total processes than the number of writers.
        class ParallelInputArchive : public BaseParallelArchive<BinaryFstreamInputArchive>, public  BaseInputArchive {
        public:
            ParallelInputArchive() {}

            /// Creates a parallel archive for input
            ParallelInputArchive(World& world, const char* filename, int nio=1) {
                open(world, filename, nio);
            }
        };


        /// Disable type info for parallel output archives
        template <class T>
        struct ArchivePrePostImpl<ParallelOutputArchive,T> {
            static void preamble_store(const ParallelOutputArchive& ar) {};
            static inline void postamble_store(const ParallelOutputArchive& ar) {};
        };

        /// Disable type info for parallel input archives
        template <class T>
        struct ArchivePrePostImpl<ParallelInputArchive,T> {
            static inline void preamble_load(const ParallelInputArchive& ar) {};
            static inline void postamble_load(const ParallelInputArchive& ar) {};
        };


        template <class T>
        struct ArchiveImpl<ParallelOutputArchive, T> {
            /// Parallel objects are forwarded to their implementation of parallel store
            template <typename Q>
            static inline
            typename madness::enable_if<is_derived_from<Q, ParallelSerializableObject>, const ParallelOutputArchive&>::type
            wrap_store(const ParallelOutputArchive& ar, const Q& t) {
                ArchiveStoreImpl<ParallelOutputArchive,T>::store(ar,t);
                return ar;
            }

            /// Serial objects write only from process 0
            template <typename Q>
            static inline
            typename madness::disable_if<is_derived_from<Q, ParallelSerializableObject>, const ParallelOutputArchive&>::type
            wrap_store(const ParallelOutputArchive& ar, const Q& t) {
                if (ar.get_world()->rank()==0) {
                    ar.local_archive() & t;
                }
                return ar;
            }
        };

        template <class T>
        struct ArchiveImpl<ParallelInputArchive, T> {
            /// Parallel objects are forwarded to their implementation of parallel load
            template <typename Q>
            static inline
            typename madness::enable_if<is_derived_from<Q, ParallelSerializableObject>, const ParallelInputArchive&>::type
            wrap_load(const ParallelInputArchive& ar, const Q& t) {
                ArchiveLoadImpl<ParallelInputArchive,T>::load(ar,const_cast<T&>(t));
                return ar;
            }

            /// Serial objects read only from process 0 and broadcast results
            template <typename Q>
            static inline
            typename madness::disable_if<is_derived_from<Q, ParallelSerializableObject>, const ParallelInputArchive&>::type
            wrap_load(const ParallelInputArchive& ar, const Q& t) {
                if (ar.get_world()->rank()==0) {
                    ar.local_archive() & t;
                }
                ar.broadcast(const_cast<T&>(t), 0);
                return ar;
            }
        };


        /// Write archive array only from process zero
        template <class T>
        struct ArchiveImpl< ParallelOutputArchive, archive_array<T> > {
            static inline const ParallelOutputArchive& wrap_store(const ParallelOutputArchive& ar, const archive_array<T>& t) {
                if (ar.get_world()->rank() == 0) ar.local_archive() & t;
                return ar;
            }
        };

        /// Read archive array and broadcast
        template <class T>
        struct ArchiveImpl< ParallelInputArchive, archive_array<T> > {
            static inline const ParallelInputArchive& wrap_load(const ParallelInputArchive& ar, const archive_array<T>& t) {
                if (ar.get_world()->rank() == 0) ar.local_archive() & t;
                ar.broadcast(t, 0);
                return ar;
            }
        };

        /// Forward fixed size array to archive_array
        template <class T, std::size_t n>
        struct ArchiveImpl<ParallelOutputArchive, T[n]> {
            static inline const ParallelOutputArchive& wrap_store(const ParallelOutputArchive& ar, const T(&t)[n]) {
                ar << wrap(&t[0],n);
                return ar;
            }
        };

        /// Forward fixed size array to archive_array
        template <class T, std::size_t n>
        struct ArchiveImpl<ParallelInputArchive, T[n]> {
            static inline const ParallelInputArchive& wrap_load(const ParallelInputArchive& ar, const T(&t)[n]) {
                ar >> wrap(&t[0],n);
                return ar;
            }
        };
    }
}

#endif // MADNESS_WORLD_PARAR_H__INCLUDED
