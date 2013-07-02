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


  $Id: worldptr.h 2217 2011-03-13 19:09:11Z justus.c79@gmail.com $
*/

#ifndef MADNESS_WORLD_WORLDPTR_H__INCLUDED
#define MADNESS_WORLD_WORLDPTR_H__INCLUDED

#include <world/worldexc.h>     // for MADNESS_ASSERT
#include <world/worldtypes.h>   // for ProcessID
#include <world/archive.h>      // for wrap_opaque
#include <algorithm>            // for std::swap
#include <iostream>             // for std::iostream

namespace madness {

    // Forward decl
    class World;

    namespace detail {

        // Forward decl
        ProcessID world_rank(const World&);
        unsigned int world_id(const World&);
        World* world_from_id(unsigned int);

        template<typename U>
        struct ptr_traits {
            typedef U & reference;
        };

        template<>
        struct ptr_traits<void> {
            typedef void reference;
        };

        /// A global pointer address, valid anywhere in the world.

        /// Stores a globally addressable pointer. It can be sent to any
        /// process in the world.
        /// \tparam T The pointer type
        template <typename T>
        class WorldPtr {
        public:
            typedef unsigned long worldidT; ///< World ID type

        private:
            World* world_;          ///< A pointer to the world
            worldidT worldid_;      ///< The world id
            ProcessID rank_;        ///< The rank of the node that the pointer belongs to
            T* pointer_;            ///< The pointer being referenced

            template<typename>
            friend class WorldPtr;

            /// Current local rank

            /// \return The rank of the current node. If the pointer is not
            /// not set, then -2.
            /// \note -2 is returned so it is not equal to the null value of -1.
            /// \throw nothing
            ProcessID local_rank() const { return (world_ != NULL ? world_rank(*world_) : -2); }

        public:

            typedef T* pointer;
            typedef typename ptr_traits<T>::reference reference;

            /// Default constructor

            /// Creates a NULL pointer. There is no owner (i.e. the owner is set
            /// to -1).
            /// \throw nothing
            WorldPtr() :
                world_(NULL),
                worldid_(0),
                rank_(-1),
                pointer_(NULL)
            { }

            /// World pointer constructor

            /// Construct a world pointer form a local pointer.
            /// \param w A reference to the local world.
            /// \param p The local pointer.
            /// \throw nothing
            WorldPtr(World& w, T* p) :
                world_(&w),
                worldid_(detail::world_id(w) + 1),
                rank_(detail::world_rank(w)),
                pointer_(p)
            { }

            /// Copy constructor

            /// \param other The world pointer to be copied
            /// \throw nothing
            WorldPtr(const WorldPtr<T>& other) :
                world_(other.world_),
                worldid_(other.worldid_),
                rank_(other.rank_),
                pointer_(other.pointer_)
            { }

            /// Copy conversion constructor

            /// Copy and convert pointer from \c U* to \c T* type.
            /// \tparam U The pointer type of the \c other pointer
            /// \param other The world pointer to be copied
            /// \throw nothing
            /// \note \c U* must be implicitly convertible to T* type.
            template <typename U>
            WorldPtr(const WorldPtr<U>& other) :
                world_(other.world_),
                worldid_(other.worldid_),
                rank_(other.rank_),
                pointer_(other.pointer_)
            { }

            /// Copy assignment operator

            /// \param other The world pointer to be copied
            /// \return A reference to this object
            /// \throw nothing
            WorldPtr<T>& operator=(const WorldPtr<T>& other) {
                world_ = other.world_;
                worldid_ = other.worldid_;
                rank_ = other.rank_;
                pointer_ = other.pointer_;

                return *this;
            }

            /// Copy conversion assignment operator

            /// Copy and convert pointer from \c U* to \c T* type.
            /// \tparam U The pointer type of the \c other pointer
            /// \param other The world pointer to be copied
            /// \return A reference to this object
            /// \throw nothing
            /// \note \c U* must be implicitly convertible to T* type.
            template <typename U>
            WorldPtr<T>& operator=(const WorldPtr<U>& other) {
                world_ = other.world_;
                worldid_ = other.worldid_;
                rank_ = other.rank_;
                pointer_ = other.pointer_;

                return *this;
            }

            /// Check that the world pointer references a local pointer.

            /// \return \c true when the pointer points to a local address, and
            /// \c false when it points to a remote address or is NULL.
            /// \throw nothing
            bool is_local() const { return local_rank() == rank_; }

            /// Check that the world pointer has an owner

            /// \return true when the pointer has a valid owner
            /// \throw nothing
            bool has_owner() const { return (rank_ != -1) && (world_ != NULL); }

            /// Pointer accessor

            /// Get the pointer of the world pointer. Note: A default initialized
            /// pointer is not considered to be local because it is not associated
            /// with a world.
            /// \return The local pointer.
            /// \throw MadnessException When the pointer references a remote
            /// address.
            pointer get() const {
                // It is not safe to access this pointer remotely unless null.
                MADNESS_ASSERT(is_local());
                return pointer_;
            }

            /// Dereference operator

            /// Dereference the local pointer.
            /// \return A reference to the local pointer
            /// \throw MadnessException When the pointer references a remote
            /// address, or when the pointer is NULL.
            reference operator*() const {
                // It is not safe to access a NULL pointer with this operator.
                MADNESS_ASSERT(pointer_ != NULL);
                // It is not safe to access this pointer remotely.
                MADNESS_ASSERT(is_local());
                return *pointer_;
            }

            /// Pointer arrow operator

            /// Access members of the pointer
            /// \return The local pointer
            /// \throw MadnessException When the pointer references a remote
            /// address, or when the pointer is NULL.
            pointer operator->() const {
                // It is not safe to access a NULL pointer with this operator.
                MADNESS_ASSERT(pointer_ != NULL);
                // It is not safe to access this pointer remotely.
                MADNESS_ASSERT(is_local());
                return pointer_;
            }

            /// Boolean conversion operator

            /// \return \c true when the pointer is non-null and \c false otherwise
            /// \throw nothing
            operator bool () const { return pointer_; }

            /// Boolean conversion operator

            /// \return \c true when the pointer is null and \c false otherwise
            /// \throw nothing
            bool operator ! () const { return !pointer_; }

            /// Equality comparison operator

            /// \tparam U Another pointer type
            /// \return \c true when the pointers refer to the same address from
            /// the same node in the same world.
            /// \throw nothing
            template <typename U>
            bool operator==(const WorldPtr<U>& other) const {
                return (pointer_ == other.pointer_) && (rank_ == other.rank_)
                    && (worldid_ == other.worldid_);
            }

            /// Inequality comparison operator

            /// \tparam U Another pointer type
            /// \return \c true when the pointers refer to different addresses or
            /// different nodes or different worlds.
            /// \throw nothing
            template <typename U>
            bool operator!=(const WorldPtr<U>& other) const {
                return (pointer_ != other.pointer_) || (rank_ != other.rank_)
                    || (worldid_ != other.worldid_);
            }

            /// Less-than comparison operator

            /// This operator does a lexicographical comparison of world ID,
            /// rank, and pointer (in that order).
            /// \tparam U Another pointer type
            /// \return \c true when the lexicographical comparison of world ID,
            /// rank, and pointer is true, and false otherwise.
            /// \throw nothing
            template <typename U>
            bool operator<(const WorldPtr<U>& other) const {
                return (worldid_ < other.worldid_) ||
                        ((worldid_ == worldid_) && ((rank_ < other.rank_) ||
                        ((rank_ == other.rank_) && (pointer_ < other.pointer_))));
            }

            /// World accessor

            /// \return A pointer to the world (may be NULL)
            /// \throw MadnessException When the pointer world has not been set
            /// (i.e. When \c has_owner()==false)
            World& get_world() const {
                MADNESS_ASSERT(world_ != NULL);
                return *world_;
            }

            /// World ID accessor

            /// \return The world ID of the world that the pointer belongs to.
            /// \throw nothing
            worldidT get_worldid() const { return worldid_ - 1; }

            /// Rank accessor

            /// If the pointer is not associated with
            /// \return The rank of the process that owns the pointer
            /// \throw nothing
            ProcessID owner() const { return rank_; }

            /// Swap the content of \c this with \c other

            /// \tparam U The other world pointer type
            /// \throw nothing
            /// \note \c U* must be implicitly convertible to T* type.
            template <typename U>
            void swap(WorldPtr<U>& other) {
                std::swap(world_, other.world_);
                std::swap(worldid_, other.worldid_);
                std::swap(rank_, other.rank_);
                std::swap(pointer_, other.pointer_);
            }

            /// Serialize/deserialize the world pointer.

            /// Serialize/deserialize the world pointer for remote communication
            /// or write to disk.
            /// \tparam Archive The archive object type.
            template <class Archive>
            inline void load_internal_(const Archive& ar) {
                ar & worldid_ & rank_ & archive::wrap_opaque(pointer_);
                world_ = (worldid_ != 0 ? detail::world_from_id(get_worldid()) : NULL);
            }

            template <class Archive>
            inline void store_internal_(const Archive& ar) const {
                ar & worldid_ & rank_ & archive::wrap_opaque(pointer_);
            }

            friend std::ostream& operator<<(std::ostream& out, const WorldPtr<T>& p) {
                out << "WorldPointer(ptr=" << p.pointer_ << ", rank=";
                if(p.rank_ >= 0)
                    out << p.rank_;
                else
                    out << "none";
                out << ", worldid=";
                if(p.worldid_ != 0)
                    out << (p.worldid_ - 1);
                else
                    out << "none";
                out << ")";
                return out;
            }
        }; // class WorldPtr

        /// Swap the content of \c l with \c r

        /// \tparam U The other world pointer type
        /// \throw nothing
        /// \note \c U* must be implicitly convertible to T* type.
        template <typename T, typename U>
        void swap(WorldPtr<T>& l, WorldPtr<U>& r) {
            l.swap(r);
        }

        template <typename T, typename U>
        bool operator>(const WorldPtr<T>& left, const WorldPtr<U>& right) {
            return right < left;
        }

        template <typename T, typename U>
        bool operator<=(const WorldPtr<T>& left, const WorldPtr<U>& right) {
            return !(right < left);
        }

        template <typename T, typename U>
        bool operator>=(const WorldPtr<T>& left, const WorldPtr<U>& right) {
            return !(left < right);
        }



    } // namespace detail

    namespace archive {
        template <typename, typename>
        struct ArchiveLoadImpl;

        template <typename, typename>
        struct ArchiveStoreImpl;

        template <typename Archive, typename T>
        struct ArchiveLoadImpl<Archive, detail::WorldPtr<T> > {
            static inline void load(const Archive& ar, detail::WorldPtr<T>& p) {
                p.load_internal_(ar);
            }
        };

        template <typename Archive, typename T>
        struct ArchiveStoreImpl<Archive, detail::WorldPtr<T> > {
            static inline void store(const Archive& ar, const detail::WorldPtr<T>& p) {
                p.store_internal_(ar);
            }
        };

    } // namespace archive
} // namespace madness
#endif // MADNESS_WORLD_WORLDPTR_H__INCLUDED
