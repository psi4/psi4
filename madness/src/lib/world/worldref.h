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


  $Id: worldref.h 2379 2011-06-17 03:01:22Z justus.c79@gmail.com $
*/


#ifndef MADNESS_WORLD_WORLDREF_H__INCLUDED
#define MADNESS_WORLD_WORLDREF_H__INCLUDED

/// \file worldref.h
/// \brief Implements RemoteReference which is for internal use

#include <world/atomicint.h>    // for AtomicInt
#include <world/sharedptr.h>    // for shared_ptr
#include <world/worldtypes.h>   // for ProcessID
#include <world/archive.h>      // for wrap_opaque
#include <world/worldam.h>      // for new_am_arg
#include <world/worldptr.h>     // for WorldPtr
#include <world/worldhashmap.h> // for ConcurrentHashMap
#include <iosfwd>               // for std::ostream

//#define MADNESS_REMOTE_REFERENCE_DEBUG
#ifdef MADNESS_REMOTE_REFERENCE_DEBUG
#include <world/print.h>        // for print
#endif

namespace madness {

    class World;
    template <typename T> class RemoteReference;
//
//    template <typename T>
//    std::ostream& operator<<(std::ostream& s, const RemoteReference<T>& ref);

    namespace archive {
        template <typename, typename>
        struct ArchiveLoadImpl;

        template <typename, typename>
        struct ArchiveStoreImpl;
    }

    namespace detail {

        /// Base class for remote counter implementation objects

        /// This class only holds an atomic counter. The use counter tracks
        /// local copies of the counter an references that have been copied as
        /// part of the communication process. This class also provides a
        /// mechanism for hiding the pointer type.
        /// \note The actual counter manipulation is handled by RemoteCounter.
        /// This class only provides the counter interface.
        /// \note This class is considered an implementation detail and may
        /// change at any time. You should not use this class directly.
        class RemoteCounterBase {
        private:
            madness::AtomicInt count_;          ///< reference count

            // Copy not allowed
            RemoteCounterBase(const RemoteCounterBase&);
            RemoteCounterBase& operator=(const RemoteCounterBase&);

        public:

            RemoteCounterBase() { count_ = 1; }
            virtual ~RemoteCounterBase() { }

            /// Counter key accessor

            /// The key is the pointer for which the remote counter is counting
            /// references.
            /// \return The pointer that is being counted.
            virtual void* key() const = 0;

            /// Remote and local counter accessor

            /// The use counter tracks local copies of the counter an references
            /// that have been copied as part othe communication process
            long use_count() const { return count_; }

            /// Increment the reference count

            /// The reference count should be incremented when a local copy of
            /// the counter is created or the when the counter is serialized as
            /// part of communication.
            /// \throw nothing
            void add_ref() {
#ifdef MADNESS_REMOTE_REFERENCE_DEBUG
                long c = count_++;
                print(">>> RemoteCounterBase(", this->key(), ") +ref count=", c + 1);
#else
                count_++;
#endif
            }

            /// Decrement the reference count

            /// \return true if the reference count has dropped to zero
            /// \throw nothing
            bool release() {
#ifdef MADNESS_REMOTE_REFERENCE_DEBUG
                long c = count_;
                print(">>> RemoteCounterBase(", this->key(), ") -ref count=", c - 1);
#endif
                return count_.dec_and_test();
            }
        }; // class RemoteCounterBase

        /// Remote counter implementation object.

        /// This class stores a shared pointer in memory to ensure that the
        /// referenced object is valid as long as there are outstanding remote
        /// references.
        /// \tparam T The type of the referenced shared_ptr object.
        /// \note This class is considered an implementation detail and may
        /// change at any time. You should not use this class directly.
        template <typename T>
        class RemoteCounterImpl : public RemoteCounterBase {
        private:
            // At some point this should probably be changed to a quick allocator
            // When that happens, also uncomment the new and delete operators
//            typedef std::allocator<RemoteCounterImpl<T> > A;

            // Keep a copy of the shared pointer to make sure it stays in memory
            // while we have outstanding remote references to it.
            std::shared_ptr<T> pointer_; ///< pointer that is remotely referenced

        public:
            explicit RemoteCounterImpl(const std::shared_ptr<T>& p) :
                RemoteCounterBase(), pointer_(p)
            { }

            virtual ~RemoteCounterImpl() { }

            /// Counter key accessor

            /// The key is the pointer for which the remote counter is counting
            /// references.
            /// \return The pointer that is being counted.
            virtual void* key() const { return static_cast<void*>(pointer_.get()); }

//            void* operator new(std::size_t) {
//                return A().allocate(1);
//            }
//
//            void operator delete(void * p) {
//                A().deallocate(static_cast<RemoteCounterImpl<T> *>(p), 1);
//            }
        }; // class RemoteCounterImpl

        /// Remote reference counter

        /// Automatically counts local and remote references to an object. The
        /// reference count is incremented when the object is copied locally or
        /// serialized as part of communication.
        class RemoteCounter {
        private:
            typedef RemoteCounterBase implT;
            typedef ConcurrentHashMap<void*, WorldPtr<implT> > pimpl_mapT;

            static pimpl_mapT pimpl_map_;   ///< A map of currently registered
                                            ///< implementation objects. The key is
                                            ///< it's referenced pointer.

            /// Pointer to the shared counter implementation object
            mutable WorldPtr<implT> pimpl_;

            /// Clean-up the implementation object

            /// Here we check that the pimpl has been initialized, and if so, we
            /// release the current reference. If the count drops to zero, then
            /// this is the last reference to the pimpl and it should be deleted.
            void destroy();

            /// Register a local shared pointer

            /// This function will first search the local pointer register for
            /// the shared pointer \c p. If found the pimpl for that pointer
            /// will be returned. Otherwise a new pimpl will be created and
            /// returned.
            /// \tparam T The shared pointer type to register
            /// \param w The world where the pointer lives
            /// \param p The shared pointer to register
            /// \return A world pointer to the pimpl
            /// \throw std::bad_alloc If pimpl allocation fails.
            /// \throw madness::MadnessException If pointer cannot be inserted
            /// into the pointer registration map.
            template <typename T>
            static WorldPtr<implT> register_ptr_(World& w, const std::shared_ptr<T>& p) {
                // Check for a null pointer
                if(p.get() == NULL)
                    return WorldPtr<implT>(w, NULL);

                pimpl_mapT::accessor acc;
                // Pointer is local and non-null
                if(pimpl_map_.insert(acc,static_cast<void*>(p.get()))) {
                    // The pointer is not registered so we need to make a
                    // new pimpl.
                    implT* pimpl = new RemoteCounterImpl<T>(p);

                    try{
                        acc->second = WorldPtr<implT>(w, pimpl);
                    } catch(...) {
                        delete pimpl;
                        throw;
                    }

#ifdef MADNESS_REMOTE_REFERENCE_DEBUG
                        print(">>> RemoteCounter::register_ptr_(new): key=", p.get(), ", pimpl=", acc->second);
#endif
                } else {
                    // The pointer is already registered, so we just need
                    // increment the counter.
#ifdef MADNESS_REMOTE_REFERENCE_DEBUG
                    print(">>> RemoteCounter::register_ptr_(existing): key=", acc->second->key(), ", pimpl=", acc->second);
#endif
                    acc->second->add_ref();
                }

                return acc->second;
            }

            /// Unregister a local shared pointer reference

            /// \param key The key of the \c RemoteReference object to be unregistered.
            /// \throw MadnessException If \c key is not found in the pointer map.
            static void unregister_ptr_(void* key) {
                std::size_t ereased = pimpl_map_.erase(key);
                MADNESS_ASSERT(ereased > 0);
            }

            RemoteCounter(const WorldPtr<implT>& p) :
                pimpl_(p)
            { }

        public:

            RemoteCounter() : pimpl_() { }

            RemoteCounter(const RemoteCounter& other) :
                pimpl_(other.pimpl_)
            {
                if(pimpl_ && pimpl_.is_local())
                    pimpl_->add_ref();
            }

            template <typename T>
            explicit RemoteCounter(World& w, const std::shared_ptr<T>& p) :
                pimpl_(register_ptr_(w, p))
            { }

            ~RemoteCounter() { destroy(); }

            RemoteCounter& operator=(const RemoteCounter& other);

            /// Counter accessor

            /// \return The number of local and remote references
            /// \throw none
            long use_count() const { return (pimpl_.is_local() ? pimpl_->use_count() : 0); }
            bool unique() const { return use_count() == 1; }
            bool empty() const { return ! pimpl_; }

            bool is_local() const { return pimpl_.is_local(); }
            bool has_owner() const { return pimpl_.has_owner(); }
            ProcessID owner() const { return pimpl_.owner(); }

            WorldPtr<implT>::worldidT
            get_worldid() const { return pimpl_.get_worldid(); }
            World& get_world() const { return pimpl_.get_world(); }
            void swap(RemoteCounter& other) {
                madness::detail::swap(pimpl_, other.pimpl_);
            }

        private:

            template <typename, typename>
            friend struct archive::ArchiveLoadImpl;

            template <typename, typename>
            friend struct archive::ArchiveStoreImpl;

            template <typename Archive>
            void load_(const Archive& ar) {
                WorldPtr<implT> p;
                ar & p;
                RemoteCounter(p).swap(*this);

#ifdef MADNESS_REMOTE_REFERENCE_DEBUG
                print(">>> RemoteCounter::load: pimpl=", pimpl_);
#endif
            }

            template <typename Archive>
            void store_(const Archive& ar) const {
                ar & pimpl_;

                if(! ar.count_only()) {
#ifdef MADNESS_REMOTE_REFERENCE_DEBUG
                    print(">>> RemoteCounter::store: pimpl=", pimpl_);
#endif
                    if(pimpl_.is_local())
                        pimpl_->add_ref();
                    else
                        pimpl_ = WorldPtr<implT>();
                }
            }


        }; // class RemoteCounter

        inline void swap(RemoteCounter& l, RemoteCounter& r) { l.swap(r); }

        std::ostream& operator<<(std::ostream& out, const RemoteCounter& counter);

    } // namespace detail


    /// Simple structure used to manage references/pointers to remote instances

    /// This class was intended only for internal use and is still rather
    /// poorly thought through, however, it seems to fill a wider need.
    /// \note Do not serialize via wrap_opaque().
    /// \note Ownership of a reference is transfered when serialized on a remote
    /// node. You should not attempt to send a remote reference to more than one
    /// node except from the owning node. If you do serialize more than once,
    /// this will cause an invalid memory access on the owning node.
    /// \note !!! It is YOUR RESPONSIBILITY to release the reference count. This
    /// can be done by sending the remote reference back to the owner or by
    /// calling reset(). If this is not done, you will have a memory leak.
    template <typename T>
    class RemoteReference {
    public:
        typedef typename detail::ptr_traits<T>::reference referenceT;
        typedef T* pointerT;

    private:
        mutable pointerT pointer_;      ///< World pointer
        detail::RemoteCounter counter_; ///< Remote reference counter

        // This is for RemoteReferences of other types, so they can still access
        // private members.
        template <typename>
        friend class RemoteReference;

        // Handles reset of a remote reference from another node.
        static void reset_handler(const AmArg& arg) {
            RemoteReference<T> r;
            arg & r;
            // r resets on scope exit.
        }

    public:

        /// Makes a non-shared (no reference count) null pointer
        RemoteReference() :
            pointer_(), counter_() {};

        /// Construct a remote reference to p.

        /// \param w The world that \c p belongs to.
        /// \param p The \c shared_ptr that is to be referenced.
        /// \note \c p must be locally addressable pointer
        RemoteReference(World& w, const std::shared_ptr<T>& p) :
            pointer_(p.get()), counter_(w, p)
        { }

        /// Copy constructor

        /// \param other The reference to be copied
        RemoteReference(const RemoteReference<T>& other) :
            pointer_(other.pointer_), counter_(other.counter_)
        { }

        /// Copy conversion constructor

        /// \tparam U The remote reference type to be copied
        /// \param other The reference to be copied
        /// \note \c U* must be implicitly convertible to \c T*
        template <typename U>
        RemoteReference(const RemoteReference<U>& other) :
            pointer_(other.pointer_), counter_(other.counter_)
        { }

        ~RemoteReference() { }

        /// Copy conversion assignment operator

        /// \param other The reference to be copied
        RemoteReference<T>& operator=(const RemoteReference<T>& other) {
            RemoteReference<T>(other).swap(*this);
            return *this;
        }

        /// Copy conversion assignment operator

        /// \tparam U The remote reference type to be copied
        /// \param other The reference to be copied
        /// \note \c U* must be implicitly convertible to \c T*
        template <typename U>
        RemoteReference<T>& operator=(const RemoteReference<U>& other) {
            RemoteReference<T>(other).swap(*this);
            return *this;
        }


        /// Release this reference

        /// This function will clear the reference and leave it in the default
        /// constructed state. If the reference is non-local, then a message is
        /// sent to the reference owner that releases the reference.
        /// \warning Only call this function for non-local references when it
        /// will not otherwise be returned to the reference owner as part of a
        /// message.
        void reset() {
            if((! (counter_.is_local())) && counter_.has_owner())
                get_world().am.send(owner(), RemoteReference<T>::reset_handler, new_am_arg(*this));
            else
                RemoteReference<T>().swap(*this);
        }

        /// Boolean conversion operator

        /// \return true when the reference is initialized to a non zero value
        /// or uninitialized, otherwise false
        operator bool() const {
            return ! counter_.empty();
        }

        /// Reference pointer accessor

        /// \return The referenced pointer
        /// \throw MadnessException If the pointer is not local
        pointerT get() const {
            MADNESS_ASSERT(counter_.is_local());
            return pointer_;
        }

        /// Reference object accessor

        /// \return A reference to the referenced object
        /// \throw MadnessException If the pointer is uninitialized
        /// \throw MadnessException If the pointer is not local
        referenceT operator*() const {
            MADNESS_ASSERT(pointer_ != NULL);
            MADNESS_ASSERT(counter_.is_local());
            return *pointer_;
        }

        /// Reference object pointer accessor

        /// \return A pointer to the referenced object
        /// \throw MadnessException If the pointer is uninitialized
        /// \throw MadnessException If the pointer is not local
        pointerT operator->() const {
            MADNESS_ASSERT(pointer_ != NULL);
            MADNESS_ASSERT(counter_.is_local());
            return pointer_;
        }

        /// Reference count accessor

        /// \return The total number of local and remote references.
        /// \throw nothing
        long use_count() const { return counter_.use_count(); }

        /// Get uniqueness

        /// \return True when the use count is equal to exactly 1.
        /// \throw nothing
        bool unique() const { return counter_.unique(); }

        /// Swap references

        /// Exchange the value of this \c RemoteReference with \c other
        /// \c RemoteReference
        /// \tparam U The type of the other remote reference.
        /// \note U* must be implicitly convertible to T*.
        template <typename U>
        void swap(RemoteReference<U>& other) {
            std::swap(pointer_, other.pointer_);
            madness::detail::swap(counter_, other.counter_);
        }

        /// Locally owned reference

        /// \return true if owner is equal to the current rank of the owning
        /// world, otherwise false
        /// \throw nothing
        inline bool is_local() const { return counter_.is_local(); }

        /// Reference owner accessor

        /// \return rank of owning process, or -1 if not initialized
        /// \throw nothing
        inline ProcessID owner() const { return counter_.owner(); }

        /// Owning world accessor

        /// \return A reference to the world that owns the pointer
        /// \throw MadnessException If the reference is uninitialized
        World& get_world() const { return counter_.get_world(); }

        /// Serialize the remote reference

        /// \tparam Archive The serialization archive type
        /// \param ar The serialization archive object.
        template <typename Archive>
        void serialize(const Archive& ar) const {
            // All of the interesting stuff happens in the counter serialization.
            ar & archive::wrap_opaque(pointer_) & counter_;
        }

    public:

        /// Add the remote reference to the given \c std::ostream, \c out.

        /// \param out The output stream to add \c ref to.
        /// \param ref The remote reference to add to the out stream
        friend std::ostream& operator<<(std::ostream& out, const RemoteReference<T>& ref) {
            out << "RemoteReference( pointer=" << ref.pointer_ << " counter=" << ref.counter_ << ")";
            return out;
        }
    }; // class RemoteReference


    /// Swap the two remote references

    /// \param l The left reference to be swapped with \c r
    /// \param r The right reference to be swapped with \c l
    /// \note T* must be implicitly convertible to U* and vis versa.
    template <typename T, typename U>
    void swap(RemoteReference<T>& l, RemoteReference<U>& r) {
        l.swap(r);
    }

    namespace archive {

        // This function is not allowed. Therefore it is not implemented so that
        // a compiler error is generated it it is called. This still does not
        // prevent remote references from being wrapped as part of another object.
        template <typename T>
        archive_array<unsigned char> wrap_opaque(const RemoteReference<T>& t);

        // Remote counter serialization

        template <typename Archive>
        struct ArchiveLoadImpl<Archive, detail::RemoteCounter > {
            static inline void load(const Archive& ar, detail::RemoteCounter& c) {
                c.load_(ar);
            }
        };

        template <typename Archive>
        struct ArchiveStoreImpl<Archive, detail::RemoteCounter > {
            static inline void store(const Archive& ar, const detail::RemoteCounter& c) {
                c.store_(ar);
            }
        };

    } // namespace archive
} // namespace madness

#endif // MADNESS_WORLD_WORLDREF_H__INCLUDED
