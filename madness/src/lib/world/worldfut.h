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


  $Id: worldfut.h 2294 2011-05-02 18:55:18Z justus.c79@gmail.com $
*/


#ifndef MADNESS_WORLD_WORLDFUT_H__INCLUDED
#define MADNESS_WORLD_WORLDFUT_H__INCLUDED

/// \file worldfut.h
/// \brief Implements Future
/// \ingroup futures

#include <vector>
#include <world/nodefaults.h>
#include <world/worlddep.h>
#include <world/array.h>
#include <world/sharedptr.h>
#include <world/worldref.h>
#include <world/typestuff.h>
#include <world/world.h>

namespace madness {

    //extern SharedCounter future_count; // For tracking memory leak


/*

\section gotchas Gotchas

\subsection futures Futures and STL vectors (e.g., \c vectors<Future<int>> )

 This to be turned back into documentation eventually

A common misconception is that STL containers initialize their
contents by \invoking the default constructor of each item in
the container since we are told that the items must be default
constructable.  But this is \em incorrect.  The items are initialized
by invoking the copy constructor for each element on a \em single
object made with the default constructor.   For futures this
is a very bad problem.  For instance,
\code
   vector< Future<double> > v(3);
\endcode
is equivalent to the following with an array of three elements
\code
   Future<double> junk;
   Future<double> v[3] = {junk,junk,junk};
\endcode
Since the Future copy constructor is by necessity shallow, each
element of \c v ends up referring to the future implementation that
underlies \c junk.  When you assign to an element of \c v, you'll also
be assigning to junk.  But since futures are single assignment
variables, you can only do that once.  Hence, when you assign a
second element of \c v you'll get a runtime exception.

The fix (other than using arrays) is to initialize STL vectors and
other containers from the special element returned by
\c Future<T>::default_initializer() which if passed into the copy
constructor will cause it to behave just like the default contructor.
Thus, the following code is what you actually need to use an STL
vector of futures
\code
   vector< Future<double> > v(3,Future<double>::default_initializer());
\endcode
which sucks, so we provide the factory function
\code
   template <typename T>
   vector< Future<T> > future_vector_factory(std::size_t n);
\endcode
which enables you to write
\code
   vector< Future<double> > v = future_vector_factory<double>(3);
\endcode
which merely blows instead of sucking.
*/

    template <typename T> class Future;

    /// Boost-type-trait-like testing of if a type is a future

    /// \ingroup futures
    template <typename T>
    struct is_future : public std::false_type { };

    /// Boost-type-trait-like testing of if a type is a future

    /// \ingroup futures
    template <typename T>
    struct is_future< Future<T> > : public std::true_type { };

    /// Boost-type-trait-like mapping of Future<T> to T

    /// \ingroup futures
    template <typename T>
    struct remove_future {
        typedef T type;
    };

    /// Boost-type-trait-like mapping of Future<T> to T

    /// \ingroup futures
    template <typename T>
    struct remove_future< Future<T> > {
        typedef T type;
    };

    /// Macro to determine type of future (by removing wrapping future template)

    /// \ingroup futures
#define REMFUTURE(T) typename remove_future< T >::type

    /// Human readable printing of future to stream

    /// \ingroup futures
    template <typename T>
    std::ostream& operator<<(std::ostream& out, const Future<T>& f);


    /// Implements the functionality of Futures

    /// \ingroup futures
    template <typename T>
    class FutureImpl : private Spinlock {
        friend class Future<T>;
        friend std::ostream& operator<< <T>(std::ostream& out, const Future<T>& f);

    private:
        static const int MAXCALLBACKS = 4;
        typedef Stack<CallbackInterface*,MAXCALLBACKS> callbackT;
        typedef Stack<std::shared_ptr< FutureImpl<T> >,MAXCALLBACKS> assignmentT;
        volatile callbackT callbacks;
        volatile assignmentT assignments;
        volatile bool assigned;
        World * world;
        RemoteReference< FutureImpl<T> > remote_ref;
        volatile T t;

        /// Private: AM handler for remote set operations
        static void set_handler(const AmArg& arg) {
            RemoteReference< FutureImpl<T> > ref;
            T t;
            arg & ref & t;
            FutureImpl<T>* f = ref.get();
            f->set(t);
            ref.reset();
        }


        /// Private:  invoked locally by set routine after assignment
        inline void set_assigned() {
            // ASSUME THAT WE HAVE THE MUTEX WHILE IN HERE (so that
            // destructor is not invoked by another thread as a result
            // of our actions)
            MADNESS_ASSERT(!assigned);
            assigned = true;

            callbackT& cb = const_cast<callbackT&>(callbacks);
            assignmentT& as = const_cast<assignmentT&>(assignments);

//             // if taking copy not reference so that destructor can still check for uninvoked cbs
//             const_cast<callbackT&>(callbacks).clear();
//             const_cast<assignmentT&>(assignments).clear();

            while (as.size()) {
                std::shared_ptr< FutureImpl<T> >& p = as.pop();
                MADNESS_ASSERT(p);
                p->set(const_cast<T&>(t));
                p.reset();
            }
            while (cb.size()) {
                CallbackInterface* p = cb.pop();
                MADNESS_ASSERT(p);
                p->notify();
            }
        }

        inline void add_to_assignments(const std::shared_ptr< FutureImpl<T> >& f) {
            // ASSUME lock is already acquired by Future<T>::set()
            if (assigned) {
                f->set(const_cast<T&>(t));
            }
            else {
                assignmentT& nvas = const_cast<assignmentT&>(assignments);
                nvas.push(f);
            }
        }


    public:

        // Local unassigned value
        FutureImpl()
                : callbacks()
                , assignments()
                , assigned(false)
                , world(0)
                , remote_ref()
                , t() {
            //print("FUTCON(a)",(void*) this);

            //future_count.inc();
        }


        // Local assigned value
        FutureImpl(const T& t)
                : callbacks()
                , assignments()
                , assigned(false)
                , world(0)
                , remote_ref()
                , t(t) {
            //print("FUTCON(b)",(void*) this);
            set_assigned();
            //future_count.inc();
        }


        // Wrapper for a remote future
        FutureImpl(const RemoteReference< FutureImpl<T> >& remote_ref)
                : callbacks()
                , assignments()
                , assigned(false)
                , world(& remote_ref.get_world())
                , remote_ref(remote_ref)
                , t() {
            //print("FUTCON(c)",(void*) this);

            //future_count.inc();
        }


        // Returns true if the value has been assigned
        inline bool probe() const {
            return assigned;
        }


        // Registers a function to be invoked when future is assigned

        // Callbacks are invoked in the order registered.  If the
        // future is already assigned the callback is immediately
        // invoked.
        inline void register_callback(CallbackInterface* callback) {
            ScopedMutex<Spinlock> fred(this);
            if (assigned) callback->notify();
            else const_cast<callbackT&>(callbacks).push(callback);
        }


        // Sets the value of the future (assignment)
        void set(const T& value) {
            ScopedMutex<Spinlock> fred(this);
            if (world) {
                if (remote_ref.owner() == world->rank()) {
                    remote_ref.get()->set(value);
                    set_assigned();
                    remote_ref.reset();
                }
                else {
                    const ProcessID owner = remote_ref.owner();
                    world->am.send(owner,
                                   FutureImpl<T>::set_handler,
                                   new_am_arg(remote_ref, value));
                    set_assigned();
                }
            }
            else {
                const_cast<T&>(t) = value;
                set_assigned();
            }
        }


        // Gets/forces the value, waiting if necessary (error if not local)
        T& get() {
            MADNESS_ASSERT(!world);  // Only for local futures
            World::await(bind_nullary_mem_fun(this,&FutureImpl<T>::probe));
            return *const_cast<T*>(&t);
        }

        // Gets/forces the value, waiting if necessary (error if not local)
        const T& get() const {
            MADNESS_ASSERT(!world);  // Only for local futures
            World::await(bind_nullary_mem_fun(this,&FutureImpl<T>::probe));
            return *const_cast<const T*>(&t);
        }

        T& operator->() {
            return get();
        }

        const T& operator->() const {
            return get();
        }

        bool is_local() const {
            return world == 0;
        }

        bool replace_with(FutureImpl<T>* f) {
            MADNESS_EXCEPTION("IS THIS WORKING? maybe now we have the mutex", 0);
//            ScopedMutex<Spinlock> fred(this);
//             MADNESS_ASSERT(!world); // was return false;
//             MADNESS_ASSERT(!assigned || f->assigned);
//             if (f->world) {
//                 world = f->world;
//                 remote_ref = f->remote_ref;
//                 f->world = 0;
//             }
//             while(f->callbacks.size()) callbacks.push(f->callbacks.pop());
//             while(f->assignments.size()) assignments.push(f->assignments.pop());
            return true;
        }

        virtual ~FutureImpl() {

            //future_count.dec_and_test();

            ScopedMutex<Spinlock> fred(this);
            //print("FUTDEL",(void*) this);
//             if (!assigned && world) {
//                 print("Future: unassigned remote future being destroyed?");
//                 //remote_ref.dec();
//                 abort();
//             }

            if (const_cast<callbackT&>(callbacks).size()) {
                print("Future: uninvoked callbacks being destroyed?", assigned);
                abort();
            }
            if (const_cast<assignmentT&>(assignments).size()) {
                print("Future: uninvoked assignment being destroyed?", assigned);
                abort();
            }
        }
    };


    /// A future is a possibly yet unevaluated value

    /// \ingroup futures
    /// Uses delegation to FutureImpl to provide desired
    /// copy/assignment semantics as well as safe reference counting
    /// for remote futures.
    ///
    /// Since we are using Futures a lot to store local values coming
    /// from containers and inside task wrappers for messages, we
    /// included in this class a value.  If a future is assigned
    /// before a copy/remote-reference is taken, the shared ptr is
    /// never made.  The point of this is to eliminate the two mallocs
    /// that must be peformed for every new shared_ptr.
    template <typename T>
    class Future {

        friend std::ostream& operator<< <T>(std::ostream& out, const Future<T>& f);

    private:
      // If f==0 ... can only happen in Future(value) ... i.e., the future is constructed
      // as assigned to value ... to optimize away the new(futureimpl) we instead set f=0
      // and copy the value

      // If f==nonzero ... there is an underlying future that may/maynot be assigned

        mutable std::shared_ptr< FutureImpl<T> > f;
        T value;
        const bool is_the_default_initializer;

        class dddd {};
        explicit Future(const dddd&)
                : f()
                , value()
                , is_the_default_initializer(true)
        {
        }

    public:
        typedef RemoteReference< FutureImpl<T> > remote_refT;

        /// Makes an unassigned future
        Future()
                : f(new FutureImpl<T>())
                , value()
                , is_the_default_initializer(false)
        {
        }

        /// Makes an assigned future
        explicit Future(const T& t)
                : f()
                , value(t)
                , is_the_default_initializer(false)
        {
        }


        /// Makes a future wrapping a remote reference
        explicit Future(const remote_refT& remote_ref)
                : f(new FutureImpl<T>(remote_ref))
                , value()
                , is_the_default_initializer(false)
        {
        }


        /// Copy constructor is shallow
        Future(const Future<T>& other)
	        : f(other.f)
                , value(other.value)
                , is_the_default_initializer(false)
        {
            if (other.is_the_default_initializer) {
                f.reset(new FutureImpl<T>());
            }
	    // Otherwise nothing to do ... done in member initializers
        }


        /// See Gotchas on the documentation mainpage about why this exists and how to use it.
        static const Future<T> default_initializer() {
            return Future<T>(dddd());
        }


        /// Assignment future = future makes a shallow copy just like copy constructor
        Future<T>& operator=(const Future<T>& other) {
            if (this != &other) {
                MADNESS_ASSERT(!probe());
                if (other.f) {
                    f = other.f;
                } else {
                    set(other);
                }
            }
            return *this;
        }

        /// A.set(B) where A & B are futures ensures A has/will-have the same value as B.

        /// An exception is thrown if A is already assigned since a
        /// Future is a single assignment variable.  We don't yet
        /// track multiple assignments from unassigned futures.
        ///
        /// If B is already assigned, this is the same as A.set(B.get())
        /// which sets A to the value of B.
        ///
        /// If B has not yet been assigned, the behavior is to ensure
        /// that when B is assigned that both A and B will be assigned
        /// and have the same value (though they may/may not refer to
        /// the same underlying copy of the data and indeed may even
        /// be in different processes).
        void set(const Future<T>& other) {
            if (this != &other) {
                MADNESS_ASSERT(!probe());
                if (other.probe()) {
                    set(other.get());     // The easy case
                }
                else {
                    // Assignment is supposed to happen just once so
                    // safe to assume that this is not being messed
                    // with ... also other might invoke the assignment
                    // callback since it could have been assigned
                    // between the test above and now (and this does
                    // happen)
                    other.f->lock();     // BEGIN CRITICAL SECTION
                    other.f->add_to_assignments(f); // Recheck of assigned is performed in here
                    other.f->unlock(); // END CRITICAL SECTION
                }
            }
        }

        /// Assigns the value ... it can only be set ONCE.
        inline void set(const T& value) {
	    MADNESS_ASSERT(f);
            f->set(value);
        }


        /// Gets the value, waiting if necessary (error if not a local future)
        inline T& get() {
	        if (f) {
                return f->get();
            } else {
                return value;
            }
        }

        /// Gets the value, waiting if necessary (error if not a local future)
        inline const T& get()  const {
            if (f) {
                return f->get();
            } else {
                return value;
            }
        }


        /// Returns true if the future has been assigned
        inline bool probe() const {
            if (f) return f->probe();
            else return true;
        }

        /// Same as get()
        inline operator T() {
            return get();
        }

        /// Same as get() const
        inline operator const T() const {
            return get();
        }


        /// Returns a structure used to pass references to another process.

        /// This is used for passing pointers/references to another
        /// process.  To make remote references completely safe, the
        /// RemoteReference increments the internal reference count of
        /// the Future.  The counter is decremented by either
        /// assigning to the remote Future or its destructor if it is
        /// never assigned.  The remote Future is ONLY useful for
        /// setting the future.  It will NOT be notified if the value
        /// is set elsewhere.
        ///
        /// If this is already a reference to a remote future, the
        /// actual remote reference is returned ... i.e., \em not a
        /// a reference to the local future.  Therefore, the local
        /// future will not be notified when the result is set
        /// (i.e., the communication is short circuited).
        inline remote_refT remote_ref(World& world) const {
            MADNESS_ASSERT(!probe());
            if (f->world)
                return f->remote_ref;
            else
                return RemoteReference< FutureImpl<T> >(world, f);
        }


        inline bool is_local() const {
            return ((!f) || f->is_local() || probe());
        }

        inline bool is_remote() const {
            return !is_local();
        }


        /// Registers an object to be called when future is assigned

        /// Callbacks are invoked in the order registered.  If the
        /// future is already assigned the callback is immediately
        /// invoked.
        inline void register_callback(CallbackInterface* callback) {
            if (probe())
                callback->notify();
            else
                f->register_callback(callback);
        }
    };


    /// A future of a future is forbidden (by private constructor)

    /// \ingroup futures
    template <typename T> class Future< Future<T> > {
        Future() {}
    };


    /// Specialization of FutureImpl<void> for internal convenience ... does nothing useful!

    /// \ingroup futures
    template <> class FutureImpl<void> {};

    /// Specialization of Future<void> for internal convenience ... does nothing useful!

    /// \ingroup futures
    template <> class Future<void> {
    public:
        typedef RemoteReference< FutureImpl<void> > remote_refT;

        remote_refT remote_ref(World&) const {
            return remote_refT();
        }

        Future() {}

        Future(const RemoteReference< FutureImpl<void> >&) {}

        Future(const Future<Void>&) {}

        inline void set(const Future<void>&) {}

        inline Future<void>& operator=(const Future<void>&) {
            return *this;
        }

        inline void set() {}
    };

    /// Specialization of FutureImpl<Void> for internal convenience ... does nothing useful!

    /// \ingroup futures
    template <> class FutureImpl<Void> {};

    /// Specialization of Future<Void> for internal convenience ... does nothing useful!

    /// \ingroup futures
    template <> class Future<Void> {
    public:
        typedef RemoteReference< FutureImpl<Void> > remote_refT;

        remote_refT remote_ref(World& /*world*/) const {
            return remote_refT();
        }

        Future() {}

        Future(const RemoteReference< FutureImpl<Void> >& /*ref*/) {}

        Future(const Future<void>& /*f*/) {}

        inline void set(const Future<Void>& /*f*/) {}

        inline Future<Void>& operator=(const Future<Void>& /*f*/) {
            return *this;
        }

        inline void set(const Void& /*f*/) {}
    };

    /// Specialization of Future for vector of Futures

    /// \ingroup futures
    /// Enables passing a vector of futures into a task and having
    /// the dependencies correctly tracked.  Does not directly
    /// support most operations that other futures do ... that is
    /// the responsiblility of the individual futures in the vector.
    template <typename T>
    class Future< std::vector< Future<T> > > : public DependencyInterface, private NO_DEFAULTS {
    private:
        typedef typename std::vector< Future<T> > vectorT;
        vectorT v;

    public:
        Future() : v() {}

        Future(const vectorT& v) : DependencyInterface(v.size()), v(v) {
            for (int i=0; i<(int)v.size(); ++i) {
                this->v[i].register_callback(this);
            }
        }
        vectorT& get() {
            return v;
        }
        const vectorT& get() const {
            return v;
        }
        operator vectorT& () {
            return get();
        }
        operator const vectorT& () const {
            return get();
        }
    };


    /// Probes a future for readiness, other types are always ready

    /// \ingroup futures
    template <typename T>
    ENABLE_IF(is_future<T>,bool) future_probe(const T& t) {
        return t.probe();
    }

    /// Probes a future for readiness, other types are always ready

    /// \ingroup futures
    template <typename T>
    DISABLE_IF(is_future<T>,bool) future_probe(const T& t) {
        return true;
    }


    /// Friendly I/O to streams for futures

    /// \ingroup futures
    template <typename T>
    std::ostream& operator<<(std::ostream& out, const Future<T>& f);

    template <>
    std::ostream& operator<<(std::ostream& out, const Future<void>& f);

    template <>
    std::ostream& operator<<(std::ostream& out, const Future<Void>& f);

    /// Factory for vectors of futures (see section Gotchas on the mainpage)

    /// \ingroup futures
    template <typename T>
    std::vector< Future<T> > future_vector_factory(std::size_t n) {
        return std::vector< Future<T> >(n, Future<T>::default_initializer());
    }


    namespace archive {
        /// Serialize an assigned future

        /// \ingroup futures
        template <class Archive, typename T>
        struct ArchiveStoreImpl< Archive, Future<T> > {
            static inline void store(const Archive& ar, const Future<T>& f) {
                MAD_ARCHIVE_DEBUG(std::cout << "serializing future" << std::endl);
                MADNESS_ASSERT(f.probe());
                ar & f.get();
            }
        };


        /// Deserialize a future into an unassigned future

        /// \ingroup futures
        template <class Archive, typename T>
        struct ArchiveLoadImpl< Archive, Future<T> > {
            static inline void load(const Archive& ar, Future<T>& f) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserializing future" << std::endl);
                MADNESS_ASSERT(!f.probe());
                T value;
                ar & value;
                f.set(value);
            }
        };


        /// Serialize an assigned future

        /// \ingroup futures
        template <class Archive>
        struct ArchiveStoreImpl< Archive, Future<void> > {
            static inline void store(const Archive& ar, const Future<void>& f) {
            }
        };


        /// Deserialize a future into an unassigned future

        /// \ingroup futures
        template <class Archive>
        struct ArchiveLoadImpl< Archive, Future<void> > {
            static inline void load(const Archive& ar, Future<void>& f) {
            }
        };

        /// Serialize an assigned future

        /// \ingroup futures
        template <class Archive>
        struct ArchiveStoreImpl< Archive, Future<Void> > {
            static inline void store(const Archive& ar, const Future<Void>& f) {
            }
        };


        /// Deserialize a future into an unassigned future

        /// \ingroup futures
        template <class Archive>
        struct ArchiveLoadImpl< Archive, Future<Void> > {
            static inline void load(const Archive& ar, Future<Void>& f) {
            }
        };
    }

    template <typename T>
    std::ostream& operator<<(std::ostream& out, const Future<T>& f) ;

    template <>
    std::ostream& operator<<(std::ostream& out, const Future<void>& f) ;

    template <>
    std::ostream& operator<<(std::ostream& out, const Future<Void>& f) ;

#ifdef WORLD_INSTANTIATE_STATIC_TEMPLATES
    template <typename T>
    std::ostream& operator<<(std::ostream& out, const Future<T>& f) {
        if (f.probe()) out << f.get();
        else if (f.is_remote()) out << f.f->remote_ref;
        else if (f.f) out << "<unassigned refcnt=" << f.f.use_count() << ">";
        else out << "<unassigned>";
        return out;
    }

#endif

}


#endif // MADNESS_WORLD_WORLDFUT_H__INCLUDED
