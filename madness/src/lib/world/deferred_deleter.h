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


  $Id: deferred_deleter.h 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/

#ifndef MADNESS_WORLD_DEFERRED_DELETER_H__INCLUDED
#define MADNESS_WORLD_DEFERRED_DELETER_H__INCLUDED

#include <world/sharedptr.h>
#include <world/typestuff.h>
#include <world/deferred_cleanup.h>

namespace madness {

    /// Deferred deleter for smart pointers.

    /// This deleter object places the shared pointer on the deferred deletion
    /// list of a world when the last reference to the pointer is destroyed.
    /// Once the pointer is placed in the deferred deletion list, it is
    /// destroyed by the world object at the next global fence of that world.
    /// You may pass any arbitrary deleter function pointer/functor to to the
    /// DeferredDeleter constructor to handle cleanup. If no deleter function
    /// pointer/functor is provided by the user, the pointer will be freed with
    /// the \c delete operator.
    /// \tparam ptrT The pointer type that will be deleted
    /// \tparam deleterT The deleter function pointer/functor type that will be use to
    /// cleanup the pointer [Default = void(*)(ptrT*) ].
    /// \note D type must be void(*)(T*) for function pointers or a functor type
    /// that includes a void D::operator() (ptrT*) function and have an accessible
    /// copy constructor and assignment operator.
    template <typename ptrT, typename deleterT = void(*)(ptrT*)>
    class DeferredDeleter {
    private:
        // Hold a shared pointer to the deferred deleter object to make sure it
        // is available for all deleter objects, even if the gop object is gone.
        std::shared_ptr<detail::DeferredCleanup> deferred_; ///< Deferred cleanup object.
        deleterT deleter_;  ///< The deleter function or functor that deletes
                            ///< the pointer

    public:
        /// Construct a default deleter for a function pointer
        template <typename D>
        static typename enable_if<std::is_same<D, void(*)(ptrT*)>, D>::type
        default_deleter() { return & detail::checked_delete<ptrT>; }

        /// Construct a default deleter for a functor
        template <typename D>
        static typename disable_if<std::is_same<D, void(*)(ptrT*)>, D>::type
        default_deleter() { return D(); }

        /// Constructs a deferred deleter object.

        /// The deleter function pointer \c d will be used to delete the pointer
        /// at a global fence of world \c w.
        /// \param w A reference to the world object, which will be responsible
        /// for pointer deletion.
        /// \param d A deleter function pointer/functor [Default = if \c D
        /// \c == \c void(*)(ptrT*) then \c d \c = \c &detail::checked_delete<ptrT>
        /// else \c d \c = \c D() ].
        DeferredDeleter(World& w, deleterT d = default_deleter<deleterT>()) :
            deferred_(detail::DeferredCleanup::get_deferred_cleanup(w)), deleter_(d)
        { }

        /// Copy constructor

        /// \param other The deleter object to be copied.
        DeferredDeleter(const DeferredDeleter<ptrT, deleterT>& other) :
            deferred_(other.deferred_), deleter_(other.deleter_)
        { }

        /// Copy assignment operator.

        /// \param other The deleter object to be copied.
        /// \return A reference to this object.
        DeferredDeleter<ptrT, deleterT>& operator=(const DeferredDeleter<ptrT, deleterT>& other) {
            deferred_ = other.deferred_;
            deleter_ = other.deleter_;
            return *this;
        }

        /// The deferred deletion function.

        /// This function is called when the last reference to the shared
        /// pointer is destroyed. It will place the pointer in the deferred
        /// cleanup list of world.
        void operator()(ptrT* p) const {
            std::shared_ptr<ptrT> temp(p, deleter_); // do not store this pointer
                                                // that is the job of deferred_
            deferred_->add(std::static_pointer_cast<void>(temp));
        }
    }; // class DeferredDeleter

    /// Make a defered deleter object for use with std::shared_ptr

    /// The following must be valid:
    /// \code
    /// T* p = new T(/* arguments here */);
    /// d(p);
    /// \endcode
    /// \tparam ptrT The pointer type that the deleter function will delete
    /// \tparam deleterT The deleter function type. This may be a functor or a
    /// function pointer type.
    /// \param w The world object that will be responsible for deleting the
    /// pointer at global sync points.
    /// \param d The deleter function pointer or functor
    /// \note A copy of the deleter is stored in the deferred deleter.
    template <typename ptrT, typename deleterT>
    inline DeferredDeleter<ptrT,deleterT> make_deferred_deleter(World& w, deleterT d) {
        return DeferredDeleter<ptrT, deleterT>(w, d);
    }

    /// Make a defered deleter object for use with std::shared_ptr

    /// The pointer passed to the defered deleter will be deleted with the
    /// \c delete operator.
    /// \tparam T The pointer type that the deleter function will delete
    /// pointer type.
    /// \param w The world object that will be responsible for deleting the
    /// pointer at global sync points.
    template <typename ptrT>
    inline DeferredDeleter<ptrT> make_deferred_deleter(World& w) {
        return DeferredDeleter<ptrT>(w);
    }

}  // namespace madness

#endif // MADNESS_WORLD_DEFERRED_DELETER_H__INCLUDED
