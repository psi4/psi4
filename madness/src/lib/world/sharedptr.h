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


 $Id: sharedptr.h 2250 2011-04-03 14:31:26Z justus.c79@gmail.com $
 */

#ifndef MADNESS_WORLD_SHAREDPTR_H__INCLUDED
#define MADNESS_WORLD_SHAREDPTR_H__INCLUDED

/// \file sharedptr.h
/// \brief Includes TR1 shared_ptr. If shared_ptr is in std::tr1 namespace, it
/// is imported into the std namespace. It also includes make_shared and
/// allocate_shared helper functions which are a part of the current C++0x
/// draft.

#include <madness_config.h>

// Select header that contains shared_ptr
#if defined(MADNESS_USE_MEMORY)
#include <memory>
#elif defined(MADNESS_USE_TR1_MEMORY)
#include <tr1/memory>
#elif defined(MADNESS_USE_BOOST_TR1_MEMORY_HPP)
#include <boost/tr1/memory.hpp>
#else
#error No acceptable memory include directive was found.
#endif // MEMORY

#ifndef MADNESS_BEGIN_NAMESPACE_TR1

#if defined(BOOST_TR1_MEMORY_INCLUDED) || defined(BOOST_TR1_MEMORY_HPP_INCLUDED)

// We are using boost
#define MADNESS_BEGIN_NAMESPACE_TR1 namespace boost {
#define MADNESS_END_NAMESPACE_TR1 } // namespace std

#elif defined(MADNESS_HAS_STD_TR1_SHARED_PTR)

// We are using TR1
#define MADNESS_BEGIN_NAMESPACE_TR1 namespace std { namespace tr1 {
#define MADNESS_END_NAMESPACE_TR1 } } // namespace std namespace tr1

#elif defined(MADNESS_HAS_STD_SHARED_PTR)

// We are using C++0x
#define MADNESS_BEGIN_NAMESPACE_TR1 namespace std {
#define MADNESS_END_NAMESPACE_TR1 } // namespace std

#else
// We do not know.
#error Unable to determine the correct namespace for TR1 fuctional.

#endif

#endif // MADNESS_BEGIN_NAMESPACE_TR1

#if defined(MADNESS_HAS_STD_TR1_SHARED_PTR) && !defined(MADNESS_HAS_STD_SHARED_PTR)
#define MADNESS_HAS_STD_SHARED_PTR 1
// shard_ptr is in std::tr1 but we want it in std namespace
namespace std {
    using ::std::tr1::bad_weak_ptr;
    using ::std::tr1::shared_ptr;
    using ::std::tr1::swap;
    using ::std::tr1::static_pointer_cast;
    using ::std::tr1::dynamic_pointer_cast;
    using ::std::tr1::const_pointer_cast;
    using ::std::tr1::get_deleter;
    using ::std::tr1::weak_ptr;
    using ::std::tr1::enable_shared_from_this;
}

#endif // defined(MADNESS_HAS_STD_TR1_SHARED_PTR) && !defined(MADNESS_HAS_STD_SHARED_PTR)
namespace madness {
    namespace detail {

        // These checked delete and deleters are copied from Boost.
        // They ensure that compilers issue warnings if T is an incomplete type.

        /// Checked pointer delete function

        /// This function ensures that the pointer is a complete type.
        /// \tparam T The pointer type (must be a complete type).
        /// \param p The pointer to be deleted.
        template <typename T>
        inline void checked_delete(T* p) {
            // intentionally complex - simplification causes regressions
            typedef char type_must_be_complete[sizeof(T) ? 1 : -1];
            (void)sizeof(type_must_be_complete);
            delete p;
        }

        /// Checked array pointer delete function

        /// This function ensures that the pointer is a complete type.
        /// \tparam T The pointer type (must be a complete type).
        /// \param a The array pointer to be deleted.
        template <typename T>
        inline void checked_array_delete(T* a) {
            typedef char type_must_be_complete[sizeof(T) ? 1 : -1];
            (void)sizeof(type_must_be_complete);
            delete[] a;
        }

        /// Function to free memory for a shared_ptr using free()

        /// Checks the pointer to make sure it is a complete type, you will get
        /// a compiler error if it is not.
        template <typename T>
        inline void checked_free(T* p) {
            typedef char type_must_be_complete[sizeof(T) ? 1 : -1];
            (void)sizeof(type_must_be_complete);
            free(p);
        }

        /// Use this function with shared_ptr to do nothing for the pointer cleanup
        template <typename T>
        inline void no_delete(T*) { }

        inline void no_delete(void*) { }

        /// Checked pointer delete functor

        /// This functor is used to delete a pointer. It ensures that the
        /// pointer is a complete type.
        /// \tparam T The pointer type (must be a complete type).
        template <typename T>
        struct CheckedDeleter {
            typedef void result_type;
            typedef T * argument_type;

            void operator()(T* p) const {
                checked_delete(p);
            }
        };

        /// Checked array pointer delete functor

        /// This functor is used to delete an array pointer. It ensures that the
        /// pointer is a complete type.
        /// \tparam T The pointer type (must be a complete type).
        template <typename T>
        struct CheckedArrayDeleter {
            typedef void result_type;
            typedef T * argument_type;

            void operator()(T* a) const {
                checked_array_delete(a);
            }
        };

        /// Deleter to free memory for a shared_ptr using free()

        /// Checks the pointer to make sure it is a complete type, you will get
        /// a compiler error if it is not.
        template <typename T>
        struct CheckedFree {
            typedef void result_type;
            typedef T * argument_type;

            void operator()(T* p) const {
                checked_free(p);
            }
        };

        /// Use this deleter with shared_ptr to do nothing for the pointer cleanup
        template <typename T>
        struct NoDeleter {
            typedef void result_type;
            typedef T * argument_type;

            void operator()(T*) const { }
        };

    } // namespace detail
} // namespace madness

// make_shared / allocate_shared
#ifndef MADNESS_HAS_STD_MAKE_SHARED
#define MADNESS_HAS_STD_MAKE_SHARED 1

#if defined(BOOST_TR1_MEMORY_INCLUDED) || defined(BOOST_TR1_MEMORY_HPP_INCLUDED)

// We are using Boost tr1 so we can use the Boost make_shared function
#include <boost/make_shared.hpp>
namespace std {
    using ::boost::make_shared;
    using ::boost::allocate_shared;
} // namespace std

#else

// We do not have make_shared/allocate_shared so we need to implement it here.
namespace std {

    /// Construct an shared pointer

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T());
    /// \endcode
    /// \tparam T The shared_ptr type
    /// \return A shared_ptr constructed with the given arguments
    template <class T>
    std::shared_ptr<T> make_shared() {
        return std::shared_ptr<T>(new T());
    }

    /// Construct an shared pointer with an allocator

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(), deleter(), a);
    /// \endcode
    /// where deleter is a deleter object supplied by this function.
    /// \tparam T The shared_ptr type
    /// \tparam A the allocator type
    /// \param a The allocator object used to allocate the shared_ptr
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A> std::shared_ptr<T> allocate_shared(A const & a) {
        return std::shared_ptr<T>(new T(), &madness::detail::checked_delete, a);
    }

    /// Construct an shared pointer

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(t1));
    /// \endcode
    /// \tparam T The shared_ptr type
    /// \tparam T1 pointer constructor argument 1 type
    /// \param t1 pointer constructor argument 1
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class T1>
    std::shared_ptr<T> make_shared(T1& t1) {
        return std::shared_ptr<T>(new T(t1));
    }

    template <class T, class T1>
    std::shared_ptr<T> make_shared(T1 const& t1) {
        return std::shared_ptr<T>(new T(t1));
    }

    /// Construct an shared pointer with an allocator

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(t1), deleter(), a);
    /// \endcode
    /// where deleter is a deleter object supplied by this function.
    /// \tparam T The shared_ptr type
    /// \tparam A the allocator type
    /// \tparam T1 pointer constructor argument 1 type
    /// \param a The allocator object used to allocate the shared_ptr
    /// \param t1 pointer constructor argument 1
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A, class T1>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1) {
        return std::shared_ptr<T>(new T(t1), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const& t1) {
        return std::shared_ptr<T>(new T(t1), &madness::detail::checked_delete, a);
    }

    /// Construct an shared pointer

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(t1, t2));
    /// \endcode
    /// \tparam T The shared_ptr type
    /// \tparam T1 pointer constructor argument 1 type
    /// \tparam T2 pointer constructor argument 2 type
    /// \param t1 pointer constructor argument 1
    /// \param t2 pointer constructor argument 2
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class T1, class T2>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2) {
        return std::shared_ptr<T>(new T(t1, t2));
    }

    template <class T, class T1, class T2>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2) {
        return std::shared_ptr<T>(new T(t1, t2));
    }

    template <class T, class T1, class T2>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2) {
        return std::shared_ptr<T>(new T(t1, t2));
    }

    template <class T, class T1, class T2>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2) {
        return std::shared_ptr<T>(new T(t1, t2));
    }

    /// Construct an shared pointer with an allocator

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(t1, t2), deleter(), a);
    /// \endcode
    /// where deleter is a deleter object supplied by this function.
    /// \tparam T The shared_ptr type
    /// \tparam A the allocator type
    /// \tparam T1 pointer constructor argument 1 type
    /// \tparam T2 pointer constructor argument 2 type
    /// \param a The allocator object used to allocate the shared_ptr
    /// \param t1 pointer constructor argument 1
    /// \param t2 pointer constructor argument 2
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A, class T1, class T2>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2) {
        return std::shared_ptr<T>(new T(t1, t2), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2) {
        return std::shared_ptr<T>(new T(t1, t2), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2) {
        return std::shared_ptr<T>(new T(t1, t2), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2) {
        return std::shared_ptr<T>(new T(t1, t2), &madness::detail::checked_delete, a);
    }

    /// Construct an shared pointer

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(t1, t2, t3));
    /// \endcode
    /// \tparam T The shared_ptr type
    /// \tparam T1 pointer constructor argument 1 type
    /// \tparam T2 pointer constructor argument 2 type
    /// \tparam T3 pointer constructor argument 3 type
    /// \param t1 pointer constructor argument 1
    /// \param t2 pointer constructor argument 2
    /// \param t3 pointer constructor argument 3
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class T1, class T2, class T3>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3& t3) {
        return std::shared_ptr<T>(new T(t1, t2, t3));
    }

    template <class T, class T1, class T2, class T3>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3& t3) {
        return std::shared_ptr<T>(new T(t1, t2, t3));
    }

    template <class T, class T1, class T2, class T3>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3& t3) {
        return std::shared_ptr<T>(new T(t1, t2, t3));
    }

    template <class T, class T1, class T2, class T3>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3& t3) {
        return std::shared_ptr<T>(new T(t1, t2, t3));
    }

    template <class T, class T1, class T2, class T3>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3 const & t3) {
        return std::shared_ptr<T>(new T(t1, t2, t3));
    }

    template <class T, class T1, class T2, class T3>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3 const & t3) {
        return std::shared_ptr<T>(new T(t1, t2, t3));
    }

    template <class T, class T1, class T2, class T3>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3 const & t3) {
        return std::shared_ptr<T>(new T(t1, t2, t3));
    }

    template <class T, class T1, class T2, class T3>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3 const & t3) {
        return std::shared_ptr<T>(new T(t1, t2, t3));
    }

    /// Construct an shared pointer with an allocator

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(t1, t2, t3), deleter(), a);
    /// \endcode
    /// where deleter is a deleter object supplied by this function.
    /// \tparam T The shared_ptr type
    /// \tparam A the allocator type
    /// \tparam T1 pointer constructor argument 1 type
    /// \tparam T2 pointer constructor argument 2 type
    /// \tparam T3 pointer constructor argument 3 type
    /// \param a The allocator object used to allocate the shared_ptr
    /// \param t1 pointer constructor argument 1
    /// \param t2 pointer constructor argument 2
    /// \param t3 pointer constructor argument 3
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A, class T1, class T2, class T3>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3& t3) {
        return std::shared_ptr<T>(new T(t1, t2, t3), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3& t3) {
        return std::shared_ptr<T>(new T(t1, t2, t3), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3& t3) {
        return std::shared_ptr<T>(new T(t1, t2, t3), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3& t3) {
        return std::shared_ptr<T>(new T(t1, t2, t3), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3 const & t3) {
        return std::shared_ptr<T>(new T(t1, t2, t3), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3 const & t3) {
        return std::shared_ptr<T>(new T(t1, t2, t3), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3 const & t3) {
        return std::shared_ptr<T>(new T(t1, t2, t3), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3 const & t3) {
        return std::shared_ptr<T>(new T(t1, t2, t3), &madness::detail::checked_delete, a);
    }

    /// Construct an shared pointer

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(t1, t2, t3, t4));
    /// \endcode
    /// \tparam T The shared_ptr type
    /// \tparam T1 pointer constructor argument 1 type
    /// \tparam T2 pointer constructor argument 2 type
    /// \tparam T3 pointer constructor argument 3 type
    /// \tparam T4 pointer constructor argument 4 type
    /// \param t1 pointer constructor argument 1
    /// \param t2 pointer constructor argument 2
    /// \param t3 pointer constructor argument 3
    /// \param t4 pointer constructor argument 4
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3& t3, T4& t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4));
    }

    template <class T, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3& t3, T4& t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4));
    }

    template <class T, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3& t3, T4& t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4));
    }

    template <class T, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3& t3, T4& t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4));
    }

    template <class T, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3 const & t3, T4& t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4));
    }

    template <class T, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3 const & t3, T4& t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4));
    }

    template <class T, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3 const & t3, T4& t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4));
    }

    template <class T, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3 const & t3, T4& t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4));
    }

    template <class T, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3& t3, T4 const & t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4));
    }

    template <class T, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3& t3, T4 const & t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4));
    }

    template <class T, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3& t3, T4 const & t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4));
    }

    template <class T, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3& t3, T4 const & t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4));
    }

    template <class T, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3 const & t3, T4 const & t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4));
    }

    template <class T, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3 const & t3, T4 const & t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4));
    }

    template <class T, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3 const & t3, T4 const & t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4));
    }

    template <class T, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4));
    }

    /// Construct an shared pointer with an allocator

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(t1, t2, t3, t4), deleter(), a);
    /// \endcode
    /// where deleter is a deleter object supplied by this function.
    /// \tparam T The shared_ptr type
    /// \tparam A the allocator type
    /// \tparam T1 pointer constructor argument 1 type
    /// \tparam T2 pointer constructor argument 2 type
    /// \tparam T3 pointer constructor argument 3 type
    /// \tparam T4 pointer constructor argument 4 type
    /// \param a The allocator object used to allocate the shared_ptr
    /// \param t1 pointer constructor argument 1
    /// \param t2 pointer constructor argument 2
    /// \param t3 pointer constructor argument 3
    /// \param t4 pointer constructor argument 4
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3& t3, T4& t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3& t3, T4& t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3& t3, T4& t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3& t3, T4& t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3 const & t3, T4& t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3 const & t3, T4& t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3 const & t3, T4& t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3 const & t3,
            T4& t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3& t3, T4 const & t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3& t3, T4 const & t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3& t3, T4 const & t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3& t3,
            T4 const & t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3 const & t3, T4 const & t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3 const & t3,
            T4 const & t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3 const & t3,
            T4 const & t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3 const & t3,
            T4 const & t4) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4), &madness::detail::checked_delete, a);
    }

    /// Construct an shared pointer

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    /// \endcode
    /// \tparam T The shared_ptr type
    /// \tparam T1 pointer constructor argument 1 type
    /// \tparam T2 pointer constructor argument 2 type
    /// \tparam T3 pointer constructor argument 3 type
    /// \tparam T4 pointer constructor argument 4 type
    /// \tparam T5 pointer constructor argument 5 type
    /// \param t1 pointer constructor argument 1
    /// \param t2 pointer constructor argument 2
    /// \param t3 pointer constructor argument 3
    /// \param t4 pointer constructor argument 4
    /// \param t5 pointer constructor argument 5
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3& t3, T4& t4, T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3& t3, T4& t4, T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3& t3, T4& t4, T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3& t3, T4& t4, T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3 const & t3, T4& t4, T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3 const & t3, T4& t4, T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3 const & t3, T4& t4, T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3 const & t3, T4& t4, T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3& t3, T4 const & t4, T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3& t3, T4 const & t4, T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3& t3, T4 const & t4, T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3& t3, T4 const & t4, T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3 const & t3, T4 const & t4, T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3 const & t3, T4 const & t4, T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3 const & t3, T4 const & t4, T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4,
            T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3& t3, T4& t4, T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3& t3, T4& t4, T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3& t3, T4& t4, T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3& t3, T4& t4, T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3 const & t3, T4& t4, T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3 const & t3, T4& t4, T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3 const & t3, T4& t4, T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3 const & t3, T4& t4,
            T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3& t3, T4 const & t4, T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3& t3, T4 const & t4, T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3& t3, T4 const & t4, T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3& t3, T4 const & t4,
            T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3 const & t3, T4 const & t4, T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3 const & t3, T4 const & t4,
            T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3 const & t3, T4 const & t4,
            T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4,
            T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5));
    }

    /// Construct an shared pointer with an allocator

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), deleter(), a);
    /// \endcode
    /// where deleter is a deleter object supplied by this function.
    /// \tparam T The shared_ptr type
    /// \tparam A the allocator type
    /// \tparam T1 pointer constructor argument 1 type
    /// \tparam T2 pointer constructor argument 2 type
    /// \tparam T3 pointer constructor argument 3 type
    /// \tparam T4 pointer constructor argument 4 type
    /// \tparam T5 pointer constructor argument 5 type
    /// \param a The allocator object used to allocate the shared_ptr
    /// \param t1 pointer constructor argument 1
    /// \param t2 pointer constructor argument 2
    /// \param t3 pointer constructor argument 3
    /// \param t4 pointer constructor argument 4
    /// \param t5 pointer constructor argument 5
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3& t3, T4& t4, T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3& t3, T4& t4, T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3& t3, T4& t4, T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3& t3, T4& t4,
            T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3 const & t3, T4& t4, T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3 const & t3, T4& t4,
            T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3 const & t3, T4& t4,
            T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3 const & t3,
            T4& t4, T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3& t3, T4 const & t4, T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3& t3, T4 const & t4,
            T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3& t3, T4 const & t4,
            T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3& t3,
            T4 const & t4, T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3 const & t3, T4 const & t4,
            T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3 const & t3,
            T4 const & t4, T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3 const & t3,
            T4 const & t4, T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3 const & t3,
            T4 const & t4, T5& t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3& t3, T4& t4, T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3& t3, T4& t4,
            T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3& t3, T4& t4,
            T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3& t3, T4& t4,
            T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3 const & t3, T4& t4,
            T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3 const & t3, T4& t4,
            T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3 const & t3, T4& t4,
            T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3 const & t3,
            T4& t4, T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3& t3, T4 const & t4,
            T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3& t3, T4 const & t4,
            T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3& t3, T4 const & t4,
            T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3& t3,
            T4 const & t4, T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3 const & t3, T4 const & t4,
            T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3 const & t3,
            T4 const & t4, T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3 const & t3,
            T4 const & t4, T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3 const & t3,
            T4 const & t4, T5 const & t5) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5), &madness::detail::checked_delete, a);
    }

    /// Construct an shared pointer

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    /// \endcode
    /// \tparam T The shared_ptr type
    /// \tparam T1 pointer constructor argument 1 type
    /// \tparam T2 pointer constructor argument 2 type
    /// \tparam T3 pointer constructor argument 3 type
    /// \tparam T4 pointer constructor argument 4 type
    /// \tparam T5 pointer constructor argument 5 type
    /// \tparam T6 pointer constructor argument 6 type
    /// \param t1 pointer constructor argument 1
    /// \param t2 pointer constructor argument 2
    /// \param t3 pointer constructor argument 3
    /// \param t4 pointer constructor argument 4
    /// \param t5 pointer constructor argument 5
    /// \param t6 pointer constructor argument 6
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3& t3, T4& t4, T5& t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3& t3, T4& t4, T5& t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3& t3, T4& t4, T5& t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3& t3, T4& t4, T5& t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3 const & t3, T4& t4, T5& t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3 const & t3, T4& t4, T5& t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3 const & t3, T4& t4, T5& t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3 const & t3, T4& t4, T5& t5,
            T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3& t3, T4 const & t4, T5& t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3& t3, T4 const & t4, T5& t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3& t3, T4 const & t4, T5& t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3& t3, T4 const & t4, T5& t5,
            T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3 const & t3, T4 const & t4, T5& t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3 const & t3, T4 const & t4, T5& t5,
            T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3 const & t3, T4 const & t4, T5& t5,
            T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4,
            T5& t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3& t3, T4& t4, T5 const & t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3& t3, T4& t4, T5 const & t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3& t3, T4& t4, T5 const & t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3& t3, T4& t4, T5 const & t5,
            T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3 const & t3, T4& t4, T5 const & t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3 const & t3, T4& t4, T5 const & t5,
            T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3 const & t3, T4& t4, T5 const & t5,
            T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3 const & t3, T4& t4,
            T5 const & t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3& t3, T4 const & t4, T5 const & t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3& t3, T4 const & t4, T5 const & t5,
            T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3& t3, T4 const & t4, T5 const & t5,
            T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3& t3, T4 const & t4,
            T5 const & t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3 const & t3, T4 const & t4, T5 const & t5,
            T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3 const & t3, T4 const & t4,
            T5 const & t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3 const & t3, T4 const & t4,
            T5 const & t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4,
            T5 const & t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3& t3, T4& t4, T5& t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3& t3, T4& t4, T5& t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3& t3, T4& t4, T5& t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3& t3, T4& t4, T5& t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3 const & t3, T4& t4, T5& t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3 const & t3, T4& t4, T5& t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3 const & t3, T4& t4, T5& t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3 const & t3, T4& t4, T5& t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3& t3, T4 const & t4, T5& t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3& t3, T4 const & t4, T5& t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3& t3, T4 const & t4, T5& t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3& t3, T4 const & t4, T5& t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3 const & t3, T4 const & t4, T5& t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3 const & t3, T4 const & t4, T5& t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3 const & t3, T4 const & t4, T5& t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4,
            T5& t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3& t3, T4& t4, T5 const & t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3& t3, T4& t4, T5 const & t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3& t3, T4& t4, T5 const & t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3& t3, T4& t4, T5 const & t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3 const & t3, T4& t4, T5 const & t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3 const & t3, T4& t4, T5 const & t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3 const & t3, T4& t4, T5 const & t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3 const & t3, T4& t4,
            T5 const & t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3& t3, T4 const & t4, T5 const & t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3& t3, T4 const & t4, T5 const & t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3& t3, T4 const & t4, T5 const & t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3& t3, T4 const & t4,
            T5 const & t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2& t2, T3 const & t3, T4 const & t4, T5 const & t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2& t2, T3 const & t3, T4 const & t4,
            T5 const & t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1& t1, T2 const & t2, T3 const & t3, T4 const & t4,
            T5 const & t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    template <class T, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4,
            T5 const & t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6));
    }

    /// Construct an shared pointer with an allocator

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), deleter(), a);
    /// \endcode
    /// where deleter is a deleter object supplied by this function.
    /// \tparam T The shared_ptr type
    /// \tparam A the allocator type
    /// \tparam T1 pointer constructor argument 1 type
    /// \tparam T2 pointer constructor argument 2 type
    /// \tparam T3 pointer constructor argument 3 type
    /// \tparam T4 pointer constructor argument 4 type
    /// \tparam T5 pointer constructor argument 5 type
    /// \tparam T6 pointer constructor argument 6 type
    /// \param a The allocator object used to allocate the shared_ptr
    /// \param t1 pointer constructor argument 1
    /// \param t2 pointer constructor argument 2
    /// \param t3 pointer constructor argument 3
    /// \param t4 pointer constructor argument 4
    /// \param t5 pointer constructor argument 5
    /// \param t6 pointer constructor argument 6
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3& t3, T4& t4, T5& t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3& t3, T4& t4, T5& t5,
            T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3& t3, T4& t4, T5& t5,
            T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3& t3, T4& t4,
            T5& t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3 const & t3, T4& t4, T5& t5,
            T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3 const & t3, T4& t4,
            T5& t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3 const & t3, T4& t4,
            T5& t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3 const & t3,
            T4& t4, T5& t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3& t3, T4 const & t4, T5& t5,
            T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3& t3, T4 const & t4,
            T5& t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3& t3, T4 const & t4,
            T5& t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3& t3,
            T4 const & t4, T5& t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3 const & t3, T4 const & t4,
            T5& t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3 const & t3,
            T4 const & t4, T5& t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3 const & t3,
            T4 const & t4, T5& t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3 const & t3,
            T4 const & t4, T5& t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3& t3, T4& t4, T5 const & t5,
            T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3& t3, T4& t4,
            T5 const & t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3& t3, T4& t4,
            T5 const & t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3& t3, T4& t4,
            T5 const & t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3 const & t3, T4& t4,
            T5 const & t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3 const & t3, T4& t4,
            T5 const & t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3 const & t3, T4& t4,
            T5 const & t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3 const & t3,
            T4& t4, T5 const & t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3& t3, T4 const & t4,
            T5 const & t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3& t3, T4 const & t4,
            T5 const & t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3& t3, T4 const & t4,
            T5 const & t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3& t3,
            T4 const & t4, T5 const & t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3 const & t3, T4 const & t4,
            T5 const & t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3 const & t3,
            T4 const & t4, T5 const & t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3 const & t3,
            T4 const & t4, T5 const & t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3 const & t3,
            T4 const & t4, T5 const & t5, T6& t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3& t3, T4& t4, T5& t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3& t3, T4& t4, T5& t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3& t3, T4& t4, T5& t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3& t3, T4& t4,
            T5& t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3 const & t3, T4& t4, T5& t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3 const & t3, T4& t4,
            T5& t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3 const & t3, T4& t4,
            T5& t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3 const & t3,
            T4& t4, T5& t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3& t3, T4 const & t4, T5& t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3& t3, T4 const & t4,
            T5& t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3& t3, T4 const & t4,
            T5& t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3& t3,
            T4 const & t4, T5& t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3 const & t3, T4 const & t4,
            T5& t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3 const & t3,
            T4 const & t4, T5& t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3 const & t3,
            T4 const & t4, T5& t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3 const & t3,
            T4 const & t4, T5& t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3& t3, T4& t4, T5 const & t5,
            T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3& t3, T4& t4,
            T5 const & t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3& t3, T4& t4,
            T5 const & t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3& t3, T4& t4,
            T5 const & t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3 const & t3, T4& t4,
            T5 const & t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3 const & t3, T4& t4,
            T5 const & t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3 const & t3, T4& t4,
            T5 const & t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3 const & t3,
            T4& t4, T5 const & t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3& t3, T4 const & t4,
            T5 const & t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3& t3, T4 const & t4,
            T5 const & t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3& t3, T4 const & t4,
            T5 const & t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3& t3,
            T4 const & t4, T5 const & t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2& t2, T3 const & t3, T4 const & t4,
            T5 const & t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2& t2, T3 const & t3,
            T4 const & t4, T5 const & t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1& t1, T2 const & t2, T3 const & t3,
            T4 const & t4, T5 const & t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3 const & t3,
            T4 const & t4, T5 const & t5, T6 const & t6) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6), &madness::detail::checked_delete,
                a);
    }

    /// Construct an shared pointer

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6, t7));
    /// \endcode
    /// \tparam T The shared_ptr type
    /// \tparam T1 pointer constructor argument 1 type
    /// \tparam T2 pointer constructor argument 2 type
    /// \tparam T3 pointer constructor argument 3 type
    /// \tparam T4 pointer constructor argument 4 type
    /// \tparam T5 pointer constructor argument 5 type
    /// \tparam T6 pointer constructor argument 6 type
    /// \tparam T7 pointer constructor argument 7 type
    /// \param t1 pointer constructor argument 1
    /// \param t2 pointer constructor argument 2
    /// \param t3 pointer constructor argument 3
    /// \param t4 pointer constructor argument 4
    /// \param t5 pointer constructor argument 5
    /// \param t6 pointer constructor argument 6
    /// \param t7 pointer constructor argument 7
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class T1, class T2, class T3, class T4, class T5, class T6, class T7>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4,
            T5 const & t5, T6 const & t6, T7 const & t7) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6, t7));
    }

    /// Construct an shared pointer with an allocator

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6, t7), deleter(), a);
    /// \endcode
    /// where deleter is a deleter object supplied by this function.
    /// \tparam T The shared_ptr type
    /// \tparam A the allocator type
    /// \tparam T1 pointer constructor argument 1 type
    /// \tparam T2 pointer constructor argument 2 type
    /// \tparam T3 pointer constructor argument 3 type
    /// \tparam T4 pointer constructor argument 4 type
    /// \tparam T5 pointer constructor argument 5 type
    /// \tparam T6 pointer constructor argument 6 type
    /// \tparam T7 pointer constructor argument 7 type
    /// \param a The allocator object used to allocate the shared_ptr
    /// \param t1 pointer constructor argument 1
    /// \param t2 pointer constructor argument 2
    /// \param t3 pointer constructor argument 3
    /// \param t4 pointer constructor argument 4
    /// \param t5 pointer constructor argument 5
    /// \param t6 pointer constructor argument 6
    /// \param t7 pointer constructor argument 7
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6,
            class T7>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3 const & t3,
            T4 const & t4, T5 const & t5, T6 const & t6, T7 const & t7) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6, t7),
                &madness::detail::checked_delete, a);
    }

    /// Construct an shared pointer

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6, t7, t8));
    /// \endcode
    /// \tparam T The shared_ptr type
    /// \tparam T1 pointer constructor argument 1 type
    /// \tparam T2 pointer constructor argument 2 type
    /// \tparam T3 pointer constructor argument 3 type
    /// \tparam T4 pointer constructor argument 4 type
    /// \tparam T5 pointer constructor argument 5 type
    /// \tparam T6 pointer constructor argument 6 type
    /// \tparam T7 pointer constructor argument 7 type
    /// \tparam T8 pointer constructor argument 8 type
    /// \param t1 pointer constructor argument 1
    /// \param t2 pointer constructor argument 2
    /// \param t3 pointer constructor argument 3
    /// \param t4 pointer constructor argument 4
    /// \param t5 pointer constructor argument 5
    /// \param t6 pointer constructor argument 6
    /// \param t7 pointer constructor argument 7
    /// \param t8 pointer constructor argument 8
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class T1, class T2, class T3, class T4, class T5, class T6, class T7,
            class T8>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4,
            T5 const & t5, T6 const & t6, T7 const & t7, T8 const & t8) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6, t7, t8));
    }

    /// Construct an shared pointer with an allocator

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6, t7, t8), deleter(), a);
    /// \endcode
    /// where deleter is a deleter object supplied by this function.
    /// \tparam T The shared_ptr type
    /// \tparam A the allocator type
    /// \tparam T1 pointer constructor argument 1 type
    /// \tparam T2 pointer constructor argument 2 type
    /// \tparam T3 pointer constructor argument 3 type
    /// \tparam T4 pointer constructor argument 4 type
    /// \tparam T5 pointer constructor argument 5 type
    /// \tparam T6 pointer constructor argument 6 type
    /// \tparam T7 pointer constructor argument 7 type
    /// \tparam T8 pointer constructor argument 8 type
    /// \param a The allocator object used to allocate the shared_ptr
    /// \param t1 pointer constructor argument 1
    /// \param t2 pointer constructor argument 2
    /// \param t3 pointer constructor argument 3
    /// \param t4 pointer constructor argument 4
    /// \param t5 pointer constructor argument 5
    /// \param t6 pointer constructor argument 6
    /// \param t7 pointer constructor argument 7
    /// \param t8 pointer constructor argument 8
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6,
            class T7, class T8>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3 const & t3,
            T4 const & t4, T5 const & t5, T6 const & t6, T7 const & t7, T8 const & t8) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6, t7, t8),
                &madness::detail::checked_delete, a);
    }

    /// Construct an shared pointer

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6, t7, t8, t9));
    /// \endcode
    /// where deleter is a deleter object supplied by this function.
    /// \tparam T The shared_ptr type
    /// \tparam T1 pointer constructor argument 1 type
    /// \tparam T2 pointer constructor argument 2 type
    /// \tparam T3 pointer constructor argument 3 type
    /// \tparam T4 pointer constructor argument 4 type
    /// \tparam T5 pointer constructor argument 5 type
    /// \tparam T6 pointer constructor argument 6 type
    /// \tparam T7 pointer constructor argument 7 type
    /// \tparam T8 pointer constructor argument 8 type
    /// \tparam T9 pointer constructor argument 9 type
    /// \param t1 pointer constructor argument 1
    /// \param t2 pointer constructor argument 2
    /// \param t3 pointer constructor argument 3
    /// \param t4 pointer constructor argument 4
    /// \param t5 pointer constructor argument 5
    /// \param t6 pointer constructor argument 6
    /// \param t7 pointer constructor argument 7
    /// \param t8 pointer constructor argument 8
    /// \param t9 pointer constructor argument 9
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class T1, class T2, class T3, class T4, class T5, class T6, class T7,
            class T8, class T9>
    std::shared_ptr<T> make_shared(T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4,
            T5 const & t5, T6 const & t6, T7 const & t7, T8 const & t8, T9 const & t9) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6, t7, t8, t9));
    }

    /// Construct an shared pointer with an allocator

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6, t7, t8, t9), deleter(), a);
    /// \endcode
    /// where deleter is a deleter object supplied by this function.
    /// \tparam T The shared_ptr type
    /// \tparam A the allocator type
    /// \tparam T1 pointer constructor argument 1 type
    /// \tparam T2 pointer constructor argument 2 type
    /// \tparam T3 pointer constructor argument 3 type
    /// \tparam T4 pointer constructor argument 4 type
    /// \tparam T5 pointer constructor argument 5 type
    /// \tparam T6 pointer constructor argument 6 type
    /// \tparam T7 pointer constructor argument 7 type
    /// \tparam T8 pointer constructor argument 8 type
    /// \tparam T9 pointer constructor argument 9 type
    /// \param a The allocator object used to allocate the shared_ptr
    /// \param t1 pointer constructor argument 1
    /// \param t2 pointer constructor argument 2
    /// \param t3 pointer constructor argument 3
    /// \param t4 pointer constructor argument 4
    /// \param t5 pointer constructor argument 5
    /// \param t6 pointer constructor argument 6
    /// \param t7 pointer constructor argument 7
    /// \param t8 pointer constructor argument 8
    /// \param t9 pointer constructor argument 9
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A, class T1, class T2, class T3, class T4, class T5, class T6,
            class T7, class T8, class T9>
    std::shared_ptr<T> allocate_shared(A const & a, T1 const & t1, T2 const & t2, T3 const & t3,
            T4 const & t4, T5 const & t5, T6 const & t6, T7 const & t7, T8 const & t8,
            T9 const & t9) {
        return std::shared_ptr<T>(new T(t1, t2, t3, t4, t5, t6, t7, t8, t9),
                &madness::detail::checked_delete, a);
    }

} // namespace std

#endif // MADNESS_USE_BOOST_TR1_MEMORY_HPP
#endif // MADNESS_HAS_STD_MAKE_SHARED
#endif // MADNESS_WORLD_SHAREDPTR_H__INCLUDED
