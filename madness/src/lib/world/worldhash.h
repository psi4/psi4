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


  $Id: worldhash.h 2249 2011-04-01 22:01:37Z justus.c79@gmail.com $
*/


#ifndef MADNESS_WORLD_WORLDHASH_H__INCLUDED
#define MADNESS_WORLD_WORLDHASH_H__INCLUDED

/*!
 \file worldhash.h
 \brief Defines hash functions for use in distributed containers

MADNESS uses hashing functions are modeled after Boost.Functional/Hash. It has
many similar function calls including hash_value, hash_combine, and hash_range.
In addition, it is also compatible with C++ TR1 hashing functors. The
\c madness::Hash functor interface is identical to both boost::hash and
\c std::hash , and any one of these may be used. By default, MADNESS hashing
functions can hash all fundamental types including integral types, floating
point types, pointer types, std::string, std::wstring, and std::array. Presently
\c hashT is typedef to \c std::size_t .

Since having a good hash is important, we are using Bob Jenkin's "lookup v3"
hash from http://www.burtleburtle.net/bob/c/lookup3.c. The preferred interface
for these function is with the \c madness::Hash<T> functor, but you can also use
hash_value or the "lookup v3" interface directly.

\b Note: Bob Jenkin's "lookup v3" hash returns a uint32_t hash value, which is
cast to std::size_t. This is done for compatibility with std::hash.

\b WARNING: While both std::hash and madness::Hash have the same interface, they
will not generate the same hash values from the same key value.

MADNESS hashing consists of one functor and three functions
 \li \c madness::Hash is the hash functor and the primary interface for hashing
 \li \c hash_value() is for hashing a single value and is the function used by the Hash functor
 \li \c hash_range() is for hashing a group of values. You can use this function to hash
   an iterator range, a C-style array, or a pointer.
 \li \c hash_combine() hashs a given value and combines it with a seed. This is
 useful for combining multiple elements into a single hash value.

There are several options for creating and using hash functions for your custom
types. The easiest method is to define a \c hash_value function for your key
type in the same namespace. This function will automatically be used by MADNESS
hashing containers. The hashing function should have the following form.
\code
namespace MyNamespace {

    class Key {
        // ...
    };

    madness::hashT hash_value(const Key& t) {
        // ...
    }

} // namespace MyNamespace
\endcode
You may also use the intrusive method and define a hash member function for your
key.
\code
class Key {
public:

    // ...

    madness::hashT hash() const {
        // ...
    }
};
\endcode
You can create a specialization of madness::Hash for your type directly.
\code
class Key {
  // ...
};

namespace madness {

    template <>
    struct Hash<Key> {
        hashT operator()(const Key& a) const {
            // ...
        }
    };

} // namespace madness
\endcode
If you use any of the above methods, MADNESS hash_combine and hash_range functions
will be able to hash your custom types.

In addition to these methods, you can use std::hash, boost::hash, or create your
own custom hashing functor that has the same form as madness::hash. However, if
you want to use a hashing functor other than madness::Hash, you will need to
provide the appropriate template parameter to the hashing container.

*/

#include <madness_config.h>
#include <world/typestuff.h>
#include <world/enable_if.h>
#include <stdint.h>
#include <cstddef>
#include <iterator>

// Select header that contains hash
#if defined(MADNESS_USE_FUNCTIONAL)
#include <functional>

#elif defined(MADNESS_USE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(MADNESS_USE_BOOST_TR1_FUNCTIONAL_HPP)
#include <boost/tr1/functional.hpp>
#else
#error No acceptable functional include directive was found.
#endif // FUNCTIONAL

#ifndef MADNESS_BEGIN_NAMESPACE_TR1

#if defined(BOOST_TR1_FUNCTIONAL_INCLUDED) || defined(BOOST_TR1_FUNCTIONAL_HPP_INCLUDED)

// We are using boost
#define MADNESS_BEGIN_NAMESPACE_TR1 namespace boost {
#define MADNESS_END_NAMESPACE_TR1 } // namespace std

#elif defined(MADNESS_HAS_STD_TR1_HASH)

// We are using TR1
#define MADNESS_BEGIN_NAMESPACE_TR1 namespace std { namespace tr1 {
#define MADNESS_END_NAMESPACE_TR1 } } // namespace std namespace tr1

#elif defined(MADNESS_HAS_STD_HASH)

// We are using C++0x
#define MADNESS_BEGIN_NAMESPACE_TR1 namespace std {
#define MADNESS_END_NAMESPACE_TR1 } // namespace std

#else
// We do not know.
#error Unable to determine the correct namespace for TR1 fuctional.

#endif

#endif // MADNESS_BEGIN_NAMESPACE_TR1


#if defined(MADNESS_HAS_STD_TR1_HASH) && !defined(MADNESS_HAS_STD_HASH)
#define MADNESS_HAS_STD_HASH 1
// hash is in std::tr1 but we want it in std namespace
namespace std {
    using ::std::tr1::hash;
}
#endif

// Bob Jenkin's "lookup v3" hash from http://www.burtleburtle.net/bob/c/lookup3.c.
extern "C" {
    uint32_t hashword(const uint32_t *k, size_t length, uint32_t initval);
    uint32_t hashlittle(const void *key, size_t length, uint32_t initval);
}

namespace madness {

    // Hash, hash_value, hash_combine, and hash_range a la Boost.

    /// The hash value type
    typedef std::size_t hashT;

    /// Hash a single fundamental object

    /// \tparam T The fundamental type
    /// \param t The object to hash
    /// \return The hashed value
    /// \note Use heavily optimized hashword when sizeof(T) is multiple
    /// of sizeof(uint32_t) and presumably correctly aligned.
    template <class T>
    inline typename enable_if_c<std::is_fundamental<T>::value &&
        ((sizeof(T)%sizeof(uint32_t)) == 0),
        hashT>::type
    hash_value(const T t) {
        hashT result = 0ul;
        result = hashword(reinterpret_cast<const uint32_t*>(&t),
            std::integral_constant<std::size_t,sizeof(T)/sizeof(uint32_t)>::value, 0u);
        return result;
    }

    /// Hash a single fundamental object

    /// \tparam T The fundamental type
    /// \param t The object to hash
    /// \return The hashed value
    template <class T>
    inline typename madness::enable_if_c<std::is_fundamental<T>::value &&
        ((sizeof(T)%sizeof(uint32_t)) != 0),
        hashT>::type
    hash_value(const T t) {
        hashT result = 0ul;
        result = hashlittle(reinterpret_cast<const void*>(&t), sizeof(T), 0u);
        return result;
    }

    /// Hash a pointer address

    /// \tparam T The pointer type
    /// \param t The pointer to be hashed
    /// \return The hashed value
    template <typename T>
    inline hashT hash_value(const T* t) {
        const unsigned long n = reinterpret_cast<unsigned long>(t);
        return hash_value(n);
    }

    /// Hash a class object

    /// \tparam T The class type
    /// \param t The object to hash
    /// \return \c t.hash()
    template <typename T>
    inline typename enable_if_c<!(std::is_fundamental<T>::value ||
        std::is_pointer<T>::value || std::is_array<T>::value), hashT>::type
    hash_value(const T& t) {
        return t.hash();
    }

    /// Hash a string

    /// \tparam T The character type
    /// \param t The string to hash
    /// \return The hashed value of the string
    template <typename T>
    inline hashT hash_value(const std::basic_string<T>& t) {
        return hash_range(t.c_str(), t.size());
    }

    /// Hash functor

    /// This hash functor calls hash_value for the given type, \c T . The
    /// namespace for hash_value function is not specified so you are free to
    /// implement your own version for your data type as follows:
    /// \code
    /// namespace MyNamespace {
    ///     class MyClass;
    ///     madness::hashT hash_value(const MyClass& t) {
    ///         // ...
    ///     }
    /// } // namespace MyNamespace
    /// \endcode
    /// or you can specialize this functor directly.
    /// \tparam T The object type to hash
    /// \param t The object to be hashed
    /// \return The hashed value
    template <typename T>
    struct Hash {
        hashT operator()(const T& t) const {
            return hash_value(t);
        }
    }; // struct Hash

    namespace detail {
        /// Internal use only
        // We don't hash anything here. It is just used for combining a hashed
        // value with a seed.
        inline void combine_hash(hashT& seed, hashT hash) {
            seed ^= hash + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
    }

    /// Combine hash values

    /// This function uses the standard hash function.
    /// \tparam T The type to hash
    /// \param[in,out] seed The initial hash seed value
    /// \param[in] v The value to be hashed
    template <class T>
    inline void hash_combine(hashT& seed, const T& v) {
        Hash<T> hasher;
        return detail::combine_hash(seed, hasher(v));
    }

    /// Combine the hash values of an iterator range

    /// \tparam It the iterator type
    /// \param[in,out] seed The initial hash seed value
    /// \param[in] first The first element of the iterator range to be hashed
    /// \param[in] last The end of the iterator range to be hashed
    template <class It>
    inline void hash_range(hashT& seed, It first, It last) {
        Hash<typename std::iterator_traits<It>::value_type> hasher;
        for(; first != last; ++first)
            detail::combine_hash(seed, hasher(*first));
    }

    /// Combine the hash values of an iterator range

    /// \tparam It the iterator type
    /// \param[in] first The first element of the iterator range to be hashed
    /// \param[in] last The end of the iterator range to be hashed
    /// \return The hashed iterator range
    template <class It>
    inline hashT hash_range(It first, It last) {
        hashT seed = 0;
        hash_range(seed, first, last);

        return seed;
    }

    /// Combine the hash values of a C-style array

    /// \tparam T The type to be hashed
    /// \tparam n The size of the C-style array
    /// \param[in] t The array to be hashed
    /// \return The hashed array value
    template <class T, std::size_t n>
    inline hashT hash_range(const T(&t)[n]) {
        return hash_range(t, n);
    }

    /// Combine the hash values of a C-style array

    /// \tparam T The type to be hashed
    /// \tparam n The size of the C-style array
    /// \param[in,out] seed The initial hash seed value
    /// \param[in] t The array to be hashed
    /// \note This function uses std::hash
    template <class T, std::size_t n>
    inline void hash_range(hashT& seed, const T(&t)[n]) {
        hash_range(seed, t, n);
    }

    /// Combine the hash values of a pointer range

    /// \tparam T The type to be hashed
    /// \tparam n The size of the C-style array
    /// \param[in,out] seed The initial hash seed value
    /// \param[in] t A pointer to the beginning of the range to be hashed
    /// \param[in] n The number of elements to hashed
    /// \note May use heavily optimized hashword when n * sizeof(T) is multiple
    /// of sizeof(uint32_t) and presumably correctly aligned.
    template <class T>
    inline typename enable_if<std::is_fundamental<T> >::type
    hash_range(hashT& seed, const T* t, std::size_t n) {
        const std::size_t bytes = n * sizeof(T);
        if((bytes % sizeof(uint32_t)) == 0)
            seed = hashword(reinterpret_cast<const uint32_t *>(t), bytes/sizeof(uint32_t), seed);
        else
            seed = hashlittle(static_cast<const void *>(t), bytes, seed);
    }

    /// Combine the hash values of a pointer range

    /// \tparam T The type to be hashed
    /// \tparam n The size of the C-style array
    /// \param[in,out] seed The initial hash seed value
    /// \param[in] t A pointer to the beginning of the range to be hashed
    /// \param[in] n The number of elements to hashed
    template <class T>
    inline typename disable_if<std::is_fundamental<T> >::type
    hash_range(hashT& seed, const T* t, std::size_t n) {
        hash_range(seed, t, t + n);
    }

    /// Combine the hash values of a pointer range

    /// \tparam T The type to be hashed
    /// \tparam n The size of the C-style array
    /// \param t A pointer to the beginning of the range to be hashed
    /// \param n The number of elements to hashed
    /// \return The hashed pointer range value
    template <class T>
    inline hashT hash_range(const T* t, std::size_t n) {
        hashT seed = 0ul;
        hash_range(seed, t, n);
        return seed;
    }

} // namespace madness


#endif // MADNESS_WORLD_WORLDHASH_H__INCLUDED
