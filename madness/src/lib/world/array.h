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


  $Id: array.h 2249 2011-04-01 22:01:37Z justus.c79@gmail.com $
*/

#ifndef MADNESS_WORLD_ARRAY_H__INCLUDED
#define MADNESS_WORLD_ARRAY_H__INCLUDED

#include <madness_config.h>
#include <world/worldexc.h>
#include <world/worldhash.h>
#include <vector>
#include <algorithm>
#include <iostream>

// Select header that contains array
#if defined(MADNESS_USE_ARRAY)
#include <array>
#elif defined(MADNESS_USE_TR1_ARRAY)
#include <tr1/array>
#elif defined(MADNESS_USE_BOOST_TR1_ARRAY_HPP)
#include <boost/tr1/array.hpp>
#else
#error No acceptable array include directive was found.
#endif // ARRAY

#ifndef MADNESS_BEGIN_NAMESPACE_TR1

#if defined(BOOST_TR1_ARRAY_INCLUDED) || defined(BOOST_TR1_ARRAY_HPP_INCLUDED)

// We are using boost
#define MADNESS_BEGIN_NAMESPACE_TR1 namespace boost {
#define MADNESS_END_NAMESPACE_TR1 } // namespace std

#elif defined(MADNESS_HAS_STD_TR1_ARRAY)

// We are using TR1
#define MADNESS_BEGIN_NAMESPACE_TR1 namespace std { namespace tr1 {
#define MADNESS_END_NAMESPACE_TR1 } } // namespace std namespace tr1

#elif defined(MADNESS_HAS_STD_ARRAY)

// We are using C++0x
#define MADNESS_BEGIN_NAMESPACE_TR1 namespace std {
#define MADNESS_END_NAMESPACE_TR1 } // namespace std

#else
// We do not know.
#error Unable to determine the correct namespace for TR1 fuctional.

#endif

#endif // MADNESS_BEGIN_NAMESPACE_TR1

// Insert the tr1 array class into the std namespace.
namespace std {

#if defined(MADNESS_HAS_STD_TR1_ARRAY) && !defined(MADNESS_HAS_STD_ARRAY)
#define MADNESS_HAS_STD_ARRAY 1

    using ::std::tr1::array;
    using ::std::tr1::swap;
    using ::std::tr1::tuple_size;
    using ::std::tr1::tuple_element;
    using ::std::tr1::tuple_size;
    using ::std::tr1::tuple_element;
    using ::std::tr1::get;

#endif
} // namespace std

MADNESS_BEGIN_NAMESPACE_TR1

    /// Output std::array to stream for human consumption
    template <typename T, std::size_t N>
    std::ostream& operator<<(std::ostream& s, const array<T,N>& a) {
        s << "[";
        for(std::size_t i=0; i<N; ++i) {
            s << a[i];
            if (i != (N-1)) s << ",";
        }
        s << "]";
        return s;
    }

    /// Hash std::array with madness::Hash
    template <typename T, std::size_t N>
    madness::hashT hash_value(const array<T,N>& a) {
        // Use this version of range for potential optimization.
        return madness::hash_range(a.data(), N);
    }

    /// Hash std::array with MADNESS Hashing
    template <typename T, std::size_t N>
    struct hash<array<T,N> > {
        std::size_t operator()(const array<T,N>& a) const {
            return hash_value(a);
        }
    };

MADNESS_END_NAMESPACE_TR1

namespace madness {

    // Serialize std::array objects
    namespace archive {

        template <class Archive, class T>
        struct ArchiveStoreImpl;
        template <class Archive, class T>
        struct ArchiveLoadImpl;

        template <class Archive, typename T, std::size_t N>
        struct ArchiveStoreImpl<Archive, std::array<T,N> > {
            static void store(const Archive& ar, const std::array<T,N>& a) {
                for(typename std::array<T,N>::const_iterator it = a.begin(); it != a.end(); ++it)
                    ar & (*it);
            }
        };

        template <class Archive, typename T, std::size_t N>
        struct ArchiveLoadImpl<Archive, std::array<T,N> > {
            static void load(const Archive& ar, std::array<T,N>& a) {
                for(typename std::array<T,N>::iterator it = a.begin(); it != a.end(); ++it)
                    ar & (*it);
            }
        };

    } // namespace archive

    /// A simple, fixed dimension Coordinate

    /// Eliminates memory allocation cost, is just POD so can be
    /// copied easily and allocated on the stack, and the known
    /// dimension permits aggressive compiler optimizations.
    template <typename T, std::size_t N>
    class Vector {
    public:
        typedef std::array<T,N> arrayT;

    private:
        arrayT data_;

    public:
        // type defs
        typedef typename arrayT::value_type             value_type;
        typedef typename arrayT::iterator               iterator;
        typedef typename arrayT::const_iterator         const_iterator;
        typedef typename arrayT::reverse_iterator       reverse_iterator;
        typedef typename arrayT::const_reverse_iterator const_reverse_iterator;
        typedef typename arrayT::reference              reference;
        typedef typename arrayT::const_reference        const_reference;
        typedef typename arrayT::size_type              size_type;
        typedef typename arrayT::difference_type        difference_type;

        /// The size of the vector
        static const size_type static_size = N;

        /// Default constructor does not initialize vector contents
        Vector() {}

        /// Initialize all elements to value t
        template <typename Q>
        explicit Vector(Q t) {
            fill(t);
        }

        /// Construct from a C-style array of the same dimension
        template <typename Q>
        explicit Vector(const Q (&t)[N]) {
            std::copy(t, t + N, data_.begin());
        }

        /// Construct from an STL vector of equal or greater length
        template <typename Q, typename A>
        explicit Vector(const std::vector<Q, A>& t) {
            operator=(t);
        }

        template <typename Q>
        explicit Vector(const std::array<Q, N>& t) {
            data_ = t;
        }

        /// Copy constructor is deep (since a coordinate is POD)
        Vector(const Vector<T,N>& other) {
            data_ = other.data_;
        }

        /// Copy constructor is deep (since a coordinate is POD)
        template <typename Q>
        Vector(const Vector<Q,N>& other) {
            data_ = other.data_;
        }

        /// Assignment is deep (since a vector is POD)
        Vector<T,N>& operator=(const Vector<T,N>& other) {
            data_ = other.data_;
            return *this;
        }

        /// Assignment is deep (since a vector is POD)
        template <typename Q>
        Vector<T,N>& operator=(const Vector<Q,N>& other) {
            data_ = other.data_;
            return *this;
        }

        /// Assignment is deep (since a vector is POD)
        template <typename Q, typename A>
        Vector<T,N>& operator=(const std::vector<Q, A>& other) {
            MADNESS_ASSERT(other.size() >= N);
            std::copy(other.begin(), other.begin() + N, data_.begin());
            return *this;
        }

        /// Fill from scalar value
        Vector<T,N>& operator=(const T& t) {
            fill(t);
            return *this;
        }

        // type conversion
        operator std::array<T,N> () { return data_; }

         // iterator support
         iterator begin() { return data_.begin(); }
         const_iterator begin() const { return data_.begin(); }
         iterator end() { return data_.end(); }
         const_iterator end() const { return data_.end(); }

         // reverse iterator support
         reverse_iterator rbegin() { return data_.rbegin(); }
         const_reverse_iterator rbegin() const { return data_.rbegin(); }
         reverse_iterator rend() { return data_.rend(); }
         const_reverse_iterator rend() const { return data_.rend(); }

         // capacity
         size_type size() const { return data_.size(); }
         bool empty() const { return data_.empty(); }
         size_type max_size() const { return data_.max_size(); }

         // element access
         reference operator[](size_type i) { return data_[i]; }
         const_reference operator[](size_type i) const { return data_[i]; }
         reference at(size_type i) { return data_.at(i); }
         const_reference at(size_type i) const { return data_.at(i); }
         reference front() { return data_.front(); }
         const_reference front() const { return data_.front(); }
         reference back() { return data_.back(); }
         const_reference back() const { return data_.back(); }
         const T* data() const { return data_.data(); }
         T* c_array() { return data_.data(); }

         // modifiers
         void swap(Vector<T, N>& other) { data_.swap(other.data_); }
         void fill(const T& t) {
#ifdef MADNESS_ARRAY_HAS_FILL
             data_.fill(t);
 #else
             data_.assign(t);
 #endif
         }

        /// In-place element-wise multiplcation by a scalar

        /// Returns a reference to this for chaining operations
        template <typename Q>
        Vector<T,N>& operator*=(Q q) {
            for(size_type i = 0; i < N; ++i)
                data_[i] *= q;
            return *this;
        }

        /// In-place element-wise addition of another vector

        /// Returns a reference to this for chaining operations
        template <typename Q>
        Vector<T,N>& operator+=(const Vector<Q,N>& q) {
            for(size_type i = 0; i < N; ++i)
                data_[i] += q[i];
            return *this;
        }

        /// In-place element-wise subtraction of another vector

        /// Returns a reference to this for chaining operations
        template <typename Q>
        Vector<T,N>& operator-=(const Vector<Q,N>& q) {
            for(size_type i = 0; i < N; ++i)
                data_[i] -= q[i];
            return *this;
        }

        /// Support for MADNESS serialization
        template <typename Archive>
        void serialize(Archive& ar) {
            ar & data_;
        }

        /// Support for MADNESS hashing
        hashT hash() const {
            return hash_value(data_);
        }

        // comparisons
        friend bool operator==(const Vector<T, N>& l, const Vector<T, N>& r) {
            return l.data_ == r.data_;
        }

        friend bool operator!=(const Vector<T, N>& l, const Vector<T, N>& r) {
            return l.data_ != r.data_;
        }

        friend bool operator<(const Vector<T, N>& l, const Vector<T, N>& r) {
            return l.data_ < r.data_;
        }

        friend bool operator>(const Vector<T, N>& l, const Vector<T, N>& r) {
            return l.data_ > r.data_;
        }

        friend bool operator<=(const Vector<T, N>& l, const Vector<T, N>& r) {
            return l.data_ <= r.data_;
        }

        friend bool operator>=(const Vector<T, N>& l, const Vector<T, N>& r) {
            return l.data_ >= r.data_;
        }

        /// Output Vector to stream for human consumption
        friend std::ostream& operator<<(std::ostream& s, const Vector<T,N>& v) {
            s << v.data_;
            return s;
        }
    }; // class Vector

    template <typename T, std::size_t N>
    void swap(Vector<T,N>& l, Vector<T,N>& r) {
        l.swap(r);
    }

    // Arithmetic operators


    /// Scale a coordinate

    /// Multiply the scalar value \c l by each \c Vector element
    /// \tparam T The left-hand \c Vector element type
    /// \tparam N The \c Vector size
    /// \tparam U The right-hand scalar type
    /// \param l The left-hand \c Vector
    /// \param r The right-hand scalar value to be multiplied by the \c Vector
    /// elements
    /// \return A new coordinate, \c c, where \c c[i]==(l[i]*r)
    template <typename T, std::size_t N, typename U>
    Vector<T,N> operator*(Vector<T,N> l, U r) {
        // coordinate passed by value to allow compiler optimization
        for (std::size_t i = 0; i < N; ++i)
            l[i] *= r;
        return l;
    }


    /// Scale a \c Vector

    /// Multiply the scalar value \c r by each \c Vector element
    /// \tparam T The left-hand \c Vector element type
    /// \tparam N The \c Vector size
    /// \tparam U The right-hand scalar type
    /// \param l The left-hand coordinate
    /// \param r The right-hand scalar value to be multiplied by the \c Vector
    /// elements
    /// \return A new coordinate, \c c, where \c c[i]==(l*r[i])
    template <typename T, typename U, std::size_t N>
    Vector<T,N> operator*(T l, Vector<U,N> r) {
        // coordinate passed by value to allow compiler optimization
        for (std::size_t i = 0; i < N; ++i)
            r[i] *= l;
        return l;
    }

    /// Multiply two \c Vector objects

    /// Do an element-wise multiplication of \c r and \c q and return the result
    /// in a new coordinate.
    /// \tparam T The left-hand \c Vector element type
    /// \tparam N The \c Vector size
    /// \tparam U The right-hand \c Vector element type
    /// \param l The left-hand \c Vector
    /// \param r The right-hand \c Vector
    /// \return A new \c Vector, \c c, where \c c[i]==(l[i]*r[i])
    template <typename T, std::size_t N, typename U>
    Vector<T,N> operator*(Vector<T,N> l, const Vector<U,N>& r) {
        // coordinate r passed by value to allow compiler optimization
        for (std::size_t i = 0; i < N; ++i)
            l[i] *= r[i];
        return l;
    }

    /// Add a scalar to a \c Vector

    /// Add the scalar value \c r to each \c Vector element
    /// \tparam T The left-hand \c Vector element type
    /// \tparam N The \c Vector size
    /// \tparam U The right-hand scalar type
    /// \param l The left-hand coordinate
    /// \param r The right-hand scalar value to be added to the \c Vector
    /// \return A new \c Vector, \c c, where \c c[i]==(l[i]+r)
    template <typename T, std::size_t N, typename U>
    Vector<T,N> operator+(Vector<T,N> l, U r) {
        // coordinate passed by value to allow compiler optimization
        for (std::size_t i = 0; i < N; ++i)
            l[i] += r;
        return l;
    }

    /// Add two \c Vector opbjects

    /// Do an element-wise addition of \c l and \c r and return the result in a
    /// new \c Vector.
    /// \tparam T The left-hand \c Vector element type
    /// \tparam N The \c Vector size
    /// \tparam U The right-hand \c Vector element type
    /// \param l The left-hand \c Vector
    /// \param r The right-hand \c Vector
    /// \return A new \c Vector, \c c, where \c c[i]==(l[i]+r[i])
    template <typename T, std::size_t N, typename U>
    Vector<T,N> operator+(Vector<T,N> l, const Vector<U,N>& r) {
        // coordinate r passed by value to allow compiler optimization
        for (std::size_t i = 0; i < N; ++i)
            l[i] += r[i];
        return l;
    }

    /// Subtract a scalar from a \c Vector

    /// Subtract the scalar value \c r from the \c Vector elements \c l[i]
    /// \tparam T The left-hand \c Vector element type
    /// \tparam N The \c Vector size
    /// \tparam U The right-hand scalar type
    /// \param l The left-hand \c Vector
    /// \param r The right-hand scalar value to be added to the \c Vector
    /// \return A new \c Vector, \c c, where \c c[i]==(l[i]-r)
    template <typename T, std::size_t N, typename U>
    Vector<T,N> operator-(Vector<T,N> l, U r) {
        // coordinate passed by value to allow compiler optimization
        for (std::size_t i = 0; i < N; ++i)
            l[i] -= r;
        return l;
    }

    /// Subtract two \c Vector

    /// Do an element-wise subtraction of \c l and \c r and return the result in
    /// a new coordinate.
    /// \tparam T The left-hand \c Vector element type
    /// \tparam N The \c Vector size
    /// \tparam U The right-hand \c Vector element type
    /// \param l The left-hand \c Vector
    /// \param r The right-hand \c Vector
    /// \return A new coordinate, \c c, where \c c[i]==(l[i]-r[i])
    template <typename T, std::size_t N, typename U>
    Vector<T,N> operator-(Vector<T,N> l, const Vector<U,N>& r) {
        // coordinate r passed by value to allow compiler optimization
        for (std::size_t i = 0; i < N; ++i)
            l[i] -= r[i];
        return l;
    }



    /// Your friendly neighborhood factory function
    template <typename T>
    Vector<T,1> vec(T x) {
        Vector<T,1> r; r[0] = x;
        return r;
    }

    /// Your friendly neighborhood factory function
    template <typename T>
    Vector<T,2> vec(T x, T y) {
        Vector<T,2> r; r[0] = x; r[1] = y;
        return r;
    }

    /// Your friendly neighborhood factory function
    template <typename T>
    Vector<T,3> vec(T x, T y, T z) {
        Vector<T,3> r; r[0] = x; r[1] = y; r[2] = z;
        return r;
    }

    /// Your friendly neighborhood factory function
    template <typename T>
    Vector<T,4> vec(T x, T y, T z, T xx) {
        Vector<T,4> r; r[0] = x; r[1] = y; r[2] = z; r[3] = xx;
        return r;
    }

    /// Your friendly neighborhood factory function
    template <typename T>
    Vector<T,5> vec(T x, T y, T z, T xx, T yy) {
        Vector<T,5> r; r[0] = x; r[1] = y; r[2] = z; r[3] = xx; r[4] = yy;
        return r;
    }

    /// Your friendly neighborhood factory function
    template <typename T>
    Vector<T,6> vec(T x, T y, T z, T xx, T yy, T zz) {
        Vector<T,6> r; r[0] = x; r[1] = y; r[2] = z; r[3] = xx; r[4] = yy; r[5] = zz;
        return r;
    }


    /// A simple, fixed-size, stack
    template <typename T, std::size_t N>
    class Stack {
    private:
        std::array<T,N> t;
        std::size_t n;

    public:
        Stack() : n(0) {}

        void push(const T& value) {
            MADNESS_ASSERT(n < N);
            t[n++] = value;
        }

        T& pop() {
            MADNESS_ASSERT(n > 0);
            return t[--n];
        }

        T& front() {
            MADNESS_ASSERT(n > 0);
            return t[n-1];
        }

        std::size_t size() const {
            return n;
        }

        void clear() {
            n = 0;
        }
    }; // class Stack


    /// Returns a Vector<T,1> initialized from the arguments
    template <typename T>
    inline std::array<T,1> array_factory(const T& v0) {
        std::array<T,1> v;
        v[0] = v0;
        return v;
    }

    /// Returns a Vector<T,2> initialized from the arguments
    template <typename T>
    inline std::array<T,2> array_factory(const T& v0, const T& v1) {
        std::array<T,2> v;
        v[0] = v0;
        v[1] = v1;
        return v;
    }

    /// Returns a Vector<T,3> initialized from the arguments
    template <typename T>
    inline std::array<T,3> array_factory(const T& v0, const T& v1,
                                     const T& v2) {
        std::array<T,3> v;
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        return v;
    }

    /// Returns a Vector<T,4> initialized from the arguments
    template <typename T>
    inline std::array<T,4> array_factory(const T& v0, const T& v1,
                                     const T& v2, const T& v3) {
        std::array<T,4> v;
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        v[3] = v3;
        return v;
    }

    /// Returns a Vector<T,5> initialized from the arguments
    template <typename T>
    inline std::array<T,5> array_factory(const T& v0, const T& v1,
                                     const T& v2, const T& v3,
                                     const T& v4) {
        std::array<T,5> v;
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        v[3] = v3;
        v[4] = v4;
        return v;
    }

    /// Returns a Vector<T,6> initialized from the arguments
    template <typename T>
    inline std::array<T,6> array_factory(const T& v0, const T& v1,
                                     const T& v2, const T& v3,
                                     const T& v4, const T& v5) {
        std::array<T,6> v;
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        v[3] = v3;
        v[4] = v4;
        v[5] = v5;
        return v;
    }
} // namespace madness

#endif // MADNESS_WORLD_ARRAY_H__INCLUDED
