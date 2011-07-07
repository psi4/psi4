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


  $Id: tensor.h 2372 2011-06-15 00:23:15Z rjharrison $
*/

#ifndef MADNESS_TENSOR_TENSOR_H__INCLUDED
#define MADNESS_TENSOR_TENSOR_H__INCLUDED

#include <madness_config.h>
#include <misc/ran.h>
#include <world/posixmem.h>
#include <world/sharedptr.h>

#include <complex>
#include <vector>
#include <cmath>
#include <cstdlib>

#include <world/archive.h>

typedef std::complex<float> float_complex;
typedef std::complex<double> double_complex;

// These probably have to be included in this order
#include <tensor/tensor_macros.h>
#include <tensor/type_data.h>
#include <tensor/slice.h>
#include <tensor/vector_factory.h>
#include <tensor/basetensor.h>
#include <tensor/tensoriter.h>
#include <tensor/tensorexcept.h>
#include <tensor/aligned.h>
#include <tensor/mxm.h>
#include <tensor/mtxmq.h>


/*!
  \file tensor.h
  \brief Defines and implements most of Tensor
  \ingroup tensor
  \addtogroup tensor

  \par Introduction

  A tensor is a multi-dimensional array and does not incorporate any concepts
  of covariance and contravariance.

  When a new tensor is created, the underlying data is also allocated.
  E.g.,
  \code
  Tensor<double> a(3,4,5)
  \endcode
  creates a new 3-dimensional tensor and allocates a contiguous
  block of 60 doubles which are initialized to zero.  The dimensions
  (numbered from the left starting at 0) are in C or row-major
  order.  Thus, for the tensor \c a , the stride between successive
  elements of the right-most dimension is 1.  For the middle
  dimension it is 5.  For the left-most dimension it is 20.  Thus,
  the loops
  \code
  for (i=0; i<3; ++i)
      for (j=0; j<4; ++j)
          for (k=0; k<5; ++k)
              a(i,j,k) = ...
  \endcode
  will go sequentially (and thus efficiently) through memory.
  If the dimensions have been reordered (e.g., with \c swapdim()
  or \c map() ), or if the tensor is actually a slice of another
  tensor, then the layout in memory may be more complex and
  may not reflect a contiguous block of memory.

  Multiple tensors may be used to provide multiple identical or
  distinct views of the same data.  E.g., in the following
  \code
  Tensor<double> a(2,3);  // A new tensor initialized to zero
  Tensor<double> b = a;
  \endcode
  \c a and \c b provide identical views of the same data, thus
  \code
  b(1,2) = 99;
  cout << a(1,2) << endl;  // Outputs 99
  cout << b(1,2) << endl;  // Outputs 99
  \endcode

  \par Shallow copy and assignment

  It is important to appreciate that the views and the data are
  quite independent.  In particular, the default copy constructor
  and assignment operations only copy the tensor (the view) and not
  the data --- <em> i.e., the copy constructor and assigment operations
  only take shallow copies</em>.  This is for both consistency and
  efficiency.  Thus, assigning one tensor to another generates another
  view of the same data, replacing any previous view and not moving
  or copying any of the data.
  E.g.,
  \code
  Tensor<double> a(2,3);   // A new tensor initialized to zero
  Tensor<double> c(3,3,3); // Another new tensor
  Tensor<double> b = a;    // b is a view of the same data as a
  a = c;                   // a is now a view of c's data
  b = c                    // b is now also a view of c's data and the
  // data allocated originally for a is freed
  \endcode
  The above example also illustrates how reference counting is used
  to keep track of the underlying data.  Once there are no views
  of the data, it is automatically freed.

  There are only two ways to actually copy the underlying data.  A
  new, complete, and contigous copy of a tensor and its data may be
  generated with the \c copy() function.  Or, to copy data from one tensor
  into the data viewed by another tensor, you must use a Slice.

  \par Indexing

  One dimensional tensors (i.e., vectors) may be indexed using
  either square brackets (e.g., \c v[i] ) or round brackets (e.g.,
  \c v(i) ).  All higher-dimensional tensors must use only round
  brackets (e.g., \c t(i,j,k) ).  This is due to C++'s restriction
  that the indexing operator (\c [] ) can only have one argument.
  The indexing operation should generate efficient code.

  For the sake of efficiency, no bounds checking is performed by
  default by most single element indexing operations.  Checking can
  be enabled at compile time by defining \c -DTENSOR_BOUNDS_CHECKING for
  application files including \c tensor.h.  The MADNESS configure script
  has the option \c --enable-tensor-bound-checking to define the macro
  in \c madness_config.h .  The general indexing
  operation that takes a \c std::vector<long> index and all slicing
  operations always perform bounds checking.  To make indexing with
  checking a bit easier, a factory function has been provided for
  vectors ... but note you need to explicitly use longs as the
  index.
  \code
  Tensor<long> a(7,7,7);
  a(3,4,5) += 1;                    // OK ... adds 1 to element (3,4,5)
  a(3,4,9) += 1;                    // BAD ... undetected out-of-bounds access
  a(vector_factory(3L,4L,9L)) += 1; // OK ... out-bounds access will
  // be detected at runtime.
  \endcode

  \par Slicing

  Slices generate sub-tensors --- i.e., views of patches of the
  data.  E.g., to refer to all but the first and last elements in
  each dimension of a matrix use
  \code
  a(Slice(1,-2),Slice(1,-2))
  \endcode
  Or to view odd elements in each dimension
  \code
  a(Slice(0,-1,2),Slice(0,-1,2))
  \endcode
  A slice or patch of a
  tensor behaves exactly like a tensor \em except for assignment.
  When a slice is assigned to, the data is copied with the
  requirement that the source and destinations agree in size and
  shape (i.e., they conform).  Thus, to copy the all of the data
  from a to b,
  \code
  Tensor<double> a(3,4), b(3,4);
  a = 1;                              // Set all elements of a to 1
  b = 2;                              // Set all elements of b to 2
  a(Slice(0,-1,1),Slice(0,-1,1)) = b; // Copy all data from b to a
  a(_,_) = b(_,_);                    // Copy all data from b to a
  a(___) = b(___);                    // Copy all data from b to a
  a(Slice(1,2),Slice(1,2) = b;        // Error, do not conform
  \endcode
  Special slice values \c _ ,\c  _reverse, and \c  ___ have
  been defined to refer to all elements in a dimension, all
  elements in a dimension but reversed, and all elements in all
  dimensions, respectively.

  \par Iteration and algorithms

  See tensor_macros.h for documentation on the easiest mechanisms for iteration over
  elements of tensors and tips for optimization.  See \c TensorIterator for
  the most general form of iteration.
*/


#ifndef HAVE_STD_ABS_LONG
#ifndef HAVE_STD_LABS
namespace std {
    static long abs(long a) {
        return a>=0 ? a : -a;
    }
}
#else
namespace std {
    static long abs(long a) {
        return std::labs(a);
    }
}
#endif
#endif


namespace madness {
#define IS_ODD(n) ((n)&0x1)
#define IS_UNALIGNED(p) (((unsigned long)(p))&0x7)


    /// For real types return value, for complex return conjugate
    template <typename Q, bool iscomplex>
    struct conditional_conj_struct {
        static Q op(const Q& coeff) {
            return coeff;
        }
    };

    /// For real types return value, for complex return conjugate
    template <typename Q>
    struct conditional_conj_struct<Q,true> {
        static Q op(const Q& coeff) {
            return conj(coeff);
        }
    };

    /// For real types return value, for complex return conjugate
    template <typename Q>
    Q conditional_conj(const Q& coeff) {
        return conditional_conj_struct<Q,TensorTypeData<Q>::iscomplex>::op(coeff);
    }

    namespace detail {
        template <typename T> T mynorm(T t) {
            return t*t;
        }

        template <typename T> T mynorm(std::complex<T> t) {
            return std::norm(t);
        }
    }

    template <class T> class SliceTensor;


    /// A tensor is a multidimension array

    /// \ingroup tensor
    template <class T> class Tensor : public BaseTensor {
        template <class U> friend class SliceTensor;

    protected:
        T* restrict _p;
        std::shared_ptr<T> _shptr;


        void allocate(long nd, const long d[], bool dozero) {
            _id = TensorTypeData<T>::id;
            if (nd < 0) {
                _p = 0;
                _shptr.reset();
                _size = 0;
                _ndim = -1;
                return;
            }

            TENSOR_ASSERT(nd>0 && nd <= TENSOR_MAXDIM,"invalid ndim in new tensor", nd, 0);
            // sanity check ... 2GB in doubles
            for (int i=0; i<nd; ++i) {
                TENSOR_ASSERT(d[i]>=0 && d[i]<268435456, "invalid dimension size in new tensor",d[i],0);
            }
            set_dims_and_size(nd, d);
            if (_size) {
                TENSOR_ASSERT(_size>=0 && _size<268435456, "invalid size in new tensor",_size,0);
                try {
#if HAVE_IBMBGP
#define TENSOR_ALIGNMENT 32
#else
#define TENSOR_ALIGNMENT 16
#endif

#ifdef WORLD_GATHER_MEM_STATS
                    _p = new T[size];
                    _shptr = std::shared_ptr<T>(_p);
#else
                    if (posix_memalign((void **) &_p, TENSOR_ALIGNMENT, sizeof(T)*_size)) throw 1;
                    _shptr.reset(_p, &::madness::detail::checked_free<T>);
#endif
                }
                catch (...) {
                    std::printf("new failed nd=%ld type=%ld size=%ld\n", nd, id(), _size);
                    std::printf("  %ld %ld %ld %ld %ld %ld\n",
                                d[0], d[1], d[2], d[3], d[4], d[5]);
                    TENSOR_EXCEPTION("new failed",_size,this);
                }
                //std::printf("allocated %p [%ld]  %ld\n", _p, size, p.use_count());
                if (dozero) {
                    //T zero = 0; for (long i=0; i<_size; ++i) _p[i] = zero;
                    // or
#ifdef HAVE_MEMSET
                    memset((void *) _p, 0, _size*sizeof(T));
#else
                    aligned_zero(_size, _p);
#endif
                }
            }
            else {
                _p = 0;
                _shptr.reset();
            }
        }

        // Free memory and restore default constructor state
        void deallocate() {
            _p = 0;
            _shptr.reset();
            _size = 0;
            _ndim = -1;
        }

    public:
        /// C++ typename of this tensor.
        typedef T type;

        /// C++ typename of the real type associated with a complex type.
        typedef typename TensorTypeData<T>::scalar_type scalar_type;

        /// C++ typename of the floating point type associated with scalar real type
        typedef typename TensorTypeData<T>::float_scalar_type float_scalar_type;

        /// Default constructor does not allocate any data and sets ndim=-1, size=0, _p=0, and id.
        Tensor() : _p(0) {
            _id = TensorTypeData<T>::id;
        }

        /// Copy constructor is shallow (same as assignment)

        /// \em Caveat \em emptor: The shallow copy constructor has many virtues but
        /// enables you to violate constness with simple code such as
        /// \code
        /// const Tensor<double> a(5);
        /// Tensor<double> b(a);
        /// b[1] = 3; // a[1] is now also 3
        /// \endcode
        Tensor(const Tensor<T>& t) {
            _id = TensorTypeData<T>::id;
            *this = t;
        }

        /// Assignment is shallow (same as copy constructor)

        /// \em Caveat \em emptor: The shallow assignment has many virtues but
        /// enables you to violate constness with simple code such as
        /// \code
        /// const Tensor<double> a(5);
        /// Tensor<double> b;
        /// b = a;
        /// b[1] = 3; // a[1] is now also 3
        /// \endcode
        Tensor<T>& operator=(const Tensor<T>& t) {
            if (this != &t) {
                _p = t._p;
                _shptr = t._shptr;
                _size = t._size;
                _ndim = t._ndim;
                for (int i=0; i<TENSOR_MAXDIM; ++i) {
                    _dim[i] = t._dim[i];
                    _stride[i] = t._stride[i];
                }
            }
            return *this;
        }


        /// Type conversion makes a deep copy
        template <class Q> operator Tensor<Q>() const { // type conv => deep copy
            Tensor<Q> result = Tensor<Q>(this->_ndim,this->_dim,false);
            BINARY_OPTIMIZED_ITERATOR(Q, result, const T, (*this), *_p0 = (Q)(*_p1));
            return result;
        }

        /// Create and zero new 1-d tensor

        /// @param[in] d0 Size of dimension 0
        explicit Tensor(int d0) : _p(0) {
            _dim[0] = d0;
            allocate(1, _dim, true);
        }

        /// Create and zero new 1-d tensor

        /// @param[in] d0 Size of dimension 0
        explicit Tensor(long d0) : _p(0) {
            _dim[0] = d0;
            allocate(1, _dim, true);
        }

        /// Create and zero new 2-d tensor

        /// @param[in] d0 Size of dimension 0
        /// @param[in] d1 Size of dimension 1
        explicit Tensor(long d0, long d1) : _p(0) {
            _dim[0] = d0; _dim[1] = d1;
            allocate(2, _dim, true);
        }

        /// Create and zero new 3-d tensor

        /// @param[in] d0 Size of dimension 0
        /// @param[in] d1 Size of dimension 1
        /// @param[in] d2 Size of dimension 2
        explicit Tensor(long d0, long d1, long d2) : _p(0) {
            _dim[0] = d0; _dim[1] = d1; _dim[2] = d2;
            allocate(3, _dim, true);
        }

        /// Create and zero new 4-d tensor

        /// @param[in] d0 Size of dimension 0
        /// @param[in] d1 Size of dimension 1
        /// @param[in] d2 Size of dimension 2
        /// @param[in] d3 Size of dimension 3
        explicit Tensor(long d0, long d1, long d2, long d3) : _p(0) {
            _dim[0] = d0; _dim[1] = d1; _dim[2] = d2; _dim[3] = d3;
            allocate(4, _dim, true);
        }

        /// Create and zero new 5-d tensor

        /// @param[in] d0 Size of dimension 0
        /// @param[in] d1 Size of dimension 1
        /// @param[in] d2 Size of dimension 2
        /// @param[in] d3 Size of dimension 3
        /// @param[in] d4 Size of dimension 4
        explicit Tensor(long d0, long d1, long d2, long d3, long d4) : _p(0) {
            _dim[0] = d0; _dim[1] = d1; _dim[2] = d2; _dim[3] = d3; _dim[4] = d4;
            allocate(5, _dim, true);
        }

        /// Create and zero new 6-d tensor

        /// @param[in] d0 Size of dimension 0
        /// @param[in] d1 Size of dimension 1
        /// @param[in] d2 Size of dimension 2
        /// @param[in] d3 Size of dimension 3
        /// @param[in] d4 Size of dimension 4
        /// @param[in] d5 Size of dimension 5
        explicit Tensor(long d0, long d1, long d2, long d3, long d4, long d5) {
            _dim[0] = d0; _dim[1] = d1; _dim[2] = d2; _dim[3] = d3; _dim[4] = d4; _dim[5] = d5;
            allocate(6, _dim, true);
        }

        /// Create and optionally zero new n-d tensor. This is the most general constructor.

        /// @param[in] d Vector containing size of each dimension, number of dimensions inferred from vcector size.
        /// @param[in] dozero If true (default) the tensor is initialized to zero
        explicit Tensor(const std::vector<long>& d, bool dozero=true) : _p(0) {
            allocate(d.size(),&(d[0]),dozero);
        }

        /// Politically incorrect general constructor.

        /// @param[in] nd Number of dimensions
        /// @param[in] d Size of each dimension
        /// @param[in] dozero If true (default) the tensor is initialized to zero
        explicit Tensor(long nd, const long d[], bool dozero=true) : _p(0) {
            allocate(nd,d,dozero);
        }

        /// Inplace fill tensor with scalar

        /// @param[in] x Value used to fill tensor via assigment
        /// @return %Reference to this tensor
        Tensor<T>& operator=(T x) {
            UNARY_OPTIMIZED_ITERATOR(T,(*this),*_p0 = x);
            return *this;
        }

        /// Inplace fill with a scalar (legacy name)

        /// @param[in] x Value used to fill tensor via assigment
        /// @return %Reference to this tensor
        Tensor<T>& fill(T x) {
            *this = x;
            return *this;
        }

        /// Inplace addition of two tensors

        /// @param[in] t Conforming tensor to be added in-place to this tensor
        /// @return %Reference to this tensor
        template <typename Q>
        Tensor<T>& operator+=(const Tensor<Q>& t) {
            BINARY_OPTIMIZED_ITERATOR(T, (*this), const T, t, *_p0 += *_p1);
            return *this;
        }

        /// Inplace subtraction of two tensors

        /// @param[in] t Conforming tensor to be subtracted in-place from this tensor
        /// @return %Reference to this tensor
        template <typename Q>
        Tensor<T>& operator-=(const Tensor<Q>& t) {
            BINARY_OPTIMIZED_ITERATOR(T, (*this), const T, t, *_p0 -= *_p1);
            return *this;
        }

        /// Addition of two tensors to produce a new tensor

        /// @param[in] t Conforming tensor to be added out-of-place to this tensor
        /// @return New tensor
        template <typename Q>
        Tensor< TENSOR_RESULT_TYPE(T,Q) > operator+(const Tensor<Q>& t) const {
            typedef TENSOR_RESULT_TYPE(T,Q) resultT;
            Tensor<resultT> result(_ndim,_dim,false);
            TERNARY_OPTIMIZED_ITERATOR(resultT, result, const T, (*this), const Q, t, *_p0 = *_p1 + *_p2);
            return result;
        }

        /// Subtraction of two tensors to produce a new tensor

        /// @param[in] t Conforming tensor to be subtracted out-of-place from this tensor
        /// @return New tensor
        template <typename Q>
        Tensor< TENSOR_RESULT_TYPE(T,Q) > operator-(const Tensor<Q>& t) const {
            typedef TENSOR_RESULT_TYPE(T,Q) resultT;
            Tensor<resultT> result(_ndim,_dim,false);
            TERNARY_OPTIMIZED_ITERATOR(resultT, result, const T, (*this), const Q, t, *_p0 = *_p1 - *_p2);
            return result;
        }

        /// Multiplication of tensor by a scalar of a supported type to produce a new tensor

        /// @param[in] x Scalar value
        /// @return New tensor
        template <typename Q>
        typename IsSupported<TensorTypeData<Q>, Tensor<TENSOR_RESULT_TYPE(T,Q)> >::type
        operator*(const Q& x) const {
            typedef TENSOR_RESULT_TYPE(T,Q) resultT;
            Tensor<resultT> result(_ndim,_dim,false);
            BINARY_OPTIMIZED_ITERATOR(resultT, result, const T, (*this), *_p0 = *_p1 * x);
            return result;
        }

        /// Divide tensor by a scalar of a supported type to produce a new tensor

        /// @param[in] x Scalar value
        /// @return New tensor
        template <typename Q>
        typename IsSupported<TensorTypeData<Q>, Tensor<TENSOR_RESULT_TYPE(T,Q)> >::type
        operator/(const Q& x) const {
            typedef TENSOR_RESULT_TYPE(T,Q) resultT;
            Tensor<resultT> result(_ndim,_dim);
            BINARY_OPTIMIZED_ITERATOR(resultT, result, const T, (*this), *_p0 = *_p1 / x);
            return result;
        }

        /// Add a scalar of the same type to all elements of a tensor producing a new tensor

        /// @param[in] x Scalar value
        /// @return New tensor
        template <typename Q>
        typename IsSupported<TensorTypeData<Q>, Tensor<TENSOR_RESULT_TYPE(T,Q)> >::type
        operator+(const Q& x) const {
            typedef TENSOR_RESULT_TYPE(T,Q) resultT;
            Tensor<resultT> result(_ndim,_dim);
            BINARY_OPTIMIZED_ITERATOR(resultT, result, const T, (*this), *_p0 = *_p1 + x);
            return result;
        }

        /// Subtract a scalar of the same type from all elements producing a new tensor

        /// @param[in] x Scalar value
        /// @return New tensor
        template <typename Q>
        typename IsSupported<TensorTypeData<Q>, Tensor<TENSOR_RESULT_TYPE(T,Q)> >::type
        operator-(const Q& x) const {
            return (*this) + (-x);
        }

        /// Unary negation producing a new tensor

        /// @return New tensor
        Tensor<T> operator-() const {
            Tensor<T> result = Tensor<T>(_ndim,_dim,false);
            BINARY_OPTIMIZED_ITERATOR(T, result, const T, (*this), *(_p0) = - (*_p1));
            return result;
        }

        /// Inplace multiplication by scalar of supported type

        /// @param[in] x Scalar value
        /// @return %Reference to this tensor
        template <typename Q>
        typename IsSupported<TensorTypeData<Q>,Tensor<T>&>::type
        operator*=(const Q& x) {
            UNARY_OPTIMIZED_ITERATOR(T, (*this), *_p0 *= x);
            return *this;
        }

        /// Inplace multiplication by scalar of supported type (legacy name)

        /// @param[in] x Scalar value
        /// @return %Reference to this tensor
        template <typename Q>
        typename IsSupported<TensorTypeData<Q>,Tensor<T>&>::type
        scale(Q x) {
            return (*this)*=x;
        }

        /// Inplace increment by scalar of supported type

        /// @param[in] x Scalar value
        /// @return %Reference to this tensor
        template <typename Q>
        typename IsSupported<TensorTypeData<Q>,Tensor<T>&>::type
        operator+=(const Q& x) {
            UNARY_OPTIMIZED_ITERATOR(T, (*this), *_p0 += x);
            return *this;
        }

        /// Inplace decrement by scalar of supported type

        /// @param[in] x Scalar value
        /// @return %Reference to this tensor
        template <typename Q>
        typename IsSupported<TensorTypeData<Q>,Tensor<T>&>::type
        operator-=(const Q& x) {
            UNARY_OPTIMIZED_ITERATOR(T, (*this), *_p0 -= x);
            return *this;
        }


        /// Inplace complex conjugate

        /// @return %Reference to this tensor
        Tensor<T>& conj() {
            UNARY_OPTIMIZED_ITERATOR(T, (*this), *_p0 = conditional_conj(*_p0));
            return *this;
        }

        /// Inplace fill with random values ( \c [0,1] for floats, \c [0,MAXSIZE] for integers)

        /// @return %Reference to this tensor
        Tensor<T>& fillrandom() {
            if (iscontiguous()) {
                madness::RandomVector<T>(size(), ptr());
            }
            else {
                UNARY_OPTIMIZED_ITERATOR(T,(*this), *_p0 = madness::RandomValue<T>());
            }
            return *this;
        }

        /// Inplace fill with the index of each element

        /// Each element is assigned it's logical index according to this loop structure
        /// \code
        /// Tensor<float> t(5,6,7,...)
        /// long index=0;
        /// for (long i=0; i<_dim[0]; ++i)
        ///    for (long j=0; j<_dim[1]; ++j)
        ///       for (long k=0; k<_dim[2]; ++k)
        ///          ...
        ///          tensor(i,j,k,...) = index++
        /// \endcode
        ///
        /// @return %Reference to this tensor
        Tensor<T>& fillindex() {
            long count = 0;
            UNARY_UNOPTIMIZED_ITERATOR(T,(*this), *_p0 = count++); // Fusedim would be OK
            return *this;
        }

        /// Inplace set elements of \c *this less than \c x in absolute magnitude to zero.

        /// @param[in] x Scalar value
        /// @return %Reference to this tensor
        Tensor<T>& screen(double x) {
            T zero = 0;
            UNARY_OPTIMIZED_ITERATOR(T,(*this), if (std::abs(*_p0)<x) *_p0=zero);
            return *this;
        }


        /// Return true if bounds checking was enabled at compile time

        /// @return True if bounds checking was enabled at compile time
        static bool bounds_checking() {
#ifdef TENSOR_BOUNDS_CHECKING
            return true;
#else
            return false;
#endif
        }

        /// 1-d indexing operation using \c [] \em without bounds checking.

        /// @param[in] i index for dimension 0
        /// @return %Reference to element
        T& operator[](long i) {
#ifdef TENSOR_BOUNDS_CHECKING
            TENSOR_ASSERT(i>=0 && i<_dim[0],"1d bounds check failed dim=0",i,this);
#endif
            return _p[i*_stride[0]];
        }

        /// 1-d indexing operation using \c [] \em without bounds checking.

        /// @param[in] i index for dimension 0
        /// @return %Reference to element
        const T& operator[](long i) const {
#ifdef TENSOR_BOUNDS_CHECKING
            TENSOR_ASSERT(i>=0 && i<_dim[0],"1d bounds check failed dim=0",i,this);
#endif
            return _p[i*_stride[0]];
        }

        /// 1-d indexing operation \em without bounds checking.

        /// @param[in] i index for dimension 0
        /// @return %Reference to element
        T& operator()(long i) {
#ifdef TENSOR_BOUNDS_CHECKING
            TENSOR_ASSERT(i>=0 && i<_dim[0],"1d bounds check failed dim=0",i,this);
#endif
            return _p[i*_stride[0]];
        }

        /// 1-d indexing operation \em without bounds checking.

        /// @param[in] i index for dimension 0
        /// @return %Reference to element
        const T& operator()(long i) const {
#ifdef TENSOR_BOUNDS_CHECKING
            TENSOR_ASSERT(i>=0 && i<_dim[0],"1d bounds check failed dim=0",i,this);
#endif
            return _p[i*_stride[0]];
        }

        /// 2-d indexing operation \em without bounds checking.

        /// @param[in] i index for dimension 0
        /// @param[in] j index for dimension 1
        /// @return %Reference to element
        T& operator()(long i, long j) {
#ifdef TENSOR_BOUNDS_CHECKING
            TENSOR_ASSERT(i>=0 && i<_dim[0],"2d bounds check failed dim=0",i,this);
            TENSOR_ASSERT(j>=0 && j<_dim[1],"2d bounds check failed dim=1",j,this);
#endif
            return _p[i*_stride[0]+j*_stride[1]];
        }

        /// 2-d indexing operation \em without bounds checking.

        /// @param[in] i index for dimension 0
        /// @param[in] j index for dimension 1
        /// @return %Reference to element
        const T& operator()(long i, long j) const {
#ifdef TENSOR_BOUNDS_CHECKING
            TENSOR_ASSERT(i>=0 && i<_dim[0],"2d bounds check failed dim=0",i,this);
            TENSOR_ASSERT(j>=0 && j<_dim[1],"2d bounds check failed dim=1",j,this);
#endif
            return _p[i*_stride[0]+j*_stride[1]];
        }

        /// 3-d indexing operation \em without bounds checking.

        /// @param[in] i index for dimension 0
        /// @param[in] j index for dimension 1
        /// @param[in] k index for dimension 2
        /// @return %Reference to element
        T& operator()(long i, long j, long k) {
#ifdef TENSOR_BOUNDS_CHECKING
            TENSOR_ASSERT(i>=0 && i<_dim[0],"3d bounds check failed dim=0",i,this);
            TENSOR_ASSERT(j>=0 && j<_dim[1],"3d bounds check failed dim=1",j,this);
            TENSOR_ASSERT(k>=0 && k<_dim[2],"3d bounds check failed dim=2",k,this);
#endif
            return _p[i*_stride[0]+j*_stride[1]+k*_stride[2]];
        }

        /// 3-d indexing operation \em without bounds checking.

        /// @param[in] i index for dimension 0
        /// @param[in] j index for dimension 1
        /// @param[in] k index for dimension 2
        /// @return %Reference to element
        const T& operator()(long i, long j, long k) const {
#ifdef TENSOR_BOUNDS_CHECKING
            TENSOR_ASSERT(i>=0 && i<_dim[0],"3d bounds check failed dim=0",i,this);
            TENSOR_ASSERT(j>=0 && j<_dim[1],"3d bounds check failed dim=1",j,this);
            TENSOR_ASSERT(k>=0 && k<_dim[2],"3d bounds check failed dim=2",k,this);
#endif
            return _p[i*_stride[0]+j*_stride[1]+k*_stride[2]];
        }

        /// 4-d indexing operation \em without bounds checking.

        /// @param[in] i index for dimension 0
        /// @param[in] j index for dimension 1
        /// @param[in] k index for dimension 2
        /// @param[in] l index for dimension 3
        /// @return %Reference to element
        T& operator()(long i, long j, long k, long l) {
#ifdef TENSOR_BOUNDS_CHECKING
            TENSOR_ASSERT(i>=0 && i<_dim[0],"4d bounds check failed dim=0",i,this);
            TENSOR_ASSERT(j>=0 && j<_dim[1],"4d bounds check failed dim=1",j,this);
            TENSOR_ASSERT(k>=0 && k<_dim[2],"4d bounds check failed dim=2",k,this);
            TENSOR_ASSERT(l>=0 && l<_dim[3],"4d bounds check failed dim=3",l,this);
#endif
            return _p[i*_stride[0]+j*_stride[1]+k*_stride[2]+
                           l*_stride[3]];
        }

        /// 4-d indexing operation \em without bounds checking.

        /// @param[in] i index for dimension 0
        /// @param[in] j index for dimension 1
        /// @param[in] k index for dimension 2
        /// @param[in] l index for dimension 3
        /// @return %Reference to element
        const T& operator()(long i, long j, long k, long l) const {
#ifdef TENSOR_BOUNDS_CHECKING
            TENSOR_ASSERT(i>=0 && i<_dim[0],"4d bounds check failed dim=0",i,this);
            TENSOR_ASSERT(j>=0 && j<_dim[1],"4d bounds check failed dim=1",j,this);
            TENSOR_ASSERT(k>=0 && k<_dim[2],"4d bounds check failed dim=2",k,this);
            TENSOR_ASSERT(l>=0 && l<_dim[3],"4d bounds check failed dim=3",l,this);
#endif
            return _p[i*_stride[0]+j*_stride[1]+k*_stride[2]+
                           l*_stride[3]];
        }

        /// 5-d indexing operation \em without bounds checking.

        /// @param[in] i index for dimension 0
        /// @param[in] j index for dimension 1
        /// @param[in] k index for dimension 2
        /// @param[in] l index for dimension 3
        /// @param[in] m index for dimension 4
        /// @return %Reference to element
        T& operator()(long i, long j, long k, long l, long m) {
#ifdef TENSOR_BOUNDS_CHECKING
            TENSOR_ASSERT(i>=0 && i<_dim[0],"5d bounds check failed dim=0",i,this);
            TENSOR_ASSERT(j>=0 && j<_dim[1],"5d bounds check failed dim=1",j,this);
            TENSOR_ASSERT(k>=0 && k<_dim[2],"5d bounds check failed dim=2",k,this);
            TENSOR_ASSERT(l>=0 && l<_dim[3],"5d bounds check failed dim=3",l,this);
            TENSOR_ASSERT(m>=0 && m<_dim[4],"5d bounds check failed dim=4",m,this);
#endif
            return _p[i*_stride[0]+j*_stride[1]+k*_stride[2]+
                           l*_stride[3]+m*_stride[4]];
        }

        /// 5-d indexing operation \em without bounds checking.

        /// @param[in] i index for dimension 0
        /// @param[in] j index for dimension 1
        /// @param[in] k index for dimension 2
        /// @param[in] l index for dimension 3
        /// @param[in] m index for dimension 4
        /// @return %Reference to element
        const T& operator()(long i, long j, long k, long l, long m) const {
#ifdef TENSOR_BOUNDS_CHECKING
            TENSOR_ASSERT(i>=0 && i<_dim[0],"5d bounds check failed dim=0",i,this);
            TENSOR_ASSERT(j>=0 && j<_dim[1],"5d bounds check failed dim=1",j,this);
            TENSOR_ASSERT(k>=0 && k<_dim[2],"5d bounds check failed dim=2",k,this);
            TENSOR_ASSERT(l>=0 && l<_dim[3],"5d bounds check failed dim=3",l,this);
            TENSOR_ASSERT(m>=0 && m<_dim[4],"5d bounds check failed dim=4",m,this);
#endif
            return _p[i*_stride[0]+j*_stride[1]+k*_stride[2]+
                           l*_stride[3]+m*_stride[4]];
        }

        /// 6-d indexing operation \em without bounds checking.

        /// @param[in] i index for dimension 0
        /// @param[in] j index for dimension 1
        /// @param[in] k index for dimension 2
        /// @param[in] l index for dimension 3
        /// @param[in] m index for dimension 4
        /// @param[in] n index for dimension 5
        /// @return %Reference to element
        T& operator()(long i, long j, long k, long l, long m, long n) {
#ifdef TENSOR_BOUNDS_CHECKING
            TENSOR_ASSERT(i>=0 && i<_dim[0],"6d bounds check failed dim=0",i,this);
            TENSOR_ASSERT(j>=0 && j<_dim[1],"6d bounds check failed dim=1",j,this);
            TENSOR_ASSERT(k>=0 && k<_dim[2],"6d bounds check failed dim=2",k,this);
            TENSOR_ASSERT(l>=0 && l<_dim[3],"6d bounds check failed dim=3",l,this);
            TENSOR_ASSERT(m>=0 && m<_dim[4],"6d bounds check failed dim=4",m,this);
            TENSOR_ASSERT(n>=0 && n<_dim[5],"6d bounds check failed dim=5",n,this);
#endif
            return _p[i*_stride[0]+j*_stride[1]+k*_stride[2]+
                           l*_stride[3]+m*_stride[4]+n*_stride[5]];
        }

        /// 6-d indexing operation \em without bounds checking.

        /// @param[in] i index for dimension 0
        /// @param[in] j index for dimension 1
        /// @param[in] k index for dimension 2
        /// @param[in] l index for dimension 3
        /// @param[in] m index for dimension 4
        /// @param[in] n index for dimension 5
        /// @return %Reference to element
        const T& operator()(long i, long j, long k, long l, long m, long n) const {
#ifdef TENSOR_BOUNDS_CHECKING
            TENSOR_ASSERT(i>=0 && i<_dim[0],"6d bounds check failed dim=0",i,this);
            TENSOR_ASSERT(j>=0 && j<_dim[1],"6d bounds check failed dim=1",j,this);
            TENSOR_ASSERT(k>=0 && k<_dim[2],"6d bounds check failed dim=2",k,this);
            TENSOR_ASSERT(l>=0 && l<_dim[3],"6d bounds check failed dim=3",l,this);
            TENSOR_ASSERT(m>=0 && m<_dim[4],"6d bounds check failed dim=4",m,this);
            TENSOR_ASSERT(n>=0 && n<_dim[5],"6d bounds check failed dim=5",n,this);
#endif
            return _p[i*_stride[0]+j*_stride[1]+k*_stride[2]+
                           l*_stride[3]+m*_stride[4]+n*_stride[5]];
        }

        /// Politically incorrect general indexing operation \em without bounds checking.

        /// @param[in] ind Array containing index for each dimension
        /// @return %Reference to element
        T& operator()(const long ind[]) {
            long offset = 0;
            for (int d=0; d<_ndim; ++d) {
                long i = ind[d];
#ifdef TENSOR_BOUNDS_CHECKING
                TENSOR_ASSERT(i>=0 && i<_dim[0],"non-PC general indexing bounds check failed dim=",d,this);
#endif
                offset += i*_stride[d];
            }
            return _p[offset];
        }

        /// Politically incorrect general indexing operation \em without bounds checking.

        /// @param[in] ind Array containing index for each dimension
        /// @return %Reference to element
        const T& operator()(const long ind[]) const {
            long offset = 0;
            for (int d=0; d<_ndim; ++d) {
                long i = ind[d];
#ifdef TENSOR_BOUNDS_CHECKING
                TENSOR_ASSERT(i>=0 && i<_dim[0],"non-PC general indexing bounds check failed dim=",d,this);
#endif
                offset += i*_stride[d];
            }
            return _p[offset];
        }

        /// General indexing operation \em with bounds checking.

        /// @param[in] ind Vector containing index for each dimension
        /// @return %Reference to element
        T& operator()(const std::vector<long> ind) {
            TENSOR_ASSERT(ind.size()>=(unsigned int) _ndim,"invalid number of dimensions",ind.size(),this);
            long index=0;
            for (long d=0; d<_ndim; ++d) {
                TENSOR_ASSERT(ind[d]>=0 && ind[d]<_dim[d],"out-of-bounds access",ind[d],this);
                index += ind[d]*_stride[d];
            }
            return _p[index];
        }

        /// General indexing operation \em with bounds checking.

        /// @param[in] ind Vector containing index for each dimension
        /// @return %Reference to element
        const T& operator()(const std::vector<long> ind) const {
            TENSOR_ASSERT(ind.size()>=(unsigned int) _ndim,"invalid number of dimensions",ind.size(),this);
            long index=0;
            for (long d=0; d<_ndim; ++d) {
                TENSOR_ASSERT(ind[d]>=0 && ind[d]<_dim[d],"out-of-bounds access",ind[d],this);
                index += ind[d]*_stride[d];
            }
            return _p[index];
        }

        /// General slicing operation

        /// @param[in] s Vector containing slice for each dimension
        /// @return SliceTensor viewing patch of original tensor
        SliceTensor<T> operator()(const std::vector<Slice>& s) {
            TENSOR_ASSERT(s.size()>=(unsigned)(this->ndim()), "invalid number of dimensions",
                          this->ndim(),this);
            return SliceTensor<T>(*this,&(s[0]));
        }

        /// General slicing operation (const)

        /// @param[in] s Vector containing slice for each dimension
        /// @return Constant Tensor viewing patch of original tensor
        const Tensor<T> operator()(const std::vector<Slice>& s) const {
            TENSOR_ASSERT(s.size()>=(unsigned)(this->ndim()), "invalid number of dimensions",
                          this->ndim(),this);
            return SliceTensor<T>(*this,&(s[0]));
        }

        /// Return a 1d SliceTensor that views the specified range of the 1d Tensor

        /// @return SliceTensor viewing patch of original tensor
        SliceTensor<T> operator()(const Slice& s0) {
            TENSOR_ASSERT(this->ndim()==1,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[1] = {s0};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 1d SliceTensor that views the specified range of the 1d Tensor

        /// @return Constant Tensor viewing patch of original tensor: \f$ R(*,*,\ldots) \rightarrow I(*,*,\ldots) \f$
        const Tensor<T> operator()(const Slice& s0) const {
            TENSOR_ASSERT(this->ndim()==1,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[1] = {s0};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 1d SliceTensor that views the specified range of the 2d Tensor

        /// @return SliceTensor viewing patch of original tensor: \f$ R(*) \rightarrow I(i,*) \f$
        SliceTensor<T> operator()(long i, const Slice& s1) {
            TENSOR_ASSERT(this->ndim()==2,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[2] = {Slice(i,i,0),s1};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 1d SliceTensor that views the specified range of the 2d Tensor

        /// @return Constant Tensor viewing patch of original tensor
        const Tensor<T> operator()(long i, const Slice& s1) const {
            TENSOR_ASSERT(this->ndim()==2,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[2] = {Slice(i,i,0),s1};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 1d SliceTensor that views the specified range of the 2d Tensor

        /// @return SliceTensor viewing patch of original tensor
        SliceTensor<T> operator()(const Slice& s0, long j) {
            TENSOR_ASSERT(this->ndim()==2,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[2] = {s0,Slice(j,j,0)};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 1d constant Tensor that views the specified range of the 2d Tensor

        /// @return Constant Tensor viewing patch of original tensor
        const Tensor<T> operator()(const Slice& s0, long j) const {
            TENSOR_ASSERT(this->ndim()==2,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[2] = {s0,Slice(j,j,0)};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 2d SliceTensor that views the specified range of the 2d Tensor

        /// @return SliceTensor viewing patch of original tensor
        SliceTensor<T> operator()(const Slice& s0, const Slice& s1) {
            TENSOR_ASSERT(this->ndim()==2,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[2] = {s0,s1};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 2d constant Tensor that views the specified range of the 2d Tensor

        /// @return Constant Tensor viewing patch of original tensor
        const Tensor<T> operator()(const Slice& s0, const Slice& s1) const {
            TENSOR_ASSERT(this->ndim()==2,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[2] = {s0,s1};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 3d SliceTensor that views the specified range of the 3d Tensor

        /// @return SliceTensor viewing patch of original tensor
        SliceTensor<T> operator()(const Slice& s0, const Slice& s1, const Slice& s2) {
            TENSOR_ASSERT(this->ndim()==3,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[3] = {s0,s1,s2};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 3d constant Tensor that views the specified range of the 3d Tensor

        /// @return Constant Tensor viewing patch of original tensor
        const Tensor<T> operator()(const Slice& s0, const Slice& s1, const Slice& s2) const {
            TENSOR_ASSERT(this->ndim()==3,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[3] = {s0,s1,s2};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 2d SliceTensor that views the specified range of the 3d Tensor

        /// @return SliceTensor viewing patch of original tensor
        SliceTensor<T> operator()(long i, const Slice& s1, const Slice& s2) {
            TENSOR_ASSERT(this->ndim()==3,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[3] = {Slice(i,i,0),s1,s2};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 2d constant Tensor that views the specified range of the 3d Tensor

        /// @return Constant Tensor viewing patch of original tensor
        const Tensor<T> operator()(long i, const Slice& s1, const Slice& s2) const {
            TENSOR_ASSERT(this->ndim()==3,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[3] = {Slice(i,i,0),s1,s2};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 2d SliceTensor that views the specified range of the 3d Tensor

        /// @return SliceTensor viewing patch of original tensor
        SliceTensor<T> operator()(const Slice& s0, long j, const Slice& s2) {
            TENSOR_ASSERT(this->ndim()==3,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[3] = {s0,Slice(j,j,0),s2};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 2d constant Tensor that views the specified range of the 3d Tensor

        /// @return Constant Tensor viewing patch of original tensor
        const Tensor<T> operator()(const Slice& s0, long j, const Slice& s2) const {
            TENSOR_ASSERT(this->ndim()==3,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[3] = {s0,Slice(j,j,0),s2};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 2d SliceTensor that views the specified range of the 3d Tensor

        /// @return SliceTensor viewing patch of original tensor
        SliceTensor<T> operator()(const Slice& s0, const Slice& s1, long k) {
            TENSOR_ASSERT(this->ndim()==3,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[3] = {s0,s1,Slice(k,k,0)};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 2d constant Tensor that views the specified range of the 3d Tensor

        /// @return Constant Tensor viewing patch of original tensor
        const Tensor<T> operator()(const Slice& s0, const Slice& s1, long k) const {
            TENSOR_ASSERT(this->ndim()==3,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[3] = {s0,s1,Slice(k,k,0)};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 1d SliceTensor that views the specified range of the 3d Tensor

        /// @return SliceTensor viewing patch of original tensor
        SliceTensor<T> operator()(long i, long j, const Slice& s2) {
            TENSOR_ASSERT(this->ndim()==3,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[3] = {Slice(i,i,0),Slice(j,j,0),s2};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 1d constant Tensor that views the specified range of the 3d Tensor

        /// @return Constant Tensor viewing patch of original tensor
        const Tensor<T> operator()(long i, long j, const Slice& s2) const {
            TENSOR_ASSERT(this->ndim()==3,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[3] = {Slice(i,i,0),Slice(j,j,0),s2};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 1d SliceTensor that views the specified range of the 3d Tensor

        /// @return SliceTensor viewing patch of original tensor
        SliceTensor<T> operator()(long i, const Slice& s1, long k) {
            TENSOR_ASSERT(this->ndim()==3,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[3] = {Slice(i,i,0),s1,Slice(k,k,0)};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 1d constant Tensor that views the specified range of the 3d Tensor

        /// @return Constant Tensor viewing patch of original tensor
        const Tensor<T> operator()(long i, const Slice& s1, long k) const {
            TENSOR_ASSERT(this->ndim()==3,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[3] = {Slice(i,i,0),s1,Slice(k,k,0)};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 1d SliceTensor that views the specified range of the 3d Tensor

        /// @return SliceTensor viewing patch of original tensor
        SliceTensor<T> operator()(const Slice& s0, long j, long k) {
            TENSOR_ASSERT(this->ndim()==3,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[3] = {s0,Slice(j,j,0),Slice(k,k,0)};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 1d constant Tensor that views the specified range of the 3d Tensor

        /// @return Constant Tensor viewing patch of original tensor
        const Tensor<T> operator()(const Slice& s0, long j, long k) const {
            TENSOR_ASSERT(this->ndim()==3,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[3] = {s0,Slice(j,j,0),Slice(k,k,0)};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 1-4d SliceTensor that views the specified range of the 4d Tensor

        /// @return SliceTensor viewing patch of original tensor
        SliceTensor<T> operator()(const Slice& s0, const Slice& s1, const Slice& s2,
                                  const Slice& s3) {
            TENSOR_ASSERT(this->ndim()==4,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[4] = {s0,s1,s2,s3};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 1-4d constant Tensor that views the specified range of the 4d Tensor

        /// @return Constant Tensor viewing patch of original tensor
        const Tensor<T> operator()(const Slice& s0, const Slice& s1, const Slice& s2,
                                  const Slice& s3) const {
            TENSOR_ASSERT(this->ndim()==4,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[4] = {s0,s1,s2,s3};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 1-5d SliceTensor that views the specified range of the 5d Tensor

        /// @return SliceTensor viewing patch of original tensor
        SliceTensor<T> operator()(const Slice& s0, const Slice& s1, const Slice& s2,
                                  const Slice& s3, const Slice& s4) {
            TENSOR_ASSERT(this->ndim()==5,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[5] = {s0,s1,s2,s3,s4};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 1-5d constant Tensor that views the specified range of the 5d Tensor

        /// @return Constant Tensor viewing patch of original tensor
        const Tensor<T> operator()(const Slice& s0, const Slice& s1, const Slice& s2,
                                  const Slice& s3, const Slice& s4) const {
            TENSOR_ASSERT(this->ndim()==5,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[5] = {s0,s1,s2,s3,s4};
            return SliceTensor<T>(*this,s);
        }

        /// Return a 1-6d SliceTensor that views the specified range of the 6d Tensor

        /// @return SliceTensor viewing patch of original tensor
        SliceTensor<T> operator()(const Slice& s0, const Slice& s1, const Slice& s2,
                                  const Slice& s3, const Slice& s4, const Slice& s5) {
            TENSOR_ASSERT(this->ndim()==6,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[6] = {s0,s1,s2,s3,s4,s5};
            return SliceTensor<T>(*this,s);
        }


        /// Return a 1-6d constant Tensor that views the specified range of the 6d Tensor

        /// @return Constant Tensor viewing patch of original tensor
        const Tensor<T> operator()(const Slice& s0, const Slice& s1, const Slice& s2,
                                  const Slice& s3, const Slice& s4, const Slice& s5) const {
            TENSOR_ASSERT(this->ndim()==6,"invalid number of dimensions",
                          this->ndim(),this);
            Slice s[6] = {s0,s1,s2,s3,s4,s5};
            return SliceTensor<T>(*this,s);
        }

        /// Returns new view/tensor reshaping size/number of dimensions to conforming tensor

        /// @param[in] ndimnew Number of dimensions in the result
        /// @param[in] d Array containing size of each new dimension
        /// @return New tensor (viewing same underlying data as the original but with different shape)
        Tensor<T> reshape(int ndimnew, const long* d) {
            Tensor<T> result(*this);
            result.reshape_inplace(ndimnew,d);
            return result;
        }

        /// Returns new view/tensor reshaping size/number of dimensions to conforming tensor

        /// @param[in] ndimnew Number of dimensions in the result
        /// @param[in] d Array containing size of each new dimension
        /// @return New tensor (viewing same underlying data as the original but with different shape)
        const Tensor<T> reshape(int ndimnew, const long* d) const {
            Tensor<T> result(*const_cast<Tensor<T>*>(this));
            result.reshape_inplace(ndimnew,d);
            return result;
        }

        /// Returns new view/tensor reshaping size/number of dimensions to conforming tensor

        /// @param[in] d Array containing size of each new dimension
        /// @return New tensor (viewing same underlying data as the original but with different shape)
        Tensor<T> reshape(const std::vector<long>& d) {
            return reshape(d.size(), &d[0]);
        }

        /// Returns new view/tensor reshaping size/number of dimensions to conforming tensor

        /// @param[in] d Array containing size of each new dimension
        /// @return New tensor (viewing same underlying data as the original but with different shape)
        const Tensor<T> reshape(const std::vector<long>& d) const {
            return reshape(d.size(), &d[0]);
        }

        /// Returns new view/tensor rehapings to conforming 1-d tensor with given dimension

        /// @param[in] dim0 Size of new dimension 0
        /// @return New tensor (viewing same underlying data as the original but with different shape)
        Tensor<T> reshape(long dim0) {
            long d[1] = {dim0};
            return reshape(1,d);
        }
        /// Returns new view/tensor rehapings to conforming 1-d tensor with given dimension

        /// @param[in] dim0 Size of new dimension 0
        /// @return New tensor (viewing same underlying data as the original but with different shape)
        const Tensor<T> reshape(long dim0) const {
            long d[1] = {dim0};
            return reshape(1,d);
        }

        /// Returns new view/tensor rehaping to conforming 2-d tensor with given dimensions

        /// @param[in] dim0 Size of new dimension 0
        /// @param[in] dim1 Size of new dimension 1
        /// @return New tensor (viewing same underlying data as the original but with different shape)
        Tensor<T> reshape(long dim0, long dim1) {
            long d[2] = {dim0,dim1};
            return reshape(2,d);
        }

        /// Returns new view/tensor rehaping to conforming 2-d tensor with given dimensions

        /// @param[in] dim0 Size of new dimension 0
        /// @param[in] dim1 Size of new dimension 1
        /// @return New tensor (viewing same underlying data as the original but with different shape)
        const Tensor<T> reshape(long dim0, long dim1) const {
            long d[2] = {dim0,dim1};
            return reshape(2,d);
        }

        /// Returns new view/tensor rehaping to conforming 3-d tensor with given dimensions

        /// @param[in] dim0 Size of new dimension 0
        /// @param[in] dim1 Size of new dimension 1
        /// @param[in] dim2 Size of new dimension 2
        /// @return New tensor (viewing same underlying data as the original but with different shape)
        Tensor<T> reshape(long dim0, long dim1, long dim2) {
            long d[3] = {dim0,dim1,dim2};
            return reshape(3,d);
        }

        /// Returns new view/tensor rehaping to conforming 3-d tensor with given dimensions

        /// @param[in] dim0 Size of new dimension 0
        /// @param[in] dim1 Size of new dimension 1
        /// @param[in] dim2 Size of new dimension 2
        /// @return New tensor (viewing same underlying data as the original but with different shape)
        const Tensor<T> reshape(long dim0, long dim1, long dim2) const {
            long d[3] = {dim0,dim1,dim2};
            return reshape(3,d);
        }

        /// Returns new view/tensor rehaping to conforming 4-d tensor with given dimensions

        /// @param[in] dim0 Size of new dimension 0
        /// @param[in] dim1 Size of new dimension 1
        /// @param[in] dim2 Size of new dimension 2
        /// @param[in] dim3 Size of new dimension 3
        /// @return New tensor (viewing same underlying data as the original but with different shape)
        Tensor<T> reshape(long dim0, long dim1, long dim2, long dim3) {
            long d[4] = {dim0,dim1,dim2,dim3};
            return reshape(4,d);
        }

        /// Returns new view/tensor rehaping to conforming 4-d tensor with given dimensions

        /// @param[in] dim0 Size of new dimension 0
        /// @param[in] dim1 Size of new dimension 1
        /// @param[in] dim2 Size of new dimension 2
        /// @param[in] dim3 Size of new dimension 3
        /// @return New tensor (viewing same underlying data as the original but with different shape)
        const Tensor<T> reshape(long dim0, long dim1, long dim2, long dim3) const {
            long d[4] = {dim0,dim1,dim2,dim3};
            return reshape(4,d);
        }

        /// Returns new view/tensor rehaping to conforming 5-d tensor with given dimensions

        /// @param[in] dim0 Size of new dimension 0
        /// @param[in] dim1 Size of new dimension 1
        /// @param[in] dim2 Size of new dimension 2
        /// @param[in] dim3 Size of new dimension 3
        /// @param[in] dim4 Size of new dimension 4
        /// @return New tensor (viewing same underlying data as the original but with different shape)
        Tensor<T> reshape(long dim0, long dim1, long dim2, long dim3, long dim4) {
            long d[5] = {dim0,dim1,dim2,dim3,dim4};
            return reshape(5,d);
        }

        /// Returns new view/tensor rehaping to conforming 5-d tensor with given dimensions

        /// @param[in] dim0 Size of new dimension 0
        /// @param[in] dim1 Size of new dimension 1
        /// @param[in] dim2 Size of new dimension 2
        /// @param[in] dim3 Size of new dimension 3
        /// @param[in] dim4 Size of new dimension 4
        /// @return New tensor (viewing same underlying data as the original but with different shape)
        const Tensor<T> reshape(long dim0, long dim1, long dim2, long dim3, long dim4) const {
            long d[5] = {dim0,dim1,dim2,dim3,dim4};
            return reshape(5,d);
        }

        /// Returns new view/tensor rehaping to conforming 6-d tensor with given dimensions

        /// @param[in] dim0 Size of new dimension 0
        /// @param[in] dim1 Size of new dimension 1
        /// @param[in] dim2 Size of new dimension 2
        /// @param[in] dim3 Size of new dimension 3
        /// @param[in] dim4 Size of new dimension 4
        /// @param[in] dim5 Size of new dimension 5
        /// @return New tensor (viewing same underlying data as the original but with different shape)
        Tensor<T> reshape(long dim0, long dim1, long dim2, long dim3, long dim4, long dim5) {
            long d[6] = {dim0,dim1,dim2,dim3,dim4,dim5};
            return reshape(6,d);
        }

        /// Returns new view/tensor rehaping to conforming 6-d tensor with given dimensions

        /// @param[in] dim0 Size of new dimension 0
        /// @param[in] dim1 Size of new dimension 1
        /// @param[in] dim2 Size of new dimension 2
        /// @param[in] dim3 Size of new dimension 3
        /// @param[in] dim4 Size of new dimension 4
        /// @param[in] dim5 Size of new dimension 5
        /// @return New tensor (viewing same underlying data as the original but with different shape)
        const Tensor<T> reshape(long dim0, long dim1, long dim2, long dim3, long dim4, long dim5) const {
            long d[6] = {dim0,dim1,dim2,dim3,dim4,dim5};
            return reshape(6,d);
        }

        /// Returns new view/tensor rehshaping to flat (1-d) tensor
        Tensor<T> flat() {
            long d[1] = {_size};
            return reshape(1,d);
        }

        /// Returns new view/tensor rehshaping to flat (1-d) tensor
        const Tensor<T> flat() const {
            long d[1] = {_size};
            return reshape(1,d);
        }

        /// Returns new view/tensor splitting dimension \c i as \c dimi0*dimi1 to produce conforming d+1 dimension tensor

        /// @return New tensor (viewing same underlying data as the original but with additional dimensions)
        Tensor<T> splitdim(long i, long dimi0, long dimi1) {
            Tensor<T> result(*this);
            result.splitdim_inplace(i, dimi0, dimi1);
            return result;
        }

        /// Returns new view/tensor splitting dimension \c i as \c dimi0*dimi1 to produce conforming d+1 dimension tensor

        /// @return New tensor (viewing same underlying data as the original but with additional dimensions)
        const Tensor<T> splitdim(long i, long dimi0, long dimi1) const {
            Tensor<T> result(*const_cast<Tensor<T>*>(this));
            result.splitdim_inplace(i, dimi0, dimi1);
            return result;
        }

        /// Returns new view/tensor fusing contiguous dimensions \c i and \c i+1

        /// @return New tensor (viewing same underlying data as the original but with fewer dimensions)
        Tensor<T> fusedim(long i) {
            Tensor<T> result(*this);
            result.fusedim_inplace(i);
            return result;
        }

        /// Returns new view/tensor fusing contiguous dimensions \c i and \c i+1

        /// @return New tensor (viewing same underlying data as the original but with fewer dimensions)
        const Tensor<T> fusedim(long i) const {
            Tensor<T> result(*const_cast<Tensor<T>*>(this));
            result.fusedim_inplace(i);
            return result;
        }

        /// Returns new view/tensor swaping dimensions \c i and \c j

        /// @return New tensor (viewing same underlying data as the original but with reordered dimensions)
        Tensor<T> swapdim(long idim, long jdim) {
            Tensor<T> result(*this);
            result.swapdim_inplace(idim, jdim);
            return result;
        }

        /// Returns new view/tensor swaping dimensions \c i and \c j

        /// @return New tensor (viewing same underlying data as the original but with reordered dimensions)
        const Tensor<T> swapdim(long idim, long jdim) const {
            Tensor<T> result(*const_cast<Tensor<T>*>(this));
            result.swapdim_inplace(idim, jdim);
            return result;
        }

        /// Returns new view/tensor permuting the dimensions

        /// @param[in] map Old dimension i becomes new dimension \c map[i]
        /// @return New tensor (viewing same underlying data as the original but with reordered dimensions)
        Tensor<T> mapdim(const std::vector<long>& map) {
            Tensor<T> result(*this);
            result.mapdim_inplace(map);
            return result;
        }

        /// Returns new view/tensor permuting the dimensions

        /// @return New tensor (viewing same underlying data as the original but with reordered dimensions)
        const Tensor<T> mapdim(const std::vector<long>& map) const {
            Tensor<T> result(*const_cast<Tensor<T>*>(this));
            result.mapdim_inplace(map);
            return result;
        }


        /// Returns new view/tensor cycling the sub-dimensions `(start,...,end)` with `shift` steps
        Tensor<T> cycledim(long nshift, long start, long end) {
            Tensor<T> result(*this);
            result.cycledim_inplace(nshift, start, end);
            return result;
        }


        /// Returns new view/tensor cycling the sub-dimensions `(start,...,end)` with `shift` steps
        const Tensor<T> cycledim(long nshift, long start, long end) const {
            Tensor<T> result(*const_cast<Tensor<T>*>(this));
            result.cycledim_inplace(nshift, start, end);
            return result;
        }


        /// Test if \c *this and \c t conform.
        template <class Q> bool conforms(const Tensor<Q>& t) const {
            return BaseTensor::conforms(&t);
        }

        /// Returns the sum of all elements of the tensor
        T sum() const {
            T result = 0;
            UNARY_OPTIMIZED_ITERATOR(const T,(*this),result += *_p0);
            return result;
        }

        /// Returns the sum of the squares of the elements
        T sumsq() const {
            T result = 0;
            UNARY_OPTIMIZED_ITERATOR(const T,(*this),result += (*_p0) * (*_p0));
            return result;
        }

        /// Return the product of all elements of the tensor
        T product() const {
            T result = 1;
            UNARY_OPTIMIZED_ITERATOR(const T,(*this),result *= *_p0);
            return result;
        }

        /// Return the minimum value (and if ind is non-null, its index) in the Tensor
        T min(long* ind=0) const {
            T result = *(this->_p);
            if (ind) {
                for (long i=0; i<_ndim; ++i) ind[i]=0;
                long nd = _ndim-1;
                UNARY_UNOPTIMIZED_ITERATOR(const T,(*this),
                                           if (result > *_p0) {
                                               result = *_p0;
                                               for (long i=0; i<nd; ++i) ind[i]=iter.ind[i];
                                               ind[nd] = _j;
                                           }
                                           );
            }
            else {
                UNARY_OPTIMIZED_ITERATOR(const T,(*this),result=std::min<T>(result,*_p0));
            }
            return result;
        }

        /// Return the maximum value (and if ind is non-null, its index) in the Tensor
        T max(long* ind=0) const {
            T result = *(this->_p);
            if (ind) {
                for (long i=0; i<_ndim; ++i) ind[i]=0;
                long nd = _ndim-1;
                UNARY_UNOPTIMIZED_ITERATOR(const T,(*this),
                                           if (result < *_p0) {
                                               result = *_p0;
                                               for (long i=0; i<nd; ++i) ind[i]=iter.ind[i];
                                               ind[nd] = _j;
                                           }
                                           );
            }
            else {
                UNARY_OPTIMIZED_ITERATOR(const T,(*this),result=std::max<T>(result,*_p0));
            }
            return result;
        }

        // For complex types, this next group returns the appropriate real type
        // For real types, the same type as T is returned (type_data.h)

        /// Returns the Frobenius norm of the tensor
        float_scalar_type normf() const {
            float_scalar_type result = 0;
            UNARY_OPTIMIZED_ITERATOR(const T,(*this),result += ::madness::detail::mynorm(*_p0));
            return (float_scalar_type) std::sqrt(result);
        }

        /// Return the absolute minimum value (and if ind is non-null, its index) in the Tensor
        scalar_type absmin(long *ind = 0) const {
            scalar_type result = std::abs(*(this->_p));
            if (ind) {
                for (long i=0; i<_ndim; ++i) ind[i]=0;
                long nd = _ndim-1;
                UNARY_UNOPTIMIZED_ITERATOR(const T,(*this),
                                           scalar_type absval = std::abs(*_p0);
                                           if (result > absval) {
                                               result = absval;
                                               for (long i=0; i<nd; ++i) ind[i]=iter.ind[i];
                                               ind[nd] = _j;
                                           }
                                           );
            }
            else {
                UNARY_OPTIMIZED_ITERATOR(const T,(*this),result=std::min<scalar_type>(result,std::abs(*_p0)));
            }
            return result;
        }

        /// Return the absolute maximum value (and if ind is non-null, its index) in the Tensor
        scalar_type absmax(long *ind = 0) const {
            scalar_type result = std::abs(*(this->_p));
            if (ind) {
                for (long i=0; i<_ndim; ++i) ind[i]=0;
                long nd = _ndim-1;
                UNARY_UNOPTIMIZED_ITERATOR(T,(*this),
                                           scalar_type absval = std::abs(*_p0);
                                           if (result < absval) {
                                               result = absval;
                                               for (long i=0; i<nd; ++i) ind[i]=iter.ind[i];
                                               ind[nd] = _j;
                                           }
                                           );
            }
            else {
                UNARY_OPTIMIZED_ITERATOR(const T,(*this),result=std::max<scalar_type>(result,std::abs(*_p0)));
            }
            return result;
        }


        /// Return the trace of two tensors (no complex conjugate invoked)
        T trace(const Tensor<T>& t) const {
            T result = 0;
            BINARY_OPTIMIZED_ITERATOR(const T,(*this),const T,t,result += (*_p0)*(*_p1));
            return result;
        }

        /// Return the trace of two tensors with complex conjugate of the leftmost (i.e., this)
        template <class Q>
        TENSOR_RESULT_TYPE(T,Q) trace_conj(const Tensor<Q>& t) const {
            TENSOR_RESULT_TYPE(T,Q) result = 0;
            BINARY_OPTIMIZED_ITERATOR(const T,(*this),const Q,t,result += conditional_conj(*_p0)*(*_p1));
            return result;
        }

        /// Inplace apply a unary function to each element of the tensor
        template <typename opT>
        Tensor<T>& unaryop(opT& op) {
            UNARY_OPTIMIZED_ITERATOR(T,(*this),*_p0=op(*_p0));
            return *this;
        }

        /// Inplace multiply by corresponding elements of argument Tensor
        Tensor<T>& emul(const Tensor<T>& t) {
            BINARY_OPTIMIZED_ITERATOR(T,(*this),const T,t,*_p0 *= *_p1);
            return *this;
        }

        /// Inplace generalized saxpy ... this = this*alpha + other*beta
        Tensor<T>& gaxpy(T alpha, const Tensor<T>& t, T beta) {
            if (iscontiguous() && t.iscontiguous()) {
                T* restrict a = ptr();
                const T* restrict b = t.ptr();
                if (alpha == T(1.0)) {
                    for (long i=0; i<_size; ++i) a[i] += b[i]*beta;
                }
                else {
                    for (long i=0; i<_size; ++i) a[i] = a[i]*alpha + b[i]*beta;
                }
            }
            else {
                //BINARYITERATOR(T,(*this),T,t, (*_p0) = alpha*(*_p0) + beta*(*_p1));
                BINARY_OPTIMIZED_ITERATOR(T,(*this),const T,t, (*_p0) = alpha*(*_p0) + beta*(*_p1));
                //ITERATOR((*this),(*this)(IND) = alpha*(*this)(IND) + beta*t(IND));
            }
            return *this;
        }

        /// Returns a pointer to the internal data
        T* ptr() {
            return _p;
        }

        /// Returns a pointer to the internal data
        const T* ptr() const {
            return _p;
        }

        /// Returns a pointer to the base class
        BaseTensor* base() {
            return static_cast<BaseTensor*>(this);
        }

        /// Returns a pointer to the base class
        const BaseTensor* base() const {
            return static_cast<BaseTensor*>(this);
        }

        /// Return iterator over single tensor
        TensorIterator<T> unary_iterator(long iterlevel=0,
                                         bool optimize=true,
                                         bool fusedim=true,
                                         long jdim=default_jdim) const {
            return TensorIterator<T>(this,(const Tensor<T>*) 0, (const Tensor<T>*) 0,
                                     iterlevel, optimize, fusedim, jdim);
        }

        /// Return iterator over two tensors
        template <class Q>
        TensorIterator<T,Q> binary_iterator(const Tensor<Q>& q,
                                            long iterlevel=0,
                                            bool optimize=true,
                                            bool fusedim=true,
                                            long jdim=default_jdim) const {
            return TensorIterator<T,Q>(this,&q,(const Tensor<T>*) 0,
                                       iterlevel, optimize, fusedim, jdim);
        }

        /// Return iterator over three tensors
        template <class Q, class R>
        TensorIterator<T,Q,R> ternary_iterator(const Tensor<Q>& q,
                                               const Tensor<R>& r,
                                               long iterlevel=0,
                                               bool optimize=true,
                                               bool fusedim=true,
                                               long jdim=default_jdim) const {
            return TensorIterator<T,Q,R>(this,&q,&r,
                                         iterlevel, optimize, fusedim, jdim);
        }

        /// End point for forward iteration
        const TensorIterator<T>& end() const {
            static TensorIterator<T> theend(0,0,0,0,0,0);
            return theend;
        }

        virtual ~Tensor() {}

        /// Frees all memory and resests to state of default constructor
        void clear() {deallocate();}
    };

    template <class T>
    std::ostream& operator << (std::ostream& out, const Tensor<T>& t);


    namespace archive {
        /// Serialize a tensor
        template <class Archive, typename T>
        struct ArchiveStoreImpl< Archive, Tensor<T> > {
            static void store(const Archive& s, const Tensor<T>& t) {
                if (t.iscontiguous()) {
                    s & t.size() & t.id();
                    if (t.size()) s & t.ndim() & wrap(t.dims(),TENSOR_MAXDIM) & wrap(t.ptr(),t.size());
                }
                else {
                    s & copy(t);
                }
            };
        };


        /// Deserialize a tensor ... existing tensor is replaced
        template <class Archive, typename T>
        struct ArchiveLoadImpl< Archive, Tensor<T> > {
            static void load(const Archive& s, Tensor<T>& t) {
                long sz, id;
                s & sz & id;
                if (id != t.id()) throw "type mismatch deserializing a tensor";
                if (sz) {
                    long _ndim, _dim[TENSOR_MAXDIM];
                    s & _ndim & wrap(_dim,TENSOR_MAXDIM);
                    t = Tensor<T>(_ndim, _dim, false);
                    if (sz != t.size()) throw "size mismatch deserializing a tensor";
                    s & wrap(t.ptr(), t.size());
                }
                else {
                    t = Tensor<T>();
                }
            };
        };

    }

    /// The class defines tensor op scalar ... here define scalar op tensor.

    /// \ingroup tensor
    template <typename T, typename Q>
    typename IsSupported < TensorTypeData<Q>, Tensor<T> >::type
    operator+(Q x, const Tensor<T>& t) {
        return t+x;
    }

    /// The class defines tensor op scalar ... here define scalar op tensor.

    /// \ingroup tensor
    template <typename T, typename Q>
    typename IsSupported < TensorTypeData<Q>, Tensor<T> >::type
    operator*(const Q& x, const Tensor<T>& t) {
        return t*x;
    }

    /// The class defines tensor op scalar ... here define scalar op tensor.

    /// \ingroup tensor
    template <typename T, typename Q>
    typename IsSupported < TensorTypeData<Q>, Tensor<T> >::type
    operator-(Q x, const Tensor<T>& t) {
        return (-t)+=x;
    }

    /// Returns a new contiguous tensor that is a deep copy of the input

    /// \ingroup tensor
    /// @result Returns a new contiguous tensor that is a deep copy of the input
    template <class T> Tensor<T> copy(const Tensor<T>& t) {
        if (t.size()) {
            Tensor<T> result = Tensor<T>(t.ndim(),t.dims(),false);
            BINARY_OPTIMIZED_ITERATOR(T, result, const T, t, *_p0 = *_p1);
            return result;
        }
        else {
            return Tensor<T>();
        }
    }

    /// Transforms one dimension of the tensor t by the matrix c, returns new contiguous tensor

    /// \ingroup tensor
    /// \code
    /// transform_dir(t,c,1) = r(i,j,k,...) = sum(j') t(i,j',k,...) * c(j',j)
    /// \endcode
    /// @param[in] t Tensor to transform (size of dimension to be transformed must match size of first dimension of \c c )
    /// @param[in] c Matrix used for the transformation
    /// @param[in] axis Dimension (or axis) to be transformed
    /// @result Returns a new, contiguous tensor
    template <class T, class Q>
    Tensor<TENSOR_RESULT_TYPE(T,Q)> transform_dir(const Tensor<T>& t, const Tensor<Q>& c, int axis) {
        if (axis == 0) {
            return inner(c,t,0,axis);
        }
        else if (axis == t.ndim()-1) {
            return inner(t,c,axis,0);
        }
        else {
            return copy(inner(t,c,axis,0).cycledim(1,axis, -1)); // Copy to make contiguous
        }
    }

    /// Returns a new deep copy of the transpose of the input tensor

    /// \ingroup tensor
    template <class T>
    Tensor<T> transpose(const Tensor<T>& t) {
        TENSOR_ASSERT(t.ndim() == 2, "transpose requires a matrix", t.ndim(), &t);
        return copy(t.swapdim(0,1));
    }

    /// Returns a new deep copy of the complex conjugate transpose of the input tensor

    /// \ingroup tensor
    template <class T>
    Tensor<T> conj_transpose(const Tensor<T>& t) {
        TENSOR_ASSERT(t.ndim() == 2, "conj_transpose requires a matrix", t.ndim(), &t);
        return conj(t.swapdim(0,1));
    }

    /// Indexing a non-constant tensor with slices returns a SliceTensor

    /// \ingroup tensor
    /// A slice tensor differs from a tensor only in that assignment
    /// causes the data to be copied rather than entire new copy
    /// generated.  You will usually not instantiate one except as a
    /// temporary produced by indexing a tensor with slice and then
    /// assigning it back to a tensor, or performing some other
    /// operation and discarding.
    template <class T> class SliceTensor : public Tensor<T> {
    private:
        SliceTensor<T>();

    public:
        SliceTensor(const Tensor<T>& t, const Slice s[])
            : Tensor<T>(const_cast<Tensor<T>&>(t)) //!!!!!!!!!!!
        {
            // C++ standard says class derived from parameterized base class cannot
            // directly access the base class elements ... must explicitly reference.

            long nd = 0, size=1;
            for (long i=0; i<t._ndim; ++i) {
                long start=s[i].start, end=s[i].end, step=s[i].step;
                //std::printf("%ld input start=%ld end=%ld step=%ld\n",
                //i, start, end, step);
                if (start < 0) start += this->_dim[i];
                if (end < 0) end += this->_dim[i];
                long len = end-start+1;
                if (step) len /= step;	// Rounds len towards zero

                // if input length is not exact multiple of step, round end towards start
                // for the same behaviour of for (i=start; i<=end; i+=step);
                end = start + (len-1)*step;

                //std::printf("%ld munged start=%ld end=%ld step=%ld len=%ld _dim=%ld\n",
                //		i, start, end, step, len, this->_dim[i]);

                TENSOR_ASSERT(start>=0 && start<this->_dim[i],"slice start invalid",start,this);
                TENSOR_ASSERT(end>=0 && end<this->_dim[i],"slice end invalid",end,this);
                TENSOR_ASSERT(len>0,"slice length must be non-zero",len,this);

                this->_p += start * t._stride[i];

                if (step) {
                    size *= len;
                    this->_dim[nd] = len;
                    this->_stride[nd] = step * t._stride[i];
                    ++nd;
                }
            }
            //For Python interface need to be able to return a scalar inside a tensor with nd=0
            //TENSOR_ASSERT(nd>0,"slicing produced a scalar, but cannot return one",nd,this);
            for (long i=nd; i<TENSOR_MAXDIM; ++i) { // So can iterate over missing dimensions
                this->_dim[i] = 1;
                this->_stride[i] = 0;
            }

            this->_ndim = nd;
            this->_size = size;
        }

        SliceTensor<T>& operator=(const SliceTensor<T>& t) {
            BINARY_OPTIMIZED_ITERATOR(T, (*this), const T, t, *_p0 = (T)(*_p1));
            return *this;
        }

        template <class Q>
        SliceTensor<T>& operator=(const SliceTensor<Q>& t) {
            BINARY_OPTIMIZED_ITERATOR(T, (*this), const Q, t, *_p0 = (T)(*_p1));
            return *this;
        }

        SliceTensor<T>& operator=(const Tensor<T>& t) {
            BINARY_OPTIMIZED_ITERATOR(T, (*this), const T, t, *_p0 = (T)(*_p1));
            return *this;
        }

        template <class Q>
        SliceTensor<T>& operator=(const Tensor<Q>& t) {
            BINARY_OPTIMIZED_ITERATOR(T, (*this), const Q, t, *_p0 = (T)(*_p1));
            return *this;
        }

        SliceTensor<T>& operator=(const T& t) {
            UNARY_OPTIMIZED_ITERATOR(T, (*this), *_p0 = t);
            return *this;
        }

        virtual ~SliceTensor() {};		// Tensor<T> destructor does enough
    };


    // Specializations for complex types
    template<> float_complex Tensor<float_complex>::min(long* ind) const ;
    template<> double_complex Tensor<double_complex>::min(long* ind) const ;
    template<> float_complex Tensor<float_complex>::max(long* ind) const ;
    template<> double_complex Tensor<double_complex>::max(long* ind) const ;

    // Stream stuff

    /// Print (for human consumption) a tensor to the stream

    /// \ingroup tensor
    template <class T>
    std::ostream& operator << (std::ostream& s, const Tensor<T>& t) {
        if (t.size() == 0) {
            s << "[empty tensor]\n";
            return s;
        }

        long maxdim = 0;
        long index_width = 0;
        for (int i = 0; i<(t.ndim()-1); ++i) {
            if (maxdim < t.dim(i)) maxdim = t.dim(i);
        }
        if (maxdim < 10)
            index_width = 1;
        else if (maxdim < 100)
            index_width = 2;
        else if (maxdim < 1000)
            index_width = 3;
        else if (maxdim < 10000)
            index_width = 4;
        else
            index_width = 6;

        std::ios::fmtflags oldflags = s.setf(std::ios::scientific);
        long oldprec = s.precision();
        long oldwidth = s.width();

        // C++ formatted IO is worse than Fortran !!
        for (TensorIterator<T> iter=t.unary_iterator(1,false,false); iter!=t.end(); ++iter) {
            const T* p = iter._p0;
            long inc = iter._s0;
            long dimj = iter.dimj;
            s.unsetf(std::ios::scientific);
            s << '[';
            for (long i=0; i<iter.ndim; ++i) {
                s.width(index_width);
                s << iter.ind[i];
                if (i != iter.ndim) s << ",";
            }
            s << "*]";
            s.setf(std::ios::scientific);
            for (long j=0; j<dimj; ++j, p+=inc) {
                s.precision(4);
                s.width(12);
                s << *p;
            }
            s.unsetf(std::ios::scientific);
            s << std::endl;
        }
        s.setf(oldflags);
        s.precision(oldprec);
        s.width(oldwidth);

        return s;
    }


    /// Outer product ... result(i,j,...,p,q,...) = left(i,k,...)*right(p,q,...)

    /// \ingroup tensor
    template <class T>
    Tensor<T> outer(const Tensor<T>& left, const Tensor<T>& right) {
        long nd = left.ndim() + right.ndim();
        TENSOR_ASSERT(nd <= TENSOR_MAXDIM,"too many dimensions in result",
                      nd,0);
        long d[TENSOR_MAXDIM];
        for (long i=0; i<left.ndim(); ++i) d[i] = left.dim(i);
        for (long i=0; i<right.ndim(); ++i) d[i+left.ndim()] = right.dim(i);
        Tensor<T> result(nd,d,false);
        T* ptr = result.ptr();

        TensorIterator<T> iter=right.unary_iterator(1,false,true);
        for (TensorIterator<T> p=left.unary_iterator(); p!=left.end(); ++p) {
            T val1 = *p;
            // Cannot reorder dimensions, but can fuse contiguous dimensions
            for (iter.reset(); iter._p0; ++iter) {
                long dimj = iter.dimj;
                T* _p0 = iter._p0;
                long Tstride = iter._s0;
                for (long _j=0; _j<dimj; ++_j, _p0+=Tstride) {
                    *ptr++ = val1 * (*_p0);
                }
            }
        }

        return result;
    }


    /// Inner product ... result(i,j,...,p,q,...) = sum(z) left(i,j,...,z)*right(z,p,q,...)

    /// \ingroup tensor
    /// By default it contracts the last dimension of the left tensor and
    /// the first dimension of the right tensor.  These defaults can be
    /// changed by specifying \c k0 and \c k1 , the index to contract in
    /// the left and right side tensors, respectively.  The defaults
    /// correspond to (\c k0=-1 and \c k1=0 ).
    template <class T, class Q>
    Tensor<TENSOR_RESULT_TYPE(T,Q)> inner(const Tensor<T>& left, const Tensor<Q>& right,
                                          long k0=-1, long k1=0) {
        if (k0 < 0) k0 += left.ndim();
        if (k1 < 0) k1 += right.ndim();
        long nd = left.ndim() + right.ndim() - 2;
        TENSOR_ASSERT(nd!=0, "result is a scalar but cannot return one ... use dot",
                      nd, &left);
        TENSOR_ASSERT(left.dim(k0) == right.dim(k1),"common index must be same length",
                      right.dim(k1), &left);

        TENSOR_ASSERT(nd > 0 && nd <= TENSOR_MAXDIM,
                      "invalid number of dimensions in the result", nd,0);

        long d[TENSOR_MAXDIM];

        long base=0;
        for (long i=0; i<k0; ++i) d[i] = left.dim(i);
        for (long i=k0+1; i<left.ndim(); ++i) d[i-1] = left.dim(i);
        base = left.ndim()-1;
        for (long i=0; i<k1; ++i) d[i+base] = right.dim(i);
        base--;
        for (long i=k1+1; i<right.ndim(); ++i) d[i+base] = right.dim(i);

        Tensor<TENSOR_RESULT_TYPE(T,Q)> result(nd,d);

        inner_result(left,right,k0,k1,result);

        return result;
    }

    /// Accumulate inner product into user provided, contiguous, correctly sized result tensor

    /// \ingroup tensor
    /// This routine may be used to optimize away the tensor constructor
    /// of the result tensor in inner loops when the result tensor may be
    /// reused or accumulated into.  If the user calls this routine
    /// directly very little checking is done since it is intended as an
    /// optimization for small tensors.  As far as the result goes, the
    /// caller is completely responsible for providing a contiguous tensor
    /// that has the correct dimensions and is appropriately initialized.
    /// The inner product is accumulated into result.
    template <class T, class Q>
    void inner_result(const Tensor<T>& left, const Tensor<Q>& right,
                      long k0, long k1, Tensor< TENSOR_RESULT_TYPE(T,Q) >& result) {

        typedef TENSOR_RESULT_TYPE(T,Q) resultT;
        // Need to include explicit optimizations for common special cases
        // E.g., contiguous, matrix-matrix, and 3d-tensor*matrix

        resultT* ptr = result.ptr();

        if (k0 < 0) k0 += left.ndim();
        if (k1 < 0) k1 += right.ndim();

        if (left.iscontiguous() && right.iscontiguous()) {
            if (k0==0 && k1==0) {
                // c[i,j] = a[k,i]*b[k,j] ... collapsing extra indices to i & j
                long dimk = left.dim(k0);
                long dimj = right.stride(0);
                long dimi = left.stride(0);
                ::mTxm(dimi,dimj,dimk,ptr,left.ptr(),right.ptr());
                return;
            }
            else if (k0==(left.ndim()-1) && k1==(right.ndim()-1)) {
                // c[i,j] = a[i,k]*b[j,k] ... collapsing extra indices to i & j
                long dimk = left.dim(k0);
                long dimi = left.size()/dimk;
                long dimj = right.size()/dimk;
                ::mxmT(dimi,dimj,dimk,ptr,left.ptr(),right.ptr());
                return;
            }
            else if (k0==0 && k1==(right.ndim()-1)) {
                // c[i,j] = a[k,i]*b[j,k] ... collapsing extra indices to i & j
                long dimk = left.dim(k0);
                long dimi = left.stride(0);
                long dimj = right.size()/dimk;
                ::mTxmT(dimi,dimj,dimk,ptr,left.ptr(),right.ptr());
                return;
            }
            else if (k0==(left.ndim()-1) && k1==0) {
                // c[i,j] = a[i,k]*b[k,j] ... collapsing extra indices to i & j
                long dimk = left.dim(k0);
                long dimi = left.size()/dimk;
                long dimj = right.stride(0);
                ::mxm(dimi,dimj,dimk,ptr,left.ptr(),right.ptr());
                return;
            }
        }

        long dimj = left.dim(k0);
        TensorIterator<Q> iter1=right.unary_iterator(1,false,false,k1);

        for (TensorIterator<T> iter0=left.unary_iterator(1,false,false,k0);
                iter0._p0; ++iter0) {
            T* restrict xp0 = iter0._p0;
            long s0 = iter0._s0;
            for (iter1.reset(); iter1._p0; ++iter1) {
                T* restrict p0 = xp0;
                Q* restrict p1 = iter1._p0;
                long s1 = iter1._s0;
                resultT sum = 0;
                for (long j=0; j<dimj; ++j,p0+=s0,p1+=s1) {
                    sum += (*p0) * (*p1);
                }
                *ptr++ += sum;
            }
        }
    }

    /// Transform all dimensions of the tensor t by the matrix c

    /// \ingroup tensor
    /// Often used to transform all dimensions from one basis to another
    /// \code
    /// result(i,j,k...) <-- sum(i',j', k',...) t(i',j',k',...) c(i',i) c(j',j) c(k',k) ...
    /// \endcode
    /// The input dimensions of \c t must all be the same and agree with
    /// the first dimension of \c c .  The dimensions of \c c may differ in
    /// size.  If the dimensions of \c c are the same, and the operation
    /// is being performed repeatedly, then you might consider calling \c
    /// fast_transform instead which enables additional optimizations and
    /// can eliminate all constructor overhead and improve cache locality.
    ///
    template <class T, class Q>
    Tensor<TENSOR_RESULT_TYPE(T,Q)> transform(const Tensor<T>& t, const Tensor<Q>& c) {
        typedef TENSOR_RESULT_TYPE(T,Q) resultT;
        TENSOR_ASSERT(c.ndim() == 2,"second argument must be a matrix",c.ndim(),&c);
        if (c.dim(0)==c.dim(1) && t.iscontiguous() && c.iscontiguous()) {
            Tensor<resultT> result(t.ndim(),t.dims(),false);
            Tensor<resultT> work(t.ndim(),t.dims(),false);
            return fast_transform(t, c, result, work);
        }
        else {
            Tensor<resultT> result = t;
            for (long i=0; i<t.ndim(); ++i) {
                result = inner(result,c,0,0);
            }
            return result;
        }
    }

    /// Transform all dimensions of the tensor t by distinct matrices c

    /// \ingroup tensor
    /// Similar to transform but each dimension is transformed with a
    /// distinct matrix.
    /// \code
    /// result(i,j,k...) <-- sum(i',j', k',...) t(i',j',k',...) c[0](i',i) c[1](j',j) c[2](k',k) ...
    /// \endcode
    /// The first dimension of the matrices c must match the corresponding
    /// dimension of t.
    template <class T, class Q>
    Tensor<TENSOR_RESULT_TYPE(T,Q)> general_transform(const Tensor<T>& t, const Tensor<Q> c[]) {
        typedef TENSOR_RESULT_TYPE(T,Q) resultT;
        Tensor<resultT> result = t;
        for (long i=0; i<t.ndim(); ++i) {
            result = inner(result,c[i],0,0);
        }
        return result;
    }

    /// Restricted but heavily optimized form of transform()

    /// \ingroup tensor
    /// Both dimensions of \c c must be the same and match all dimensions
    /// of the input tensor \c t.  All tensors must be contiguous.
    ///
    /// Performs the same operation as \c transform but it requires
    /// that the caller pass in workspace and a preallocated result,
    /// hoping that that both can be reused.  If the result and
    /// workspace are reused between calls, then no tensor
    /// constructors need be called and cache locality should be
    /// improved.  By passing in the workspace, this routine is kept
    /// thread safe.
    ///
    /// The input, result and workspace tensors must be distinct.
    ///
    /// All input tensors must be contiguous and fastest execution
    /// will result if all dimensions are even and data is aligned on
    /// 16-byte boundaries.  The workspace and the result must be of
    /// the same size as the input \c t .  The result tensor need not
    /// be initialized before calling fast_transform.
    ///
    /// \code
    ///     result(i,j,k,...) <-- sum(i',j', k',...) t(i',j',k',...)  c(i',i) c(j',j) c(k',k) ...
    /// \endcode
    ///
    /// The input dimensions of \c t must all be the same .
    template <class T, class Q>
    Tensor< TENSOR_RESULT_TYPE(T,Q) >& fast_transform(const Tensor<T>& t, const Tensor<Q>& c,  Tensor< TENSOR_RESULT_TYPE(T,Q) >& result,
            Tensor< TENSOR_RESULT_TYPE(T,Q) >& workspace) {
        typedef  TENSOR_RESULT_TYPE(T,Q) resultT;
        const Q *pc=c.ptr();
        resultT *t0=workspace.ptr(), *t1=result.ptr();
        if (t.ndim()&1) {
            t0 = result.ptr();
            t1 = workspace.ptr();
        }

        long dimj = c.dim(1);
        long dimi = 1;
        for (int n=1; n<t.ndim(); ++n) dimi *= dimj;
        long nij = dimi*dimj;


        if (IS_ODD(dimi) || IS_ODD(dimj) ||
                IS_UNALIGNED(pc) || IS_UNALIGNED(t0) || IS_UNALIGNED(t1)) {
            for (long i=0; i<nij; ++i) t0[i] = 0.0;
            mTxm(dimi, dimj, dimj, t0, t.ptr(), pc);
            for (int n=1; n<t.ndim(); ++n) {
                for (long i=0; i<nij; ++i) t1[i] = 0.0;
                mTxm(dimi, dimj, dimj, t1, t0, pc);
                std::swap(t0,t1);
            }
        }
        else {
            mTxmq(dimi, dimj, dimj, t0, t.ptr(), pc);
            for (int n=1; n<t.ndim(); ++n) {
                mTxmq(dimi, dimj, dimj, t1, t0, pc);
                std::swap(t0,t1);
            }
        }

        return result;
    }

    /// Return a new tensor holding the absolute value of each element of t

    /// \ingroup tensor
    template <class T>
    Tensor< typename Tensor<T>::scalar_type > abs(const Tensor<T>& t) {
        typedef typename Tensor<T>::scalar_type scalar_type;
        Tensor<scalar_type> result(t.ndim(),t.dims(),false);
        BINARY_OPTIMIZED_ITERATOR(scalar_type,result,const T,t,*_p0 = std::abs(*_p1));
        return result;
    }

    /// Return a new tensor holding the argument of each element of t (complex types only)

    /// \ingroup tensor
    template <class T>
    Tensor< typename Tensor<T>::scalar_type > arg(const Tensor<T>& t) {
        typedef typename Tensor<T>::scalar_type scalar_type;
        Tensor<scalar_type> result(t.ndim(),t.dims(),false);
        BINARY_OPTIMIZED_ITERATOR(scalar_type,result,T,t,*_p0 = std::arg(*_p1));
        return result;
    }

    /// Return a new tensor holding the real part of each element of t (complex types only)

    /// \ingroup tensor
    template <class T>
    Tensor< typename Tensor<T>::scalar_type > real(const Tensor<T>& t) {
        typedef typename Tensor<T>::scalar_type scalar_type;
        Tensor<scalar_type> result(t.ndim(),t.dims(),false);
        BINARY_OPTIMIZED_ITERATOR(scalar_type,result,const T,t,*_p0 = std::real(*_p1));
        return result;
    }

    /// Return a new tensor holding the imaginary part of each element of t (complex types only)

    /// \ingroup tensor
    template <class T>
    Tensor< typename Tensor<T>::scalar_type > imag(const Tensor<T>& t) {
        typedef typename Tensor<T>::scalar_type scalar_type;
        Tensor<scalar_type> result(t.ndim(),t.dims(),false);
        BINARY_OPTIMIZED_ITERATOR(scalar_type,result,const T,t,*_p0 = std::imag(*_p1));
        return result;
    }

    /// Returns a new deep copy of the complex conjugate of the input tensor (complex types only)

    /// \ingroup tensor
    template <class T>
    Tensor<T> conj(const Tensor<T>& t) {
        Tensor<T> result(t.ndim(),t.dims(),false);
        BINARY_OPTIMIZED_ITERATOR(T,result,const T,t,*_p0 = std::conj(*_p1));
        return result;
    }
}

#endif // MADNESS_TENSOR_TENSOR_H__INCLUDED
