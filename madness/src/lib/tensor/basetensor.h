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


  $Id: basetensor.h 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/


#ifndef MADNESS_TENSOR_BASETENSOR_H__INCLUDED
#define MADNESS_TENSOR_BASETENSOR_H__INCLUDED

/// \file basetensor.h
/// \brief Declares BaseTensor

#include <madness_config.h>
#include <tensor/slice.h>
#include <tensor/tensor_macros.h>

#ifdef TENSOR_INSTANCE_COUNT
#include <world/atomicint.h>
#endif

namespace madness {
    /*!
      \ingroup tensor
      \brief The base class for tensors defines generic capabilities.

      The base class manages the size, dimension and
      stride information, and provides operations to manipulate
      them.

      It also provides methods for type-safe operation on tensors using
      just the base class pointers. This interface is primarily useful
      only to the interface to Python, since Python is largely neutral
      to (and ignorant of) the type.  These are still being
      re-implemented after the big clean up.

      Since the base tensor class is virtual, you cannot have an
      instance of it.  Thus, in addition to methods that return information
      or perform checks, there are two types of base tensor
      operations.
        - Inplace operations change \c *this , and return \c void .
        - Operations that must return a new tensor return a pointer to a tensor
          allocated with \c new on the heap.  The caller is responsible for
          eventually freeing the memory using \c delete .
    */
    class BaseTensor {
    private:
#ifdef TENSOR_INSTANCE_COUNT
        static madness::AtomicInt instance_count; ///< For debug, count total# instances
#endif

    protected:

        long _size;			///< Number of elements in the tensor
        long _ndim;			///< Number of dimensions (-1=invalid; 0=scalar; >0=tensor)
        long _id; 			///< Id from TensorTypeData<T> in type_data.h
        long _dim[TENSOR_MAXDIM];	///< Size of each dimension
        long _stride[TENSOR_MAXDIM];     ///< Increment between elements in each dimension

        void set_dims_and_size(long nd, const long d[]) {
            _ndim = nd;
            _size = 1;
            if (_ndim < 0) _size=0;
            for (long i=_ndim-1; i>=0; --i) {
                _dim[i] = d[i];
                _stride[i] = _size;
                _size *= d[i];
            }
            for (long i=std::max(_ndim,0L); i<TENSOR_MAXDIM; ++i) { // So can iterate over missing dimensions
                _dim[i] = 1;
                _stride[i] = 0;
            }
        }

    public:

        BaseTensor() : _size(0), _ndim(-1) {
#ifdef TENSOR_INSTANCE_COUNT
            instance_count++;
#endif
        }

        virtual ~BaseTensor() {
#ifdef TENSOR_INSTANCE_COUNT
            instance_count--;
#endif
        }

        /// Returns the count of all current instances of tensors & slice tensors of all types.
        static inline int get_instance_count() {
#ifdef TENSOR_INSTANCE_COUNT
            return instance_count;
#else
            return 0;
#endif
        }

        /// Returns the number of elements in the tensor
        long size() const {return _size;}

        /// Returns the typeid of the tensor (c.f., \c TensorTypeData<T> )
        long id() const {return _id;}

        /// Returns the number of dimensions in the tensor
        long ndim() const {return _ndim;}

        /// Returns the size of dmension \c i
        long dim(int i) const {return _dim[i];}

        /// Returns the stride associated with dimension \c i
        long stride(int i) const {return _stride[i];}

        /// Returns the array of tensor dimensions
        const long* dims() const {return _dim;}

        /// Returns the array of tensor strides
        const long* strides() const {return _stride;}

        /// Returns true if this and *t are the same shape and size
        bool conforms(const BaseTensor *t) const {
            if (_ndim != t->_ndim) return false;
            for (long i=0; i<_ndim; ++i) {
                if (_dim[i] != t->_dim[i]) return false;
            }
            return true;
        }

        /// Returns true if the tensor refers to contiguous memory locations.
        bool iscontiguous() const {
            if (_size <= 0) return true;
            long sz = 1;
            for (long i=_ndim-1; i>=0; --i) {
                if (_stride[i] != sz) return false;
                sz *= _dim[i];
            }
            return true;
        }

    protected:

        /// Reshapes the tensor inplace
        void reshape_inplace(const std::vector<long>& d);

        /// Reshapes the tensor inplace
        void reshape_inplace(int ndimnew, const long* d);

        /// Reshapes the tensor inplace into 1D
        void flat_inplace();

        /// Splits dimension \c i
        void splitdim_inplace(long i, long dimi0, long dimi1);

        /// Fuses dimensions \c i and \c i+1
        void fusedim_inplace(long i);

        /// Swaps the dimensions
        void swapdim_inplace(long i, long j);

        /// Cyclic shift of dimensions
        void cycledim_inplace(long shift, long start, long end);

        /// General permutation of dimensions
        void mapdim_inplace(const std::vector<long>& map);
    };

}

#endif // MADNESS_TENSOR_BASETENSOR_H__INCLUDED
