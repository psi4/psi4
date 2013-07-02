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


  $Id: basetensor.cc 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/


#include <iostream>

#include <algorithm>

#include "basetensor.h"

#include "tensorexcept.h"

/// \file basetensor.cc
/// \brief Implements BaseTensor

namespace madness {

#ifdef TENSOR_INSTANCE_COUNT
    MADATOMIC_INT BaseTensor::instance_count;
#endif


    /// Reshape the size and number of dimensions.

    /// Modifies the current tensor to have the number and size of
    /// dimensions as described in the \c vector \c d .  The total number
    /// of elements must be the same before and after, and the current
    /// tensor must be contiguous.
    void BaseTensor::reshape_inplace(int nd, const long* d) {
        TENSOR_ASSERT(iscontiguous(),
                      "cannot reshape non-contiguous tensor ... consider fuse/splitdim",
                      0,this);
        long newsize=1;
        for (long i=0; i<nd; ++i) newsize *= d[i];
        TENSOR_ASSERT(_size == newsize,"old and new sizes do not match",_size,this);
        set_dims_and_size(nd,&(d[0]));
    }

    /// Reshape the size and number of dimensions.

    /// Modifies the current tensor to have the number and size of
    /// dimensions as described in the \c vector \c d .  The total number
    /// of elements must be the same before and after, and the current
    /// tensor must be contiguous.
    void BaseTensor::reshape_inplace(const std::vector<long>& d) {
        reshape_inplace(d.size(), &d[0]);
    }

    /// Reshape the current tensor to be the same size and 1-d. It must be contiguous.
    void BaseTensor::flat_inplace() {
        TENSOR_ASSERT(iscontiguous(),"not contiguous",0,this);
        long d[] = {_size};
        set_dims_and_size(1,d);
    }

    /// Split dimension i in two ... the product of the new dimensions must match the old.
    void BaseTensor::splitdim_inplace(long i, long dimi0, long dimi1) {
        if (i < 0) i += _ndim;
        TENSOR_ASSERT(i>=0 && i<_ndim, "invalid dimension", i, this);
        TENSOR_ASSERT(dimi0*dimi1 == _dim[i], "before & after sizes do not match",
                      _dim[i], this);
        TENSOR_ASSERT(_ndim+1 <= TENSOR_MAXDIM, "resulting tensor has too many dimensions",
                      _ndim+1, this);
        for (long j=_ndim-1; j>i; --j) {
            _dim[j+1] = _dim[j];
            _stride[j+1] = _stride[j];
        }
        _dim[i+1] = dimi1;
        _stride[i+1] = _stride[i];
        _dim[i] = dimi0;
        _stride[i] *= dimi1;
        ++_ndim;
    }

    /// Fuse the contiguous dimensions i and i+1.
    void BaseTensor::fusedim_inplace(long i) { // fuse i,i+1 --> i
        if (i < 0) i += _ndim;
        TENSOR_ASSERT(i>=0 && i<(_ndim-1) && _ndim>1,"invalid dimension",i,this);
        TENSOR_ASSERT(_stride[i] == _dim[i+1]*_stride[i+1],"dimensions are not contiguous",
                      i, this);
        _dim[i] *= _dim[i+1];
        _stride[i] = _stride[i+1];
        for (long j=i+1; j<=_ndim-1; ++j) {
            _dim[j] = _dim[j+1];
            _stride[j] = _stride[j+1];
        }
        _ndim--;
        _dim[_ndim] = 1;		// So can iterate over missing dimensions
        _stride[_ndim] = 0;
    }

    /// Swap the dimensions i and j.
    void BaseTensor::swapdim_inplace(long i, long j) {
        if (i < 0) i += _ndim;
        if (j < 0) j += _ndim;
        TENSOR_ASSERT(i>=0 && i<_ndim,"invalid dimension i",i,this);
        TENSOR_ASSERT(j>=0 && j<_ndim,"invalid dimension j",j,this);
        std::swap<long>(_dim[i],_dim[j]);
        std::swap<long>(_stride[i],_stride[j]);
    }

    /// Cyclic shift by nshift places of the inclusive range of dimensions [start,....,end]
    void BaseTensor::cycledim_inplace(long nshift, long start, long end) {
        long ndshift, dimtmp[TENSOR_MAXDIM], _stridetmp[TENSOR_MAXDIM];
        if (start < 0) start += _ndim; // Same convention for -ve as in Slice
        if (end < 0) end += _ndim;
        TENSOR_ASSERT(start>=0 && start<_ndim,"invalid start dimension",start,this);
        TENSOR_ASSERT(end>=0 && end>=start,"invalid end dimension",end,this);

        ndshift = end - start + 1;
        for (long i=start; i<=end; ++i) {
            dimtmp[i] = _dim[i];
            _stridetmp[i] = _stride[i];
        }
        for (long i=end; i>=start; --i) {
            long j = i + nshift;
            while (j > end) j -= ndshift;
            while (j < start) j += ndshift;
            _dim[j] = dimtmp[i];
            _stride[j] = _stridetmp[i];
        }
    }

    /// General permuation of the dimensions
    void BaseTensor::mapdim_inplace(const std::vector<long>& map) {
        TENSOR_ASSERT(_ndim == (int) map.size(),"map[] must include all dimensions",map.size(),this);
        long tmpd[TENSOR_MAXDIM], tmps[TENSOR_MAXDIM];
        for (long i=0; i<_ndim; ++i) {
            tmpd[map[i]] = _dim[i];
            tmps[map[i]] = _stride[i];
        }
        for (long i=0; i<_ndim; ++i) {
            _dim[i] = tmpd[i];
            _stride[i] = tmps[i];
        }
    }
}
