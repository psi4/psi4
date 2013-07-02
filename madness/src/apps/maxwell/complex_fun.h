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
  
  $Id: complex_fun.h 1602 2009-12-27 19:53:06Z rjharrison $
*/
#ifndef __complex_fun__
#define __complex_fun__

template <typename Q, int NDIM>
struct real_op
{
  typedef typename TensorTypeData<Q>::scalar_type resultT;
  Tensor<resultT> operator()(const Key<NDIM>& key, const Tensor<Q>& t) const
  {
    Tensor<resultT> result(t.ndim, t.dim);
    BINARY_OPTIMIZED_ITERATOR(Q, t, resultT, result, *_p1 = real(*_p0););
    return result;
  }
  template <typename Archive>
  void serialize(Archive& ar) {}
};

template<typename Q, int NDIM>
Function<typename TensorTypeData<Q>::scalar_type,NDIM> real(const Function<Q,NDIM>& func)
{
  return unary_op_coeffs(func, real_op<Q,NDIM>());
}

inline void ln(const Key<3> &key, Tensor<std::complex<double> > &t) {
	UNARY_OPTIMIZED_ITERATOR(std::complex<double>, t,
		*_p0 = log(*_p0));
}

#endif
