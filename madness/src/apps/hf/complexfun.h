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

  $Id: complexfun.h 2329 2011-05-27 16:48:59Z wsttiger@gmail.com $
*/
/*
 * complexfun.h
 *
 *  Created on: Jan 28, 2009
 *      Author: eh7
 */

#ifndef COMPLEXFUN_H_
#define COMPLEXFUN_H_

#include <mra/mra.h>

namespace madness {

//***************************************************************************
double abs(double x) {return fabs(x);}
//***************************************************************************

//***************************************************************************
double real(double x) {return x;}
//***************************************************************************

//***************************************************************************
double imag(double x) {return 0.0;}
//***************************************************************************

////***************************************************************************
//double conj(double x) {return 0.0;}
////***************************************************************************

//***************************************************************************
template <typename Q>
Tensor< std::complex<Q> > tensor_real2complex(const Tensor<Q>& r)
{
  Tensor< std::complex<Q> > c(r.ndim(), r.dims());
  BINARY_OPTIMIZED_ITERATOR(const Q, r, std::complex<Q>, c, *_p1 = std::complex<Q>(*_p0,0.0););
  return c;
}
//***************************************************************************

//***************************************************************************
template <typename Q>
Tensor<Q> tensor_real(const Tensor< std::complex<Q> >& c)
{
  Tensor<Q> r(c.ndim(), c.dims());
  BINARY_OPTIMIZED_ITERATOR(Q, r, const std::complex<Q>, c, *_p0 = real(*_p1););
  return r;
}
//***************************************************************************

//***************************************************************************
template <typename Q>
Tensor<Q> tensor_imag(const Tensor< std::complex<Q> >& c)
{
  Tensor<Q> r(c.ndim(), c.dims());
  BINARY_OPTIMIZED_ITERATOR(Q, r, const std::complex<Q>, c, *_p0 = imag(*_p1););
  return r;
}
//***************************************************************************

//***************************************************************************
template <typename Q>
Tensor<Q> tensor_abs(const Tensor< std::complex<Q> >& c)
{
  Tensor<Q> r(c.ndim(), c.dims());
  BINARY_OPTIMIZED_ITERATOR(Q, r, const std::complex<Q>, c, *_p0 = abs(*_p1););
  return r;
}
//***************************************************************************

//***************************************************************************
template <typename Q, int NDIM>
struct abs_square_op
{
  typedef typename TensorTypeData<Q>::scalar_type resultT;
  Tensor<resultT> operator()(const Key<NDIM>& key, const Tensor<Q>& t) const
  {
    Tensor<resultT> result(t.ndim(), t.dims());
    BINARY_OPTIMIZED_ITERATOR(const Q, t, resultT, result, resultT d = abs(*_p0); *_p1 = d*d);
    return result;
  }
  template <typename Archive>
  void serialize(Archive& ar) {}
};
//***************************************************************************

//***************************************************************************
template<typename Q, int NDIM>
Function<typename TensorTypeData<Q>::scalar_type,NDIM> abs_square(const Function<Q,NDIM>& func)
{
  return unary_op(func, abs_square_op<Q,NDIM>());
}
//***************************************************************************

//***************************************************************************
template <typename Q, int NDIM>
struct real_op
{
  typedef typename TensorTypeData<Q>::scalar_type resultT;
  Tensor<resultT> operator()(const Key<NDIM>& key, const Tensor<Q>& t) const
  {
    Tensor<resultT> result(t.ndim(), t.dims());
    BINARY_OPTIMIZED_ITERATOR(const Q, t, resultT, result, *_p1 = real(*_p0););
    return result;
  }
  template <typename Archive>
  void serialize(Archive& ar) {}
};
//***************************************************************************

//***************************************************************************
template<typename Q, int NDIM>
Function<typename TensorTypeData<Q>::scalar_type,NDIM> real(const Function<Q,NDIM>& func)
{
  return unary_op_coeffs(func, real_op<Q,NDIM>());
}
//***************************************************************************

//***************************************************************************
template <typename Q, int NDIM>
struct imag_op
{
  typedef typename TensorTypeData<Q>::scalar_type resultT;
  Tensor<resultT> operator()(const Key<NDIM>& key, const Tensor<Q>& t) const
  {
    Tensor<resultT> result(t.ndim(), t.dims());
    BINARY_OPTIMIZED_ITERATOR(const Q, t, resultT, result, *_p1 = imag(*_p0););
    return result;
  }
  template <typename Archive>
  void serialize(Archive& ar) {}
};
//***************************************************************************

//***************************************************************************
template<typename Q, int NDIM>
Function<typename TensorTypeData<Q>::scalar_type,NDIM> imag(const Function<Q,NDIM>& func)
{
  return unary_op_coeffs(func, imag_op<Q,NDIM>());
}
//***************************************************************************

//***************************************************************************
template <typename Q, int NDIM>
struct abs_op
{
  typedef typename TensorTypeData<Q>::scalar_type resultT;
  Tensor<resultT> operator()(const Key<NDIM>& key, const Tensor<Q>& t) const
  {
    Tensor<resultT> result(t.ndim(), t.dims());
    BINARY_OPTIMIZED_ITERATOR(const Q, t, resultT, result, *_p1 = abs(*_p0););
    return result;
  }
  template <typename Archive>
  void serialize(Archive& ar) {}
};
//***************************************************************************

//***************************************************************************
template<typename Q, int NDIM>
Function<typename TensorTypeData<Q>::scalar_type,NDIM> abs(const Function<Q,NDIM>& func)
{
  return unary_op_coeffs(func, abs_op<Q,NDIM>());
}
//***************************************************************************

//***************************************************************************
template <typename Q, int NDIM>
struct conj_op
{
  typedef Q resultT;
  Tensor<resultT> operator()(const Key<NDIM>& key, const Tensor<Q>& t) const
  {
    Tensor<resultT> result(t.ndim(), t.dims());
    BINARY_OPTIMIZED_ITERATOR(const Q, t, resultT, result, *_p1 = conj(*_p0););
    return result;
  }
  template <typename Archive>
  void serialize(Archive& ar) {}
};
//***************************************************************************

////***************************************************************************
//template<typename Q, int NDIM>
//Function<Q,NDIM> conj(const Function<Q,NDIM>& func)
//{
//  return unary_op_coeffs(func, conj_op<Q,NDIM>());
//}
////***************************************************************************

//***************************************************************************
template <typename Q, int NDIM>
struct function_real2complex_op
{
  typedef std::complex<Q> resultT;
  Tensor<resultT> operator()(const Key<NDIM>& key, const Tensor<Q>& t) const
  {
    Tensor<resultT> result(t.ndim(), t.dims());
    BINARY_OPTIMIZED_ITERATOR(const Q, t, resultT, result, *_p1 = resultT(*_p0,0.0););
    return result;
  }
  template <typename Archive>
  void serialize(Archive& ar) {}
};
//***************************************************************************

//***************************************************************************
template <typename Q, int NDIM>
Function<std::complex<Q>,NDIM> function_real2complex(const Function<Q,NDIM>& r)
{
  return unary_op_coeffs(r, function_real2complex_op<Q,NDIM>());
}
//***************************************************************************


////***************************************************************************
//template <typename Q, int NDIM>
//bool is_real(const Function<std::complex<Q>,NDIM>& f)
//{
//  Function<Q,NDIM> fim = imag(f);
//  return (fim.norm2() < 1e-8) ? true : false;
//}
////***************************************************************************
//
////***************************************************************************
//template <typename Q, int NDIM>
//bool is_imag(const Function<std::complex<Q>,NDIM>& f)
//{
//  Function<Q,NDIM> fre = real(f);
//  return (fre.norm2() < 1e-8) ? true : false;
//}
////***************************************************************************

}
#endif /* COMPLEXFUN_H_ */
