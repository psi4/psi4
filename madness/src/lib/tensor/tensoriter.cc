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


  $Id: tensoriter.cc 1602 2009-12-27 19:53:06Z rjharrison $
*/


#include <tensor/tensor.h>

#include <iostream>
#include <algorithm>
#include <complex>

#include <cmath>


/// \file tensoriter.cc
/// \brief Implements TensorIterator

// These here only to instantiate specializations
typedef std::complex<float> float_complex;
typedef std::complex<double> double_complex;


namespace madness {


    //#include "tensoriter_spec.h"

    /*
    template class TensorIterator<double,double,double>;
    template class TensorIterator<double,double,double_complex>;
    template class TensorIterator<double,double,float>;
    template class TensorIterator<double,double,float_complex>;
    template class TensorIterator<double,double,long>;
    template class TensorIterator<double,double_complex,double>;
    template class TensorIterator<double,double_complex,double_complex>;
    template class TensorIterator<double,double_complex,float>;
    template class TensorIterator<double,double_complex,float_complex>;
    template class TensorIterator<double,double_complex,long>;
    template class TensorIterator<double,float,double>;
    template class TensorIterator<double,float,double_complex>;
    template class TensorIterator<double,float,float>;
    template class TensorIterator<double,float,float_complex>;
    template class TensorIterator<double,float,long>;
    template class TensorIterator<double,float_complex,double>;
    template class TensorIterator<double,float_complex,double_complex>;
    template class TensorIterator<double,float_complex,float>;
    template class TensorIterator<double,float_complex,float_complex>;
    template class TensorIterator<double,float_complex,long>;
    template class TensorIterator<double,long,double>;
    template class TensorIterator<double,long,double_complex>;
    template class TensorIterator<double,long,float>;
    template class TensorIterator<double,long,float_complex>;
    template class TensorIterator<double,long,long>;
    template class TensorIterator<double_complex,double,double>;
    template class TensorIterator<double_complex,double,double_complex>;
    template class TensorIterator<double_complex,double,float>;
    template class TensorIterator<double_complex,double,float_complex>;
    template class TensorIterator<double_complex,double,long>;
    template class TensorIterator<double_complex,double_complex,double>;
    template class TensorIterator<double_complex,double_complex,double_complex>;
    template class TensorIterator<double_complex,double_complex,float>;
    template class TensorIterator<double_complex,double_complex,float_complex>;
    template class TensorIterator<double_complex,double_complex,long>;
    template class TensorIterator<double_complex,float,double>;
    template class TensorIterator<double_complex,float,double_complex>;
    template class TensorIterator<double_complex,float,float>;
    template class TensorIterator<double_complex,float,float_complex>;
    template class TensorIterator<double_complex,float,long>;
    template class TensorIterator<double_complex,float_complex,double>;
    template class TensorIterator<double_complex,float_complex,double_complex>;
    template class TensorIterator<double_complex,float_complex,float>;
    template class TensorIterator<double_complex,float_complex,float_complex>;
    template class TensorIterator<double_complex,float_complex,long>;
    template class TensorIterator<double_complex,long,double>;
    template class TensorIterator<double_complex,long,double_complex>;
    template class TensorIterator<double_complex,long,float>;
    template class TensorIterator<double_complex,long,float_complex>;
    template class TensorIterator<double_complex,long,long>;
    template class TensorIterator<float,double,double>;
    template class TensorIterator<float,double,double_complex>;
    template class TensorIterator<float,double,float>;
    template class TensorIterator<float,double,float_complex>;
    template class TensorIterator<float,double,long>;
    template class TensorIterator<float,double_complex,double>;
    template class TensorIterator<float,double_complex,double_complex>;
    template class TensorIterator<float,double_complex,float>;
    template class TensorIterator<float,double_complex,float_complex>;
    template class TensorIterator<float,double_complex,long>;
    template class TensorIterator<float,float,double>;
    template class TensorIterator<float,float,double_complex>;
    template class TensorIterator<float,float,float>;
    template class TensorIterator<float,float,float_complex>;
    template class TensorIterator<float,float,long>;
    template class TensorIterator<float,float_complex,double>;
    template class TensorIterator<float,float_complex,double_complex>;
    template class TensorIterator<float,float_complex,float>;
    template class TensorIterator<float,float_complex,float_complex>;
    template class TensorIterator<float,float_complex,long>;
    template class TensorIterator<float,long,double>;
    template class TensorIterator<float,long,double_complex>;
    template class TensorIterator<float,long,float>;
    template class TensorIterator<float,long,float_complex>;
    template class TensorIterator<float,long,long>;
    template class TensorIterator<float_complex,double,double>;
    template class TensorIterator<float_complex,double,double_complex>;
    template class TensorIterator<float_complex,double,float>;
    template class TensorIterator<float_complex,double,float_complex>;
    template class TensorIterator<float_complex,double,long>;
    template class TensorIterator<float_complex,double_complex,double>;
    template class TensorIterator<float_complex,double_complex,double_complex>;
    template class TensorIterator<float_complex,double_complex,float>;
    template class TensorIterator<float_complex,double_complex,float_complex>;
    template class TensorIterator<float_complex,double_complex,long>;
    template class TensorIterator<float_complex,float,double>;
    template class TensorIterator<float_complex,float,double_complex>;
    template class TensorIterator<float_complex,float,float>;
    template class TensorIterator<float_complex,float,float_complex>;
    template class TensorIterator<float_complex,float,long>;
    template class TensorIterator<float_complex,float_complex,double>;
    template class TensorIterator<float_complex,float_complex,double_complex>;
    template class TensorIterator<float_complex,float_complex,float>;
    template class TensorIterator<float_complex,float_complex,float_complex>;
    template class TensorIterator<float_complex,float_complex,long>;
    template class TensorIterator<float_complex,long,double>;
    template class TensorIterator<float_complex,long,double_complex>;
    template class TensorIterator<float_complex,long,float>;
    template class TensorIterator<float_complex,long,float_complex>;
    template class TensorIterator<float_complex,long,long>;
    template class TensorIterator<long,double,double>;
    template class TensorIterator<long,double,double_complex>;
    template class TensorIterator<long,double,float>;
    template class TensorIterator<long,double,float_complex>;
    template class TensorIterator<long,double,long>;
    template class TensorIterator<long,double_complex,double>;
    template class TensorIterator<long,double_complex,double_complex>;
    template class TensorIterator<long,double_complex,float>;
    template class TensorIterator<long,double_complex,float_complex>;
    template class TensorIterator<long,double_complex,long>;
    template class TensorIterator<long,float,double>;
    template class TensorIterator<long,float,double_complex>;
    template class TensorIterator<long,float,float>;
    template class TensorIterator<long,float,float_complex>;
    template class TensorIterator<long,float,long>;
    template class TensorIterator<long,float_complex,double>;
    template class TensorIterator<long,float_complex,double_complex>;
    template class TensorIterator<long,float_complex,float>;
    template class TensorIterator<long,float_complex,float_complex>;
    template class TensorIterator<long,float_complex,long>;
    template class TensorIterator<long,long,double>;
    template class TensorIterator<long,long,double_complex>;
    template class TensorIterator<long,long,float>;
    template class TensorIterator<long,long,float_complex>;
    template class TensorIterator<long,long,long>;
    */
}
