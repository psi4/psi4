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


  $Id: twoscale.h 1602 2009-12-27 19:53:06Z rjharrison $
*/


#ifndef MADNESS_MRA_TWOSCALE_H__INCLUDED
#define MADNESS_MRA_TWOSCALE_H__INCLUDED

#include <madness_config.h>
#include <tensor/tensor.h>
#include <world/world.h>

namespace madness {
    extern void load_coeffs(World& world, const char* dir);
    extern bool two_scale_coefficients(int k,
                                           Tensor<double>* h0, Tensor<double>* h1,
                                           Tensor<double>* g0, Tensor<double>* g1);
    extern bool two_scale_hg(int k, Tensor<double>* hg);
    extern bool test_two_scale_coefficients();

    extern bool autoc(int k, Tensor<double>* c);
    extern bool test_autoc();
}

#endif // MADNESS_MRA_TWOSCALE_H__INCLUDED
