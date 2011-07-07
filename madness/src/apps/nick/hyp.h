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
  
  $Id$
*/
//By: Robert Harrison
#ifndef HYP_H
#define HYP_H

#include <iostream>
#include <complex>
#include <algorithm>
#include <cstdio> //NEEDED
#include <cmath>
#include <nick/mpreal.h>


typedef mpfr::mpreal extended_real;
typedef std::complex<extended_real> extended_complex;
typedef std::complex<double> complexd;


/// Computes 1F1(a,b,z) internally using extended precision

/// If result is larger than 1.0, result should be accurate to
/// full double precision, otherwise result is accurate to
/// somewhat better than 1e-17.  However, the termination
/// test is not very smart so there may be failures.
complexd conhyp(const complexd& a_arg,
                const complexd& b_arg,
                const complexd& z_arg);
#endif
