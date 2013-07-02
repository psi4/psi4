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
#define WORLD_INSTANTIATE_STATIC_TEMPLATES  

/*!
  \file examples/sininteg.cc
  \brief Compute the integral sin(x) x=0..10
  \defgroup sininteg First example from getting started guide
  \ingroup examples

  Computes the integral 
  \f[
     \int_0^{10} sin(x) dx
  \f]
  by projecting \f$ sin(x) \f$ into the discontinuous spectral element
  basis and using the \c trace() method.

 */


#include <mra/mra.h>

using namespace madness;

double myf(const coord_1d& r) {
    return std::sin(r[0]);
}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);

    startup(world,argc,argv);

    FunctionDefaults<1>::set_cubic_cell(0,10);

    real_function_1d f = real_factory_1d(world).f(myf);

    double integral = f.trace();

    if (world.rank() == 0) print("The result is", integral);

    finalize();
    return 0;
}
