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
/*!
  \file functionio.cc
  \brief Example of I/O with functions
  \defgroup functionioeg Example of function I/O from getting started guide
  \ingroup examples

  \par Points of interest
  - Moving MADNESS functions to/from disk

  \par Background
  MADNESS functions are parallel objects distributed across the whole
  machine.  If you try to store a function into a sequential archive,
  only the part local to the calling process will be stored.
  I.e., you need to use a parallel archive.  The parallel archive
  has an adjustable number of proceses that actually perform I/O
  for effiency.  If your functions are few and small, the default of one
  writer is OK.  For large (or many) functions and I/O is a bottleneck,
  then increase the number of I/O nodes.

  \par Implementation
  Three different functions are created, written to disk, read back,
  and finally compared with the originals.

 */

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <constants.h>
using namespace madness;
using namespace std;

static const double L = 10;     // Half box size

/// Normalized gaussian
static double function(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    static const double alpha = 1.9; // Exponent of gaussian
    static const double fac = pow(constants::pi/alpha,-1.5);
    return exp(-alpha*(x*x+y*y+z*z))*fac;
}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    startup(world, argc, argv);
    FunctionDefaults<3>::set_cubic_cell(-L, L);
    const char filename[] = "mi_casa_es_su_casa";

    if (world.rank() == 0) print("Making test functions");
    real_function_3d f = real_factory_3d(world).f(::function);
    real_function_3d g = f*f;
    real_function_3d h = g - f;

    if (world.rank() == 0) print("Writing test functions -->", filename);
    {
        archive::ParallelOutputArchive out(world, filename);
        out & f & g & h;
    }

    if (world.rank() == 0) print("Reading test functions <--", filename);
    real_function_3d ftest, gtest, htest;
    {
        archive::ParallelInputArchive input(world, filename);
        input & ftest & gtest & htest;
    }

    double err = (f-ftest).norm2() + (g-gtest).norm2() + (h-htest).norm2();
    if (world.rank() == 0) print("The difference is", err);

    finalize();
}
