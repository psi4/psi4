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
#include <mra/mra.h>
#include <constants.h>


using namespace std;
using namespace madness;

static const std::size_t NDIM=3;
static const double L=2.0; // Simulate in [0,L]^NDIM

typedef Vector<double,NDIM> coordT;
typedef Function<double,NDIM> functionT;
typedef FunctionFactory<double,NDIM> factoryT;

double func(const coordT& r) {
    // Test is product of cos(2*Pi*(d+1)/L + 0.2)
    double result = 1.0;
    for (std::size_t d=0; d<NDIM; ++d) {
        double fac = 2.0*constants::pi/L;
        result *= cos(fac*r[d] + 0.2);
    }
    return result;
}

std::size_t axis; // Axis we are diffing

double dfunc(const coordT& r) {
    // Test is product of cos(2*Pi*(d+1)/L + 0.2)
    double result = 1.0;
    for (std::size_t d=0; d<NDIM; ++d) {
        double fac = 2.0*constants::pi/L;
        if (d == axis) {
            result *= -fac*sin(fac*r[d] + 0.2);
        }
        else {
            result *= cos(fac*r[d] + 0.2);
        }
    }
    return result;
}

struct FunctorInterfaceWrapper : public FunctionFunctorInterface<double,NDIM> {
    double (*f)(const coordT&);

    FunctorInterfaceWrapper(double (*f)(const coordT&)) : f(f) {}

    double operator()(const coordT& x) const {return f(x);}
};

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);
    FunctionDefaults<NDIM>::set_cubic_cell(0.0,L);

    Tensor<int> bc_periodic(NDIM,2);
    bc_periodic.fill(1);

    BoundaryConditions<NDIM> bc(BC_PERIODIC);

    functionT f = factoryT(world).f(func);

    for (axis=0; axis<NDIM; ++axis) {
        Derivative<double,NDIM> dx(world, axis, bc);
        functionT df = dx(f) ;
        double err = df.err(FunctorInterfaceWrapper(dfunc));
         if (world.rank() == 0) print(axis,"error",err);
    }

    finalize();
    return 0;
}
