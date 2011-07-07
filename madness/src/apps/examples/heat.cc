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
  \file examples/heat.cc
  \brief Example Green function for the 3D heat equation
  \defgroup exampleheat Solves heat equation using the Green's function
  \ingroup examples

  The source is <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local/trunk/src/apps/examples/heat.cc>here</a>.

  \par Points of interest
  - use of a functor to compute the solution at an arbitrary future time
  - convolution with the Green's function

  \par Background

  Solves the 3D time-dependent heat equation
  \f[
  \frac{\partial u}{\partial t} = c \nabla^2 u(r,t)
  \f]
  by direct convolution with the Green's function,
  \f[
  \frac{1}{\sqrt{4 \pi c t}}  \exp \frac{-x^2}{4 c t}
  \f]

*/

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/operator.h>
#include <constants.h>

using namespace madness;

typedef Vector<double,3> coordT;
typedef std::shared_ptr< FunctionFunctorInterface<double,3> > functorT;
typedef Function<double,3> functionT;
typedef FunctionFactory<double,3> factoryT;
typedef SeparatedConvolution<double,3> operatorT;
typedef Tensor<double> tensorT;

static const double L = 10;     // Half box size
static const long k = 8;        // wavelet order
static const double thresh = 1e-6; // precision
static const double c = 2.0;       //
static const double tstep = 0.333;
static const double alpha = 1.9; // Exponent

// Initial Gaussian with exponent alpha
static double uinitial(const coordT& r) {
    const double x=r[0], y=r[1], z=r[2];
    return exp(-alpha*(x*x+y*y+z*z))*pow(constants::pi/alpha,-1.5);
}


// Exact solution at time t
class uexact : public FunctionFunctorInterface<double,3> {
    double t;
public:
    uexact(double t) : t(t) {}

    double operator()(const coordT& r) const {
        const double x=r[0], y=r[1], z=r[2];
        double rsq = (x*x+y*y+z*z);

        return exp(-rsq*alpha/(1.0+4.0*alpha*t*c)) * pow(alpha/((1+4*alpha*t*c)*constants::pi),1.5);
    }
};


int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);

    startup(world, argc, argv);
    std::cout.precision(6);

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_autorefine(false);
    FunctionDefaults<3>::set_cubic_cell(-L, L);

    functionT u0 = factoryT(world).f(uinitial);
    u0.truncate();

    double u0_norm = u0.norm2();
    double u0_trace = u0.trace();

    if (world.rank() == 0) print("Initial norm", u0_norm,"trace", u0_trace);

    world.gop.fence();


    // du/dt = c del^2 u
    //
    // Time evolution operator is 1/sqrt(4 pi c t)  exp(-x^2 / 4 c t)

    tensorT expnt(1), coeff(1);
    expnt[0] = 1.0/(4.0*c*tstep);
    coeff[0] = pow(4.0*constants::pi*c*tstep,-1.5);

    operatorT G(world, coeff, expnt);

    functionT ut = G(u0);

    double ut_norm = ut.norm2();
    double ut_trace = ut.trace();
    double err = ut.err(uexact(tstep));

    if (world.rank() == 0) print("Final norm", ut_norm,"trace", ut_trace,"err",err);

    finalize();
    return 0;
}

