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
  \file nonlinschro.cc
  \brief Solves 1D nonlinear Schr&ouml;dinger equation
  \defgroup examplenonlinsc Solves a 1D nonlinear Schr&ouml;dinger equation
  \ingroup examples

  The source is <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local/trunk/src/apps/examples/nonlinschro.cc>here</a>.

  \par Points of interest
  - Convolution with the negative energy (bound state) Helmholtz Green's function
  - Iterative solution of the integral form of the equation
  - Smooth truncation of the density to manage numerical noise amplified by nearly singular potential
  - Plotting of the solution and potential

  \par Background

  This illustrates solution of a non-linear Schr&ouml;dinger motivated by
  exploring problems associated with equations of the same form from
  nuclear physics.

  We seek the lowest eigenfunction of
  \f[
  -\nabla^2 \psi(x) + V(x) \psi(x) = E \psi(x)
  \f]
  where the potential is
  \f[
  V(x) = -a \exp(-b x^2) + \frac{c}{(n(x)+\eta)^{1/3}} - d n(x)^{5/3} + V_{\mbox{shift}}
  \f]
  The parameters \f$ a \f$, \f$ b \f$, \f$ c \f$, \f$ d \f$, and  \f$ V_{\mbox{shift}} \f$
  are given in the code.  The density is given by
  \f[
  n(x) = \psi(x)^2
  \f]
  There would normally be multiple states occupied but for simplicity we are employing
  just one.

  The first term in the potential is weak and seems to be there to
  stabilize the solution.  The second term seems to act as a confining
  potential since it becomes large and positive when the density is
  small.  The third term represents short-range attraction between
  particles, and the fourth adjusts the zero of energy.

  [These notes were written by a chemist ... if you are a nuclear physicist could you
  please clean them up?].

  \par Implementation

  The integral form of the equation is
  \f[
     \psi = - 2 G_{\mu} * \left ( V \psi \right)
  \f]
  where \f$ G_{\mu}\f$ is the Green's function for the Helmholtz equation
  \f[
     \left( - \frac{d^2}{dx^2} + \mu^2 \right) G(x,x^{\prime}) = \delta(x-x^{\prime})
  \f]
  where \f$\mu = \sqrt{-2 E}\f$.

  We employ a simple fixed-point iteration to the self-consistent
  solution, but strong damping or step restriction is necessary to
  ensure convergence.  This is due to the nearly singluar potential,
  the problem being exacerbated by small \f$ \eta \f$.  A reliable
  solution scheme seems to be to first solve with a large value of
  \f$ \eta \f$ and then to reduce it in several steps to its final value.
  A much more efficient scheme would involve use of a non-linear
  equation solver instead of simple iteration.

  The density is analytically everywhere positive, but numeric noise
  can introduce regions where it is zero or even slightly negative.
  To avoid very non-physical results, we smoothly switch values less
  than some threshold to a minimum value.  This is preferable to
  employing a sharp cutoff.

*/


#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/operator.h>
#include <algorithm>

using namespace madness;

static const double L = 10.0;           // box size
static const long k = 10;                // wavelet order
static const double thresh = 1e-8;      // precision
static const double ntol = thresh*10.0; // Cutoff for numerically small density
static const double a = 0.125;
static const double b = 1.0;
static const double c = 1.0;
static const double d = 1.0;
static const double eta_end = 1e-5; // Target value for eta
static double eta = eta_end*1024;   // Initial value foe eta


static const double vasymp = c/std::pow(eta_end,1.0/3.0);
static const double vshift = -3.0; // Only impacts BSH operator

// For x<xmax, smoothly restricts x to be greater than or equal to xmin>0
double munge(double x, double xmin, double xmax) {
    // A cubic polyn that smoothly interpolates between
    //
    // x=0    where it has value xmin and zero slope, and
    //
    // x=xmax where it has value xmax and unit slope

    if (x > xmax) {
        return x;
    }
    else if (x <= 0.0) {
        return xmin;
    }
    else {
        double x1 = x/xmax;
        double x2 = x1*x1;
        double x3 = x1*x2;
        return xmin+(2.0*xmax-3.0*xmin)*x2 - (-2.0*xmin+xmax)*x3;
    }
}

// This invoked to compute the part of the potential that is a function of the density
template <typename T>
static void Vdynamic(const Key<1> & key, Tensor<T>& t)
{
    static const double onethird = 1.0/3.0;
    static const double fivethird = 5.0/3.0;
    UNARY_OPTIMIZED_ITERATOR(T, t,
                             double n=munge(*_p0, ntol, 10.0*ntol);
                             *_p0 = c/pow(n+eta,onethird) - d*pow(n+eta,fivethird));
}

// The part of the potential that does not depend on the density
static double Vstatic(const coord_1d& r) {
    const double x=r[0];
    return -a*exp(-b*x*x) + vshift;
}

static double guess(const coord_1d& r) {
    const double x=r[0];
    return exp(-0.2*x*x);
}

real_function_1d make_potential(World& world, const real_function_1d& rho) {
    real_function_1d vstatic = real_factory_1d(world).f(Vstatic);
    real_function_1d vdynamic = copy(rho);
    vdynamic.unaryop(&Vdynamic<double>);
    return vstatic + vdynamic;
}

void iterate(World& world, real_function_1d& psi) {
    // Compute density and full potential
    real_function_1d rho = psi*psi;
    real_function_1d v = make_potential(world, rho);

    // Compute energy components and print
    real_function_1d dpsi = Derivative<double,1>(world,0)(psi);
    dpsi.verify_tree();
    double kinetic_energy = inner(dpsi,dpsi);
    double potential_energy = inner(rho, v);
    double energy = potential_energy + kinetic_energy;

    // Update the wave function
    real_function_1d Vpsi = v*psi;
    Vpsi.scale(-1.0).truncate();
    real_convolution_1d op = BSHOperator<1>(world, sqrt(-energy), 0.001, 1e-6);
    real_function_1d tmp = op(Vpsi).truncate();
    real_function_1d r = tmp-psi;
    double rnorm = r.norm2();
    //double step = std::min(1.0,0.05/rnorm);
    double step = 0.025;
    if (rnorm > 10) step = 0.1/rnorm;
    psi = psi + step*r;
    psi.scale(1.0/psi.norm2());

    print("KE =", kinetic_energy,"  PE =",potential_energy,"  E =",energy, "  err(psi) =", rnorm, "  step =", step);
}


int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);

    startup(world,argc,argv);
    std::cout.precision(6);

    FunctionDefaults<1>::set_k(k);
    FunctionDefaults<1>::set_thresh(thresh);
    FunctionDefaults<1>::set_truncate_mode(1);
    FunctionDefaults<1>::set_cubic_cell(-L,L);

    real_function_1d psi = real_factory_1d(world).f(guess);
    print(psi.norm2());
    psi.scale(1.0/psi.norm2());
    print(psi.norm2());

    while (1) {
        print("\nSOLVING WITH ETA", eta);
        for (int iter=0; iter<30; iter++)
            iterate(world, psi);
        if (eta <= eta_end) break;
        eta *= 0.5;
    }

    real_function_1d rho = psi*psi;
    real_function_1d v = make_potential(world, rho);

    coord_1d lo(-L), hi(L);
    double scale = vasymp/psi(coord_1d(0.0));
    plot_line("psi.txt", 201, lo, hi, psi*scale, v);

    world.gop.fence();

    finalize();
    return 0;
}
