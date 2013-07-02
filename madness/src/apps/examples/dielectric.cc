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
  \file examples/dielectric.cc
  \brief Example solution of Poisson's equation in a dielectric (polarizable) medium
  \defgroup exampledielectric Poisson's equation in a dielectric medium
  \ingroup examples
  
  The source is <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local/trunk/src/apps/examples/dielectric.cc>here</a>.

  \par Points of interest
  - use of iterative equation solver 
  - convolution with the Green's function
  - unary operator to compute reciprocal of a function
  - use of diffuse domain approximation to represent interior surfaces
  - line and volume plots

  \par Background

  We wish to solve Poisson's equation within a non-homogeneous medium, i.e., 
  in which the permittivity is not constant. 
  \f[
  \nabla . \left( \epsilon(r) \nabla u(r)  \right) = - 4 \pi \rho(r)
  \f]
  Expanding and rearranging yields,
  \f[
     \nabla^2 u(r) = - \frac{1}{\epsilon(r)} \left( 4 \pi \rho(r) + \nabla \epsilon(r) .  \nabla u(r) \right)
  \f]
  We can interpret \f$\rho / \epsilon\f$ as the effective free charge density
  and \f$\nabla \epsilon . \nabla u / \epsilon\f$ as the induced surface charge density.
  Assuming free-space boundary conditions at long range we can invert the Laplacian
  by convolution with the free-space Green's function that is \f$-1 / 4 \pi |r-s|\f$.
  Note that MADNESS provides convolution with \f$G(r,s) = 1/|r-s|\f$ so you have to 
  keep track of the \f$-1/4\pi\f$ yourself.

  Thus, our equation becomes (deliberately written in the form of a
  fixed point iteration --- given a guess for \f$u\f$ we can compute
  the r.h.s. and directly obtain a new value for \f$u\f$ that hopefully
  is closer to the solution)
  \f[
     u = G * \left(\rho_{\mbox{vol}} + \rho_{\mbox{surf}} \right)
  \f]
  where 
  \f[
     \rho_{\mbox{vol}} = \frac{\rho}{\epsilon}
  \f]
  and
  \f[
     \rho_{\mbox{surf}} = \frac{\nabla \epsilon .  \nabla u}{4 \pi \epsilon}
  \f]
  
  Let's solve a problem to which we know the exact answer --- a point
  charge at the center of a sphere \f$R=2\f$.  Inside the sphere the permittivity
  is \f$\epsilon_1 = 1\f$ and outside it is \f$\epsilon_2 = 10\f$.  The exact
  solution is 
  \f[
     u(r) = 
            \left \lbrace 
               \begin{array}{cc}
                    \frac{1}{\epsilon_2 |r|} & |r| > R \\
                    \frac{1}{\epsilon_1 |r|} + \left( \frac{1}{\epsilon_2} - \frac{1}{\epsilon_1} \right) \frac{1}{R}  & |r| < R
               \end{array}
            \right .
            
  \f]
  The surface charge density integrated over the suface has the value
  \f[
     q_{\mbox{surface}} = -\frac{\epsilon_2 - \epsilon_1}{\epsilon_2} = -0.9
  \f]

  To implement the problem in MADNESS we want the permittivity defined
  as a function over all space, so we define a characteristic function
  \f$C(r)\f$ (or mask) that is defined to be one inside the sphere and
  zero outside
  \f[ 
    C(r) = H(R-|r|) 
  \f] 
  where \f$H(x)\f$ is the Heaviside step function.  Hence, we have
  \f[
     \epsilon(r) = \epsilon_1 C(r) + \epsilon_2 \left( 1 - C(r) \right)
  \f]
  To smooth the discontinuity we replace the Heaviside step function
  with 
  \f[
     H(x,\sigma) = \frac{1}{2} \left( 1 + \mathop{\mathrm{erf}} \frac{x}{\sigma} \right)
  \f]
  where \f$\sigma\f$ is the effective width of the step (0.2 in the code).  Similarly, the
  point charge is replaced by a suitably normalized Gaussian with a an
  exponent sufficiently large as to appear as a point charge (delta
  function) on the scale of interest (we pick \f$\xi=100\f$)
  \f[
     \rho(r) = \pi^{-3/2} e^{-\xi |r|^2}
  \f]

  Starting from an initial guess we could try to iterate
  our equation for \f$u(r)\f$ but, sadly, this does not converge reliably.
  The simplest approach is to introduce some damping or step restriction.
  With \f$m\f$ indicating the iteration number, we write
  \f[
     u^{(m+1)} = \alpha u^{(m)} + (1-\alpha) G * \left(\rho^{(m)}_{\mbox{vol}} + \rho^{(m)}_{\mbox{surf}} \right)
  \f]
  with \f$ \alpha \in [0,1]\f$ .  This works (for sufficiently small \f$\alpha\f$)
  and is implemented in the code.
  
  A more robust approach is to use a solver that exploits information
  from previous iterations (the Krylov subspace) to estimate the
  optimal direction and length of the step to take.  Such a 
  solver is provided by nonlinsol.h.  Each iteration we pass
  the current solution and corresponding residual to the solver
  and it provides the next trial vector.

  One final point is how to compute the reciprocal of the
  permittivity.  This operation is not provided explicity by MADNESS,
  in part because even functions that are
  analytically never zero, might be (nearly) zero due to numerical truncation.
  Handling this issue correctly is problem specific.  More generally,
  we want to compute a function-of-a-function.  If \f$f(r)\f$ 
  is a MADNESS function that takes a d-dimensional vector as its
  argument (in this examaple \f$d=3\f$) and \f$F(x)\f$ is computable
  function that takes a scalar argument, we want to compute the
  new MADNESS function
  \f[
      g(r) = F(f(r))
  \f]
  There are various ways to do this.  The simplest is employed here ---
  we define a function (\c reciprocal() ) that is passed into the
  \c unaryop() method that modifies the function \em in-place.
  This simple approach is justified here since we know our 
  input function is never zero even due to numerical noise,
  and the output function is about as smooth as the input.
  If it were not, we might have to refine the input function
  to obtain the desired precision.
  
*/


#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/operator.h>
#include <mra/funcplot.h>
#include <linalg/solvers.h>
#include <examples/molecularmask.h>
#include <examples/nonlinsol.h>
#include <constants.h>
#include <vector>

using namespace madness;
using namespace std;

const int k = 6; // wavelet order
const double thresh = 1e-4; // truncation threshold
const double L = 5; // box is [-L,L]
const double sigma = 0.2; // Surface width

const double epsilon_0 = 1.0; // Interior dielectric
const double epsilon_1 =10.0; // Exterior dielectric
const double R = 2.0; // Radius of cavity

// Crude macro for timing
double XXstart;
#define TIME(MSG,X) XXstart=wall_time();         \
                    X; \
                    if (world.rank() == 0) print("timer:",MSG,"used",wall_time()-XXstart) \


double reciprocal(double x) {
    return 1.0/x;
    }

double charge_function(const coord_3d& r) {
    const double expnt = 100.0;
    const double coeff = pow(1.0/constants::pi*expnt,0.5*3);
    return coeff*exp(-expnt*(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]));
}

double exact_function(const coord_3d& x) {
    const double expnt = 100.0;
    double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

    if (r > R) {
        return 1.0/(epsilon_1*r);
    }
    else {
        return erf(sqrt(expnt)*r)/(epsilon_0*r)
            + (1.0/epsilon_1 - 1.0/epsilon_0)/R;
    }
}

int main(int argc, char **argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);

    coord_3d lo(0.0), hi(0.0); // Range for line plotting
    lo[0] =-5.0;
    hi[0] = 5.0;

    // Function defaults
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_cubic_cell(-L, L);
    FunctionDefaults<3>::set_initial_level(3);
    FunctionDefaults<3>::set_truncate_on_project(true);
    //FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<3>::set_bc(BC_FREE);

    // The Coulomb operator (this is just 1/r ... whereas the notes are -1/4pir)
    real_convolution_3d op = CoulombOperator(world, sigma*0.001, thresh*0.1);

    // Derivative operators
    real_derivative_3d Dx = free_space_derivative<double,3>(world, 0);
    real_derivative_3d Dy = free_space_derivative<double,3>(world, 1);
    real_derivative_3d Dz = free_space_derivative<double,3>(world, 2);

    // We will have one sphere of radius R centered at the origin
    vector<double> atomic_radii(1,R);
    vector<coord_3d> atomic_coords(1,coord_3d(0.0));
    print("k     ", k);
    print("thresh", thresh);
    print("L     ", L);
    print("sigma ", sigma);
    print("eps0  ", epsilon_0, "  eps1  ", epsilon_1);
    print("radii ", atomic_radii);
    print("coords", atomic_coords);

    // Functors for mask related quantities
    real_functor_3d volume_functor(new MolecularVolumeMask(sigma, atomic_radii, atomic_coords));
    real_functor_3d gradx_functor(new MolecularVolumeMaskGrad(sigma, atomic_radii, atomic_coords, 0));
    real_functor_3d grady_functor(new MolecularVolumeMaskGrad(sigma, atomic_radii, atomic_coords, 1));
    real_functor_3d gradz_functor(new MolecularVolumeMaskGrad(sigma, atomic_radii, atomic_coords, 2));
    real_functor_3d surface_functor(new MolecularSurface(sigma, atomic_radii, atomic_coords));

    // Make the actual functions
    TIME("make volume ", real_function_3d volume = real_factory_3d(world).functor(volume_functor));
    TIME("make gradx  ", real_function_3d gradx = real_factory_3d(world).functor(gradx_functor));
    TIME("make grady  ", real_function_3d grady = real_factory_3d(world).functor(grady_functor));
    TIME("make gradz  ", real_function_3d gradz = real_factory_3d(world).functor(gradz_functor));
    TIME("make surface", real_function_3d surface = real_factory_3d(world).functor(surface_functor));
    TIME("make charge ", real_function_3d charge = real_factory_3d(world).f(charge_function));
    TIME("make exact  ", real_function_3d exact = real_factory_3d(world).f(exact_function));
    
    // Reciprocal of the dielectric function
    real_function_3d rdielectric = epsilon_0*volume + epsilon_1*(1.0-volume);
    rdielectric.unaryop(reciprocal);

    // Gradient of the dielectric function
    real_function_3d di_gradx = (epsilon_0-epsilon_1)*gradx;
    real_function_3d di_grady = (epsilon_0-epsilon_1)*grady;
    real_function_3d di_gradz = (epsilon_0-epsilon_1)*gradz;

    // Print some values for sanity checking
    print("the volume is", volume.trace());
    print("the area   is", surface.trace());
    print("the charge is", charge.trace());

    // Free up stuff we are not using any more to save memory
    volume.clear();
    surface.clear();
    gradx.clear();
    grady.clear();
    gradz.clear();

    const double rfourpi = 1.0/(4.0*constants::pi);
    charge = rdielectric*charge;

    // Initial guess is constant dielectric
    real_function_3d u = op(charge).truncate();
    double unorm = u.norm2();

    const bool USE_SOLVER = true;

    if (USE_SOLVER) {
        // This section employs a non-linear equation solver from solvers.h
        // and nonlinsol.h as described in
        // http://onlinelibrary.wiley.com/doi/10.1002/jcc.10108/abstract
        NonlinearSolver solver;
        for (int iter=0; iter<20; iter++) {
            real_function_3d surface_charge = 
                rfourpi*rdielectric*(di_gradx*Dx(u) + di_grady*Dy(u) + di_gradz*Dz(u));
            real_function_3d r = (u - op(charge + surface_charge)).truncate();

            real_function_3d unew = solver.update(u, r);

            double change = (unew-u).norm2();
            double err = (u-exact).norm2();
            print("iter", iter, "change", change, "err", err, 
                  "exact(3.0)", exact(coord_3d(3.0)), "soln(3.0)", u(coord_3d(3.0)),
                  "surface charge", surface_charge.trace());

            if (change > 0.3*unorm) 
                u = 0.5*unew + 0.5*u;
            else 
                u = unew;

            if (change < 10.0*thresh) break;
        }
    }
    else {
        // This section employs a simple iteration with damping (step restriction)
        for (int iter=0; iter<20; iter++) {
            real_function_3d u_prev = u;
            real_function_3d surface_charge = 
                rfourpi*rdielectric*(di_gradx*Dx(u) + di_grady*Dy(u) + di_gradz*Dz(u));
            u = op(charge + surface_charge).truncate();
            
            double change = (u-u_prev).norm2();
            double err = (u-exact).norm2();
            print("iteration", iter, change, err, exact(coord_3d(3.0)), u(coord_3d(3.0)), surface_charge.trace());
            if (change > 0.3*unorm) u = 0.5*u + 0.5*u_prev;
            if (change < 10.0*thresh) break;
        }
    }


    plot_line("testpot.dat", 301, lo, hi, u, exact);
    real_tensor cell(3,2);
    cell(_,0) = -4.0;
    cell(_,1) =  4.0;

    plotdx(u, "testpot.dx", cell);
    plotdx(exact, "exact.dx", cell);
    plotdx(u-exact, "err.dx", cell);

    finalize();
    return 0;
}
