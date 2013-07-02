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
const double L = 10; // box is [-L,L]
const double sigma = 0.3; // Surface width

const double epsilon_0 = 1.0; // Interior dielectric
const double epsilon_1 = 10.0; // Exterior dielectric
const double R = 2.456644; // Radius of cavity

// Crude macro for timing
double XXstart;
#define TIME(MSG,X) XXstart=wall_time();         \
                    X; \
                    if (world.rank() == 0) print("timer:",MSG,"used",wall_time()-XXstart) \

double reciprocal(double x) {
    return 1.0/x;
}

double nuclear_charge_function(const coord_3d& r) {
    const double expnt = 100.0;
    const double coeff = pow(1.0/constants::pi*expnt,0.5*3);
    return coeff*exp(-expnt*(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]));
}

double electronic_charge_function(const coord_3d& r) {
    const double coeff = 1.0/constants::pi;
    return coeff*exp(-2.0*sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]));
}

double charge_function(const coord_3d& r) {
    return nuclear_charge_function(r) - electronic_charge_function(r);
}

class SVPESolver {
    const double thresh;
    const double minlen;
    const double sigma;                 //< Width of surface layer
    const double epsilon_0;             //< Interior dielectric
    const double epsilon_1;             //< Exterior dielectric
    real_convolution_3d op;             //< Coulomb operator (1/r ... no 4pi)
    //std::vector<real_convolution_3d_ptr> gop; //< Gradient of the Coulomb operator
    vector_real_function_3d dlog; //< Log-derivative of the dielectric
    real_function_3d rdielectric; //< Reciprocal of the dielectric

    // /// Apply Gbar 
    // real_function_3d gbar(const real_function_3d& f) const {
    //     real_function_3d r = (*gop[0])(f).truncate()*dlog[0] + (*gop[1])(f).truncate()*dlog[1] + (*gop[1])(f).truncate()*dlog[1];
    //     return r.truncate();
    // }

    
    // /// Solve for the full Coulomb potential using the other formulation
    // real_function_3d solve2(const real_function_3d& rho, 
    //                         const real_function_3d& uguess = real_function_3d(), 
    //                         bool printing = true) const {
        
    //     real_function_3d charge = rdielectric*rho;
    //     charge.truncate();
    //     real_function_3d sig0 = gbar(charge);
    //     real_function_3d sc = uguess.is_initialized() ? make_surface_charge(uguess) : sig0;
    //     sc = sig0 + gbar(sc);
    //     NonlinearSolver sol2;
    //     for (int iter=0; iter<10; iter++) {
    //         double start = wall_time();
    //         real_function_3d r = sc - (sig0 + gbar(sc));
    //         sc = sol2.update(sc,r);
    //         double rnorm = r.norm2();
    //         if (printing) print(iter, sc.trace(), rnorm, wall_time()-start);
    //         if (rnorm < 10.0*thresh) break;
    //     }
    //     return op(rho + sc).truncate();
    // }

public:
    SVPESolver(World& world,
               double sigma, double epsilon_0, double epsilon_1, 
               const vector_real& atomic_radii, const vector_coord_3d& atomic_coords,
               const double minlen)
        : thresh(FunctionDefaults<3>::get_thresh())
        , minlen(minlen)
        , sigma(sigma)
        , epsilon_0(epsilon_0)
        , epsilon_1(epsilon_1)
        , op(CoulombOperator(world, minlen, thresh))
          //, gop(GradCoulombOperator(world, minlen, thresh))
        , dlog(3)
    {
        // Functors for mask related quantities
        real_functor_3d rdielectric_functor(new MolecularVolumeExponentialSwitchReciprocal(sigma, epsilon_0, epsilon_1, atomic_radii, atomic_coords));
        real_functor_3d gradx_functor(new MolecularVolumeExponentialSwitchLogGrad(sigma, epsilon_0, epsilon_1, atomic_radii, atomic_coords,0));
        real_functor_3d grady_functor(new MolecularVolumeExponentialSwitchLogGrad(sigma, epsilon_0, epsilon_1, atomic_radii, atomic_coords,1));
        real_functor_3d gradz_functor(new MolecularVolumeExponentialSwitchLogGrad(sigma, epsilon_0, epsilon_1, atomic_radii, atomic_coords,2));
        
        // Make the actual functions
        const double rfourpi = 1.0/(4.0*constants::pi);
        rdielectric = real_factory_3d(world).functor(rdielectric_functor).nofence();
        dlog[0] = real_factory_3d(world).functor(gradx_functor).nofence();
        dlog[1] = real_factory_3d(world).functor(grady_functor).nofence();
        dlog[2] = real_factory_3d(world).functor(gradz_functor); // FENCE
        scale(world, dlog, rfourpi);
        rdielectric.truncate(false);
        truncate(world, dlog);
    }

    /// Given the full Coulomb potential computes the surface charge
    real_function_3d make_surface_charge(const real_function_3d& u) const {
        real_derivative_3d Dx = free_space_derivative<double,3>(u.world(), 0);
        real_derivative_3d Dy = free_space_derivative<double,3>(u.world(), 1);
        real_derivative_3d Dz = free_space_derivative<double,3>(u.world(), 2);
        return (dlog[0]*Dx(u) + dlog[1]*Dy(u) + dlog[2]*Dz(u)).truncate();
    }

    /// Solve for the full Coulomb potential using the free-particle GF
    real_function_3d solve(const real_function_3d& rho, 
                           const real_function_3d uguess = real_function_3d(), 
                           bool printing = true) const {
        real_function_3d charge = rdielectric*rho;
        charge.truncate();

        // Initial guess is constant dielectric        
        real_function_3d u0 = op(charge).truncate();
        real_function_3d u = uguess.is_initialized() ? uguess : u0;
        double unorm = u.norm2();
        NonlinearSolver solver;
        for (int iter=0; iter<20; iter++) {
            double start = wall_time();
            real_function_3d surface_charge = make_surface_charge(u);
            real_function_3d r = (u - u0 - op(surface_charge)).truncate();
            double sigtot = surface_charge.trace();
            surface_charge.clear();
            
            real_function_3d unew = solver.update(u, r);
            
            double change = (unew-u).norm2();
            if (printing) 
                print("iter", iter, "change", change,
                      "soln(3.0)", u(coord_3d(3.0)),
                      "surface charge", sigtot,"used",wall_time()-start);
            
            // Step restriction 
            if (change > 0.3*unorm) 
                u = 0.5*unew + 0.5*u;
            else 
                u = unew;
            
            if (change < std::max(1e-3,10.0*thresh)) break;
        }
        return u;
    }
};

int main(int argc, char **argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);

    // Function defaults
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_cubic_cell(-L, L);
    FunctionDefaults<3>::set_initial_level(4);
    FunctionDefaults<3>::set_truncate_on_project(true);
    FunctionDefaults<3>::set_bc(BC_FREE);

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

    TIME("make charge ", real_function_3d charge   = real_factory_3d(world).f(charge_function)); 
    charge.truncate();
    SVPESolver solver(world, sigma, epsilon_0, epsilon_1, atomic_radii, atomic_coords, min(1e-3,sigma*0.1));

    real_function_3d u = solver.solve(charge);
    print("Solving again to verify that the initial guess works");
    u = solver.solve(charge,u);

    // For comparison with Chipman make the reaction potential
    real_convolution_3d op = CoulombOperator(world, min(1e-3,sigma*0.1), thresh);
    
    TIME("make ufree  ", real_function_3d ufree = op(charge).truncate());
    print("<rhotot|ureact>", 0.5*charge.inner((u-ufree)));

    finalize();
    return 0;
}
