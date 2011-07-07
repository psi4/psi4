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
  \file h2.cc
  \brief Solves the Hartree-Fock equations for the hydrogen molecule
  \defgroup examplesh2hf Hartree-Fock equations for the hydrogen molecule
  \ingroup examples

  The Hartree-Fock wave function is computed for the hydrogen molecule
  in three dimensions without using symmetry.

  Since all of the details except for the nuclear potential are the
  same, please refer to the \ref examplehehf helium atom HF example.

*/


#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/operator.h>

using namespace madness;

static const double R = 1.4;    // bond length
static const double L = 64.0*R; // box size
static const long k = 6;        // wavelet order
static const double thresh = 1e-4; // precision

static double guess(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return (exp(-sqrt(x*x+y*y+(z-R/2)*(z-R/2)+1e-8))+
            exp(-sqrt(x*x+y*y+(z+R/2)*(z+R/2)+1e-8)));
}

static double V(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return -1.0/sqrt(x*x+y*y+(z-R/2)*(z-R/2)+1e-8)+
           -1.0/sqrt(x*x+y*y+(z+R/2)*(z+R/2)+1e-8);
}

void iterate(World& world, real_function_3d& V, real_function_3d& psi, double& eps) {
    real_convolution_3d op = BSHOperator3D(world, sqrt(-2*eps), 0.001, 1e-6);
    real_function_3d Vpsi = (V*psi);
    Vpsi.scale(-2.0).truncate();
    real_function_3d tmp = op(Vpsi).truncate();
    double norm = tmp.norm2();
    real_function_3d r = tmp-psi;
    double rnorm = r.norm2();
    double eps_new = eps - 0.5*inner(Vpsi,r)/(norm*norm);
    if (world.rank() == 0) {
        print("norm=",norm," eps=",eps," err(psi)=",rnorm," err(eps)=",eps_new-eps);
    }
    psi = tmp.scale(1.0/norm);
    eps = eps_new;
}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    
    startup(world,argc,argv);
    std::cout.precision(6);

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_initial_level(5);
    FunctionDefaults<3>::set_truncate_mode(1);  
    FunctionDefaults<3>::set_cubic_cell(-L/2, L/2);
    // for (int i=0; i<3; i++) {
    // FunctionDefaults<3>::cell(i,0) = -L/2;
    // FunctionDefaults<3>::cell(i,1) =  L/2;
    // }
    
    real_function_3d Vnuc = real_factory_3d(world).f(V);
    real_function_3d psi  = real_factory_3d(world).f(guess);
    psi.truncate();
    psi.scale(1.0/psi.norm2());

    real_convolution_3d op = CoulombOperator(world, 0.001, 1e-6);

    double eps = -0.6;
    for (int iter=0; iter<10; iter++) {
        real_function_3d rho = square(psi).truncate();
        real_function_3d potential = Vnuc + op(rho).truncate();
        iterate(world, potential, psi, eps);
    }

    double kinetic_energy = 0.0;
    for (int axis=0; axis<3; axis++) {
        real_derivative_3d D = free_space_derivative<double,3>(world, axis);
        real_function_3d dpsi = D(psi);
        kinetic_energy += inner(dpsi,dpsi);
    }

    real_function_3d rho = square(psi);
    double two_electron_energy = inner(op(rho),rho);
    double nuclear_attraction_energy = 2.0*inner(Vnuc,rho);
    double nuclear_repulsion_energy = 1.0/R;
    double total_energy = kinetic_energy + two_electron_energy + 
        nuclear_attraction_energy + nuclear_repulsion_energy;
    double virial = (nuclear_attraction_energy + two_electron_energy + nuclear_repulsion_energy) / kinetic_energy;

    if (world.rank() == 0) {
        print("            Kinetic energy ", kinetic_energy);
        print(" Nuclear attraction energy ", nuclear_attraction_energy);
        print("       Two-electron energy ", two_electron_energy);
        print(" Nuclear  repulsion energy ", nuclear_repulsion_energy);
        print("              Total energy ", total_energy);
        print("                    Virial ", virial);
    }

    world.gop.fence();

    finalize();
    return 0;
}
