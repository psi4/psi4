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
  \file hatom_1d.cc
  \brief Solves the Schrodinger equation for the 1-d hydrogen atom
  \ingroup examples

  The Hartree-Fock wave function is computed for the hydrogen atom
  in one dimension without using symmetry.

*/


#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/operator.h>

using namespace madness;

static const double L = 100;    // box size
static const long k = 6;        // wavelet order
static const double thresh = 1e-4; // precision

static double guess(const coord_1d& r) {
    const double x=r[0];
    // exact, up to normalization
    return x*exp(-fabs(x));
}

static double V(const coord_1d& r) {
    const double x=r[0];
    return - 1.0/(fabs(x) + 1e-8);
}

void iterate(World& world, real_function_1d& V, real_function_1d& psi, double& eps) {
    real_convolution_1d op = BSHOperator<1>(world, sqrt(-2.0*eps), 0.001, 1e-6);
    real_function_1d Vpsi = (V*psi);
    Vpsi.scale(-1.0).truncate();
    real_function_1d tmp = op(Vpsi).truncate();
    tmp.scale(2.0).truncate();
    double norm = tmp.norm2();
    real_function_1d r = tmp-psi;
    double rnorm = r.norm2();
    double eps_new = eps - inner(Vpsi,r)/(norm*norm);
    if (world.rank() == 0) {
        print("norm=",norm," eps=",eps," ||dpsi||=",rnorm," deps=",eps_new-eps);
    }
    psi = tmp.scale(1.0/norm);
    eps = eps_new;
}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    
    startup(world,argc,argv);
    std::cout.precision(6);

    FunctionDefaults<1>::set_k(k);
    FunctionDefaults<1>::set_thresh(thresh);
    FunctionDefaults<1>::set_refine(true);
    FunctionDefaults<1>::set_initial_level(5);
    FunctionDefaults<1>::set_truncate_mode(1);
    FunctionDefaults<1>::set_cubic_cell(-L/2, L/2);
    //FunctionDefaults<1>::set_max_refine_level(2);
    // for (int i=0; i<1; i++) {
    // FunctionDefaults<1>::cell(i,0) = -L/2;
    // FunctionDefaults<1>::cell(i,1) =  L/2;
    // }
    
    real_function_1d Vnuc = real_factory_1d(world).f(V);
    real_function_1d psi  = real_factory_1d(world).f(guess);
    psi.truncate();
    psi.scale(1.0/psi.norm2());

    double eps = -0.5;
    for (int iter=0; iter<50; iter++) {
      std::cout << "Iteration " << iter << std::endl;
      //psi.print_tree();

      // plot psi
#if 0
      std::ostringstream oss;
      oss << "psi_" << iter << ".dat";
      const Vector<double,1> lo(-L/2);
      const Vector<double,1> hi(+L/2);
      plot_line(oss.str().c_str(), 10001, lo, hi, psi);
#endif

      real_function_1d potential = Vnuc;
      iterate(world, potential, psi, eps);
    }

    real_derivative_1d D = free_space_derivative<double,1>(world, 0);
    real_function_1d dpsi = D(psi);
    const double kinetic_energy = 0.5*inner(dpsi,dpsi);

    real_function_1d rho = square(psi);
    double potential_energy = inner(Vnuc,rho);

    double total_energy = kinetic_energy + potential_energy;
    double virial = potential_energy / kinetic_energy;

    if (world.rank() == 0) {
        print("            Kinetic energy ", kinetic_energy);
        print("          Potential energy ", potential_energy);
        print("              Total energy ", total_energy);
        print("                    Virial ", virial);
    }

    world.gop.fence();

    finalize();
    return 0;
}
