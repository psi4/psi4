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
  \file examples/hatom_energy.cc
  \brief Compute the energy of the hydrogen atom ground state
  \defgroup hatom_energy Energy of the hydrogen atom ground state
  \ingroup examples

  This example computes the energy of the ground state of 
  the hydrogen atom as the expectation value of the
  non-relativistic Schr&ouml;dinger Hamiltonian within the 
  Born-Oppenheimer approximation.  Explicitly, 
  \f[
      E = \frac{\langle \psi \mid - \frac{1}{2}  \nabla^2  - \frac{1}{r} \mid \psi \rangle}{\langle \psi \mid \psi \rangle}
  \f]
  where the unnormalized wave function is
  \f[
     \psi(r) = e^{-r}
  \f]

  \par Implementation

  Two functions are required - the wave function and the potential. Note the small 
  constants introduced to eliminate the cusp in the wave function and the singularity
  in the potential. Due to the volume element this smoothing has negligible effect on
  the result (perturbation theory can make this rigorous).

  The wave function is exponentially decaying and has the value 2e-9 at r=20, so
  we pick this as a box size that is effectively infinite.

  Using integration by parts and noting that \f$ \psi(\infty)=0 \f$
  \f[
       \langle \psi \mid \nabla^2 \mid \psi \rangle = - \langle \nabla \psi \mid \nabla \psi \rangle
  \f]
  Finally, due to the spherical symmetry we only need compute one component of the Laplacian
  in Cartesian coordinates.  

  The exact answer is \f$-\frac{1}{2}\f$ in atomic units.
*/

#include <mra/mra.h>

using namespace madness;

double psi(const Vector<double,3>& r) {
  return exp(-sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]+1e-6));
}

double V(const Vector<double,3>& r) {
  return -1.0/sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]+1e-6);
}

int main(int argc, char**argv) {
  // Initialize the parallel programming environment
  initialize(argc,argv);
  World world(MPI::COMM_WORLD);

  // Load info for MADNESS numerical routines
  startup(world,argc,argv);
  std::cout.precision(8);

  // Use defaults for numerical functions except for user simulation volume
  FunctionDefaults<3>::set_cubic_cell(-20,20);

  Function<double,3> u = FunctionFactory<double,3>(world).f(psi);
  Function<double,3> v = FunctionFactory<double,3>(world).f(V);
  Function<double,3> vu = v*u;
  real_derivative_3d Dx = free_space_derivative<double,3>(world,0);
  Function<double,3> du = Dx(u);
  double KE = 3*0.5*(du.inner(du));
  double PE = vu.inner(u);
  double S = u.inner(u);

  if (world.mpi.Get_rank() == 0) {
    print("the overlap integral is",S);
    print("the kinetic energy integral",KE);
    print("the potential energy integral",PE);
    print("the total energy",(KE+PE)/S);
  }

  finalize();

  return 0;
}

