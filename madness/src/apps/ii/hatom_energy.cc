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

/// \file ii/hatom_energy.cc
/// \brief Compute the energy of the hydrogen atom ground state

#include <mra/mra.h>

using namespace madness;


double psi(const Vector<double,3>& r) {
  return exp(-sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]));
}

double V(const Vector<double,3>& r) {
  return -1.0/sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]+1e-8);
}

int main(int argc, char**argv) {
  // Initialize the parallel programming environment
  MPI::Init(argc, argv);
  World world(MPI::COMM_WORLD);
  
  // Load info for MADNESS numerical routines
  startup(world,argc,argv);
  
  // Setup defaults for numerical functions
  FunctionDefaults<3>::set_k(7);                 // Wavelet order
  FunctionDefaults<3>::set_thresh(1e-5);         // Accuracy
  FunctionDefaults<3>::set_refine(true);         // Enable adaptive refinement
  FunctionDefaults<3>::set_initial_level(2);     // Initial projection level
//  for (int i=0; i<3; i++) {
//    FunctionDefaults<3>::cell(i,0) = -20.0;   // User-defined volume
//    FunctionDefaults<3>::cell(i,1) =  20.0;
//  }
  FunctionDefaults<3>::set_cubic_cell(-20.0, 20.0);
  Function<double,3> u = FunctionFactory<double,3>(world).f(psi);
  Function<double,3> v = FunctionFactory<double,3>(world).f(V);
  Function<double,3> vu = v*u;
  Function<double,3> du = diff(u,0);
  double KE = 3*0.5*(du.inner(du));
  double PE = vu.inner(u);
  double S = u.inner(u);

  print("the overlap integral is",S);
  print("the kinetic energy integral",KE);
  print("the potential energy integral",PE);
  print("the total energy",(KE+PE)/S);
  
  MPI::Finalize();
  
  return 0;
}
