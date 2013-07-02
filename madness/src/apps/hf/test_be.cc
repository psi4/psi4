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
#include <iostream>

#include "dft.h"
#include "hartreefock.h"

using std::cout;
using std::endl;

using namespace madness;

const double PI = 3.1415926535897932384;

typedef Vector<double,3> coordT;

/// Returns radius for smoothing nuclear potential with energy precision eprec
//*****************************************************************************
static double smoothing_parameter(double Z, double eprec) {
    // The min is since asymptotic form not so good at low acc.
    // The 2 is from two electrons in 1s closed shell.
    if (Z == 0.0) return 1.0;
    double Z5 = Z*Z*Z*Z*Z;
    double c = pow(std::min(1e-3,eprec)/2.0/0.00435/Z5,1.0/3.0);
    return c;
}
//*****************************************************************************


/// Regularized 1/r potential.

/// Invoke as \c u(r/c)/c where \c c is the radius of the
/// smoothed volume.
//*****************************************************************************
static double smoothed_potential(double r) {
    const double THREE_SQRTPI = 5.31736155271654808184;
    double r2 = r*r, pot;
    if (r > 6.5){
        pot = 1.0/r;
    } else if (r > 1e-8){
        pot = erf(r)/r + (exp(-r2) + 16.0*exp(-4.0*r2))/(THREE_SQRTPI);
    } else{
        pot = (2.0 + 17.0/3.0)/sqrt(PI);
    }

    return pot;
}
//*****************************************************************************

//*****************************************************************************
static double psi_func_be1(const coordT& rr)
{
  const double x=rr[0], y=rr[1], z=rr[2];
  double r = sqrt(x*x+y*y+z*z);
  return exp(-4.0*r+1e-4);
}
//*****************************************************************************

//*****************************************************************************
static double psi_func_be2(const coordT& rr)
{
  const double x=rr[0], y=rr[1], z=rr[2];
  double r = sqrt(x*x+y*y+z*z);
  return (1.0 - 2.0*r*exp(-2.0*r));
}
//*****************************************************************************


//*****************************************************************************
static double V_func_be(const coordT& r)
{
  const double x=r[0], y=r[1], z=r[2];
  double rr = sqrt(x*x + y*y + z*z);
  double c = smoothing_parameter(4.0, 1e-7);
  return -4.0 * smoothed_potential(rr/c) / c;
}
//*****************************************************************************

//*****************************************************************************
static double rho_func_be(const coordT& rr)
{
  const double x=rr[0], y=rr[1], z=rr[2];
//  double e1 = 20.0;
//  double coeff = pow(e1/PI, 1.5);
//  return -4.0 * coeff * exp(-e1 * (x*x + y*y + z*z));
  double c = 0.1;
  double r = sqrt(x*x + y*y + z*z);
  r = r / c;
  const double RPITO1P5 = 0.1795871221251665617; // 1.0/Pi^1.5
  return 4.0 * ((-3.0/2.0+(1.0/3.0)*r*r)*exp(-r*r)+(-32.0+(256.0/3.0)*r*r)*exp(-4.0*r*r))*RPITO1P5/c/c/c;
}
//*****************************************************************************

//*****************************************************************************
void test_hf_be(World& world)
{
  //if (world.rank() == 0) cout << "Running test application HartreeFock ..." << endl;

  typedef Vector<double,3> coordT;
  typedef std::shared_ptr< FunctionFunctorInterface<double,3> > functorT;

  // Function defaults
  double bsize = 42.4;
  int funck = 6;
  double thresh = 1e-4;
  FunctionDefaults<3>::set_k(funck);
  FunctionDefaults<3>::set_thresh(thresh);
  FunctionDefaults<3>::set_autorefine(false);
  FunctionDefaults<3>::set_cubic_cell(-bsize, bsize);

  // Nuclear potential (Be)
  const coordT origin(0.0);
  if (world.rank() == 0) madness::print("Creating Function object for nuclear charge density  ...");
  Function<double,3> rhon = FunctionFactory<double,3>(world).f(rho_func_be).thresh(thresh).initial_level(4);
  Function<double,3> vnuc = FunctionFactory<double,3>(world).f(V_func_be).thresh(thresh);

  double rhontrace = rhon.trace();
  if (world.rank() == 0)
    printf("trace of rhon = %.8f\n\n", rhontrace);

  if (world.rank() == 0)
    cout << "Operating on nuclear charge density ..." << endl;
  SeparatedConvolution<double, 3> op = CoulombOperator<double, 3> (world,
    FunctionDefaults<3>::get_k(), 1e-8, thresh);
  Function<double, 3> V_from_rho_nuc = apply(op, rhon);
  if (world.rank() == 0)  printf("\n");
  double L = 2.0 * bsize;
  double bstep = L / 100.0;
  vnuc.reconstruct();
  V_from_rho_nuc.reconstruct();
  for (int i = 0; i < 101; i++)
  {
    coordT p(-L / 2 + i * bstep);
    double error = fabs(vnuc(p) - V_from_rho_nuc(p));
    if (world.rank() == 0)
      printf("%.2f\t\t%.8f\t%.8f\t%.8f\t%.8f\n", p[0], vnuc(p), V_from_rho_nuc(p),
          error, error / vnuc(p));
  }
  if (world.rank() == 0) printf("\n");


  // Guess for the wavefunctions
  if (world.rank() == 0) madness::print("Creating wavefunction's ...");
  Function<double,3> psi1 = FunctionFactory<double,3>(world).f(psi_func_be1);
  psi1.scale(1.0/psi1.norm2());
  Function<double,3> psi2 = FunctionFactory<double,3>(world).f(psi_func_be2);
  psi2.scale(1.0/psi2.norm2());

  // Create list of wavefunctions
  std::vector<Function<double,3> > phis;
  phis.push_back(psi1);
  phis.push_back(psi2);

  // Creat list of eigenvalues
  std::vector<double> eigs;
  eigs.push_back(-5.0);
  eigs.push_back(-0.5);

  // Create DFT object
  if (world.rank() == 0) madness::print("Creating DFT object...");
  //HartreeFock hf(world, Vnuc, phis, eigs, true, true, thresh);
  ElectronicStructureParams params;
  params.periodic = false;
  DFT<double,3> dftcalc(world, rhon, phis, eigs, params);
  if (world.rank() == 0) madness::print("Running DFT object...");
  dftcalc.solve(51);

}
//*****************************************************************************

#define TO_STRING(s) TO_STRING2(s)
#define TO_STRING2(s) #s

//*****************************************************************************
int main(int argc, char** argv)
{
  MPI::Init(argc, argv);
  World world(MPI::COMM_WORLD);
  if (world.rank() == 0)
  {
    print("");
    print("--------------------------------------------");
    print("   MADNESS", " multiresolution testsuite");
    print("--------------------------------------------");
    print("");
    print("   number of processors ...", world.size());
    print("    processor frequency ...", cpu_frequency());
    print("            host system ...", TO_STRING(HOST_SYSTEM));
    print("             byte order ...", TO_STRING(MADNESS_BYTE_ORDER));
    print("          configured by ...", MADNESS_CONFIGURATION_USER);
    print("          configured on ...", MADNESS_CONFIGURATION_HOST);
    print("          configured at ...", MADNESS_CONFIGURATION_DATE);
    print("                    CXX ...", MADNESS_CONFIGURATION_CXX);
    print("               CXXFLAGS ...", MADNESS_CONFIGURATION_CXXFLAGS);
#ifdef WORLD_WATCHDOG
    print("               watchdog ...", WATCHDOG_BARK_INTERVAL,
        WATCHDOG_TIMEOUT);
#endif
#ifdef OPTERON_TUNE
    print("             tuning for ...", "opteron");
#elif defined(CORE_DUO_TUNE)
    print("             tuning for ...", "core duo");
#else
    print("             tuning for ...", "core2");
#endif
#ifdef BOUNDS_CHECKING
    print(" tensor bounds checking ...", "enabled");
#endif
#ifdef TENSOR_INSTANCE_COUNT
    print("  tensor instance count ...", "enabled");
#endif
    print(" ");
  }

  try
  {
    printf("WSTHORNTON: Starting up the world ... \n");

    startup(world,argc,argv);
    if (world.rank() == 0) print("Initial tensor instance count", BaseTensor::get_instance_count());
    test_hf_be(world);
  }
  catch (const MPI::Exception& e)
  {
    //        print(e);
    error("caught an MPI exception");
  }
  catch (const madness::MadnessException& e)
  {
    print(e);
    error("caught a MADNESS exception");
  }
  catch (const madness::TensorException& e)
  {
    print(e);
    error("caught a Tensor exception");
  }
  catch (const char* s)
  {
    print(s);
    error("caught a string exception");
  }
  catch (const std::string& s)
  {
    print(s);
    error("caught a string (class) exception");
  }
  catch (const std::exception& e)
  {
    print(e.what());
    error("caught an STL exception");
  }
  catch (...)
  {
    error("caught unhandled exception");
  }

  if (world.rank() == 0)
    print("entering final fence");
  world.gop.fence();
  if (world.rank() == 0)
    print("done with final fence");
  if (world.rank() == 0)
    print("Final tensor instance count", BaseTensor::get_instance_count());
  MPI::Finalize();

  return 0;
}
//*****************************************************************************

