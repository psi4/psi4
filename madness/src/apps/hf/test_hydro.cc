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
//#include "hartreefock.h"

using std::cout;
using std::endl;

using namespace madness;

const double PI = 3.1415926535897932384;

typedef Vector<double,3> coordT;

  //***************************************************************************
  template <typename T, int NDIM>
  class NuclearChargeDensityOp : public EigSolverOp<T,NDIM>
  {
  public:
    typedef Function<T,NDIM> funcT;
    //*************************************************************************
    // Constructor
    NuclearChargeDensityOp(World& world, funcT rhon, double coeff, double thresh)
    : EigSolverOp<T,NDIM>(world, coeff, thresh)
    {
      _cop =
        CoulombOperatorPtr<T,NDIM>(const_cast<World&>(world),
            FunctionDefaults<3>::get_k(), 1e-4, thresh);
      _rhon = rhon;
      // Apply operator to get potential
      _Vnuc = apply(*_cop, rhon);
    }
    //*************************************************************************

    //*************************************************************************
    ~NuclearChargeDensityOp()
    {
      delete _cop;
    }
    //*************************************************************************

    //*************************************************************************
    // Is there an orbitally-dependent term?
    virtual bool is_od() {return false;}
    //*************************************************************************

    //*************************************************************************
    // Is there a density-dependent term?
    virtual bool is_rd() {return true;}
    //*************************************************************************

    //*************************************************************************
    void prepare_op(Function<double,NDIM> rho) {}
    //*************************************************************************

    //*************************************************************************
    virtual funcT op_r(const funcT& rho, const funcT& psi)
    {
      return _Vnuc * psi;
    }
    //*************************************************************************

  private:
    //*************************************************************************
    funcT _rhon;
    //*************************************************************************

    //*************************************************************************
    SeparatedConvolution<T,NDIM>* _cop;
    //*************************************************************************

    //*************************************************************************
    funcT _Vnuc;
    //*************************************************************************
  };
  //***************************************************************************

//*****************************************************************************
static double rho_func_hydro(const coordT& rr)
{
  const double x=rr[0], y=rr[1], z=rr[2];
//  double e1 = 100.0;
//  double coeff = pow(e1/PI, 1.5);
//  return -1.0 * coeff * exp(-e1 * (x*x + y*y + z*z));
  double c = 0.1;
  double r = sqrt(x*x + y*y + z*z);
  r = r / c;
  const double RPITO1P5 = 0.1795871221251665617; // 1.0/Pi^1.5
  return ((-3.0/2.0+(1.0/3.0)*r*r)*exp(-r*r)+(-32.0+(256.0/3.0)*r*r)*exp(-4.0*r*r))*RPITO1P5/c/c/c;
}
//*****************************************************************************

//*****************************************************************************
double psi_func_hydro(const Vector<double,3>& r)
{
  double x = r[0];
  double y = r[1];
  double z = r[2];
  return exp(-(x*x + y*y + z*z));
}
//*****************************************************************************

//*****************************************************************************
void test_hydro(World& world)
{
  // Box size
  double L = 40.0;

  // Function defaults
  int funck = 8;
  double thresh = 1e-6;
  FunctionDefaults<3>::set_k(funck);
  FunctionDefaults<3>::set_thresh(thresh);
  FunctionDefaults<3>::set_cubic_cell(-L/2, L/2);

  // Nuclear potential (Be)
  //const coordT origin(0.0);
  if (world.rank() == 0) madness::print("Creating Function object for nuclear charge density ...\n\n");
  Function<double,3> rhon = FunctionFactory<double,3>(world).f(rho_func_hydro).thresh(thresh).initial_level(4);
  rhon.truncate();
  double rhontrace = rhon.trace();
  if (world.rank() == 0) madness::print("trace of rho is ", rhontrace);

  // Guess for the wavefunctions
  if (world.rank() == 0) madness::print("Creating wavefunction's ...\n\n");
  Function<double,3> psi = FunctionFactory<double,3>(world).f(psi_func_hydro);
  psi.truncate();
  psi.scale(1.0/psi.norm2());
  std::vector<Function<double,3> > phis;
  phis.push_back(psi);
  // Create list of eigenvalues
  std::vector<double> eigs;
  eigs.push_back(-0.9);
  // Create ops
  std::vector<EigSolverOp<double,3>*> ops;
  ops.push_back(new DFTNuclearChargeDensityOp<double,3>(world, rhon, 1.0, thresh, true));
  // Create eigensolver
  if (world.rank() == 0) madness::print("Creating Eigensolver object...\n\n");
  // Constructor for non-periodic system
  EigSolver<double,3> solver(world, rhon, phis, eigs, ops, thresh, true);
  if (world.rank() == 0) madness::print("Diagonalizing Hamiltonian ...\n\n");
  solver.solve(15);

  double eval = solver.get_eig(0);
  Function<double,3> func = solver.get_phi(0);
  if (world.rank() == 0) printf("reconstructing func ...\n\n");

  MPI::Finalize();
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
    startup(world,argc,argv);
    test_hydro(world);
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
