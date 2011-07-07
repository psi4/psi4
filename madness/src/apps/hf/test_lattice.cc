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
#include "poperator.h"

using namespace madness;

typedef Vector<double,3> coordT3d;
typedef Vector<double,1> coordT1d;

const double L = 15.0;
const double N = 2.0;

//*****************************************************************************
// Test function for the periodic BSH operator.
static double rho_bsh_func3d(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
  double twopi = 2 * WST_PI;
  return  cos(twopi*x/L) * cos(twopi*y/L) * cos(twopi*z/L);
}
//*****************************************************************************

//*****************************************************************************
static double func_ones(const coordT3d& r)
{
  return 1.0;
}
//*****************************************************************************

//*****************************************************************************
// Test charge density for the periodic coulomb operator.
static double_complex rho_coulomb_func3d(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
  double npi = N * WST_PI;
  return  cos(npi*x/L) * cos(npi*y/L) * cos(npi*z/L);
}
//*****************************************************************************

//*****************************************************************************
// complex exponential to test for george
static std::complex<double> complex_exp_1d(const coordT1d r)
{
  const double x=r[0];
  double npi = N * WST_PI;
  return  exp(std::complex<double>(0.0,npi*x/L));
}
//*****************************************************************************

//*****************************************************************************
// Test function used to check the consistency of the coulomb operator and the
// periodic coulomb operator. We must insure that the total charge density is zero.
static double rho_gaussian_func3d(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
  const double expnt1 = 1.0;
  const double expnt2 = 2.0;
  const double coeff1 = pow(expnt1/WST_PI, 1.5);
  const double coeff2 = pow(expnt2/WST_PI, 1.5);
  double piece1 = coeff1 * exp(-expnt1 * (x*x + y*y + z*z));
  double piece2 = -coeff2 * exp(-expnt2 * (x*x + y*y + z*z));

  //  return piece1-(1/L/L/L);
    return piece1;
}
//*****************************************************************************

//*****************************************************************************
// Test function for the periodic coulomb operator.
static double_complex phi_coulomb_func3d(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
  double npi = N * WST_PI;
  double threepi = 3 * WST_PI;
  return  (4.0*L*L/(threepi*N*N)) * cos(npi*x/L) * cos(npi*y/L) * cos(npi*z/L);
}
//*****************************************************************************

////*****************************************************************************
//static double phi_func3d(const coordT3d& r)
//{
//  const double x=r[0], y=r[1], z=r[2];
//  double twopi = 2 * WST_PI;
//  double threepi = 3 * WST_PI;
//  return  (1/threepi) * cos(twopi*x) * cos(twopi*y) * cos(twopi*z);
//}
////*****************************************************************************

//*****************************************************************************
// Test function for periodic BSH operator.
static double phi_bsh_func3d(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
  double twopi = 2 * WST_PI;
  double sixteenpisquared = 16.0 * WST_PI * WST_PI;
  return  (L*L/sixteenpisquared) * cos(twopi*x/L) * cos(twopi*y/L) * cos(twopi*z/L);
}
//*****************************************************************************

//*****************************************************************************
template <typename Q, int NDIM>
Q laplacian(const Q& f, bool periodic = false)
{
  // Check for periodic boundary conditions
  Tensor<int> oldbc = FunctionDefaults<NDIM>::get_bc();
  if (periodic)
  {
    Tensor<int> bc(NDIM,2);
    bc(___) = 1;
    FunctionDefaults<NDIM>::set_bc(bc);
  }
  else
  {
    Tensor<int> bc(NDIM,2);
    bc(___) = 0;
    FunctionDefaults<NDIM>::set_bc(bc);
  }
  // Do calculation
  Q lapf = diff(diff(f,0),0);
  for (int i=1; i<NDIM; ++i) lapf += diff(diff(f,i),i);
  // Restore previous boundary conditions
  FunctionDefaults<NDIM>::set_bc(oldbc);
  return lapf;
};
//*****************************************************************************

//*****************************************************************************
template <typename Q> class wstFunctor
{
public:
  int kmax;
  double coeff, expnt;
  wstFunctor(int kmax, double coeff, double expnt) :
    kmax(kmax), coeff(coeff), expnt(expnt)
  {
  }
  Q operator()(double x) const
  { //x \in [-1,1] as checked
    Q sx0 = 0.0;
    for (int ki = -kmax; ki <= kmax; ki++)
    {
      double k = (double) ki;
      sx0 += exp(-expnt*(x-k)*(x-k));
    }
    return sx0*coeff;
  }
};
//*****************************************************************************

////*****************************************************************************
//struct PeriodicConditionalRefineTest
//{
//  bool operator() (const Key<3>& key, const Tensor<double>& t) const
//  {
//    // Is the box above the level from where we want to refine?
//    int n = key.level();
//    if (n >= 6) return false;
//    // Are we on the boundary?
//    Translation l1 = key.translation()[0];
//    Translation l2 = key.translation()[1];
//    Translation l3 = key.translation()[2];
//    Translation maxl = (1ul<<n) - 1;
//    if ((l1 == 0 || l1 == maxl) || (l2 == 0 || l2 == maxl) || (l3 == 0 || l3 == maxl))
//    {
//      print("Refining ...\n", key);
//      return true;
//    }
//    return false;
//  }
//
//  template <typename Archive>
//  void serialize(const Archive& arch) {}
//};
////*****************************************************************************

////*****************************************************************************
//// This code tests the convolution of a periodic charge density with just one gaussian
//// with periodic boundary conditions.
//void testPeriodicGaussian(World& world, double coeff, double expnt,
//    int lmax, int k, double thresh, double* data)
//{
//  // Function defaults
//  FunctionDefaults<3>::set_k(k);
//  FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);
//  FunctionDefaults<3>::set_thresh(thresh);
//
//  // Test function
//  Function<double,3> rho = FunctionFactory<double,3>(world).f(rho_coulomb_func3d);
//
//  // Create operator
//  std::vector< std::shared_ptr< Convolution1D<double> > > ops(1);
//  ops[0] = std::shared_ptr< Convolution1D<double> >(new PeriodicGaussianConvolution1D<double>(k, 10, coeff*L, expnt*L*L));
//  SeparatedConvolution<double,3> op(world, k, ops);
//
//  // Apply operator
//  Function<double,3> phi = apply(op, rho);
//
//  for (int i=0; i<21; i++)
//  {
//    double start = L/2;
//    double step = L / 20;
//    coordT3d p(-start + i*step);
//    printf("%.2f\t\t%.8f\t%.8f\t%.8f\n", p[0], phi(p), data[i], phi(p) / data[i]);
//    //printf("%.2f\t\t%.8f\t%.8f\t%.8f\n", p[0], phi(p), data[i], fabs(phi(p) - data[i]));
//  }
//}
////*****************************************************************************

////*****************************************************************************
//void testSinglePeriodicGaussians(int argc, char** argv)
//{
//  MPI::Init(argc, argv);
//  World world(MPI::COMM_WORLD);
//  startup(world,argc,argv);
//
//  double maple_data_2500[21] =
//    {
//      -44.02214685,
//      -37.86955441,
//      -23.31010082,
//      -8.939789115,
//      -1.299027398,
//      0.0,
//      1.299027398,
//      8.939789115,
//      23.31010082,
//      37.86955441,
//      44.02214685,
//      37.86955441,
//      23.31010082,
//      8.939789115,
//      1.299027398,
//      0.0,
//      -1.299027398,
//      -8.939789115,
//      -23.31010082,
//      -37.86955441,
//      -44.02214685
//    };
//
//  double maple_data_25000[21] =
//    {
//      -1.407020542,
//      -1.210373524,
//      -0.7450293332,
//      -0.2857304296,
//      -0.04151906172,
//      0.0,
//      0.04151906172,
//      0.2857304296,
//      0.7450293332,
//      1.210373524,
//      1.407020542,
//      1.210373524,
//      0.7450293332,
//      0.2857304296,
//      0.04151906172,
//      0.0,
//      -0.04151906172,
//      -0.2857304296,
//      -0.7450293332,
//      -1.210373524,
//      -1.407020542
//    };
//
//  double maple_data_55000[21] =
//    {
//      -0.4314663948,
//      -0.3711640906,
//      -0.2284651223,
//      -0.08761995618,
//      -0.01273192489,
//      0.0,
//      0.01273192489,
//      0.08761995618,
//      0.2284651223,
//      0.3711640906,
//      0.4314663948,
//      0.3711640906,
//      0.2284651223,
//      0.08761995618,
//      0.01273192489,
//      0.0,
//      -0.01273192489,
//      -0.08761995618,
//      -0.2284651223,
//      -0.3711640906,
//      -0.4314663948
//    };
//
//  double maple_data_95000[21] =
//    {
//      -0.1901095982,
//      -0.1635396336,
//      -0.1006646476,
//      -0.03860647055,
//      -0.005609848544,
//      0.0,
//      0.005609848544,
//      0.03860647055,
//      0.1006646476,
//      0.1635396336,
//      0.1901095982,
//      0.1635396336,
//      0.1006646476,
//      0.03860647055,
//      0.005609848544,
//      0.0,
//      -0.005609848544,
//      -0.03860647055,
//      -0.1006646476,
//      -0.1635396336,
//      -0.1901095982
//    };
//
//  int k = 8;
//  double thresh = 1e-6;
//  printf("\nTesting with exponent = 2500\n\n");
//  testPeriodicGaussian(world, 100, 2500, 16, k, thresh, &maple_data_2500[0]);
//  printf("\nTesting with exponent = 25000\n\n");
//  testPeriodicGaussian(world, 100, 25000, 16, k, thresh, &maple_data_25000[0]);
//  printf("\nTesting with exponent = 55000\n\n");
//  testPeriodicGaussian(world, 100, 55000, 16, k, thresh, &maple_data_55000[0]);
//  printf("\nTesting with exponent = 95000\n\n");
//  testPeriodicGaussian(world, 100, 95000, 16, k, thresh, &maple_data_95000[0]);
//  MPI::Finalize();
//}
////*****************************************************************************

////*****************************************************************************
//void testSinglePeriodicGaussians_L10(int argc, char** argv)
//{
//  MPI::Init(argc, argv);
//  World world(MPI::COMM_WORLD);
//  startup(world,argc,argv);
//
//  double maple_data_2500[21] =
//    {
//      -1.407020542,
//      -1.210373524,
//      -0.7450293332,
//      -0.2857304296,
//      -0.04151906172,
//      0.0,
//      0.04151906172,
//      0.2857304296,
//      0.7450293332,
//       1.210373524,
//       1.407020542,
//       1.210373524,
//      0.7450293332,
//      0.2857304296,
//      0.04151906172,
//      0.0,
//      -0.04151906172,
//      -0.2857304296,
//      -0.7450293332,
//      -1.210373524,
//      -1.407020542
//    };
//
//  double maple_data_25000[21] =
//    {
//      -1.407020542,
//      -1.210373524,
//      -0.7450293332,
//      -0.2857304296,
//      -0.04151906172,
//      0.0,
//      0.04151906172,
//      0.2857304296,
//      0.7450293332,
//      1.210373524,
//      1.407020542,
//      1.210373524,
//      0.7450293332,
//      0.2857304296,
//      0.04151906172,
//      0.0,
//      -0.04151906172,
//      -0.2857304296,
//      -0.7450293332,
//      -1.210373524,
//      -1.407020542
//    };
//
//  double maple_data_55000[21] =
//    {
//      -0.4314663948,
//      -0.3711640906,
//      -0.2284651223,
//      -0.08761995618,
//      -0.01273192489,
//      0.0,
//      0.01273192489,
//      0.08761995618,
//      0.2284651223,
//      0.3711640906,
//      0.4314663948,
//      0.3711640906,
//      0.2284651223,
//      0.08761995618,
//      0.01273192489,
//      0.0,
//      -0.01273192489,
//      -0.08761995618,
//      -0.2284651223,
//      -0.3711640906,
//      -0.4314663948
//    };
//
//  double maple_data_95000[21] =
//    {
//      -0.1901095982,
//      -0.1635396336,
//      -0.1006646476,
//      -0.03860647055,
//      -0.005609848544,
//      0.0,
//      0.005609848544,
//      0.03860647055,
//      0.1006646476,
//      0.1635396336,
//      0.1901095982,
//      0.1635396336,
//      0.1006646476,
//      0.03860647055,
//      0.005609848544,
//      0.0,
//      -0.005609848544,
//      -0.03860647055,
//      -0.1006646476,
//      -0.1635396336,
//      -0.1901095982
//    };
//
//  int k = 8;
//  double thresh = 1e-6;
//  printf("\nTesting with exponent = 2500\n\n");
//  testPeriodicGaussian(world, 100, 2500, 16, k, thresh, &maple_data_2500[0]);
////  printf("\nTesting with exponent = 25000\n\n");
////  testPeriodicGaussian(world, 100, 25000, 16, k, thresh, &maple_data_25000[0]);
////  printf("\nTesting with exponent = 55000\n\n");
////  testPeriodicGaussian(world, 100, 55000, 16, k, thresh, &maple_data_55000[0]);
////  printf("\nTesting with exponent = 95000\n\n");
////  testPeriodicGaussian(world, 100, 95000, 16, k, thresh, &maple_data_95000[0]);
//  MPI::Finalize();
//}
////*****************************************************************************

//*****************************************************************************
void testPeriodicCoulomb3d(int argc, char**argv)
{
  MPI::Init(argc, argv);
  World world(MPI::COMM_WORLD);
  startup(world,argc,argv);

  // Function defaults
  int k = 10;
  double thresh = 1e-8;
  FunctionDefaults<3>::set_k(k);
  FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);
  FunctionDefaults<3>::set_thresh(thresh);

  // Create test charge density and the exact solution to Poisson's equation
  // with said charge density
  printf("building rho ...\n\n");
  Function<double_complex,3> rho = FunctionFactory<double_complex,3>(world).f(rho_coulomb_func3d);
  rho.truncate();
  printf("building phi_exact ...\n\n");
  Function<double_complex,3> phi_exact = FunctionFactory<double_complex,3>(world).f(phi_coulomb_func3d);

  // Create operator and apply
  Tensor<double> cellsize = FunctionDefaults<3>::get_cell_width();
  SeparatedConvolution<double_complex,3> op = PeriodicCoulombOp<double_complex,3>(world, k, 1e-8, thresh, cellsize);
  printf("applying operator ...\n\n");
  Function<double_complex,3> phi_test = apply(op, rho);

  // Apply the laplacian operator the phi_test and see if we get rho back
//  printf("applying laplacian operator to phi_test ...\n\n");
//  Function<double,3> rho_test = laplacian<Function<double,3>, 3>(phi_test, false);
//  rho_test.scale(-1.0/4/WST_PI);
//  Function<double,3> rho_diff = rho - rho_test;

//  rho.reconstruct();
//  rho_test.reconstruct();
//  rho_diff.reconstruct();
  phi_exact.reconstruct();
  phi_test.reconstruct();

  double bstep = L / 100.0;
  for (int i=0; i<101; i++)
  {
    coordT3d p(-L/2 + i*bstep);
    double error = abs(phi_exact(p) - phi_test(p));
    printf("%.2f\t\t%.8f\t%.8f\t%.8f\t%.8f\n", p[0], abs(phi_exact(p)), abs(phi_test(p)), error, error / abs(phi_exact(p)));
    //printf("%.2f\t\t%.8f\t%.8f\t%.8f\t%.8f\n", p[0], rho(p), rho_test(p), rho_diff(p), 0.0);
  }

  // Plot to OpenDX
  std::vector<long> npt(3,101);
  Function<double_complex,3> phi_diff = phi_exact - phi_test;
  plotdx(phi_test, "phitest.dx", FunctionDefaults<3>::get_cell(), npt);
  plotdx(phi_exact, "phiexact.dx", FunctionDefaults<3>::get_cell(), npt);
  plotdx(phi_diff, "phidiff.dx", FunctionDefaults<3>::get_cell(), npt);

  MPI::Finalize();
}
//*****************************************************************************

////*****************************************************************************
//// This method is just a simple test of the coulomb operator without periodic
//// boundary conditions.
//void testNonPeriodicCoulomb3d(int argc, char**argv)
//{
//  MPI::Init(argc, argv);
//  World world(MPI::COMM_WORLD);
//  startup(world,argc,argv);
//
//  // Function defaults
//  int k = 14;
//  double thresh = 1e-12;
//  double eps = 1e-12;
//  FunctionDefaults<3>::set_k(k);
//  FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);
//  FunctionDefaults<3>::set_thresh(thresh);
//
//  // Create test charge density and the exact solution to Poisson's equation
//  // with said charge density
//  printf("building rho ...\n\n");
//  Function<double,3> rho = FunctionFactory<double,3>(world).f(rho_gaussian_func3d);
//
//  // Create operator and apply
//  SeparatedConvolution<double,3> op = CoulombOperator<double>(world, k, 1e-6, eps);
//  printf("applying operator ...\n\n");
//  Function<double,3> phi_test = apply(op, rho);
//
//  // Apply the laplacian operator the phi_test and see if we get rho back
//  printf("applying laplacian operator to phi_test ...\n\n");
//  Function<double,3> rho_test = laplacian<Function<double,3>, 3>(phi_test, false);
//  rho_test.scale(-1.0/4/WST_PI);
//  Function<double,3> rho_diff = rho - rho_test;
//
//  rho.reconstruct();
//  rho_test.reconstruct();
//  rho_diff.reconstruct();
//
//  double bstep = L / 100.0;
//  for (int i=0; i<101; i++)
//  {
//    coordT3d p(-L/2 + i*bstep);
//    printf("%.2f\t\t%.8f\t%.8f\t%.8f\t%.8f\n", p[0], rho(p), rho_test(p), rho_diff(p), rho_test(p) / rho(p));
//  }
//
//  MPI::Finalize();
//}
////*****************************************************************************

////*****************************************************************************
//void testPeriodicBSH3d(int argc, char**argv)
//{
//  MPI::Init(argc, argv);
//  World world(MPI::COMM_WORLD);
//  startup(world,argc,argv);
//
//  // Function defaults
//  int k = 8;
//  double thresh = 1e-6;
//  double eps = 1e-6;
//  FunctionDefaults<3>::set_k(k);
//  FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);
//  FunctionDefaults<3>::set_thresh(thresh);
//
//  // Create test charge density and the exact solution to Poisson's equation
//  // with said charge density
//  Function<double,3> rho = FunctionFactory<double,3>(world).f(rho_bsh_func3d);
//  Function<double,3> phi_exact = FunctionFactory<double,3>(world).f(phi_bsh_func3d);
//
//  // Create operator and apply
//  double twopi = 2 * WST_PI;
//  Tensor<double> cellsize = FunctionDefaults<3>::get_cell_width();
//  SeparatedConvolution<double,3> op = PeriodicBSHOp<double,3>(world, twopi/L, k, 1e-8, eps, cellsize);
//  Function<double,3> phi_test = apply(op, rho);
//
//  double bstep = L / 100.0;
//  for (int i=0; i<101; i++)
//  {
//    coordT3d p(-L/2 + i*bstep);
//    printf("%.2f\t\t%.8f\t%.8f\t%.8f\t%.8f\n", p[0], phi_exact(p), phi_test(p), phi_exact(p)/phi_test(p), fabs(phi_exact(p)-phi_test(p)));
//  }
//
//  // Plot to OpenDX
////  printf("plotting to openDX ...\n\n");
////  vector<long> npt(3,101);
////  Function<double,3> phi_diff = phi_exact - phi_test;
////  plotdx(phi_test, "phitest.dx", FunctionDefaults<3>::get_cell(), npt);
////  plotdx(phi_exact, "phiexact.dx", FunctionDefaults<3>::get_cell(), npt);
////  plotdx(phi_diff, "phidiff.dx", FunctionDefaults<3>::get_cell(), npt);
////  printf("finished plotting to openDX ...\n\n");
//
//  MPI::Finalize();
//  printf("done!\n\n");
//}
////*****************************************************************************

////*****************************************************************************
//// This function test both the periodic and non-periodic versions of the Coulomb
//// operator. In order to make this test valid set L to a high value so that
//// charge distribution should not be able to see its neighbor.
//void testPeriodicCoulomb3d_gauss(int argc, char**argv)
//{
//  MPI::Init(argc, argv);
//  World world(MPI::COMM_WORLD);
//  startup(world,argc,argv);
//
//  // Function defaults
//  int k = 10;
//  double thresh = 1e-8;
//  double eps = 1e-8;
//  FunctionDefaults<3>::set_k(k);
//  FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);
//  FunctionDefaults<3>::set_thresh(thresh);
//
//  // Test for a noncubic system
//  //  Tensor<double> csize(3,3);
////  double sz = L/2;
////  csize(0,0) = -1.2 * sz; csize(0,1) = 1.2 * sz;
////  csize(1,0) = -1.5 * sz; csize(1,1) = 1.5 * sz;
////  csize(2,0) = -1.0 * sz; csize(2,1) = 1.0 * sz;
////  FunctionDefaults<3>::set_cell(csize);
//
//  // Create test charge density and the exact solution to Poisson's equation
//  // with said charge density
//  printf("building gaussian charge distribution ...\n\n");
//  Function<double,3> rho = FunctionFactory<double,3>(world).f(rho_gaussian_func3d);
//  rho.truncate();
//
//  // Average value of the test function
//  Function<double,3> ones = FunctionFactory<double,3>(world).f(func_ones);
//  //double avgval = inner(rho,ones);
//  double avgval = rho.trace();
//  printf("Average value of rho is %.8f\n\n", avgval);
//
//  // Create operator and apply
//  Tensor<double> cellsize = FunctionDefaults<3>::get_cell_width();
//  SeparatedConvolution<double,3> pop = PeriodicCoulombOp<double,3>(world, k,1e-10, eps, cellsize);
//  printf("applying periodic operator ...\n\n");
//  Function<double,3> phi_periodic = apply(pop, rho);
//  SeparatedConvolution<double,3> op = CoulombOperator<double>(world, FunctionDefaults<3>::get_k(),
//      1e-10, eps);
//  printf("applying non-periodic operator ...\n\n");
//  Function<double,3> phi_nonperiodic = apply(op, rho);
//  phi_nonperiodic.truncate();
//  phi_periodic.truncate();
//
//  phi_periodic.reconstruct();
//  phi_nonperiodic.reconstruct();
//
//  double bstep = L / 100.0;
//  for (int i=0; i<101; i++)
//  {
//    coordT3d p(-L/2 + i*bstep);
//    double error = fabs(phi_periodic(p) - phi_nonperiodic(p));
//    printf("%.2f\t\t%.8f\t%.8f\t%.8f\t%.8f\n", p[0], phi_periodic(p), phi_nonperiodic(p), error, phi_periodic(p) / phi_nonperiodic(p));
//  }
//
//  // Plot to OpenDX
////  vector<long> npt(3,101);
////  Function<double,3> phi_diff = phi_periodic - phi_nonperiodic;
////  plotdx(phi_periodic, "phiperiodic.dx", FunctionDefaults<3>::get_cell(), npt);
////  plotdx(phi_nonperiodic, "phinonperiodic.dx", FunctionDefaults<3>::get_cell(), npt);
////  plotdx(phi_diff, "phidiff.dx", FunctionDefaults<3>::get_cell(), npt);
//
//  MPI::Finalize();
//}
////*****************************************************************************

////*****************************************************************************
//// This function test both the periodic and non-periodic versions of the BSH
//// operator. In order to make this test valid set L to a high value so that
//// charge distribution should not be able to see its neighbor.
//void testPeriodicBSH3d_gauss(int argc, char**argv)
//{
//  MPI::Init(argc, argv);
//  World world(MPI::COMM_WORLD);
//  startup(world,argc,argv);
//
//  // Function defaults
//  int k = 10;
//  double thresh = 1e-8;
//  double eps = 1e-8;
//  FunctionDefaults<3>::set_k(k);
//  FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);
//  FunctionDefaults<3>::set_thresh(thresh);
//
//  // Test for a noncubic system
//  //  Tensor<double> csize(3,3);
////  double sz = L/2;
////  csize(0,0) = -1.2 * sz; csize(0,1) = 1.2 * sz;
////  csize(1,0) = -1.5 * sz; csize(1,1) = 1.5 * sz;
////  csize(2,0) = -1.0 * sz; csize(2,1) = 1.0 * sz;
////  FunctionDefaults<3>::set_cell(csize);
//
//  // Create test charge density and the exact solution to Poisson's equation
//  // with said charge density
//  printf("building gaussian charge distribution ...\n\n");
//  Function<double,3> rho = FunctionFactory<double,3>(world).f(rho_gaussian_func3d);
//  rho.truncate();
//
//  // Average value of the test function
//  Function<double,3> ones = FunctionFactory<double,3>(world).f(func_ones);
//  double avgval = inner(rho,ones);
//  printf("Average value of rho is %.8f\n\n", avgval);
//
//  // Create operator and apply
//  Tensor<double> cellsize = FunctionDefaults<3>::get_cell_width();
//  SeparatedConvolution<double,3> pop = PeriodicBSHOp<double,3>(world, -3.0, k,1e-6, eps, cellsize);
//  printf("applying periodic operator ...\n\n");
//  Function<double,3> phi_periodic = apply(pop, rho);
//  SeparatedConvolution<double,3> op = BSHOperator3D<double>(world, -3.0, FunctionDefaults<3>::get_k(),
//      1e-6, thresh);
//  printf("applying non-periodic operator ...\n\n");
//  Function<double,3> phi_nonperiodic = apply(op, rho);
//
//  double bstep = L / 100.0;
//  for (int i=0; i<101; i++)
//  {
//    coordT3d p(-L/2 + i*bstep);
//    double error = fabs(phi_periodic(p) - phi_nonperiodic(p));
//    printf("%.2f\t\t%.8f\t%.8f\t%.8f\t%.8f\n", p[0], phi_periodic(p), phi_nonperiodic(p), error, phi_periodic(p) / phi_nonperiodic(p));
//  }
//
////  // Plot to OpenDX
////  vector<long> npt(3,101);
////  Function<double,3> phi_diff = phi_periodic - phi_nonperiodic;
////  plotdx(phi_periodic, "phiperiodic.dx", FunctionDefaults<3>::get_cell(), npt);
////  plotdx(phi_nonperiodic, "phinonperiodic.dx", FunctionDefaults<3>::get_cell(), npt);
////  plotdx(phi_diff, "phidiff.dx", FunctionDefaults<3>::get_cell(), npt);
//
//  MPI::Finalize();
//}
////*****************************************************************************

//*****************************************************************************
int main(int argc, char**argv)
{
  BoundaryConditions<3> bc(BC_PERIODIC);
  FunctionDefaults<3>::set_bc(bc);
  //testPeriodicBSH3d(argc, argv);
  //testPeriodicCoulomb3d_gauss(argc, argv);
  //testNonPeriodicCoulomb3d(argc, argv);
  testPeriodicCoulomb3d(argc, argv);
  //testPeriodicBSH3d_gauss(argc, argv);
  //testSinglePeriodicGaussians(argc,argv);
  return 0;
}
//*****************************************************************************

