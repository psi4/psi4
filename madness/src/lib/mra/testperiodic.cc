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

/// \file testperiodic.cc
/// \brief test periodic convolutiosn

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/operator.h>
#include <constants.h>

using namespace madness;

typedef Vector<double,3> coordT;

const char* status[2] = {"FAIL !!!!!","PASS"};

const double L = 3.0;  // [-L,L]
const int nwave = 2;

double source(const coordT& r) {
    return cos(nwave*constants::pi*r[0]/L)*cos(nwave*constants::pi*r[1]/L)*cos(nwave*constants::pi*r[2]/L);
}

double source1(const coord_1d& r) {
    return cos(nwave*constants::pi*r[0]/L);
}

double u(const double x, const double expnt) {
    double fac = nwave*constants::pi/(2.0*L);
    return exp(-fac*fac/expnt) * cos(nwave*constants::pi*x/L);
}

// The electrostatic potential due to cos(n*pi*x/L)*(ditto z)*(ditto y)
// is ((4*L*L)/(3*n*n*pi))*cos(n*pi*x/L)
double potential(const coordT& r) {
    //const double fac = 1.0/(3*nwave*nwave*constants::pi);
    double fac = 4*L*L/(3*nwave*nwave*constants::pi);
    return source(r)*fac;
}

void test_periodic(World& world) {

    // the convolution of cos(n*pi*x/L) with sqrt(a/pi)*exp(-a*x^2)
    // is exp(-n^2*pi^2/(4*a*L^2)) * cos(n*pi*x/L)

    const long k = 14;
    const double thresh = 1e-12;
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_cubic_cell(-L,L);
    FunctionDefaults<3>::set_thresh(thresh);
    //FunctionDefaults<3>::set_initial_level(3);

    Function<double,3> f = FunctionFactory<double,3>(world).f(source);
    //f.truncate();

    std::vector< std::shared_ptr< Convolution1D<double> > > ops(1);

    std::cout.precision(10);
    double width = 2*L;
    for (int i=-4; i<=20; ++i) {
        double expnt = pow(2.0,double(i));
        double expnt_sim = expnt*width*width;
        double coeff = sqrt(expnt/constants::pi);
        double coeff_sim = coeff*width;
        ops[0].reset(new GaussianConvolution1D<double>(k, coeff_sim, expnt_sim, 0, true));

        SeparatedConvolution<double,3> op(world, ops);

        Function<double,3> opf = op(f);

        coordT r0(0);
        coordT r1(L-0.1);

        opf.reconstruct();
        f.reconstruct();

        double exact0 = u(r0[0],expnt); exact0 = exact0*exact0*exact0;
        double exact1 = u(r1[0],expnt); exact1 = exact1*exact1*exact1;
        double err0 = fabs(opf(r0)-exact0);
        double err1 = fabs(opf(r1)-exact1);

        print("exponent", expnt, err0, err1, status[err0<1e-10 && err1<1e-10]);

        // print(i, expnt, r0, f(r0), source(r0), opf(r0), exact0, );
        // print(i, expnt, r1, f(r1), source(r1), opf(r1), exact1, opf(r1)-exact1);
        // char fname[256];
        // sprintf(fname,"plot-%d.dat",i+1);
        // coordT lo(0.0); lo[0]=-L;
        // coordT hi(0.0); hi[0]= L;
        // plot_line(fname, 1001, lo, hi, f, opf);
    }

    world.gop.fence();

}


void test_periodic1(World& world) {

    // the convolution of cos(n*pi*x/L) with sqrt(a/pi)*exp(-a*x^2)
    // is exp(-n^2*pi^2/(4*a*L^2)) * cos(n*pi*x/L)

    const long k = 14;
    const double thresh = 1e-12;
    FunctionDefaults<1>::set_k(k);
    FunctionDefaults<1>::set_cubic_cell(-L,L);
    FunctionDefaults<1>::set_thresh(thresh);
    //FunctionDefaults<1>::set_initial_level(2);

    Function<double,1> f = FunctionFactory<double,1>(world).f(source1); //.norefine();
    //f.truncate();

    std::vector< std::shared_ptr< Convolution1D<double> > > ops(1);

    std::cout.precision(10);
    double width = 2*L;
    for (int i=-4; i<=20; ++i) {
        double expnt = pow(2.0,double(i));
        double expnt_sim = expnt*width*width;
        double coeff = sqrt(expnt/constants::pi);
        double coeff_sim = coeff*width;
        ops[0].reset(new GaussianConvolution1D<double>(k, coeff_sim, expnt_sim, 0, true));

        SeparatedConvolution<double,1> op(world, ops);

        Function<double,1> opf = op(f);

        coord_1d r0(0);
        coord_1d r1(L-0.1);

        opf.reconstruct();
        f.reconstruct();

        double exact0 = u(r0[0],expnt);
        double exact1 = u(r1[0],expnt);
        double err0 = fabs(opf(r0)-exact0);
        double err1 = fabs(opf(r1)-exact1);

        print("exponent", expnt, err0, err1, status[err0<6e-10 && err1<6e-10]);

        // print(i, expnt, r0, f(r0), source1(r0), opf(r0), exact0, opf(r0)-exact0);
        // print(i, expnt, r1, f(r1), source1(r1), opf(r1), exact1, opf(r1)-exact1);
        // char fname[256];
        // sprintf(fname,"plot-%d.dat",i+1);
        // coord_1d lo(0.0); lo[0]=-L;
        // coord_1d hi(0.0); hi[0]= L;
        // plot_line(fname, 1001, lo, hi, f, opf);
    }

    world.gop.fence();

}

void test_periodic2(World& world) {
    const long k = 10;
    const double thresh = 1e-8;
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_cubic_cell(-L,L);
    FunctionDefaults<3>::set_thresh(thresh);

    Function<double,3> f = FunctionFactory<double,3>(world).f(source);
    f.truncate();

    SeparatedConvolution<double,3> op = CoulombOperator(world, 1e-6, thresh);
    std::cout.precision(10);

    Function<double,3> opf = op(f);
    opf.reconstruct();

    //print("i,value,exact,relerr");
    for (int i=0; i<101; ++i) {
        coordT r = coordT(-L + i*2*L/100.0);
        double value = opf(r);
        double exact = potential(r);
        double relerr = fabs((value-exact)/exact);
        print(i,value,exact,relerr,status[relerr<3e-7]);
    }

    world.gop.fence();

}



int main(int argc, char**argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);

    FunctionDefaults<1>::set_bc(BC_PERIODIC);
    FunctionDefaults<3>::set_bc(BC_PERIODIC);

    try {

        print("1D gaussians");
        test_periodic1(world);
        print("\n3D gaussians");
        test_periodic(world);
        print("\n3D coulomb");
        test_periodic2(world);

    }
    catch (const MPI::Exception& e) {
        //        print(e);
        error("caught an MPI exception");
    }
    catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    }
    catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
    }
    catch (char* s) {
        print(s);
        error("caught a c-string exception");
    }
    catch (const char* s) {
        print(s);
        error("caught a c-string exception");
    }
    catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    }
    catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    }
    catch (...) {
        error("caught unhandled exception");
    }

    world.gop.fence();
    MPI::Finalize();

    return 0;
}
