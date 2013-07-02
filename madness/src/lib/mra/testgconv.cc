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

/// \file mra/testgconv.cc
/// \brief Test convolution with Gaussian * polyn

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/operator.h>
#include <constants.h>

using namespace madness;

static const int k = 10;
static const double thresh = 1e-8;
static const double L = 17;

// exp(-r^2) / sqrt(pi) = normalized gaussian
double g(const coord_1d& r) {
    static const double fac = 1.0/sqrt(constants::pi);
    return exp(-r[0]*r[0]) * fac;
}

// g() conv g()
double gconvg(const coord_1d& r) {
    static const double fac = 1.0/sqrt(2.0*constants::pi);
    return exp(-r[0]*r[0]*0.5) * fac;
}

// sqrt(8)*x*exp(-x^2)
double h(const coord_1d& r) {
    static const double fac = sqrt(8.0);
    return exp(-r[0]*r[0]) * fac * r[0];
}

// g() conv h() == h() conv g()
double gconvh(const coord_1d& r) {
    return exp(-0.5*r[0]*r[0]) * r[0];
}


void test_gconv(World& world) {
    coord_1d origin(0.0), lo(-L), hi(L);
    double width = 2.0*L;

    if (world.rank() == 0) print("Test gconv operation");

    real_function_1d f = real_factory_1d(world).f(g);
    print("error in integral(g) ", f.trace()-1.0);

    std::vector< std::shared_ptr< Convolution1D<double> > > ops(1);
    ops[0].reset(new GaussianConvolution1D<double>(k, width/sqrt(constants::pi),
            width*width, 0, false));
    real_convolution_1d op(world, ops);

    real_function_1d opf = op(f);
    print("error in integral(op(g)) ", opf.trace()-1.0);

    real_function_1d exact = real_factory_1d(world).f(gconvg);
    print("norm2(g conv g - exact)", (opf-exact).norm2());

    real_function_1d q = real_factory_1d(world).f(h);
    print("error in integral(h) ", q.trace());
    print("error in norm2(h)", q.norm2() - sqrt(sqrt(2.0*constants::pi)));

    real_function_1d opq = op(q);
    exact = real_factory_1d(world).f(gconvh);
    print("norm2(g conv h - exact)", (opq-exact).norm2());

    ops[0].reset(new GaussianConvolution1D<double>(k, width*width*sqrt(8.0),
            width*width, 1, false));
    real_convolution_1d oph(world, ops);

    opq = oph(f);
    print("norm2(h conv g - exact)", (opq-exact).norm2());

    plot_line("opf.dat", 1001, lo, hi, q, opq, exact);

    world.gop.fence();
}


int main(int argc, char**argv) {
    initialize(argc,argv);
    World world(MPI::COMM_WORLD);

    try {
        startup(world,argc,argv);

        FunctionDefaults<1>::set_cubic_cell(-L,L);
        FunctionDefaults<1>::set_k(k);
        FunctionDefaults<1>::set_thresh(thresh);
        FunctionDefaults<1>::set_initial_level(5);
        test_gconv(world);

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
    finalize();

    return 0;
}

