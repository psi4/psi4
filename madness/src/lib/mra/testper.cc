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

/// \file testper.cc
/// \brief test the periodic convolution operator

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/operator.h>
#include <constants.h>

using namespace madness;

typedef Vector<double,1> coordT;

double constant(const coordT& x) {
    return 1.0;
}

double q(double x, double a) {
    return 0.5*(erf(sqrt(a)*x)-erf(-sqrt(a)+sqrt(a)*x));
}

void test_per(World& world) {
    int k = 10;
    double thresh = 1e-10;

    FunctionDefaults<1>::set_cubic_cell(0.0,1.0);
    FunctionDefaults<1>::set_k(k);
    FunctionDefaults<1>::set_thresh(thresh);
    FunctionDefaults<1>::set_truncate_mode(0);
    FunctionDefaults<1>::set_truncate_on_project(false);

    Tensor<double> coeff(1L), expnt(1L);
    expnt[0] = 10000.0;
    coeff[0] = sqrt(expnt[0]/constants::pi);
    print(coeff,expnt);
    SeparatedConvolution<double,1> op(world, k, coeff, expnt);

    Function<double,1> f = FunctionFactory<double,1>(world).f(constant).initial_level(3).norefine();

    Function<double,1> opf = apply(op,f);

    f.reconstruct();

    for (int i=0; i<101; ++i) {
        double x = i*0.01;
        printf("%.2f %.8f %.8f %10.2e\n", x, f(x), opf(x), opf(x)-q(x,expnt[0]));
    }
}


int main(int argc, char**argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);

    try {
        startup(world,argc,argv);

        test_per(world);

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
    catch (const char* s) {
        print(s);
        error("caught a c-string exception");
    }
    catch (char* s) {
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

