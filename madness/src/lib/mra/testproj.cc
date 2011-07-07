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

/// \file testproj.cc
/// \brief test box size dependence of projection etc.

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <constants.h>
#include <vector>

using namespace madness;

template <typename T, std::size_t NDIM>
class Gaussian : public FunctionFunctorInterface<T,NDIM> {
public:
    typedef Vector<double,NDIM> coordT;
    const coordT center;
    const double exponent;
    const T coefficient;

    Gaussian(const coordT& center, double exponent, T coefficient)
            : center(center), exponent(exponent), coefficient(coefficient) {}

    T operator()(const coordT& x) const {
        double sum = 0.0;
        for (std::size_t i=0; i<NDIM; ++i) {
            double xx = center[i]-x[i];
            sum += xx*xx;
        };
        return coefficient*exp(-exponent*sum);
    }

    std::vector<coordT> special_points() const {
    	return std::vector<coordT>(1,center);
    }
};



template <typename T>
void test_proj(World& world) {
    typedef Vector<double,3> coordT;
    typedef std::shared_ptr< FunctionFunctorInterface<T,3> > functorT;

    FunctionDefaults<3>::set_k(6);
    FunctionDefaults<3>::set_thresh(1e-8);
    FunctionDefaults<3>::set_initial_level(5);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_autorefine(true);
    FunctionDefaults<3>::set_truncate_mode(0);
    FunctionDefaults<3>::set_truncate_on_project(false);

    const double expnt = 100.0;
    const double coeff = pow(expnt/constants::pi,1.5);
    double x = 1.0/11.0; // A non-dyadic point
    //double x = 0.0;
    const coordT origin(x);

    for (int i=7; i<=7; ++i) {
        double L = pow(2.0,double(i));
        FunctionDefaults<3>::set_cubic_cell(-L,L);
        print("I think the cell volume is", FunctionDefaults<3>::get_cell_volume());

        Function<T,3> f = FunctionFactory<T,3>(world).functor(functorT(new Gaussian<T,3>(origin, expnt, coeff)));
        f.truncate();
        f.reconstruct();
        print("L",L,f.trace(),f.norm2(),f.size()/6/6/6,f.max_depth());

        //f.print_tree();

    }
}


int main(int argc, char**argv) {
    initialize(argc,argv);
    World world(MPI::COMM_WORLD);

    try {
        startup(world,argc,argv);

        std::cout.precision(8);

        test_proj<double>(world);

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

