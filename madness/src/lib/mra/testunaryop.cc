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

/// \file testunaryop.cc
/// \brief test a unary op

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <constants.h>

using namespace madness;

template <typename T, std::size_t NDIM>
class Gaussian : public FunctionFunctorInterface<T,NDIM> {
public:
    typedef Vector<double,NDIM> coordT;
    const coordT center;
    const double exponent;
    const T coefficient;

    Gaussian(const coordT& center, double exponent, T coefficient)
            : center(center), exponent(exponent), coefficient(coefficient) {};

    T operator()(const coordT& x) const {
        double sum = 0.0;
        for (int i=0; i<NDIM; ++i) {
            double xx = center[i]-x[i];
            sum += xx*xx;
        };
        return coefficient*exp(-exponent*sum);
    };
};

template <typename T, std::size_t NDIM>
void squareit(const Key<NDIM>& key, Tensor<T>& t) {
    print("squareit", key,t);
    UNARY_OPTIMIZED_ITERATOR(T, t, *_p0 *= *_p0);
}

template <typename T, std::size_t NDIM>
void test_unaryop(World& world) {
    typedef Vector<double,NDIM> coordT;
    typedef std::shared_ptr< FunctionFunctorInterface<T,NDIM> > functorT;

    if (world.rank() == 0)
        print("Test unary operation (pointwise function-of-a-function), type =",
              archive::get_type_name<T>(),", ndim =",NDIM);

    FunctionDefaults<NDIM>::set_cubic_cell(-10,10);
    FunctionDefaults<NDIM>::set_k(7);
    FunctionDefaults<NDIM>::set_thresh(1e-5);
    FunctionDefaults<NDIM>::set_autorefine(false);

    const coordT origin(0.0);
    const double expnt = 1.0;
    const double coeff = 1.0;// pow(2.0/constants::pi,0.25*NDIM);

    Function<T,NDIM> f = FunctionFactory<T,NDIM>(world).functor(functorT(new Gaussian<T,NDIM>(origin, expnt, coeff)));
    double norm = f.norm2();
    print("norm of initial function", norm);

    Function<T,NDIM> g = copy(f);
    g.unaryop(&squareit<T,NDIM>);
    double gnorm = g.norm2();
    print("norm of the squared function", gnorm);

    double err = (f*f - g).norm2();
    print("norm of the error", err);

    err = g.err(Gaussian<T,NDIM>(origin, expnt*2.0, coeff*coeff));
    print("norm of the error", err);

    world.gop.fence();

}

int main(int argc, char**argv) {
    MPI::Init(argc, argv);
    World world(MPI::COMM_WORLD);

    try {
        startup(world,argc,argv);

        test_unaryop<double,1>(world);

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

