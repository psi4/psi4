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


  $Id: test.cc 257 2007-06-25 19:09:38Z HartmanBaker $
*/

/// \file testsuite.cc
/// \brief The QA/test suite for Function

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <unistd.h>
#include <cstdio>
#include <constants.h>
#include <mra/qmprop.h>

#include <tensor/random.h>

const double PI = 3.1415926535897932384;

using namespace madness;

typedef Vector<double,1> coordT;
typedef std::shared_ptr< FunctionFunctorInterface<double_complex,1> > functorT;
typedef Function<double_complex,1> functionT;
typedef FunctionFactory<double_complex,1> factoryT;
typedef SeparatedConvolution<double_complex,1> operatorT;

double_complex psi0(const coordT& r) {
    return double_complex(exp(-r[0]*r[0]*0.5),0.0);
}

double_complex V(const coordT& r) {
    return double_complex(0.5*r[0]*r[0],0.0);
}

double_complex dVsq(const coordT& r) {
    return double_complex(r[0]*r[0],0.0);
}

double_complex energy(const functionT& v, const functionT& psi) {
    functionT du = diff(psi,0);
    double_complex ke = 0.5*du.inner(du);
    double_complex pe = psi.inner(v*psi);

    print("  ke", real(ke), "pe", real(pe), "total", real(ke+pe), "norm", psi.norm2());

    return ke+pe;
}

functionT chin_chen(const functionT& expV,
                    const functionT& expVtilde,
                    const operatorT& G,
                    const functionT& psi0) {

    // psi(t) = exp(-i*V*t/6) exp(-i*T*t/2) exp(-i*2*Vtilde*t/3) exp(-i*T*t/2) exp(-i*V*t/6)

    functionT psi1;

    psi1 = expV*psi0;
    psi1.truncate();
    psi1 = apply(G,psi1);
    psi1.truncate();
    psi1 = expVtilde*psi1;
    psi1.truncate();
    psi1 = apply(G,psi1);
    psi1.truncate();
    psi1 = expV*psi1;
    psi1.truncate();

    return psi1;
}

functionT trotter(const functionT& expV,
                  const operatorT& G,
                  const functionT& psi0) {
    //    psi(t) = exp(-i*T*t/2) exp(-i*V*t) exp(-i*T*t/2) psi(0)

    functionT psi1 = apply(G,psi0);
    psi1.truncate();

    psi1 = expV*psi1;

    psi1 = apply(G,psi1);
    psi1.truncate();

    return psi1;
}

template<typename T, std::size_t NDIM>
struct unaryexp {
    void operator()(const Key<NDIM>& key, Tensor<T>& t) const {
        UNARY_OPTIMIZED_ITERATOR(T, t, *_p0 = exp(*_p0););
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};


/// Returns exp(-I*t*V) with truncation
functionT make_exp(double t, const functionT& v) {
    v.reconstruct();
    functionT expV = double_complex(0.0,-t)*v;
    expV.unaryop(unaryexp<double_complex,1>());
    expV.truncate();
    return expV;
}


void test_trotter(World& world) {

    const double L = 20.0;
    const int k = 16;
    const double thresh = 1e-12;
    FunctionDefaults<1>::set_cubic_cell(-L,L);
    FunctionDefaults<1>::set_k(k);
    FunctionDefaults<1>::set_thresh(thresh);

    cout.precision(8);

    functionT psi = factoryT(world).f(psi0);
    psi.scale(1.0/psi.norm2());
    functionT psi0 = copy(psi);

    functionT v = factoryT(world).f(V);

    double_complex e = energy(v,psi);

    double c = 10.0*sqrt(0.5) * 1.86;
    double tcrit = 2*PI/(c*c);
    double tstep = tcrit * 0.5;

    print("The time step is", tstep, "\n");

    operatorT G = qm_free_particle_propagator<1>(world, k, c, tstep*0.5, 2*L);
    functionT expV = make_exp(tstep,v);

    for (int step=0; step<100; ++step) {
        double time = step * tstep;
        double_complex phase = psi0.inner(psi);
        double radius = abs(phase);
        double theta = arg(phase);
        double theta_exact = -time*0.5;
        while (theta_exact > PI) theta_exact -= 2.0*PI;

        print("step", step, "time", time, "radius", radius, "arg", theta, "exact", theta_exact, "phase err", theta_exact-theta);
        //energy(v, psi);

        psi = trotter(expV, G, psi);
    }

}

void test_chin_chen(World& world) {

    const double L = 20.0;
    const int k = 16;
    const double thresh = 1e-12;
    double c = 10.0*sqrt(0.5) * 1.86;
    double tcrit = 2*PI/(c*c);
    double tstep = tcrit;
    print("The time step is", tstep, "\n");

    FunctionDefaults<1>::set_cubic_cell(-L,L);
    FunctionDefaults<1>::set_k(k);
    FunctionDefaults<1>::set_thresh(thresh);

    cout.precision(8);

    functionT psi = factoryT(world).f(psi0);
    psi.scale(1.0/psi.norm2());
    functionT psi0 = copy(psi);

    functionT v = factoryT(world).f(V);
    functionT dvsq = factoryT(world).f(dVsq);

    functionT vtilde = v - dvsq*(tstep*tstep/48.0);

    operatorT G = qm_free_particle_propagator<1>(world, k, c, tstep*0.5, 2*L);
    //G.doleaves = true;
    functionT expV = make_exp(tstep/6.0, v);
    functionT expVtilde = make_exp(2.0*tstep/3.0, vtilde);

    for (int step=0; step<10000; ++step) {
        double time = step * tstep;
        double_complex phase = psi0.inner(psi);
        double radius = abs(phase);
        double theta = arg(phase);
        double theta_exact = -time*0.5;
        while (theta_exact > PI) theta_exact -= 2.0*PI;
        while (theta_exact < -PI) theta_exact += 2.0*PI;

        print("step", step, "time", time, "radius", radius, "arg", theta, "exact", theta_exact, "phase err", theta_exact-theta);
        //energy(v, psi);

        psi = chin_chen(expV, expVtilde, G, psi);
    }

}

int main(int argc, char**argv) {
    MPI::Init(argc, argv);
    World world(MPI::COMM_WORLD);
    try {
        startup(world,argc,argv);

        //test_trotter(world);
        test_chin_chen(world);

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

    MPI::Finalize();

    return 0;
}

