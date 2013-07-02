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

/// \file testbsh.cc
/// \brief test the bsh operator

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/operator.h>
#include <linalg/solvers.h>
#include <constants.h>

using namespace madness;

typedef std::shared_ptr< WorldDCPmapInterface< Key<3> > > pmapT;
typedef Vector<double,3> coordT;
typedef std::shared_ptr< FunctionFunctorInterface<double,3> > functorT;
typedef Function<double,3> functionT;
typedef std::vector<functionT> vecfuncT;
typedef std::pair<vecfuncT,vecfuncT> pairvecfuncT;
typedef std::vector<pairvecfuncT> subspaceT;
typedef Tensor<double> tensorT;
typedef FunctionFactory<double,3> factoryT;
typedef SeparatedConvolution<double,3> operatorT;
typedef std::shared_ptr<operatorT> poperatorT;

// The delta fn approx is (1/sqrt(2*pi*s^2))*exp(-x^2 / (2*s^2))

class Sphere : public FunctionFunctorInterface<double,3> {
public:
    const coordT center;
    const double radius;
    const double sigma;
    const double expnt;
    const double norm;

    Sphere(const coordT& center, double radius, double sigma)
        : center(center)
        , radius(radius)
        , sigma(sigma)
        , expnt(0.5/(sigma*sigma))
        , norm(1.0/(sigma*sqrt(2.0*constants::pi)))
    {}

    double operator()(const coordT& x) const {
        double xx = (x[0]-center[0]);
        double yy = (x[1]-center[1]);
        double zz = (x[2]-center[2]);
        double r = sqrt(xx*xx + yy*yy + zz*zz);
        double dist = r - radius; // Normal distance from boundary
        double arg = dist*dist*expnt; // Argument to exponential
        if (arg > 35.0)
            return 0.0;
        else
            return norm*exp(-arg);
    }
};


class DSphere : public FunctionFunctorInterface<double,3> {
public:
    const coordT center;
    const double radius;
    const double sigma;
    const double expnt;
    const double norm;

    DSphere(const coordT& center, double radius, double sigma)
        : center(center)
        , radius(radius)
        , sigma(sigma)
        , expnt(0.5/(sigma*sigma))
        , norm(1.0/(sigma*sqrt(2.0*constants::pi)))
    {}

    double operator()(const coordT& x) const {
        double xx = (x[0]-center[0]);
        double yy = (x[1]-center[1]);
        double zz = (x[2]-center[2]);
        double r = sqrt(xx*xx + yy*yy + zz*zz);
        double dist = r - radius; // Normal distance from boundary
        double arg = dist*dist*expnt; // Argument to exponential
        if (arg > 35.0)
            return 0.0;
        else
            return -2.0*dist*expnt*norm*exp(-arg);
    }
};


int main(int argc, char**argv) {
    initialize(argc,argv);
    World world(MPI::COMM_WORLD);
    try {
        startup(world,argc,argv);
        std::cout.precision(6);

        // MADNESS simulation parameters
        const int k = 6;
        const double thresh = 1e-5;
        const double L = 4.0;

        // Boundary width
        double sigma = 0.1;

        // Plotting info
        int npt = 3001;
        coordT lo(0.0), hi(0.0);
        lo[0] = -L;
        hi[0] =  L;

        FunctionDefaults<3>::set_cubic_cell(-L,L);
        FunctionDefaults<3>::set_k(k);
        FunctionDefaults<3>::set_thresh(thresh);
        FunctionDefaults<3>::set_initial_level(3);
        FunctionDefaults<3>::set_refine(true);
        FunctionDefaults<3>::set_autorefine(false);
        FunctionDefaults<3>::set_truncate_mode(1);
        FunctionDefaults<3>::set_truncate_on_project(true);

        // Inner and outer shells
        functionT spha = factoryT(world).functor(functorT(new Sphere(coordT(), 1.0, sigma)));
        functionT sphb = factoryT(world).functor(functorT(new Sphere(coordT(), 3.0, sigma)));

        //functionT dspha = factoryT(world).functor(functorT(new DSphere(coordT(), 1.0, sigma)));
        //functionT dsphb = factoryT(world).functor(functorT(new DSphere(coordT(), 3.0, sigma)));

        print("inner shell norm", spha.norm2());
        print("outer shell norm", sphb.norm2());

        functionT phi = 3.0*spha + sphb;      // Delta function at the boundary
        functionT u0 = -3.0*spha + sphb;      // The desired boundary conditions * phi
        //functionT du0= -1.0*dspha + dsphb;  // The desired boundary conditions * phi
        //dspha.clear(); dsphb.clear();
        spha.clear(); sphb.clear();

        phi.truncate(thresh*0.1);
        u0.truncate(thresh*0.1);
        //du0.truncate(thresh*0.1);

        plot_line("phi.dat", npt, lo, hi, phi);
        plot_line("u0.dat", npt, lo, hi, u0);

        operatorT op = CoulombOperator(world, sigma*0.1, thresh);
        // Starting values
        double mu = 0.25;
        functionT f = u0; //
        //du0.clear();
        //        functionT lambda = f;  // lambda is actually lambdaphi

        functionT u; // This will be current solution
        functionT c; // This will be the value of the constraint

        for (int muiter=0; muiter<5; ++muiter) {
            //            for (int lamiter=0; lamiter<20; ++lamiter) {
                vecfuncT rvec;
                vecfuncT fvec;
                std::vector<double> rnorms;
                tensorT Q;
                for (int m=0; m<20; ++m) {
                    f.truncate();
                    fvec.push_back(f);
                    u = op(fvec[m]);
                    u.scale(0.25/constants::pi);
                    u.truncate(); // Current soln

                    c = u*phi - u0; // Constraint violation
                    double cnorm = c.norm2();

                    //                    rvec.push_back(fvec[m] - lambda + (phi*c)*(1.0/mu)); // Residual AugLag
                    rvec.push_back(fvec[m] + phi*c*(1.0/mu)); // Residual QuadPen
                    rvec.back().truncate();
                    double rnorm = rvec[m].norm2();
                    rnorms.push_back(rnorm);

                    // Update subspace matrix
                    tensorT Qnew(m+1,m+1);
                    if (m>0) {
                        Qnew(Slice(0,m-1),Slice(0,m-1)) = Q;
                    }
                    Q = Qnew;
                    for (int i=0; i<=m; ++i) {
                        Q(i,m) = fvec[i].inner(rvec[m]);
                        if (m != i) Q(m,i) = fvec[m].inner(rvec[i]);
                    }
//                     print("Subspace matrix");
//                     print(Q);

                    tensorT coeff = KAIN(Q);

//                     print("Subspace solution");
//                     print(coeff);

                    f = coeff[0]*(fvec[0]-rvec[0]);
                    for (int i=1; i<=m; ++i) {
                        f.gaxpy(1.0,fvec[i]-rvec[i],coeff[i]);
                    }

                    // Step restriction
                    double snorm = (f - fvec[m]).norm2();
                    double fnorm = fvec[m].norm2();

                    if (snorm > fnorm*0.3) {
                        double damp = fnorm*0.3/snorm;
                        print("    damping", damp);
                        f.gaxpy(damp, fvec[m], 1.0-damp);
                    }
                    f.truncate();

                    print("    iteration", m, "mu", mu, "fnorm", fnorm, "snorm", snorm, "rnorm", rnorm, "cnorm", cnorm);

                    if (snorm < 0.1*cnorm) {
                        // Current solution will be in f
                        break;
                    }
                }

                char fname[256];
                int lamiter = 0;
                sprintf(fname, "c-%2.2d-%2.2d.dat", lamiter, muiter);
                plot_line(fname, npt, lo, hi, c);
                sprintf(fname, "u-%2.2d-%2.2d.dat", lamiter, muiter);
                plot_line(fname, npt, lo, hi, u);
                sprintf(fname, "f-%2.2d-%2.2d.dat", lamiter, muiter);
                plot_line(fname, npt, lo, hi, f);
//                 sprintf(fname, "l-%2.2d-%2.2d.dat", lamiter, muiter);
//                 plot_line(fname, npt, lo, hi, lambda);

//                 double damp = 1.0;
//                 functionT lambdanew = lambda - phi* c * (damp/mu);

//                 print("\nupdating lambda", lamiter, (lambdanew-lambda).norm2());
//                 lambda = lambdanew;
//                 lambda.truncate();
//             }
            mu *= 0.5;
//             lambda = f + phi * c * (1.0/mu);
            print("\nupdating mu", muiter, mu);
        }

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

