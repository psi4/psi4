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


  $Id: testss.cc 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/


#include <mra/mra.h>
//#include <mra/loadbal.h>

using namespace madness;

double myfun(const double x[]) {
    double r2=0.0;
    for (int i=0; i < 3; ++i)
        r2 += x[i]*x[i];
    return r2;
}

const double PI = 3.1415926535897932384;
const double myg_expnt = 120.0;

double myg(const Vector<double,3>& r) {
    /* A square-normalized gaussian */
    double fac = pow(2.0*myg_expnt/PI,0.75);
    double x = r[0]-0.7;
    double y = r[1]-0.5;
    double z = r[2]-0.3;
    return fac*exp(-myg_expnt*(x*x + y*y + z*z));
};

void vector_myg(long npt, const double *x, const double *y,
                const double *z, double* restrict f) {
    const double fac = pow(2.0*myg_expnt/PI,0.75);
    for (int i=0; i<npt; ++i) {
        double xx = x[i] - 0.5;
        double yy = y[i] - 0.5;
        double zz = z[i] - 0.5;
        f[i] = fac*exp(-myg_expnt*(xx*xx + yy*yy + zz*zz));
    }
};

// double dmygdx(const double r[3]) {
//     /* Derivative of myg w.r.t. x */
//     return -2.0*myg_expnt*(r[0]-0.5)*myg(r);
// };

// double dmygdy(const double r[3]) {
//     /* Derivative of myg w.r.t. y */
//     return -2.0*myg_expnt*(r[1]-0.5)*myg(r);
// };

// double dmygdz(const double r[3]) {
//     /* Derivative of myg w.r.t. z */
//     return -2.0*myg_expnt*(r[2]-0.5)*myg(r);
// };


int main(int argc, char**argv) {
    MPI::Init(argc, argv);
    World world(MPI::COMM_WORLD);

    try {
        startup(world,argc,argv);
        Function<double,3> f = FunctionFactory<double,3>(world).f(myg).k(7).thresh(1e-5).nocompress();

        print("The tree after projection");
        f.print_tree();

        Vector<double,3> x = vec(0.5,0.5,0.5);
        print("the result is",f.eval(x).get()-myg(x));

        print("entering fence after eval");
        world.gop.fence();

        f.compress(false);
        print("entering fence after compress");
        world.gop.fence();
        print("The tree after compress");
        f.print_tree();

        f.reconstruct(false);
        print("entering fence after reconstruct");
        world.gop.fence();
        print("The tree after reconstruct");
        f.print_tree();

    }
    catch (const MPI::Exception& e) {
        print(e);
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
        error("caught a string exception");
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

    print("entering final fence");
    world.gop.fence();
    print("done with final fence");
    MPI::Finalize();

    return 0;
}
