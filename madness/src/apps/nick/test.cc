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
/*********************
 * Trying to isolate a memory error
 *********************/
#include <string>
#include <fstream>
using std::ofstream;
using std::ofstream;
#include <stdlib.h>
#include <iomanip>
#include <time.h>
#include <math.h>
#include "wavef.h"
#define PRINT(str) if(world.rank()==0) std::cout << str
#define PRINTLINE(str) if(world.rank()==0) std::cout << str << std::endl

using namespace madness;

typedef std::complex<double> complexd;
typedef Vector<double,NDIM> vector3D;
typedef Function<complexd,NDIM> complex_functionT;
typedef Function<double,NDIM> functionT;
typedef FunctionFactory<complexd,NDIM> complex_factoryT;
typedef FunctionFactory<double,NDIM> factoryT;
typedef std::shared_ptr< WorldDCPmapInterface< Key<3> > > pmapT;
const char* wave_function_filename(int step);
bool wave_function_exists(World& world, int step);
void wave_function_store(World& world, int step, const complex_functionT& psi);
complex_functionT wave_function_load(World& world, int step);

struct PhiKAdaptor : public FunctionFunctorInterface<std::complex<double>,3> {
    PhiK& phik;
    PhiKAdaptor(PhiK& phik) : phik(phik) {}

    std::complex<double> operator()(const vector3D& x) {
        return phik(x);
    }
};

int main(int argc, char** argv) {
    initialize(argc,argv);
    World world(MPI::COMM_WORLD);
    startup(world, argc, argv);
    // Setup defaults for numerical functions
    int    k = 12;
    double L = M_PI;
    //    loadParameters(world, k, L, Z, cutoff);
    FunctionDefaults<NDIM>::set_k(k);               // Wavelet order
    FunctionDefaults<NDIM>::set_thresh(1e-6);       // Accuracy
    FunctionDefaults<NDIM>::set_cubic_cell(-L, L);
    FunctionDefaults<NDIM>::set_initial_level(3);
    FunctionDefaults<NDIM>::set_apply_randomize(false);
    FunctionDefaults<NDIM>::set_autorefine(false);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_truncate_mode(0);
    //FunctionDefaults<NDIM>::set_pmap(pmapT(new LevelPmap(world)));
    FunctionDefaults<NDIM>::set_truncate_on_project(true);
    try{
        const int n = 9;
        std::vector<int> a(n);
        for(int i=0; i<n; i++) {
            a[i] = 0;
        }
        for(int i=world.rank(); i<n; i+=world.size()) {
            a[i] = 2*i + 1;
        }
        world.gop.sum(&a[0], n);
        world.gop.fence();
        if( world.rank()==0 ){
            int total = 0;
            for(int i=0; i<n; i++) {
                total += a[i];
            }
            std::cout << total << std::endl;
        }
    } catch (const MPI::Exception& e) {
        //print(e);
        error("caught an MPI exception");
    } catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    } catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
    } catch (const char* s) {
        print(s);
        error("caught a c-string exception");
    } catch (char* s) {
        print(s);
        error("caught a c-string exception");
    } catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    } catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    } catch (...) {
        error("caught unhandled exception");
    }
    world.gop.fence();
    finalize();				//FLAG
    return 0;
}
