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

#include "electronicstructureapp.h"
#include "solver.h"

using namespace madness;

int main(int argc, char** argv)
{
    initialize(argc, argv);

    World world(MPI::COMM_WORLD);

    try {
        // Load info for MADNESS numerical routines
        startup(world,argc,argv);
        std::cout.precision(6);
        FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap(world)));


//        // Process 0 reads input information and broadcasts
//        ElectronicStructureApp app(world, "input");
//
//        // Warm and fuzzy for the user
//        if (world.rank() == 0) {
//            print("\n\n");
//            print(" MADNESS Hartree-Fock and Density Functional Theory Program");
//            print(" ----------------------------------------------------------\n");
//            print("\n");
//            app.entity().print();
//            print("\n");
//            //app.params().print(world);
//        }
//
//        app.make_nuclear_potential();
//        app.initial_guess();
//        ElectronicStructureParams params = app.params();
//        Function<double,3> vnucrhon = app.vnucrhon();
//        vecfuncT orbs = app.orbitals();
//        std::vector<double> eigs;
//        std::vector< Function< std::complex<double>,3> > phis;
//        std::vector<double> tmpe = app.eigs();
//        print(tmpe);
//        int neps = eigs.size();
//        for (int i = 0; i < neps; i++)
//        {
//          phis.push_back(orbs[i]);
//          eigs.push_back(tmpe[i]);
//        }
//
//        Solver<double,3> dftcalc(world, vnucrhon, app.orbitals(), app.eigs(),
//                                 app.kpoints(), app.occs(), app.params(),
//                                 app.entity());

        Solver<double,3> dftcalc(world, "input");
        dftcalc.solve();
        world.gop.fence();

    } catch (const MPI::Exception& e) {
        //        print(e);
        error("caught an MPI exception");
    } catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    } catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
    } catch (const char* s) {
        print(s);
        error("caught a string exception");
    } catch (char* s) {
        print(s);
        error("caught a string exception");
    } catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    } catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    } catch (...) {
        error("caught unhandled exception");
    }

    finalize();

    return 0;
}
