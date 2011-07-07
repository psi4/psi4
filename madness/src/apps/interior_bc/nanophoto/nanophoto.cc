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

/** \file nanophoto.cc
    \brief

*/

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <linalg/gmres.h>
#include "basisfunction.h"
#include "atom.h"
#include "density.h"

void vtk_output(World &world, const char *funcname,
    const real_function_3d &func);
void scaled_plotvtk_begin(World &world, const char *filename,
    const Vector<double, 3> &plotlo, const Vector<double, 3> &plothi,
    const Vector<long, 3> &npt, bool binary = false);
int mol_geom(std::vector<Atom*> &atoms);
void read_states(int nstate, int nbasis, Tensor<double> &coeffs);

using namespace madness;

int main(int argc, char **argv) {
    double eps;
    int j, k, n;
    double thresh, phi, d, penalty;
    unsigned int i;
    char funcname[15];

    initialize(argc,argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);

    if (world.rank() == 0) {
        if(argc < 6) {
            print("usage: ./nanophoto k thresh epsilon potential-difference " \
                "tip-surface\n");
            print("potential difference is in mV, and tip-surface distance " \
                "in nm\n");
            error("bad number of arguments");
        }

        // read in and validate the command-line arguments
        k = atoi(argv[1]);
        if(k < 4) error("cheapskate");

        thresh = atof(argv[2]);
        if(thresh > 1.0e-4) error("use some real thresholds...");

        eps = atof(argv[3]) / 0.052918; // convert to a.u.
        if(eps <= 0.0) error("eps must be positive, and hopefully small");

        phi = atof(argv[4]) * 3.6749324e-5; // convert to a.u.

        d = atof(argv[5]) / 0.052918; // convert to a.u.
    }
    world.gop.broadcast(phi);
    world.gop.broadcast(eps);
    world.gop.broadcast(thresh);
    world.gop.broadcast(k);
    world.gop.broadcast(d);

    // box size
    Tensor<double> cell(3, 2);
    cell(0,0) = cell(1,0) = -150.0 / 0.052918;
    cell(0,1) = cell(1,1) = 150.0 / 0.052918;
    cell(2,0) = 10.0 * (constants::pi - 8.0) / 0.052918;
    cell(2,1) = (10.0 * (constants::pi - 8.0) + 300.0) / 0.052918;
    FunctionDefaults<3>::set_cell(cell);

    // make the basis functions to get the density
    std::vector<Atom*> atoms(0);
    std::vector<BasisFunc> basis(0);
    int nstate = mol_geom(atoms);
    int nbasis = 0;

    // make the set of basis functions
    for(i = 0; i < atoms.size(); ++i) {
        j = atoms[i]->dimBasis();
        nbasis += j;

        for(n = 0; n < j; ++n)
            basis.push_back(atoms[i]->getBasisFunc(n));
    }

    // read in the coefficients of the occupied states
    Tensor<double> coeffs(nstate, nbasis);
    if(world.rank() == 0)
        read_states(nstate, nbasis, coeffs);
    world.gop.broadcast_serializable(coeffs, 0);

    penalty = 1.75 / eps;

    // the key data structure: sets up the problem details and projects
    // the initial functions
    std::shared_ptr<TipMolecule> tpm(new TipMolecule(eps, penalty, coeffs, atoms,
        basis, phi, d));

    if(world.rank() == 0) {
        // print out the arguments
        printf("Tip-Surface Distance: %.6e nm\nPotential Difference: %.6e " \
            "mV\nk: %d\nthresh: %.6e\neps: %.6e nm\n", d * 0.052918,
            phi * 27211.385, k, thresh, eps * 0.052918);

        // project the surface function
        printf("Load Balancing\n   --- Projecting domain mask to low order.\n");
        fflush(stdout);
    }

    real_function_3d dmask, dens;

    // low order defaults
    FunctionDefaults<3>::set_k(6);
    FunctionDefaults<3>::set_thresh(1.0e-4);

    // domain mask
    tpm->fop = DOMAIN_MASK;
    dmask = real_factory_3d(world).functor(tpm);

    // density
    if(world.rank() == 0) {
        printf("   --- Projecting density to low order.\n");
        fflush(stdout);
    }
    tpm->fop = DENSITY;
    dens = real_factory_3d(world).functor(tpm);

    // do the load balancing
    LoadBalanceDeux<3> lb(world);
    lb.add_tree(dmask, DirichletLBCost<3>(1.0, 1.0));
    lb.add_tree(dens, DirichletLBCost<3>(1.0, 1.0));
    // set this map as default and redistribute the existing functions
    FunctionDefaults<3>::redistribute(world, lb.load_balance(2.0, false));
    dmask.clear();
    dens.clear();

    // set the defaults to the real deal
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);

#if 0
    // get the domain mask -- this segment can be commented out
    if(world.rank() == 0) {
        printf("Projecting the domain mask\n");
        fflush(stdout);
    }
    sprintf(funcname, "domainmask");
    tpm->fop = DOMAIN_MASK;
    dmask = real_factory_3d(world).functor(tpm);
    vtk_output(world, funcname, dmask);
    dmask.clear();

    // get the electron density -- this segment can be commented out
    if(world.rank() == 0) {
        printf("Projecting the electron density\n");
        fflush(stdout);
    }
    tpm->fop = ELECTRON_DENSITY;
    dens = real_factory_3d(world).functor(tpm);
    double denstrace = dens.trace();
    if(world.rank() == 0) {
        // this should be close to -2
        printf("   trace = %.6e\n", denstrace);
        fflush(stdout);
    }
    sprintf(funcname, "elecdens");
    vtk_output(world, funcname, dens);
    dens.clear();

    // get the total density -- this segment can be commented out
    if(world.rank() == 0) {
        printf("Projecting the total density\n");
        fflush(stdout);
    }
    tpm->fop = DENSITY;
    dens = real_factory_3d(world).functor(tpm);
    denstrace = dens.trace();
    if(world.rank() == 0) {
        // this should be close to 0
        printf("   trace = %.6e\n", denstrace);
        fflush(stdout);
    }
    sprintf(funcname, "density");
    vtk_output(world, funcname, dens);
    dens.clear();
#endif

    // do we already have a solution, or do we need to calculate it?
    real_function_3d usol;
    char arname[50];
    sprintf(arname, "solution/solution");
    if(ParallelInputArchive::exists(world, arname)) {
        // read it in
        if(world.rank() == 0) {
            printf("Reading solution from archive\n");
            fflush(stdout);
        }

        ParallelInputArchive input(world, arname, 10);
        input & usol;
        input.close();
    }
    else {
        // project the surface and rhs functions
        if(world.rank() == 0) {
            printf("Projecting the surface function\n");
            fflush(stdout);
        }
        tpm->fop = SURFACE;
        real_function_3d surf = real_factory_3d(world).functor(tpm);
        sprintf(funcname, "surface");
        //vtk_output(world, funcname, surf);

        if(world.rank() == 0) {
            printf("Projecting the rhs function\n");
            fflush(stdout);
        }
        tpm->fop = DIRICHLET_RHS;
        real_function_3d rhs = real_factory_3d(world).functor(tpm);
        sprintf(funcname, "rhs");
        //vtk_output(world, funcname, rhs);

        // green's function
        // note that this is really -G...
        real_convolution_3d G = BSHOperator<3>(world, 0.0, eps*0.1, thresh);

        // project the r.h.s. function
        // and then convolute with G
        real_function_3d grhs;
        grhs = G(rhs);
        grhs.truncate();
        rhs.clear();

        // make an initial guess:
        // uguess = rhs / penalty_prefact
        // the rescaling will make operator(uguess) close to rhs in magnitude
        //     for starting in GMRES
        usol = copy(grhs);
        usol.scale(penalty);
        usol.compress();
        DirichletCondIntOp<3> dcio(G, surf);

        // make the operators and prepare GMRES
        FunctionSpace<double, 3> space(world);
        double resid_thresh = 1.0e-3;
        double update_thresh = 1.0e-3;
        int maxiter = 50;
        GMRES(space, dcio, grhs, usol, maxiter, resid_thresh, update_thresh,
            true);

        if(world.rank() == 0) {
            printf("Writing solution to archive\n");
            fflush(stdout);
        }

        // scale back to mV from a.u.
        usol.scale(27211.385);

        // write the solution to archive
        ParallelOutputArchive output(world, arname, 10);
        output & usol;
        output.close();
    }

    // print out the solution function to a vtk file
    sprintf(funcname, "solution-depth");
    vtk_output(world, funcname, usol);

    finalize();

    return 0;
}

/** \brief Plots a function in the total domain, and also close to the
           tip/surface junction.

    The filename is funcname.vts
*/
void vtk_output(World &world, const char *funcname,
    const real_function_3d &func) {

    char filename[80];

    // print out the function on the total domain
    sprintf(filename, "%s-coarse.vts", funcname);
    const Tensor<double> cell = FunctionDefaults<3>::get_cell();
    Vector<double, 3> plotlo, plothi;
    Vector<long, 3> npts;
    plotlo[0] = 0.0 / 0.052918; plothi[0] = 0.0 / 0.052918; npts[0] = 1;
    for(int i = 1; i < 3; ++i) {
        plotlo[i] = cell(i, 0);
        plothi[i] = cell(i, 1);
        npts[i] = 251;
    }
    scaled_plotvtk_begin(world, filename, plotlo, plothi, npts);
    //plotvtk_data(func, funcname, world, filename, plotlo, plothi, npts);
    plotvtk_data(func, funcname, world, filename, plotlo, plothi, npts, false, true);
    plotvtk_end<3>(world, filename);

    // print out the solution function near the area of interest
    sprintf(filename, "%s-medium.vts", funcname);
    plotlo[0] = 0.0 / 0.052918; plothi[0] = 0.0 / 0.052918; npts[0] = 1;
    plotlo[1] = -20.0 / 0.052918; plothi[1] = 20.0 / 0.052918; npts[1] = 251;
    plotlo[2] = -10.0 / 0.052918; plothi[2] = 30.0 / 0.052918; npts[2] = 251;
    scaled_plotvtk_begin(world, filename, plotlo, plothi, npts);
    //plotvtk_data(func, funcname, world, filename, plotlo, plothi, npts);
    plotvtk_data(func, funcname, world, filename, plotlo, plothi, npts, false, true);
    plotvtk_end<3>(world, filename);

    // print out the solution function near the area of interest
    sprintf(filename, "%s-fine.vts", funcname);
    plotlo[0] = 0.0 / 0.052918; plothi[0] = 0.0 / 0.052918; npts[0] = 1;
    plotlo[1] = -0.25 / 0.052918; plothi[1] = 0.25 / 0.052918; npts[1] = 251;
    plotlo[2] = 4.75 / 0.052918; plothi[2] = 5.25 / 0.052918; npts[2] = 251;
    scaled_plotvtk_begin(world, filename, plotlo, plothi, npts);
    //plotvtk_data(func, funcname, world, filename, plotlo, plothi, npts);
    plotvtk_data(func, funcname, world, filename, plotlo, plothi, npts, false, true);
    plotvtk_end<3>(world, filename);

    // print out the solution function near the area of interest
    sprintf(filename, "%s-hyperfine.vts", funcname);
    plotlo[0] = 0.0 / 0.052918; plothi[0] = 0.0 / 0.052918; npts[0] = 1;
    plotlo[1] = -0.002 / 0.052918; plothi[1] = 0.002 / 0.052918; npts[1] = 251;
    plotlo[2] = 5.0355 / 0.052918; plothi[2] = 5.0395 / 0.052918; npts[2] = 251;
    scaled_plotvtk_begin(world, filename, plotlo, plothi, npts);
    //plotvtk_data(func, funcname, world, filename, plotlo, plothi, npts);
    plotvtk_data(func, funcname, world, filename, plotlo, plothi, npts, false, true);
    plotvtk_end<3>(world, filename);
}

/** \brief Same as plotvtk_begin, but scales the coordinates back to nm */
void scaled_plotvtk_begin(World &world, const char *filename,
    const Vector<double, 3> &plotlo, const Vector<double, 3> &plothi,
    const Vector<long, 3> &npt, bool binary) {

    PROFILE_FUNC;

    Tensor<double> cell(3, 2);
    int i;
    for(i = 0; i < 3; ++i) {
        cell(i, 0) = plotlo[i];
        cell(i, 1) = plothi[i];
    }

    FILE *f=0;
    if(world.rank() == 0) {
        f = fopen(filename, "w");
        if(!f)
            MADNESS_EXCEPTION("plotvtk: failed to open the plot file", 0);

        fprintf(f, "<VTKFile type=\"StructuredGrid\" version=\"0.1\"" \
            " byte_order=\"LittleEndian\" compressor=\"" \
            "vtkZLibDataCompressor\">\n");
        fprintf(f, "  <StructuredGrid WholeExtent=\"");
        for(i = 0; i < 3; ++i)
            fprintf(f, "0 %ld ", npt[i]-1);
        for(; i < 3; ++i)
            fprintf(f, "0 0 ");
        fprintf(f, "\">\n");
        fprintf(f, "    <Piece Extent=\"");
        for(i = 0; i < 3; ++i)
            fprintf(f, "0 %ld ", npt[i]-1);
        for(; i < 3; ++i)
            fprintf(f, "0 0 ");
        fprintf(f, "\">\n");
        fprintf(f, "      <Points>\n");
        fprintf(f, "        <DataArray NumberOfComponents=\"3\" " \
            "type=\"Float32\" format=\"ascii\">\n");

        Vector<double, 3> space;
        for(i = 0; i < 3; ++i) {
            if(npt[i] == 1)
                space[i] = 0.0;
            else
                space[i] = (cell(i, 1) - cell(i, 0)) / (npt[i] - 1);
        }

        // go through the grid
        for(LowDimIndexIterator it(npt); it; ++it) {
            for(i = 0; i < 3; ++i)
                fprintf(f, "%f ", (plotlo[i] + it[i]*space[i]) * 0.052918);
            for(; i < 3; ++i)
                fprintf(f, "0.0 ");
            fprintf(f, "\n");
        }

        fprintf(f, "        </DataArray>\n");
        fprintf(f, "      </Points>\n");
        fprintf(f, "      <PointData>\n");
        fclose(f);
    }
    world.gop.fence();
}

/** \brief Make the molecular geometry.

    \return The number of (occupied) states to use in the density. */
int mol_geom(std::vector<Atom*> &atoms) {
    Vector<double, 3> center;

    center[0] = 0.0;
    center[1] = 0.0;
    center[2] = 0.7018442318 + 5.0 / 0.052918;
    atoms.push_back(new Hydrogen(center));

    center[0] = 0.0;
    center[1] = 0.0;
    center[2] = -0.7018442318 + 5.0 / 0.052918;
    atoms.push_back(new Hydrogen(center));

    return 1;
}

/** \brief Read in the occupied states from file. */
void read_states(int nstate, int nbasis, Tensor<double> &coeffs) {
    // open file dens.dat
    std::ifstream file;
    file.open("dens.dat");
    if(!file.good()) {
        error("Error opening dens.dat");
    }

    int curstate = 0;
    int i;
    char cbuffer[100];
    double dbuffer;

    // coefficients are output in groups of 5 in the file
    for(curstate = 0; curstate < nstate; curstate += 5) {
        // blow past three useless lines
        file.getline(cbuffer, 100);
        file.getline(cbuffer, 100);
        file.getline(cbuffer, 100);

        // start reading real lines
        for(i = 0; i < nbasis; ++i) {
            // some useless stuff at the front of each line
            file.read(cbuffer, 14);

            file >> dbuffer;
            coeffs(curstate, i) = dbuffer;

            if(curstate + 1 < nstate) {
                file >> dbuffer;
                coeffs(curstate + 1, i) = dbuffer;
            }

            if(curstate + 2 < nstate) {
                file >> dbuffer;
                coeffs(curstate + 2, i) = dbuffer;
            }

            if(curstate + 3 < nstate) {
                file >> dbuffer;
                coeffs(curstate + 3, i) = dbuffer;
            }

            if(curstate + 4 < nstate) {
                file >> dbuffer;
                coeffs(curstate + 4, i) = dbuffer;
            }

            // flush the rest of the line
            file.getline(cbuffer, 100);
        }

        // there's a blank line in between each block of five
        file.getline(cbuffer, 100);
    }

    // close the file
    file.close();
}
