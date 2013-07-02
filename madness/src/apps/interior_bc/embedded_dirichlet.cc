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

/** \file embedded_dirichlet.cc
    \brief Provides 2-D & 3-D test problems for examining the convergence of
           embedded boundary conditions.

    The auxiliary PDE being solved is
    \f[ \nabla^2 u - p(\varepsilon) S (u-g) = \varphi f, \f]
    where
       - \f$u\f$ is the solution function
       - \f$\varepsilon\f$ is the thickness of the boundary layer
       - \f$p(\varepsilon)\f$ is the penalty prefactor, \f$2/\varepsilon\f$
         seems to work well.
       - \f$S\f$ is the surface function
       - \f$g\f$ is the Dirichlet condition to be enforced on the surface
       - \f$\varphi\f$ is the domain mask (1 inside, 0 outside, blurry on the
         border)
       - \f$f\f$ is the inhomogeneity.

    The three test problems are
       -# A sphere of radius \f$R\f$ with \f$g = Y_0^0\f$, homogeneous
          (CONSTANT)
       -# A sphere of radius \f$R\f$ with \f$g = y_1^0\f$, homogeneous
          (COSTHETA)
       -# An ellipsoid of radii \f$(a=0.5, b=1.0, c=1.5)\f$ with \f$g = 2\f$,
          inhomogeneous: \f$f = 4 (a^{-2} + b^{-2} + c^{-2})\f$.

    This program allows testing of various parameters,
       -# The surface thickness
       -# The penalty prefactor
       -# The type of domain masking (LLRV or Gaussian)
       -# The curvature / shape of the domain
       .
    for their effect on convergence of the solution. */

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <linalg/gmres.h>
#include <muParser/muParser.h>
#include "test_problems.h"

using namespace madness;

int main(int argc, char **argv) {
    double eps, penalty_prefact;
    int k, prob, dim;
    double thresh, radius;
    Mask mask;

    initialize(argc,argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);

    // the structures for the problem, unfortunately different DIMs need
    // different structures...
    std::shared_ptr<EmbeddedDirichlet<2> > functor2;
    std::shared_ptr<EmbeddedDirichlet<3> > functor3;

    if (world.rank() == 0) {
        if(argc < 6) {
            std::cerr << "Usage error: ./app_name k thresh prob eps penalty" \
                " mask [radius, prob = CONSTANT or COSTHETA]" << std::endl;
            std::cerr << "    Where prob = 1 for constant sphere,\n" \
                         "                 2 for cosine theta sphere,\n" \
                         "                 3 for ellipsoid,\n" \
                         "                 4 for unit circle\n" \
                         "                 5 for Y20 sphere\n" \
                         "                 6 for inhomogeneous constant sphere\n" << std::endl;
            std::cerr << "    Where mask = 1 for LLRV, 2 for Gaussian\n"
                << std::endl;
            std::cerr << "    Where penalty is the penalty_prefact, " \
                "specified as a function\n    of eps, i.e. 2/eps" << std::endl;
            error("bad number of arguments");
        }

        // read in and validate the command-line arguments
        k = atoi(argv[1]);
        if(k < 4) error("cheapskate");

        thresh = atof(argv[2]);
        if(thresh > 1.0e-4) error("use some real thresholds...");

        prob = atoi(argv[3]);
        if(prob < 1 || prob > 6) error("bad problem number");

        eps = atof(argv[4]);
        if(eps <= 0.0) error("eps must be positive, and hopefully small");

        mu::Parser parser;
        try {
            parser.DefineVar("eps", &eps);
            parser.SetExpr(std::string(argv[5]));
            penalty_prefact = parser.Eval();
        }
        catch(mu::Parser::exception_type &e) {
            error(e.GetMsg().c_str());
        }
        if(penalty_prefact <= 0.0) error("penalty_prefact must be positive");

        switch(atoi(argv[6])) {
        case 1:
            mask = LLRV;
            break;
        case 2:
            mask = Gaussian;
            break;
        default:
            error("unknown domain mask type, should be 1 or 2");
            break;
        }

        if(prob == 1 || prob == 2 || prob == 5 || prob == 6) {
            if(argc > 7) {
                radius = atof(argv[7]);
                if(radius <= 0.0) error("radius must be positive");
            }
            else
                radius = 1.0;
        }
    }
    world.gop.broadcast(prob);
    world.gop.broadcast(eps);
    world.gop.broadcast(thresh);
    world.gop.broadcast(k);
    world.gop.broadcast(mask);
    world.gop.broadcast(penalty_prefact);
    world.gop.broadcast(radius);

    // do the final problem setup
    switch(prob) {
    case 1:
        functor3.reset(new ConstantSphere(k,
                       thresh, eps, std::string(argv[5]), penalty_prefact,
                       radius, mask));
        dim = 3;
        break;
    case 2:
        functor3.reset(new CosineSphere(k, thresh,
                       eps, std::string(argv[5]), penalty_prefact, radius,
                       mask));
        dim = 3;
        break;
    case 3:
        functor3.reset(new Ellipsoid(k, thresh,
                       eps, std::string(argv[5]), penalty_prefact, mask));
        dim = 3;
        break;
    case 4:
        functor2.reset(new LLRVCircle(k, thresh,
                       eps, std::string(argv[5]), penalty_prefact, mask));
        dim = 2;
        break;
    case 5:
        functor3.reset(new Y20Sphere(k, thresh,
                       eps, std::string(argv[5]), penalty_prefact, radius,
                       mask));
        dim = 3;
        break;
    case 6:
        functor3.reset(new InhomoConstantSphere(k,
                       thresh, eps, std::string(argv[5]), penalty_prefact,
                       radius, mask));
        dim = 3;
        break;
    default:
        dim = 0;
        error("shouldn't be here");
        break;
    }

    if(world.rank() == 0) {
        // print out the arguments
        switch(dim) {
        case 2:
            functor2->printout();
            break;
        case 3:
            functor3->printout();
            break;
        }
    }

    // project the surface function
    if(world.rank() == 0) {
        printf("Projecting the surface function (to low order)\n");
        fflush(stdout);
    }
    real_function_2d surf2;
    real_function_3d surf3;

    switch(dim) {
    case 2:
        functor2->fop = SURFACE;
        surf2 = real_factory_2d(world).k(6).thresh(1.0e-4).functor(functor2);
        break;
    case 3:
        functor3->fop = SURFACE;
        surf3 = real_factory_3d(world).k(6).thresh(1.0e-4).functor(functor3);
        break;
    }

    if(world.rank() == 0) {
        printf("Performing load balancing\n");
        fflush(stdout);
    }
    switch(dim) {
    case 2:
        break;
    case 3:
        functor3->load_balance(world, surf3);
        break;
    }

    // reproject the surface function to the requested threshold / k
    if(k > 6 || thresh < 1.0e-4) {
        if(world.rank() == 0) {
            printf("Reprojecting the surface function to requested order\n");
            fflush(stdout);
        }
        switch(dim) {
        case 2:
            surf2.clear();
            surf2 = real_factory_2d(world).functor(functor2);
            break;
        case 3:
            surf3.clear();
            surf3 = real_factory_3d(world).functor(functor3);
            break;
        }
    }

    // project the domain mask
    real_function_2d phi2;
    real_function_3d phi3;
    if(world.rank() == 0) {
        printf("Projecting the domain mask\n");
        fflush(stdout);
    }
    switch(dim) {
    case 2:
        functor2->fop = DOMAIN_MASK;
        phi2 = real_factory_2d(world).functor(functor2);
        break;
    case 3:
        functor3->fop = DOMAIN_MASK;
        phi3 = real_factory_3d(world).functor(functor3);
        break;
    }

    // print out the errors in perimeter (2-D) and surface area (3-D)
    // these are checks of the diffuse domain approximation
    double surf_integral, anals;
    double vol_integral, analv;
    surf_integral = 0.0;
    anals = 0.0;
    vol_integral = 0.0;
    analv = 0.0;

    switch(dim) {
    case 2:
        surf_integral = surf2.trace();
        anals = functor2->SurfaceIntegral();
        vol_integral = phi2.trace();
        analv = functor2->VolumeIntegral();
        break;
    case 3:
        surf_integral = surf3.trace();
        anals = functor3->SurfaceIntegral();
        vol_integral = phi3.trace();
        analv = functor3->VolumeIntegral();
        break;
    }
    if(world.rank() == 0) {
        printf("Error in Surface Integral: %.6e\n",
            fabs(surf_integral/penalty_prefact - anals));
        printf("Error in Volume Integral: %.6e\n",
            fabs(vol_integral - analv));
    }

    // green's function
    // note that this is really -G...
    real_convolution_2d G2 = BSHOperator<2>(world, 0.0, eps*0.1, thresh);
    real_convolution_3d G3 = BSHOperator<3>(world, 0.0, eps*0.1, thresh);
    //G3.broaden();

    // project the r.h.s. function (phi*f - penalty*S*g)
    // and then convolute with G
    real_function_2d usol2, rhs2;
    real_function_3d usol3, rhs3;
    if(world.rank() == 0) {
        printf("Projecting the r.h.s. function\n");
        fflush(stdout);
    }

    switch(dim) {
    case 2:
        functor2->fop = DIRICHLET_RHS;
        usol2 = real_factory_2d(world).functor(functor2);
        rhs2 = G2(usol2);
        rhs2.truncate();
        usol2.clear();
        break;
    case 3:
        functor3->fop = DIRICHLET_RHS;
        usol3 = real_factory_3d(world).functor(functor3);
        rhs3 = G3(usol3);
        rhs3.truncate();
        usol3.clear();
        break;
    }

    // load balance using the domain mask, the surface function, and the rhs
    if(dim == 3) {
        if(world.rank() == 0){ 
            printf("Load Balancing again...\n");
            fflush(stdout);
        }

        LoadBalanceDeux<3> lb(world);
        lb.add_tree(phi3, DirichletLBCost<3>(1.0, 1.0));
        lb.add_tree(surf3, DirichletLBCost<3>(1.0, 1.0));
        lb.add_tree(rhs3, DirichletLBCost<3>(1.0, 1.0));
        FunctionDefaults<3>::redistribute(world, lb.load_balance(2.0,
            false));
    }

    // make an initial guess:
    // uguess = rhs / penalty_prefact
    // the rescaling will make operator(uguess) close to rhs in magnitude for
    //     starting in GMRES
    DirichletCondIntOp<2> *dcio2 = NULL;
    DirichletCondIntOp<3> *dcio3 = NULL;
    switch(dim) {
    case 2:
        usol2 = copy(rhs2);
        usol2.scale(1.0 / penalty_prefact);
        usol2.compress();
        dcio2 = new DirichletCondIntOp<2>(G2, surf2);
        break;
    case 3:
        usol3 = copy(rhs3);
        usol3.scale(1.0 / penalty_prefact);
        usol3.compress();
        dcio3 = new DirichletCondIntOp<3>(G3, surf3);
        break;
    }

    // make the operators and prepare GMRES
    FunctionSpace<double, 2> space2(world);
    FunctionSpace<double, 3> space3(world);
    double resid_thresh = 1.0e-5;
    double update_thresh = 1.0e-5;
    int maxiter = 30;
    switch(dim) {
    case 2:
        GMRES(space2, *dcio2, rhs2, usol2, maxiter, resid_thresh, update_thresh,
              true);
        break;
    case 3:
        GMRES(space3, *dcio3, rhs3, usol3, maxiter, resid_thresh, update_thresh,
              true);
        break;
    }

    // compare to the exact solution
    real_function_2d uexact2, uerror2;
    real_function_3d uexact3, uerror3;
    double error, ratio, exactval;
    error = 0.0;
    ratio = 0.0;
    exactval = 0.0;

    std::vector<Vector<double, 2> > check_pts2;
    std::vector<Vector<double, 3> > check_pts3;

    switch(dim) {
    case 2:
        functor2->fop = EXACT;
        uexact2 = real_factory_2d(world).functor(functor2);
        uerror2 = (usol2 - uexact2)*phi2; // only use interior solution
        error = uerror2.norm2();
        break;
    case 3:
        functor3->fop = EXACT;
        uexact3 = real_factory_3d(world).functor(functor3);
        uerror3 = (usol3 - uexact3)*phi3; // only use interior solution
        error = uerror3.norm2();
        break;
    }

    if(world.rank() == 0) {
        printf("\nu interior error: %.10e\n", error);
        fflush(stdout);
    }

    // check the points prescribed by the problem
    switch(dim) {
    case 2:
        check_pts2 = functor2->check_pts();
        for(std::vector<Vector<double, 2> >::iterator iter =
            check_pts2.begin(); iter != check_pts2.end(); ++iter) {

            ratio = usol2(*iter);
            exactval = functor2->ExactSol(*iter);

            if(world.rank() == 0) {
                printf("u/uexact ratio at (%.2f, %.2f): %.10e\n",
                    (*iter)[0], (*iter)[1], ratio/exactval);
            }
        }
        break;
    case 3:
        check_pts3 = functor3->check_pts();
        for(std::vector<Vector<double, 3> >::iterator iter =
            check_pts3.begin(); iter != check_pts3.end(); ++iter) {

            ratio = usol3(*iter);
            exactval = functor3->ExactSol(*iter);

            if(world.rank() == 0) {
                printf("u/uexact ratio at (%.2f, %.2f, %.2f): %.10e\n",
                    (*iter)[0], (*iter)[1], (*iter)[2], ratio/exactval);
            }
        }
        break;
    }

    // uncomment these lines for various plots
    if(dim == 2) {
        // make a line plot along the positive z axis
        if(world.rank() == 0)
            printf("\n\n");
        double xmin = 0.0;
        double xmax = 2.0;
        int nx = 201;
        coord_2d pt2;
        pt2[1] = 0.0;
        double dx = (xmax - xmin) / (nx - 1);
        for(int i = 0; i < nx; ++i) {
            pt2[0] = xmin + i * dx;
            double uval = usol2(pt2);

            if(world.rank() == 0) {
                printf("%.4e %.4e\n", pt2[0], uval);
            }
        }
    }

    // make a line plot along the positive z axis
    if(dim == 3) {
        if(world.rank() == 0)
            printf("\n\n");
        /*double zmin = 0.0;
        double zmax = 2.0;
        int nz = 201;*/
        double zmin = 1.0 - 10.0*eps;
        double zmax = 1.0 + 10.0*eps;
        int nz = 51;
        coord_3d pt;
        pt[0] = pt[1] = 0.0;
        double dz = (zmax - zmin) / (nz - 1);
        for(int i = 0; i < nz; ++i) {
            pt[2] = zmin + i * dz;
            double uval = usol3(pt);

            if(world.rank() == 0) {
                printf("%.4e %.4e\n", pt[2], uval);
            }
        }
    }

    // print out the solution function
    /*char filename[100];
    sprintf(filename, "spheresol.vts");
    Vector<double, 3> plotlo, plothi;
    Vector<long, 3> npts;
    for(int i = 0; i < 3; ++i) {
        plotlo[i] = -2.0;
        plothi[i] = 2.0;
        npts[i] = 101;
    }
    plotvtk_begin(world, filename, plotlo, plothi, npts);
    plotvtk_data(usol, "usol", world, filename, plotlo, plothi, npts);
    plotvtk_end<3>(world, filename);*/

    finalize();

    return 0;
}
