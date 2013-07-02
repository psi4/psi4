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

/** \file interior_dirichlet.cc
    \brief This file demonstrates solving a problem with interior (embedded)
    Dirichlet conditions.

    The signed_distance_functions shapes (mra/sdf_domainmask.h and
    mra/sdf_shape_3D.h) are used to specify a sphere of radius 1 centered
    at the origin.  A Dirichlet condition (dir_cond) is imposed on this sphere.

    After constructing the mask and the imposed boundary condition, the
    following routine is used to solve the equation (see the Li et al. paper).

    Suppose \f$\varphi\f$ is the mask function (1 on the inside, 0 on the
    outside, blurry on the boundary), \f$u\f$ is the desired function, \f$d\f$
    is the imposed Dirichlet condition on the boundary, \f$f\f$ is the
    inhomogeneity, \f$\mathcal{L}\f$ is the differential operator, and \f$G\f$
    is its free-space Green's function.

    The DE is \f$ \mathcal{L} u = f\f$ in the domain, and
        \f[ \mathcal{L}u - \varepsilon^{-2} b(\varphi) (u - d) = \varphi f, \f]
    where \f$b(\varphi) = 36 \varepsilon^{-1} \varphi^2 (1 - \varphi)^2\f$ and
    \f$\varepsilon\f$ is the thickness of the surface layer.

    Applying the Green's function:
        \f[ u - \varepsilon^{-2} G*( b(\varphi) u) == G*(\varphi f) -
            \varepsilon^{-2} G*( b(\varphi) d). \f]
    Thus, solving this problem involves a linear solve, as provided by GMRES.

    In this example, \f$\mathcal{L} = \nabla^2\f$, \f$ f(\vec{x}) = 0 \f$,
    and \f$ d(\theta, \phi) = Y_1^0(\theta, \phi) \f$ on the boundary
    (the unit sphere).

    The analytical solution \b inside the sphere is \f$u(\vec{x}) = z\f$.
*/

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <linalg/gmres.h>
#include <mra/sdf_shape_3D.h>

using namespace madness;

/** \brief The Dirichlet condition on the sphere.

    @param pt The point at which to evaluate; \f$|pt|=1\f$.
    @return The Dirichlet condition. */
static double dir_cond(const coord_3d &pt) {
   // Y_1^0
   const double r = sqrt(pt[0]*pt[0] + pt[1]*pt[1] + pt[2]*pt[2]);

   return pt[2] / r;
}

/** \brief The exact solution, for comparison (only valid inside the sphere).

    \param pt The point at which to evaluate.
    \return The exact solution. */
static double exact_sol(const coord_3d & pt) {
   const double r = sqrt(pt[0]*pt[0] + pt[1]*pt[1] + pt[2]*pt[2]);

   if(r < 1.0e-3)
      return 0.0;
   else
      return pt[2];
}

/// \brief The operator needed for solving for \f$u\f$ with GMRES
class DirichletCondIntOp : public Operator<real_function_3d> {
   protected:
      /// \brief The Green's function
      const SeparatedConvolution<double,3> &G;
      /// \brief The surface function, \f$b\f$
      const real_function_3d &b;

      /** \brief Applies the operator to \c invec

          \param[in] invec The input vector
          \param[out] outvec The action of the Green's function on \c invec */
      void action(const real_function_3d &invec, real_function_3d &outvec)
         const {

         outvec = invec - G(b*invec);
         outvec.truncate();
      }

   public:
      DirichletCondIntOp(const SeparatedConvolution<double, 3> &gin,
         const real_function_3d &bin)
         : G(gin), b(bin) {}
};


int main(int argc, char **argv) {
    double eps, inveps;

    initialize(argc,argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);

    eps = 0.1;
    inveps = 1.0 / eps;

    // Function defaults
    int k = 6;
    double thresh = 1.0e-4;
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_cubic_cell(-2.0, 2.0);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_max_refine_level(6);

    // create the Dirichlet boundary condition, expanded throughout the domain
    real_function_3d d = real_factory_3d(world).f(dir_cond);
    d.truncate();

    // create the domain mask, phi, and the surface function, b
    coord_3d pt(0.0); // Origin
    std::shared_ptr< SignedDFInterface<3> > sphere(new SDFSphere(1.0, pt));

    // use LLRV domain masking
    std::shared_ptr< DomainMaskInterface > llrvmask(new LLRVDomainMask(eps));

    // make the functor
    std::shared_ptr< DomainMaskSDFFunctor<3> > functor(new DomainMaskSDFFunctor<3>(llrvmask, sphere));

    real_function_3d phi = real_factory_3d(world).functor(functor);

    // create the surface function
    functor->setMaskFunction(DomainMaskSDFFunctor<3>::SURFACE);
    real_function_3d b = real_factory_3d(world).functor(functor);

    // check the surface area of the sphere
    double surfarea = b.trace();
    if(world.rank() == 0)
        printf("Error in surface area: %.6e\n",
            fabs(surfarea - 4.0*constants::pi));

    // scale the surface by \f$-\varepsilon^{-2}\f$
    // The two powers of \f$\varepsilon\f$ are from the auxiliary DE in
    // LLRV: the surface function \c b always appears with this factor.
    // The -1 is from the fact that MADNESS BSHOperator gives -G.
    b.scale(-inveps * inveps);
    b.truncate();

    // setup the Green's function
    // NOTE that CoulombOperator essentially makes the BSH w/ k == 0.0,
    // and then rescales by 4 pi.  This is more consistent.
    real_convolution_3d G = BSHOperator<3>(world, 0.0, eps*0.1, thresh);

    // compute the inhomogeneous portion
    real_function_3d usol = b*d; // should be -b*d, but b accounts for -G
    real_function_3d uinhomog = G(usol).truncate();
    uinhomog.scale(-1.0); // add the -1 from the Green's function
    uinhomog.truncate();
    world.gop.fence();
    usol.clear();

    // solve the linear system
    // make an initial guess -- make its norm, after one application of
    //                          dcio, close to uinhomog
    usol = copy(uinhomog);
    usol.scale(eps*eps);
    DirichletCondIntOp dcio(G, b);
    FunctionSpace<double, 3> space(world);
    int maxiters = 5;
    double resid_thresh = 1.0e-3;
    double update_thresh = 1.0e-3;
    GMRES(space, dcio, uinhomog, usol, maxiters, resid_thresh, update_thresh,
        true);

    real_function_3d exact = real_factory_3d(world).f(exact_sol);
    double error = ((usol - exact)*phi).norm2();
    if(world.rank() == 0)
        printf("   u error: %.10e\n", error);

    // set up file output
    char filename[100];
    sprintf(filename, "interior.vts");
    Vector<double, 3> plotlo, plothi;
    Vector<long, 3> npts;
    for(int i = 0; i < 3; ++i) {
        plotlo[i] = -1.1;
        plothi[i] = 1.1;
        npts[i] = 71;
    }
    plotvtk_begin(world, filename, plotlo, plothi, npts);
    plotvtk_data(usol, "usol", world, filename, plotlo, plothi, npts);
    plotvtk_data(exact, "exact", world, filename, plotlo, plothi, npts);
    plotvtk_end<3>(world, filename);

    finalize();

    return 0;
}
