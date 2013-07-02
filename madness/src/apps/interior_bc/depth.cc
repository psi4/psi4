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

/** \file depth.cc
    \brief Provides a simple driver for examining the diffuse domain
           approximation.

    Specifically, this file uses 1-D functions (the domains is left and right)
    to look at convergence properties and computational resource requirements
    (particularly inside MADNESS) for varying surface thicknesses and locations
    (dyadic vs. non-dyadic points). */

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/sdf_domainmask.h>

using namespace madness;

// sdf interface for a point in 1-D
class Division : public SignedDFInterface<1> {
    private:
    // the dividing point
    double pt;

    public:
    Division(double pt) : pt(pt) {}

    double sdf(const coord_1d &x) const {
        return x[0] - pt;
    }

    coord_1d grad_sdf(const coord_1d &x) const {
        coord_1d ret;
        ret[0] = 1.0;
        return ret;
    }
};

int main(int argc, char **argv) {
    double eps;
    int k, maxrefine;
    double thresh;
    coord_1d x;

    initialize(argc,argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);

    if (world.rank() == 0) {
        if(argc < 6) {
            std::cerr << "Usage error: ./app_name k thresh max_refine eps" \
                " pt" << std::endl;
            error("bad number of arguments");
        }

        eps = atof(argv[4]);
        if(eps <= 0.0) error("eps must be positive, and hopefully small");
        maxrefine = atoi(argv[3]);
        if(maxrefine < 6) error("cheapskate");
        thresh = atof(argv[2]);
        if(thresh > 1.0e-4) error("use some real thresholds...");
        k = atoi(argv[1]);
        if(k < 4) error("cheapskate");
        x[0] = atof(argv[5]);
        if(fabs(x[0]) >= 1.0) error("out of bounds");
    }
    world.gop.broadcast(eps);
    world.gop.broadcast(thresh);
    world.gop.broadcast(maxrefine);
    world.gop.broadcast(k);
    world.gop.broadcast(x);

    // Function defaults
    FunctionDefaults<1>::set_k(k);
    FunctionDefaults<1>::set_cubic_cell(-1.0, 1.0);
    FunctionDefaults<1>::set_thresh(thresh);
    FunctionDefaults<1>::set_max_refine_level(maxrefine);

    std::shared_ptr<SignedDFInterface<1> > pt(new Division(x[0]));

    std::shared_ptr<DomainMaskInterface> llrv(new LLRVDomainMask(eps));

    std::shared_ptr<DomainMaskSDFFunctor<1> > functor(new DomainMaskSDFFunctor<1>(llrv, pt));

    real_function_1d f = real_factory_1d(world).functor(functor).
        truncate_on_project();
    functor->setMaskFunction(DomainMaskSDFFunctor<1>::SURFACE);
    real_function_1d fs = real_factory_1d(world).functor(functor).
        truncate_on_project();

    // look at the depth of the MADNESS tree as a function of position
    /*for(int i = 0; i < 101; ++i) {
        x[0] = (0.52*i) / 100 - 0.26;

        print(x[0], f.depthpt(x), fs.depthpt(x));
    }*/

    // maximum depths of the MADNESS trees for masking and surface functions
    print(f.max_depth(), fs.max_depth());

    finalize();

    return 0;
}
