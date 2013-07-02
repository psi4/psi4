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

/*!
  \file examples/dataloadbal.cc
  \brief Illustrates how to use static data/load balancing of functions
  \defgroup loadbaleg Data and load balancing 
  \ingroup examples

  The source is <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local/trunk/src/apps/examples/dataloadbal.cc>here</a>.

  This is one of the more computationally demanding examples - either run it on 
  50 or more nodes of jaguar or reduce the number of functions employed
  (value of \c NFUNC in the source).

  \par Points of interest
  - Using functors to incorporate state into functions to be compressed
  - Using special points in a functor to provide hints to adaptive refinement algorithm
  - Operating on multiple functions in a vector
  - Heuristic analysis of function trees to generate a new mapping of data to processor
  - Installation of the new process map with automatic redistribution of data
  - Timing to estimate benefit of load balancing

  \par Background
  
  Poor distribution of work (load imbalance) is the largest reason for
  inefficient parallel execution within MADNESS.  Poor data distribution
  (data imbalance) contributes to load imbalance and also leads to
  out-of-memory problems due to one or more processes having too much
  data.  Thus, we are interested in uniform distribution of both
  work and data.
  
  Many operations in MADNESS are entirely data driven (i.e.,
  computation occurs in the processor that "owns" the data) since
  there is insufficient work to justify moving data between processes
  (e.g., computing the inner product between functions).  However, a
  few expensive operations can have work shipped to other processors.

  There are presently three load balancing mechanisms within MADNESS
  - static and driven by the distribution of data,
  - dynamic via random assignment of work, and
  - dynamic via work stealing (currently only in prototype).
  
  Until the work stealing becomes production quality we must exploit
  the first two forms.  The random work assignment is controlled by
  options in the FunctionDefaults class.
  - FunctionDefaults::set_apply_randomize(bool) controls the use of
  randomization in applying integral (convolution) operators. It is
  typically beneficial when computing to medium/high precision.
  - FunctionDefaults::set_project_randomize(bool) constrols the use
  of randomization in projecting from an analytic form (i.e., C++) 
  into the discontinuous spectral element basis.  It is typically 
  beneficial unless there is already a good static data distribution.
  Since these options are straightforward to enable, this example
  focusses on static data redistribution.

  The process map (an instance of WorldDCPmapInterface) controls
  mapping of data to processors and it is actually quite easy to write
  your own (e.g., see WorldDCDefaultPmap or LevelPmap) that ensure
  uniform data distribution.  However, you also seek to incorporate
  estimates of the computational cost into the distribution.  The
  class LBDeuxPmap (deux since it is the second such class) in
  mra/lbdeux.h does this by examining the functions you request and
  using provided weights to estimate the computational cost.
  Communication costs are proportional to the number of broken links
  in the tree.  Since some operations work in the scaling function
  basis, some in the wavelet basis, and some in non-standard form,
  there is an element of empiricism in getting best performance from
  most algorithms.

  \par Implementation

  The example times how long it takes to perform the following operations
  - constructing 4000 Gaussian functions with random exponent and origin,
  - truncating them,
  - differentiating them, and
  - applying the Coulomb Green's function to them.

  The process map (data distribution) is then modified using the LBDeux
  heuristic and the operations repeated.

  \par Results

  \verbatim
  /tmp/work/harrison % aprun -n 50 -d 12 -cc none ./dataloadbal 
  Runtime initialized with 10 threads in the pool and affinity 1 0 2
  Before load balancing
  project 2.57 truncate 9.18 differentiate 11.80 convolve 12.57 balance 0.00
  project 2.74 truncate 9.08 differentiate 11.33 convolve 13.07 balance 0.00
  project 2.80 truncate 8.84 differentiate 10.76 convolve 12.66 balance 2.19
  After load balancing
  project 2.96 truncate 2.56 differentiate 3.28 convolve 9.81 balance 0.00
  project 3.05 truncate 2.61 differentiate 3.46 convolve 9.24 balance 0.00
  project 3.01 truncate 2.60 differentiate 3.25 convolve 9.30 balance 0.00
  \endverbatim

*/


#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/operator.h>
#include <mra/vmra.h>
#include <mra/lbdeux.h>
#include <constants.h>
using namespace madness;

static const int NFUNC = 4000;

// A class that behaves like a function to compute a Gaussian of given origin and exponent
class Gaussian : public FunctionFunctorInterface<double,3> {
public:
    const coord_3d center;
    const double exponent;
    const double coefficient;
    std::vector<coord_3d> specialpt;

    Gaussian(const coord_3d& center, double exponent, double coefficient)
        : center(center), exponent(exponent), coefficient(coefficient), specialpt(1)
    {
        specialpt[0][0] = center[0]; 
        specialpt[0][1] = center[1]; 
        specialpt[0][2] = center[2]; 
    }

    // MADNESS will call this interface
    double operator()(const coord_3d& x) const {
        double sum = 0.0;
        for (int i=0; i<3; i++) {
            double xx = center[i]-x[i];
            sum += xx*xx;
        };
        return coefficient*exp(-exponent*sum);
    }
    
    // By default, adaptive projection into the spectral element basis
    // starts uniformly distributed at the initial level.  However, if
    // a function is "spiky" it may be necessary to project at a finer
    // level but doing this uniformly is expensive.  This method
    // enables us to tell MADNESS about points/areas needing deep
    // refinement (the default is no special points).
    std::vector<coord_3d> special_points() const {
        return specialpt;
    }
};

// Makes a new square-normalized Gaussian functor with random origin and exponent
real_functor_3d random_gaussian() {
    const double expntmin=1e-1;
    const double expntmax=1e4;
    const real_tensor& cell = FunctionDefaults<3>::get_cell();
    coord_3d origin;
    for (int i=0; i<3; i++) {
        origin[i] = RandomValue<double>()*(cell(i,1)-cell(i,0)) + cell(i,0);
    }
    double lo = log(expntmin);
    double hi = log(expntmax);
    double expnt = exp(RandomValue<double>()*(hi-lo) + lo);
    double coeff = pow(2.0*expnt/constants::pi,0.75);
    return real_functor_3d(new Gaussian(origin,expnt,coeff));
}

// This structure is used to estimate the cost of computing on a block of coefficients
// The constructor saves an estimate of the relative cost of computing on
// leaf (scaling function) or interior (wavelet) coefficients.  Unless you 
// are doing nearly exclusively just one type of operation the speed is
// not that sensitive to the choice, so the default is usually OK.
//
// The operator() method is invoked for each block of coefficients
// (function node) to estimate the cost. Since pretty much everything
// involves work at the top of the tree we give levels 0 and 1 a 100x
// multiplier to avoid them being a bottleneck.
struct LBCost {
    double leaf_value;
    double parent_value;
    LBCost(double leaf_value=1.0, double parent_value=1.0) 
        : leaf_value(leaf_value)
        , parent_value(parent_value) 
    {}

    double operator()(const Key<3>& key, const FunctionNode<double,3>& node) const {
        if (key.level() <= 1) {
            return 100.0*(leaf_value+parent_value);
        }
        else if (node.is_leaf()) {
            return leaf_value;
        }
        else {
            return parent_value;
        }
    }
};

void test(World& world, bool doloadbal=false) {
    double start;
    vector_real_function_3d f(NFUNC);

    default_random_generator.setstate(99); // Ensure all processes have the same state

    // By default a sychronization (fence) is performed after projecting a new function
    // but if you are projecting many functions this is inefficient.  Here we turn
    // off the fence when projecting, and do it manually just once.  
    start = wall_time();
    for (int i=0; i<NFUNC; i++) 
        f[i] = real_factory_3d(world).functor(random_gaussian()).nofence();
    world.gop.fence();
    double projection = wall_time() - start;

    start = wall_time();
    truncate(world, f);
    double truncation = wall_time() - start;

    start = wall_time();
    Derivative<double,3> Dx(world,0);
    apply(world, Dx, f); // Computes vector of derivatives and discards result
    double differentiation = wall_time() - start;

    // Make the Coulomb operator
    real_convolution_3d op = CoulombOperator(world, 1e-2, 1e-6);

    start = wall_time();
    apply(world, op, f); // Applies Coulomb GF and discards result
    double convolution = wall_time() - start;

    start = wall_time();
    if (doloadbal) {
        LoadBalanceDeux<3> lb(world);
        for (int i=0; i<NFUNC; i++)
            lb.add_tree(f[i], LBCost(2.0,1.0));

        // Calling redistribute() installs the new process map and
        // redistributes all functions using the old process map.
        // This is almost always what you want.  However, since we
        // know that we are about to throw away our functions (f) we
        // could simply call set_pmap() that installs the new map but
        // does not redistribute.
        FunctionDefaults<3>::redistribute(world, lb.load_balance(2.0,false));
    }
    double loadbal = wall_time() - start;

    if (world.rank() == 0) printf("project %.2f truncate %.2f differentiate %.2f convolve %.2f balance %.2f\n",
                                  projection, truncation, differentiation, convolution, loadbal);
}

int main(int argc, char** argv) {
  // Initialize the parallel programming environment
  initialize(argc,argv);
  World world(MPI::COMM_WORLD);

  // Load info for MADNESS numerical routines
  startup(world,argc,argv);

  // Set/override defaults
  FunctionDefaults<3>::set_cubic_cell(-20,20);
  FunctionDefaults<3>::set_apply_randomize(false);
  FunctionDefaults<3>::set_project_randomize(false);
  FunctionDefaults<3>::set_truncate_on_project(true);

  // First three without data redistribution
  if (world.rank() == 0) print("Before load balancing");
  test(world);
  test(world);
  test(world, true);

  // At end of last test data was redistributed, repeat again three times
  if (world.rank() == 0) print("After load balancing");
  test(world);
  test(world);
  test(world);

  finalize();

  return 0;
}
