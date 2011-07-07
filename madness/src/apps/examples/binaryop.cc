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
  \file examples/binaryop.cc
  \brief Illustrates general composition of two functions
  \defgroup examplebinop Illustrates general composition of two functions
  \ingroup examples
  
  The source is <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local/trunk/src/apps/examples/binaryop.cc>here</a>.

  \par Points of interest
  - use of a binary operation to apply a complex operation to two functions
  - use of asymptotic analysis to ensure good behavior in presence of numerical noise

  \par Background

  In nuclear physics density functional theory it is necessary to compute
  functions of the form

  \f[
  U(r) = \frac{\Delta^2 (r)}{\rho^{2/3} (r)}
  \f]

  The functions \f$ \Delta \f$ and \f$ \rho \f$ are both expected to go to zero
  at large \f$ r \f$ as is the ratio (i.e., \f$ \Delta \f$ goes to zero faster
  than \f$ \rho \f$).  Moreover, \f$ \rho \f$ should everywhere be positive 
  and is not expected to be zero in the interior region (for ground states only?).

  \par Implementation

  The first problem is how to compose this operation inside MADNESS.
  One could square \f$ Delta \f$, use \c unaryop() to compute the
  negative fractional power of \f$ \rho() \f$, and then multiply the
  two.  With care (see below) this should work.  Easier, faster, and
  more accurate is to do all of the above at once. This is accomplished
  with a binary operation that acts upon the input function values and
  returns the result.

  The second and most significant problem is numerical noise in the
  functions that can lead to \f$ \rho \f$ being zero and even negative
  while \f$ Delta \f$ is non-zero.  However, the expected asymptotics
  tell us that if either \f$ Delta \f$ or \f$ \rho \f$ are so small
  that noise is dominating, that the result of the binary operation
  should be zero.  [Aside. To accomplish the same using a unary
  operation the operation that computes \f$ rho^{-2/3} \f$ should
  return zero if \f$ rho \f$ is small. But this precludes us from
  simultaneously using information about the size of \f$ \Delta \f$
  and does not ensure that both are computed at the same level of
  refinement.]

  Analysis is necessary.  The threshold at which to screen values to
  their asymptotic form depends on the problem, the accuracy of
  computation, and possibly the box size.  In this problem we choose
  \f[
  \Delta(r) = exp(- | r | )
  \f]
  and 
  \f[
  \rho(r) = exp(- 2 | r | ) = \Delta^2(r)
  \f]
  Thus, the exact result is
  \f[
  U(r) = \frac{\Delta^2 (r)}{\rho^{2/3} (r)} = exp( - 2 | r | / 3)
  \f]

  Note that the result has a smaller exponent than the two input
  functions and is therefore significant at a longer range.  Since we
  cannot generate information we do not have, once the input functions
  degenerate into numerical noise we must expect that the ratio is
  also just noise.  In the physical application, the potential \f$
  U(r) \f$ is applied to another function that is also decaying
  expoentially, which makes \em small noise at long range not
  significant.  By screening to the physically expected value of zero
  we therefore ensure correct physics.
  
*/

#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>

using namespace madness;

static const double L = 30;     // Half box size
static const long k = 8;        // wavelet order
static const double thresh = 1e-6; // precision

static const double small = thresh*1e-4;

double delta(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return exp(-sqrt(x*x+y*y+z*z+1e-6));
}

double uexact(const coord_3d& r) {
    return pow(delta(r),2.0/3.0);
}

// This functor is used to perform the binary operation
struct Uop {
    void operator()(const Key<3>& key, 
                    real_tensor U,
                    const real_tensor& Delta, 
                    const real_tensor& rho) const {
        ITERATOR(U,
                 double d = Delta(IND);
                 double p = rho(IND);
                 if (p<small || d<small) 
                     U(IND) = 0.0;
                 else 
                     U(IND) = d*d/pow(p,2.0/3.0);
                 );
    }

    template <typename Archive> 
    void serialize(Archive& ar) {}
};

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    
    startup(world, argc, argv);
    std::cout.precision(6);

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_cubic_cell(-L, L);
    FunctionDefaults<3>::set_initial_level(4);
    
    real_function_3d Delta = real_factory_3d(world).f(delta);
    Delta.truncate(); // Deliberately truncate to introduce numerical noise

    real_function_3d rho = Delta*Delta;
    rho.truncate(); // Deliberately truncate to introduce numerical noise

    real_function_3d U = binary_op(Delta, rho, Uop());

    double err = U.err(uexact);
    if (world.rank() == 0) print("Estimated error norm is ", err);

    // Make the exact result just for plotting
    real_function_3d exact = real_factory_3d(world).f(uexact);

    // Make a line plot from the origin along the x axis to examine the functions in detail
    coord_3d lo(0), hi(0);
    hi[0] = L;
    plot_line("binaryop.dat", 1001, lo, hi, Delta, rho, U, exact);

    finalize();
    return 0;
}

