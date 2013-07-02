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

  $Id: nonlinsol.h 2154 2011-02-20 15:41:52Z rjharrison $
*/

#ifndef MADNESS_EXAMPLES_NONLINSOL_H__INCLUDED
#define MADNESS_EXAMPLES_NONLINSOL_H__INCLUDED

/*!
  \file examples/nonlinsol.h
  \brief Example implementation of Krylov-subspace nonlinear equation solver 
  \defgroup nonlinearsolve Simple Krylov-subspace nonlinear equation solver 
  \ingroup examples

  This class implements the solver described in 
  \verbatim
   R. J. Harrison, Krylov subspace accelerated inexact newton method for linear
   and nonlinear equations, J. Comput. Chem. 25 (2004), no. 3, 328–334.
  \endverbatim
 */

namespace madness {
    /// A simple Krylov-subspace nonlinear equation solver 

    /// \ingroup nonlinearsolve
    class NonlinearSolver {
        const unsigned int maxsub; //< Maximum size of subspace dimension
	vector_real_function_3d ulist, rlist; ///< Subspace information
	real_tensor Q;
    public:
	NonlinearSolver(unsigned int maxsub = 10) : maxsub(maxsub) {}

	/// Computes next trial solution vector

	/// You are responsible for performing step restriction or line search
	/// (not necessary for linear problems).
	///
	/// @param u Current solution vector
	/// @param r Corresponding residual
	/// @return Next trial solution vector
	real_function_3d update(const real_function_3d& u, const real_function_3d& r) {
	    int iter = ulist.size();
	    ulist.push_back(u);
	    rlist.push_back(r);

	    // Solve subspace equations
	    real_tensor Qnew(iter+1,iter+1);
	    if (iter>0) Qnew(Slice(0,-2),Slice(0,-2)) = Q;
	    for (int i=0; i<=iter; i++) {
		Qnew(i,iter) = inner(ulist[i],rlist[iter]);
		Qnew(iter,i) = inner(ulist[iter],rlist[i]);
	    }
	    Q = Qnew;
	    real_tensor c = KAIN(Q);

	    // Form new solution in u
	    real_function_3d unew = real_factory_3d(u.world());
	    unew.compress();
	    for (int i=0; i<=iter; i++) {
		unew.gaxpy(1.0,ulist[i], c[i]); 
		unew.gaxpy(1.0,rlist[i],-c[i]); 
	    }
	    unew.truncate();

            if (ulist.size() == maxsub) {
                ulist.erase(ulist.begin());
                rlist.erase(rlist.begin());
                Q = copy(Q(Slice(1,-1),Slice(1,-1)));
            }
	    return unew;
	}
    };
}
#endif
