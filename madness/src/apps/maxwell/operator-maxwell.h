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
  
  $Id: operator-maxwell.h 1602 2009-12-27 19:53:06Z rjharrison $
*/
#ifndef __operator_maxwell__
#define __operator_maxwell__

#include <string>
#include <constants.h>
#include <mra/mra.h>
#include <linalg/gmres.h>

using namespace std;
using namespace madness;

typedef double_complex complexd;
typedef Function<complexd, 3> function;
typedef Vector<function, 3> vecfunc;

/// A class for storing all the pertinent data for an enveloped incident
/// pulse.
///
/// NOTE that all units are in atomic units, unless specified otherwise.
class EFieldOperator : public Operator<vecfunc> {
protected:
	const function &epshat;
	const complexd prefact;
	const vecfunc &gradlnepshat;
	const Function<double, 3> &box_mask; // mask to make the boundary 0
	const Function<double, 3> &grad_mask; // mask to damp out derivatives
	const SeparatedConvolution<double, 3> &G;

	void action(const vecfunc &in, vecfunc &out) const {
		function dotp = (gradlnepshat[0]*in[0] + gradlnepshat[1]*in[1]
			+ gradlnepshat[2]*in[2]) * box_mask;
		vecfunc preop;
		int i;

		for(i = 0; i < 3; ++i)
			preop[i] = diff(dotp, i) * grad_mask;
		dotp.clear();

		for(i = 0; i < 3; ++i) {
			dotp = epshat*in[i];
			dotp.compress();
			preop[i].compress();
			preop[i].gaxpy(complexd(1.0, 0.0), dotp, prefact);
			out[i] = in[i] - apply(G, preop[i]);
			out[i].truncate();
			dotp.clear();
			preop[i].clear();
		}
	}

public:
	/// Needs the epsilon-hat complex dielectric, the frequency omega,
	/// mu0, eps0, the gradient of ln(epshat), the box_mask function to make
	/// the function 0 at the boundary, the grad_mask function to smooth
	/// out derivatives near the boundaries, and the Poisson Green's function
	EFieldOperator(const function &_epshat, const double _omega,
		const double _mu0, const double _eps0, const vecfunc &_gradlnepshat,
		const Function<double, 3> &_box_mask,
		const Function<double, 3> &_grad_mask,
		const SeparatedConvolution<double, 3> &_G) : epshat(_epshat),
		prefact(complexd(0.0, -_omega*_eps0*_mu0)), gradlnepshat(_gradlnepshat),
		box_mask(_box_mask), grad_mask(_grad_mask), G(_G) {}
};

#endif
