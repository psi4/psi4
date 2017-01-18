/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup ccresponse
    \brief Handle the three tensors needed for Raman Optical Activity Scattering.

    ROA Scattering requires the following polarizability tensors for displaced geometries:
      (1) electric-dipole/electric-dipole;
      (2) electric-dipole/electric-quadrupole; and
      (3) electric-dipole/magnetic-dipole.

  -TDC, August 2009
*/
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"
#include "psi4/physconst.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

/*
 *  DELETE THIS FILE IF YOU WANT.  IT SHOULD HAVE BEEN  PUSHED PREVIOUSLY, BUT THAT WOULD
 *  NOT APPEAR TO BE THE CASE SINCE WHEN I PULLED CHANGES, IT SAID I HAD AN UNTRACKED
 *  FILE THAT WAS THE SCATTER.CC FILE.  SO HERE IT IS AS AN ALTERNATE PRESENTATION OF THE FILE
 *  THAT IS COMPLETELY UNECESSARY, BUT MAYBE HELPFUL.
 *
 *  -HM
 */


namespace psi { namespace ccresponse {

/* Prototype the tensor derivative function */
std::vector < SharedMatrix > compute_tensor_deriv(std::vector < SharedMatrix > tensor_list, const double disp_size);

/* Computes all the ROA data using all the tensors */
void scatter2(void)
{
	printf("Scattering FUN-ction.\n");

	//* Put Python Lists of Lists (the various Tensors) into Vectors of Matrices
	//* --> This means that the scatter function will need to take python lists as arguments.

	//* Compute Derivatives of Tensors (perhaps using the function defined below, but it's
	//*  								really just a suggestion).

	//* Adapt TDC's roa.c Code
	 /*
	  *  This will involve aquiring other pieces of needed data (like the reference geomertry
	  *  number of atoms, etc) using PSI4 style objects, probably.
      *
      *  Also, the Hessian matrix and dipole moment derivatives are needed as well,
	  *  which for now can just be read in from files in the top level roa job directory,
	  *  either by python then fed into the "scatter" function, or simply just read in by
      *  the "scatter" function.
	  *
	  */

}

/* Computes tensor derivatives */
std::vector < SharedMatrix > compute_tensor_deriv(std::vector < SharedMatrix > tensor_list, const double disp_size)
{
	/* Set up the Vector of Atom_Coord Derivative Tensors */
	std::vector < SharedMatrix > der_tensors;

	/* Set up a temporary matrix. This convolutedly grabs the
 	*  appropriate size.  Polarizability and Opt. Rotation tensors are 3x3,
 	*  but the Dipole/Quad. is 3x9, I think.
 	*/
	Matrix temp;
	//temp.copy(tensor_list.get(0));
	temp = tensor_list.get(0);

	/* Take the derivatives --> ("atom_x_+" - "atom_x_-")/(2.0*disp_size) */
	int ntensors = tensor_list.size();
	double gfactor = 1/(2.0 * disp_size);

	for(int i=0; i < ntensors; ++i)  {
	  int p = i*2;
	  int m = i*2 + 1;
	  //printf("p=%d, m=%d\n",pgeom,mgeom);
	  //temp.copy((tensor_list.get(p)->subtract(tensor_list.get(m)));
	  temp.zero();
	  temp.copy(tensor_list.get(p))
	  temp.subtract(tensor_list.get(m));
	  temp.scale(gfactor);
	  der_tensors.push_back(temp);
	}

	return der_tensors;
}

}} // namespace psi::ccresponse
