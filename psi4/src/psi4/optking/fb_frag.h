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

#ifndef _opt_fb_frag_h_
#define _opt_fb_frag_h_

/* This class will define a fragment for use in situations where a fragment
is defined by space-fixed coordinates: 3 center-of-mass and 3 Euler angle
coordinates.

This class is intended to be used for fixed-body EFP optimizations in QChem.
In QChem, the gradient is already computed in these space-fixed coordinates.
Thus, it is not necessary to write the code for the B-matrix to do an
optimization.  Eventually, a B-matrix will be necessary in order to, for
example, transform a cartesian analytic Hessian into internal coordinates for
these optimizations.

Since the angular coodinates could be described as linear combinations of
positions of atoms. We might be able to use the principle axes to define the
coordinates (derived by WDA), but we would still have to relate these to the
space-fixed coordinates of the gradient computed by QChem.
*/

#include "package.h"

#include "frag.h"
#include "mem.h"

namespace opt {

class FB_FRAG : public FRAG {

  // values and forces are provided by QChem, not computed by optking
  // they will be of dimension 6
  double *values;
  double *forces;

  public:
  // we will build a dummy fragment with no atoms
  FB_FRAG() : FRAG(0) {
    values = init_array(6);
    forces = init_array(6);
  }

  ~FB_FRAG() {
    free_array(values);
    free_array(forces);
   }

  void set_values(double * values_in);
  void set_forces(double * forces_in);

  double * get_values_pointer(void) const { return values; }
  double * get_forces_pointer(void) const { return forces; }

  // we don't have a valid B matrix for these
  // add 6 bogus stretches
  void add_dummy_coords(int ndummy);

// We will assign a value of 0 to these on the first iteration; subsequently, the
// values will be calculated as total Delta(q) from the start of the optimization
  void print_intcos(std::string psi_fp, FILE *qc_fp);

/* Add function to return a string definition of FB fragment, if needed
  // return string of intco definition
  std::string get_coord_definition(int coord_index, int atom_offset_A=0, int atom_offset_B=0);
*/

  double **H_guess(void);

  // Tell QChem to update rotation matrix and com for FB fragment
  void displace (int fb_frag_index, double *dq);

};

}

#endif