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

/*! \file combo_coordinates.h
    \ingroup optking
    \brief Class for a collection of linear combinations of simple internal or cartesian coordinates.
*/

#ifndef _opt_combo_coordinates_h_
#define _opt_combo_coordinates_h_

#include <vector>
#include <string>

using std::vector;
using std::string;

namespace opt {

class COMBO_COORDINATES {

  private:

    vector<SIMPLE_COORDINATE *>  simples;  // collection of simple and/or cartesian coordinates
    vector<vector<int> >     index;  // collection of coordinate index for linear combination
    vector<vector<double> >  coeff;  // collection of coefficients for linear combination

  public:

  friend class FRAG;       // these can all add coordinates
  friend class FB_FRAG;
  friend class INTERFRAG; 

  // Get all values.
  double *values(GeomType geom) const;

  // Get one value.
  double value(GeomType geom, int lookup) const;

  // Fills in a provided B matrix row for one coordinate.
  // If the desired cartesian column indices/dimension spans the molecule, i.e.,
  // possibly more than just one fragment, then provide the atom offset.
  bool DqDx(GeomType geom, int lookup, double *dqdx, int frag_atom_offset=0) const;

  // Fills in a B' derivative matrix for one coordinate.
  // If the desired cartesian indices/dimension spans the molecule, i.e.,
  // possibly more than just one fragment, then provide the atom offset.
  bool Dq2Dx2(GeomType geom, int lookup, double **dq2dx2, int frag_atom_offset=0) const;

  // For error message reporting, describe coordinate briefly
  string get_coord_definition(int lookup, int frag_atom_offset) const;

  // Clear current combination coordinates. Leave simples in place.
  void clear_combos(void);

  // Remove a particular combination coordinate by index.
  void erase_combo(int cc);

  // Print s vectors.
  void print_s(std::string psi_fp, FILE *qc_fp, GeomType geom) const;

  // Print one coordinate
  void print(std::string psi_fp, FILE *qc_fp, int cc, GeomType geom, int off) const;

  // Function to print the change in a combo coordinate
  void print_disp(std::string psi_fp, FILE *qc_fp, int lookup, const double q_orig, const double f_q,
    const double dq, const double new_q, int atom_offset) const;

  // Transform vector from simples to combination via linear combination
  double * transform_simples_to_combo(double *arr_simples) const;

  // Transform vector from simples to combination via linear combination
  double ** transform_simples_to_combo(double **mat_simples) const;

  int Nsimples(void) const { return simples.size(); }

};

}

#endif