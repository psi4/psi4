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

/*! \file combo_coordinates.cc
    \ingroup optking
    \brief Functions for class handling linear combinations of simple internal or cartesian coordinates.
*/

#include "coordinates.h"
#include "psi4/optking/physconst.h"
#include <sstream>
#include "print.h"

using std::ostringstream;

namespace opt {

// Get value of one coordinate.
double COMBO_COORDINATES::value(GeomType geom, int lookup) const {
  double tval = 0.0;

  for (std::size_t s=0; s<index.at(lookup).size(); ++s)
    tval += coeff.at(lookup).at(s) * simples.at(index[lookup][s])->value(geom);

  return tval;
}

// Get all values.
double * COMBO_COORDINATES::values(GeomType geom) const {
  double *q = init_array(index.size());

  for (std::size_t cc=0; cc<index.size(); ++cc)
    q[cc] = value(geom, cc);

  return q;
}

// Fills in a provided B matrix row for one coordinate.
// offset specifies first atom in this fragment; allows direct use by molecule with additional atoms
// Calling routine must 0 memory.
bool COMBO_COORDINATES::DqDx(GeomType geom, int lookup, double *dqdx, int atom_offset) const {
  for (std::size_t s=0; s<index.at(lookup).size(); ++s) {          // loop over simples in combo
    //oprintf_out("simple %d, \n", index[lookup][s]);
    double **dqdx_simple = simples.at(index[lookup][s])->DqDx(geom);

    for (int j=0; j < simples[ index[lookup][s] ]->g_natom(); ++j) { // loop over atoms in s vector
      int atom = atom_offset + simples[ index[lookup][s] ]->g_atom(j);
      // oprintf_out("atom %d, \n", atom);

      for (int xyz=0; xyz<3; ++xyz) {
        //oprintf_out("coeff %10.5f, dqdx_simple %10.5f", coeff[lookup][s], dqdx_simple[j][xyz]);
        dqdx[3*atom + xyz] += coeff.at(lookup).at(s) * dqdx_simple[j][xyz];
      }
    }

    free_matrix(dqdx_simple);
  }
  return true;
}

// Fills in a B' derivative matrix for one coordinate.
// If the desired cartesian indices/dimension spans more than just one fragment, provide the atom offset.

// For a combination coordinate, shift of any atom within the fragment could change the value,
// so our dimension is (3*frag atoms) X (3*frag atoms)
// Memory is provided by calling function.  Must be zero'ed beforehand.
bool COMBO_COORDINATES::Dq2Dx2(GeomType geom, int lookup, double **dq2dx2, int atom_offset) const {

  for (std::size_t s=0; s<index.at(lookup).size(); ++s) {  // loop over simples
    // Dimension of dq2dx2_simple is 6x6, 9x9, or 12x12
    double **dq2dx2_simple = simples[ index[lookup].at(s) ]->Dq2Dx2(geom);

    int natom_simple = simples[ index[lookup][s] ]->g_natom();
    int atom_a, atom_b;

    for (int a=0; a<natom_simple; ++a) {
      atom_a = atom_offset + simples[ index[lookup][s] ]->g_atom(a);

      for (int b=0; b<natom_simple; ++b) {
        atom_b = atom_offset + simples[ index[lookup][s] ]->g_atom(b);

        for (int xyz_a=0; xyz_a<3; ++xyz_a)
          for (int xyz_b=0; xyz_b<3; ++xyz_b)
            dq2dx2[3*atom_a + xyz_a][3*atom_b + xyz_b] += coeff.at(lookup).at(s)
                                           * dq2dx2_simple[3*a + xyz_a][3*b + xyz_b];
      }
    }
  }
  return true;
}

// Clear current combinations.  Leaves simples in place
void COMBO_COORDINATES::clear_combos(void) {
  for (std::size_t i=0; i<index.size(); ++i)
    index[i].clear();
  for (std::size_t i=0; i<coeff.size(); ++i)
    coeff[i].clear();
  index.clear();
  coeff.clear();
  return;
}

void COMBO_COORDINATES::erase_combo(int cc) {
  index[cc].clear();
  coeff[cc].clear();
  index.erase(index.begin() + cc);
  coeff.erase(coeff.begin() + cc);
}

// Print out all s vectors
void COMBO_COORDINATES::print_s(std::string psi_fp, FILE *qc_fp, GeomType geom) const {
  oprintf(psi_fp, qc_fp, "\t---S vectors for internals---\n");
  for(std::size_t cc=0; cc<index.size(); ++cc) {
    oprintf_out("Coordinate %d\n", cc+1);
    for(std::size_t s=0; s<index[cc].size(); ++s) {
      oprintf_out("\tCoeff %15.10lf\n", coeff.at(cc).at(s));
      simples[ index[cc][s] ]->print_s(psi_fp, qc_fp, geom);
    }
  }
}

// Print out one coordinate.
void COMBO_COORDINATES::print(std::string psi_fp, FILE *qc_fp, int cc, GeomType geom, int off) const {
  if (index[cc].size() == 1)
    simples[ index[cc][0] ]->print(psi_fp, qc_fp, geom, off);
  else {
    for(std::size_t s=0; s<index[cc].size(); ++s) {
      //oprintf_out(" %d (%10.5f)", index[cc][s]+1, coeff.at(cc).at(s));
      oprintf_out("\t(%10.5f)\n", coeff.at(cc).at(s));
      simples[ index[cc][s] ]->print(psi_fp, qc_fp, geom, off);
    }
  }
}

// Function to print the change in a combo coordinate
void COMBO_COORDINATES::print_disp(std::string psi_fp, FILE *qc_fp, int lookup, const double q_orig, const double f_q,
    const double dq, const double new_q, int atom_offset) const {

  if (index[lookup].size() == 1) // Simple printout style
    simples[ index[lookup][0] ]->print_disp(psi_fp, qc_fp, q_orig, f_q, dq, new_q, atom_offset);
  else {
    ostringstream iss(ostringstream::out);
    iss << "CC(" << lookup+1 << ")" << std::flush;
    oprintf(psi_fp, qc_fp, "%-15s = %13.6lf%13.6lf%13.6lf%13.6lf\n",
      iss.str().c_str(), q_orig*_bohr2angstroms, f_q*_hartree2aJ/_bohr2angstroms,
      dq*_bohr2angstroms, new_q*_bohr2angstroms);
  }
  return;
}

// For error message reporting, describe coordinate briefly
string COMBO_COORDINATES::get_coord_definition(int lookup, int atom_offset) const {
  ostringstream iss(ostringstream::out);

  if (index.at(lookup).size() == 1) // keep it short if possible
    iss << simples.at(index[lookup][0])->get_definition_string(atom_offset);
  else {
    for (std::size_t s=0; s<index.at(lookup).size(); ++s) {
      iss << index[lookup][s] + 1 << ":" << coeff.at(lookup).at(s) << ":";
      iss << simples.at(index[lookup][s])->get_definition_string(atom_offset);
    }
  }
  return iss.str();
}

// Transform a vector whose dimension is number of simples to number of coordinates.
double * COMBO_COORDINATES::transform_simples_to_combo(double *arr_simples) const {
  double *arr_combos = init_array(index.size());

  for(std::size_t cc=0; cc<index.size(); ++cc)
    for(std::size_t s=0; s<index[cc].size(); ++s)
      arr_combos[cc] += coeff.at(cc).at(s) * arr_simples[ index[cc][s] ];
  return arr_combos;
}

// Transform a matrix from simples by simples to combinations by combinations
double ** COMBO_COORDINATES::transform_simples_to_combo(double **mat_simples) const {

  std::size_t Ns = simples.size(); // num simples
  std::size_t Nc = index.size();   // num combos

  double **Y = init_matrix(Ns, Nc);

  for (std::size_t s1=0; s1<Ns; ++s1)
    for(std::size_t cc=0; cc<Nc; ++cc)
      for(std::size_t s2=0; s2<index[cc].size(); ++s2)
        Y[s1][cc] += mat_simples[s1][ index[cc][s2] ] *  coeff[cc][s2];

  double **mat_combo = init_matrix(Nc, Nc);

  for (std::size_t c1=0; c1<Nc; ++c1)
    for (std::size_t c2=0; c2<Nc; ++c2)
      for(std::size_t s1=0; s1<index[c1].size(); ++s1) // on s1 in c1 contributes
        mat_combo[c1][c2] += Y[ index[c1][s1] ][ c2 ] *  coeff[c1][s1];

  free_matrix(Y);
  return mat_combo;
}

}
