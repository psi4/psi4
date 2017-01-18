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

/*!
   \ingroup optking
   \file frag.h : header for fragment class
*/
#ifndef _opt_frag_h_
#define _opt_frag_h_

#include <cstdio>
#include <vector>
#include <string>
#include "psi4/libparallel/ParallelPrinter.h"
#include "print.h"
#include "coordinates.h"
#include "psi4/psi4-dec.h"

namespace opt {

class INTERFRAG;
using std::vector;

/*!
  \ingroup optking
  \class   FRAG
  \brief   A group of atoms, its geometry, and its internal coordinates.
  */
class FRAG {

 protected:              // private, except the EFP derived class can access
  int natom;             ///< number of atoms in fragment
  double *Z;             ///< atomic numbers
  double **geom;         ///< cartesian coordinates
  double **grad;         ///< cartesian coordinates
  double *mass;          ///< nuclear masses
  bool   **connectivity; ///< connectivity matrix
  bool frozen;           ///< whether to optimize
  COMBO_COORDINATES coords; ///< simple or linear combinations of simple coordinates

 public:
  friend class INTERFRAG;

  /** Constructor to create fragment.

    There is one constructor to insert the atoms, and a different constructor available
    just to allocate the memory.
    \param natom_in number of atoms.
    \param Z_in     atomic charges
    \param geom_in  Cartesian coordinates (natom_in x 3)
    */
  FRAG(int natom_in, double *Z_in, double **geom_in); // use if Z and geom are known
  FRAG(int natom_in);

  ~FRAG(); // free memory

  /** Description
    \returns number of atoms.
    */
  int g_natom(void) const { return natom; };

  void set_masses(void);

  void print_geom(std::string psi_fp, FILE *qc_fp, const int id, bool print_mass = false);
  void print_geom_grad(std::string psi_fp, FILE *qc_fp, const int id, bool print_mass = false);
  void print_simples(std::string psi_fp, FILE *qc_fp, int atom_offset=0) const;
  void print_coords(std::string psi_fp, FILE *qc_fp, int atom_offset=0) const;
  void print_combinations(std::string psi_fp, FILE *qc_fp) const;
  void print_intco_dat(std::string psi_fp, FILE *qc_fp, int atom_offset=0) const;

  std::string get_coord_definition(int coord_index, int atom_offset=0);
  std::string get_simple_definition(int simple_index, int atom_offset=0);

  INTCO_TYPE get_simple_type(int simple_index);

  void update_connectivity_by_distances(void);
  void update_connectivity_by_bonds(void);

  void print_connectivity(std::string psi_fp, FILE *qc_fp, const int id, const int offset = 0) const ;

  // add simple internals based on connectivity; return number added
  int add_stre_by_connectivity(void);
  int add_bend_by_connectivity(void);
  int add_tors_by_connectivity(void);
  int add_cartesians(void);

  int form_trivial_coord_combinations(void);
  void add_trivial_coord_combination(int simple_id);
  int form_delocalized_coord_combinations(void);
  int form_natural_coord_combinations(void);

  int add_simples_by_connectivity(void) {
    int n;
    n  = add_stre_by_connectivity();
    n += add_bend_by_connectivity();
    n += add_tors_by_connectivity(); // but check bond angles
    return n;
  }

  int add_auxiliary_bonds(void);

  // add connectivity between two atoms; atom numbers are for fragment
  void connect(int i, int j) {
    connectivity[i][j] = connectivity[j][i] = true;
  }

  // Compute B matrix for only this fragment
  double ** compute_B(void) const ;

  // Compute B matrix. Use prevously allocated memory.  Offsets are ideal for molecule.
  void compute_B(double **B_in, int coord_offset, int atom_offset) const ;

  // Compute B only for the simple coordinates.
  //void compute_B_simples(double **B, int coord_offset, int atom_offset) const;

  void compute_G(double **, bool use_masses=false) const;

  // Compute B' matrix for one coordinate in given memory with given geometry.
  void compute_derivative_B(GeomType g, int coord_index, double **Bprime, int atom_offset) const;

  // Compute B' matrix for one coordinate in given memory with present geometry.
  void compute_derivative_B(int coord_index, double **Bprime, int atom_offset) const;

  // Compute B' matrix for one coordinate for fragment.
  double ** compute_derivative_B(int coord_index) const;

  // compute and print B matrix (for debugging)
  void print_B(std::string psi_fp, FILE *qc_fp) const ;

  // check nearness to 180 and save value
  void fix_tors_near_180(void);

  // check nearness to 180 and save value
  void fix_oofp_near_180(void);

  // Fix bend axes for consistency during displacments
  void fix_bend_axes(void);
  void unfix_bend_axes(void);

  // check if interior angles of torsion are near 0 or linear
  //bool check_tors_for_bad_angles(void) const;

  // return number of intrafragment coordinates
  int Ncoord(void) const { return coords.index.size(); };

  // The following 2 functions are only used by the B matrix testing routines.
  // return natom in definition of coord # coord_index
  int g_simple_natom(const int coord_index) const {
    return coords.simples.at(coord_index)->g_natom();
  }

  // return atom i in definition of coord # coord_index
  int g_simple_atom(const int coord_index, const int atom) const {
    return coords.simples.at(coord_index)->g_atom(atom);
  }

  // return array of atomic numbers
  double *g_Z(void) const;
  double *g_Z_pointer(void) { return Z; }

  //print s vectors to output file
  void print_s(std::string psi_fp, FILE *qc_fp, GeomType geom) const {
    coords.print_s(psi_fp, qc_fp, geom);
    return;
  }

  // Get all values.
  double *coord_values(void) const;
  double *coord_values(GeomType geom) const;

  // Get one value.
  double coord_value(int lookup) const;
  double coord_value(GeomType geom, int lookup) const;

  // Identify wayward angles
  std::vector<int>  validate_angles(double const * const dq, int atom_offset);

  // is simple one already present?
  bool present(const SIMPLE_COORDINATE *one) const;

  int find(const SIMPLE_COORDINATE *one) const;

  // displace fragment by dq ; forces and offset are provided for printing
  void displace(double *dq, double *fq, int atom_offset=0);
  // utility used by displace
  bool displace_util(double *dq, bool focus_on_constraints);

  double ** g_geom_pointer(void) { return geom; };           // returns pointer
  double ** g_geom(void) const;                              // returns a copy
  GeomType g_geom_const_pointer(void) const { return geom;}; // returns const pointer

  double ** g_grad(void);
  double ** g_grad_pointer(void) {return grad;};

  double * g_geom_array(void);
  double * g_grad_array(void);
  void set_geom_array(double * geom_array_in);
  void set_geom(double ** geom_in);
  void set_grad(double **grad_in);

  void print_geom(std::string psi_fp, FILE *qc_fp); // write cartesian geometry out for next step
  void print_geom_irc(std::string psi_fp, FILE *qc_fp); // write cartesian geometry out for next step

  double ** H_guess(void);
  // function to help with Lindh guess hessian
  double Lindh_rho(int A, int B, double RAB) const;
  // function to help with Lindh guess hessian - original constants
  double **Lindh_guess(void);
  bool **g_connectivity(void) const;
  const bool * const * g_connectivity_pointer(void) const;

  bool read_coord(std::vector<std::string> & tokens, int first_atom_offset);

  // return matrix of constraints on coordinates
  double ** compute_constraints(void) const;

  // add any missing hydrogen bonds within the fragment
  // return number added
  int add_hbonds(void);

  void simple_add(SIMPLE_COORDINATE * i) {
    coords.simples.push_back(i);
  }

  double g_mass(int i) { return mass[i]; }

  // relating to frozen fragments
  bool is_frozen(void) const { return frozen; }
  void freeze(void)   { frozen = true; }
  void unfreeze(void) { frozen = false; }

  // freeze coords within fragments
  void freeze_coords(void);

  /**
   * Compute center of mass of given geometry
   */
  double *com(GeomType in_geom);
  /**
   * Compute center of mass of fragment geometry
   */
  double *com(void) { return (com(geom)); }
  /**
   * Compute intertia tensor of given geometry
   */
  double **inertia_tensor (GeomType in_geom);
  /**
   * Compute intertia tensor of fragment geometry
   */
  double **inertia_tensor(void) { return (inertia_tensor(geom)); }
  /**
   * Compute principal axes of given geometry
   */
  int principal_axes(GeomType geom, double **axes, double *evals);
  /**
   * Compute principal axes of fragment geometry
   * @param evals are moments returned in ascending order; zero evals (and evects) are removed.
   * @param axes rows are principal axes.
   * @returns The number of non-zero principal axes.
   *
   * We may have to canonically order the degenerate evals later.
   */
  int principal_axes(double **axes, double *evals) {
    return (principal_axes(geom, axes, evals));
  }

  // These functions are needed for the forces() function
  // to apply user-defined equilibrium extra forces.
  bool coord_has_fixed_eq_val(int coord_index) const {
    return coords.simples[coord_index]->has_fixed_eq_val();
  }
  double coord_fixed_eq_val(int coord_index) const {
    return coords.simples[coord_index]->fixed_eq_val();
  }
  void add_combination_coord(vector<int> ids, vector<double> coeffs); // for molecule_read_coord

  /**
   * @param R_string string of atom pairs for frozen distances
   * @param B_string string of atom triples for frozen bends
   * @param D_string string of atom quartets for frozen dihedrals
   * @param C_string string with lists of atom and xyz specification for frozen cartesians
   * @returns True if any constraints are present.
  */
  bool apply_frozen_constraints(std::string R_string, std::string B_string, std::string D_string, std::string C_string);

  /**
   * @param R_string string of atom pairs + equilibrium value for fixed distances
   * @param B_string string of atom triplets + equilibrium value for fixed bends
   * @param D_string string of atom quartets + equilibrium value for fixed dihedrals
   * @returns True if any constraints are present.
  */
  bool apply_fixed_constraints(std::string R_string, std::string B_string, std::string D_string);

  void erase_combo_coord(int index) { coords.erase_combo(index); } ;

  /* Are any coordinates present that are not cartesians? */
  bool is_noncart_present(void) const;

};

}

#endif
