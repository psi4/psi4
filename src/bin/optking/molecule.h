/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

/*!
  \ingroup optking
  \file molecule.h : header for MOLECULE class
*/
#ifndef _opt_molecule_h_
#define _opt_molecule_h_

#include "package.h"

#include "frag.h"
#include "interfrag.h"
#include "fb_frag.h"
#include "print.h"
#include "libparallel/ParallelPrinter.h"
#include <fstream>

namespace opt {

class MOLECULE {

  vector<FRAG *> fragments;           // fragments with intrafragment coordinates
  vector<INTERFRAG *> interfragments; // interfragment coordinates
  vector<FB_FRAG *> fb_fragments;     // EFP fixed-body fragment - does not actually
                                      // store atoms (for now at least)
  double energy;

 public:

  MOLECULE(int num_atoms); // allocate molecule with one fragment of this size

  ~MOLECULE() {
    for (int i=0; i<fragments.size(); ++i)
      delete fragments[i];
    fragments.clear();
    for (int i=0; i<interfragments.size(); ++i)
      delete interfragments[i];
    interfragments.clear();
    for (int i=0; i<fb_fragments.size(); ++i)
      delete fb_fragments[i];
    fb_fragments.clear();
  }

  // if you have one fragment - make sure all atoms are bonded
  // if not, break up into multiple fragments
  void fragmentize(void);

  void add_interfragment(void);

  // Determine trivial, simple coordinate combinations.
  int form_trivial_coord_combinations(void);

  // Determine initial delocalized coordinate coefficients.
  int form_delocalized_coord_combinations(void);

  // Determine Pulay natural coordinate combinations.
  int form_natural_coord_combinations(void);

  int add_cartesians(void);

  int g_nfragment(void) const { return fragments.size(); };

  int g_nfb_fragment(void) const { return fb_fragments.size(); };

  int g_natom(void) const { // excludes atoms in fb fragments
    int n = 0;
    for (int f=0; f<fragments.size(); ++f)
      n += fragments[f]->g_natom();
    return n;
  }

  int Ncoord(void) const {
    int n=0;
    for (int f=0; f<fragments.size(); ++f)
      n += fragments[f]->Ncoord();
    for (int i=0; i<interfragments.size(); ++i)
      n += interfragments[i]->Ncoord();
    for (int e=0; e<fb_fragments.size(); ++e)
      n += fb_fragments[e]->Ncoord();
    return n;
  }

  int Ncoord_intrafragment(void) const {
    int n=0;
    for (int f=0; f<fragments.size(); ++f)
      n += fragments[f]->Ncoord();
    return n;
  }

  int Ncoord_interfragment(void) const {
    int n=0;
    for (int f=0; f<interfragments.size(); ++f)
      n += interfragments[f]->Ncoord();
    return n;
  }

  int Ncoord_fb_fragment(void) const {
    int n=0;
    for (int f=0; f<fb_fragments.size(); ++f)
      n += fb_fragments[f]->Ncoord();
    return n;
  }

  // given fragment index returns first atom in that fragment
  int g_atom_offset(int index) const {
    int n = 0;
    for (int f=1; f<=index; ++f)
      n += fragments[f-1]->g_natom();
    return n;
  }

  // given fragment index, returns the absolute index
  // for the first internal coordinate on that fragment
  int g_coord_offset(int index) const {
    int n = 0;
    for (int f=0; f<index; ++f)
      n += fragments[f]->Ncoord();
    return n;
  }

  // given interfragment index, returns the absolute index for the first 
  // internal coordinate of that interfragment set
  int g_interfragment_coord_offset(int index) const {
    int n = Ncoord_intrafragment();
    for (int f=0; f<index; ++f)
      n += interfragments[f]->Ncoord();
    return n;
  }

  // given fb_fragment index, returns the absolute index for the first
  // internal coordinate of that fb set
  int g_fb_fragment_coord_offset(int index) const {
    int n = Ncoord_intrafragment() + Ncoord_interfragment();
    for (int f=0; f<index; ++f)
      n += fb_fragments[f]->Ncoord();
    return n;
  }

  double g_energy(void) const { return energy; }

  void update_connectivity_by_distances(void) {
    for (int i=0; i<fragments.size(); ++i)
      fragments[i]->update_connectivity_by_distances();
  }

  void update_connectivity_by_bonds(void) {
    for (int i=0; i<fragments.size(); ++i)
      fragments[i]->update_connectivity_by_bonds();
  }

  void print_connectivity(std::string psi_fp, FILE *qc_fp) const {
    for (int i=0; i<fragments.size(); ++i)
      fragments[i]->print_connectivity(psi_fp, qc_fp, i, g_atom_offset(i));
  }

  void print_geom(std::string psi_fp, FILE *qc_fp, bool print_mass = false) {
    for (int i=0; i<fragments.size(); ++i)
      fragments[i]->print_geom(psi_fp, qc_fp, i, print_mass);
  }

  void print_xyz(int iter_shift = 0);

  void print_geom_grad(std::string psi_fp, FILE *qc_fp, bool print_mass = false) {
    for (int i=0; i<fragments.size(); ++i)
      fragments[i]->print_geom_grad(psi_fp, qc_fp, i, print_mass);
  }

  void print_coords(std::string psi_fp, FILE *qc_fp) const;
  void print_simples(std::string psi_fp, FILE *qc_fp) const;

  // print definition of an internal coordinate from global index 
  std::string get_coord_definition_from_global_index(int coord_index) const;

  void update_fb_values(void);

  void print_intco_dat(std::string psi_fp_coord, FILE *qc_fp);

  int add_intrafragment_simples_by_connectivity(void) {
    int n=0;
    for (int i=0; i<fragments.size(); ++i)
      n += fragments[i]->add_simples_by_connectivity();
    return n;
  }

  int add_intrafragment_auxiliary_bonds(void) {
    int n=0;
    for (int i=0; i<fragments.size(); ++i)
      n += fragments[i]->add_auxiliary_bonds();
    return n;
  }

  // compute coord values from frag member geometries
  double * coord_values(void) const {
    GeomType x = g_geom_2D();
    double *q = coord_values(x);
    return q;
  }

  // compute coord values from given geometry ; empty space for EFP coordinates included
  double * coord_values(GeomType new_geom) const {
    double *q, *q_frag, *q_IF;
    q = init_array(Ncoord());

    for (int f=0; f<fragments.size(); ++f) {
      q_frag = fragments[f]->coord_values( &(new_geom[g_atom_offset(f)]) );

      for (int i=0; i<fragments[f]->Ncoord(); ++i)
        q[g_coord_offset(f)+i]  = q_frag[i];

      free_array(q_frag);
    }

    for (int I=0; I<interfragments.size(); ++I) {
      int A_index = interfragments[I]->g_A_index();
      int B_index = interfragments[I]->g_B_index();
      
      q_IF = interfragments[I]->coord_values( &(new_geom[g_atom_offset(A_index)]),
        &(new_geom[g_atom_offset(B_index)]) );

      for (int i=0; i<interfragments[I]->Ncoord(); ++i)
        q[g_interfragment_coord_offset(I)+i]  = q_IF[i];

      free_array(q_IF);
    }

    return q;
  }

  void write_geom(void);
  void symmetrize_geom(void);
  void print_geom_out(void);

  double ** compute_B(void) const;
  double ** compute_derivative_B(int coord_index) const ;

  double ** compute_G(bool use_masses=false) const;

  double * g_grad_array(void) const {
    int f, i;
    double *g, *g_frag;

    g = init_array(3*g_natom());
    for (f=0; f<fragments.size(); ++f) {
      g_frag = fragments[f]->g_grad_array();
      for (i=0; i<3*fragments[f]->g_natom(); ++i)
        g[3*g_atom_offset(f)+i] = g_frag[i];
      free_array(g_frag);
    }
    return g;
  }

  double * g_masses(void) const;
  double * g_Z(void) const;
  double * g_u_vector(void) const; // reciprocal masses in vector

  double * g_geom_array(void) {
    double *g, *g_frag;

    g = init_array(3*g_natom());
    for (int f=0; f<fragments.size(); ++f) {
      g_frag = fragments[f]->g_geom_array();
      for (int i=0; i<3*fragments[f]->g_natom(); ++i)
        g[3*g_atom_offset(f)+i] = g_frag[i];
      free_array(g_frag);
    }
    return g;
  }

  double ** g_geom_2D(void) const {
    double **g_frag;
    double **g = init_matrix(g_natom(),3);

    for (int f=0; f<fragments.size(); ++f) {
      g_frag = fragments[f]->g_geom();
      for (int i=0; i<fragments[f]->g_natom(); ++i)
        for (int xyz=0; xyz<3; ++xyz)
          g[g_atom_offset(f)+i][xyz] = g_frag[i][xyz];
      free_matrix(g_frag);
    }
    return g;
  }

  double ** g_grad_2D(void) const {
    double **g, *g_frag;

    g = init_matrix(g_natom(),3);
    for (int f=0; f<fragments.size(); ++f) {
      g_frag = fragments[f]->g_grad_array();
      int cnt=0;
      for (int i=0; i<fragments[f]->g_natom(); ++i)
        for (int xyz=0; xyz<3; ++xyz)
          g[g_atom_offset(f)+i][xyz] = g_frag[cnt++];
      free_array(g_frag);
    }
    return g;
  }

  void H_guess(void) const;
  double **Lindh_guess(void) const;

  void forces(void);
  void apply_constraint_forces(void);
  bool has_fixed_eq_vals(void);
  void project_f_and_H(void);
  void project_dq(double *);
  void irc_step(void);
  void nr_step(void);
  void rfo_step(void);
  void prfo_step(void);
  void backstep(void);
  void sd_step(void);
  void sd_step_cartesians(void);
  void linesearch_step(void);

  void apply_intrafragment_step_limit(double * & dq);
  void check_intrafragment_zero_angles(double const * const dq);

  void set_geom_array(double * array_in) {
    for (int f=0; f<fragments.size(); ++f)
      fragments[f]->set_geom_array( &(array_in[3*g_atom_offset(f)]) );
  }

  void fix_tors_near_180(void) {
    for (int f=0; f<fragments.size(); ++f)
      fragments[f]->fix_tors_near_180();
    for (int I=0; I<interfragments.size(); ++I)
      interfragments[I]->fix_tors_near_180();
  }

/*
  void check_tors_for_bad_angles(void) {
    for (int f=0; f<fragments.size(); ++f)
      fragments[f]->check_tors_for_bad_angles();
    for (int I=0; I<interfragments.size(); ++I)
      interfragments[I]->check_tors_for_bad_angles();
  }
*/

  void test_B(void);
  void test_derivative_B(void);

  bool cartesian_H_to_internals(double **H_cart) const;

  void set_masses(void) {
    for (int f=0; f<fragments.size(); ++f)
      fragments[f]->set_masses();
  }

  bool read_coords(std::ifstream & fin);
  // function to obtain geometry and gradient
  void read_geom_grad(void);

  // Compute constraint matrix
  double ** compute_constraints(void);

  void add_fb_fragments(void);

  // freeze interfragment modes that break symmetry
  void freeze_interfragment_asymm(void);

  // determine whether a linear combination of coords breaks symmetry
  bool coord_combo_is_symmetric(double *coord_combo, int dim);

  // Apply string list of user-specified internals to be frozen.
  bool apply_input_constraints();

};

}

#endif

