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
#include "efp_frag.h"
#include "print.h"

#include <fstream>

namespace opt {

class MOLECULE {

  vector<FRAG *> fragments;           // fragments with intrafragment coordinates
  vector<INTERFRAG *> interfragments; // interfragment coordinates
//****AVC****//
//#if defined(OPTKING_PACKAGE_QCHEM)
  vector<EFP_FRAG *> efp_fragments;   // EFP fixed-body fragment - does not actually
                                      // store atoms (for now at least)
//#endif
//****AVC****//
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
//****AVC****//
//#if defined(OPTKING_PACKAGE_QCHEM)
    for (int i=0; i<efp_fragments.size(); ++i)
      delete efp_fragments[i];
    efp_fragments.clear();
//#endif
//****AVC****//
  }

  // if you have one fragment - make sure all atoms are bonded
  // if not, break up into multiple fragments
  void fragmentize(void);

  void add_interfragment(void);

  int g_nfragment(void) const { return fragments.size(); };

  int g_nefp_fragment(void) const { return efp_fragments.size(); };

  int g_natom(void) const { // excludes atoms in efp fragments
    int n = 0;
    for (int f=0; f<fragments.size(); ++f)
      n += fragments[f]->g_natom();

//****AVC****//
    n += 3*efp_fragments.size();
//****AVC****//

    return n;
  }

  int g_nintco(void) const {
    int n=0;
    for (int f=0; f<fragments.size(); ++f)
      n += fragments[f]->g_nintco();
    for (int i=0; i<interfragments.size(); ++i)
      n += interfragments[i]->g_nintco();
    for (int e=0; e<efp_fragments.size(); ++e)
      n += efp_fragments[e]->g_nintco();
    return n;
  }

  int g_nintco_intrafragment(void) const {
    int n=0;
    for (int f=0; f<fragments.size(); ++f)
      n += fragments[f]->g_nintco();
    return n;
  }

  int g_nintco_interfragment(void) const {
    int n=0;
    for (int f=0; f<interfragments.size(); ++f)
      n += interfragments[f]->g_nintco();
    return n;
  }

  int g_nintco_efp_fragment(void) const {
    int n=0;
    for (int f=0; f<efp_fragments.size(); ++f)
      n += efp_fragments[f]->g_nintco();
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
  int g_intco_offset(int index) const {
    int n = 0;
    for (int f=0; f<index; ++f)
      n += fragments[f]->g_nintco();
    return n;
  }

  // given interfragment index, returns the absolute index for the first 
  // internal coordinate of that interfragment set
  int g_interfragment_intco_offset(int index) const {
    int n = g_nintco_intrafragment();
    for (int f=0; f<index; ++f)
      n += interfragments[f]->g_nintco();
    return n;
  }

  // given efp_fragment index, returns the absolute index for the first
  // internal coordinate of that efp set
  int g_efp_fragment_intco_offset(int index) const {
    int n = g_nintco_intrafragment() + g_nintco_interfragment();
    for (int f=0; f<index; ++f)
      n += efp_fragments[f]->g_nintco();
    return n;
  }

  double g_energy(void) const { return energy; }

  void update_connectivity_by_distances(void) {
printf("\nupdate_connectivity_by_distances(), fragments.size() = %d\n", fragments.size());
    for (int i=0; i<fragments.size(); ++i)
      fragments[i]->update_connectivity_by_distances();
printf("\nupdate_connectivity_by_distances(), fragments.size() = %d\n", fragments.size());
  }

  void update_connectivity_by_bonds(void) {
    for (int i=0; i<fragments.size(); ++i)
      fragments[i]->update_connectivity_by_bonds();
  }

  void print_connectivity(FILE *fout) const {
    for (int i=0; i<fragments.size(); ++i)
      fragments[i]->print_connectivity(fout, i, g_atom_offset(i));
  }

  void print_geom(FILE *fout, bool print_mass = false) {
    for (int i=0; i<fragments.size(); ++i)
      fragments[i]->print_geom(fout, i, print_mass);
  }

  void print_geom_grad(FILE *fout, bool print_mass = false) {
    for (int i=0; i<fragments.size(); ++i)
      fragments[i]->print_geom_grad(fout, i, print_mass);
  }

  void print_intcos(FILE *fout) {
    int a,b;
    for (int i=0; i<fragments.size(); ++i) {
      fprintf(fout,"\t---Fragment %d Intrafragment Coordinates---\n", i+1);
      fragments[i]->print_intcos(fout, g_atom_offset(i));
    }
    for (int i=0; i<interfragments.size(); ++i) {
      a = interfragments[i]->g_A_index();
      b = interfragments[i]->g_B_index();
      interfragments[i]->print_intcos(fout, g_atom_offset(a), g_atom_offset(b));
    }

    for (int i=0; i<efp_fragments.size(); ++i) {
      fprintf(fout,"\t---Fragment %d EFP fragment Coordinates---\n", i+1);
      efp_fragments[i]->print_intcos(fout);
    }
  }

  // print definition of an internal coordinate from global index 
  std::string get_intco_definition_from_global_index(int coord_index) const;

  void update_efp_values(void);

  void print_intco_dat(FILE *fp_intco);

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

  // compute intco values from frag member geometries
  double * intco_values(void) const {
    GeomType x = g_geom_2D();
    double *q = intco_values(x);
    return q;
  }

  // compute intco values from given geometry ; empty space for EFP coordinates included
  double * intco_values(GeomType new_geom) const {
    double *q, *q_frag, *q_IF;
    q = init_array(g_nintco());

    for (int f=0; f<fragments.size(); ++f) {
      q_frag = fragments[f]->intco_values( &(new_geom[g_atom_offset(f)]) );

      for (int i=0; i<fragments[f]->g_nintco(); ++i)
        q[g_intco_offset(f)+i]  = q_frag[i];

      free_array(q_frag);
    }

    for (int I=0; I<interfragments.size(); ++I) {
      int A_index = interfragments[I]->g_A_index();
      int B_index = interfragments[I]->g_B_index();
      
      q_IF = interfragments[I]->intco_values( &(new_geom[g_atom_offset(A_index)]),
        &(new_geom[g_atom_offset(B_index)]) );

      for (int i=0; i<interfragments[I]->g_nintco(); ++i)
        q[g_interfragment_intco_offset(I)+i]  = q_IF[i];

      free_array(q_IF);
    }

    return q;
  }

  void write_geom(void);
  void symmetrize_geom(void);
  void print_geom(void);

  double ** compute_B(void) const;
  double ** compute_derivative_B(int intco_index) const ;

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

//****AVC****//
    for (int f=0; f<efp_fragments.size(); f++)
    {
      g_frag = efp_fragments[f]->get_geom_array();
      for (int i=0; i<3*3; i++)
        g[3*3*f+i] = g_frag[i];
    }
//****AVC****//

    return g;
  }

  double ** g_geom_2D(void) const {
    double **g_frag;
    double **g = init_matrix(g_natom()+3*g_nefp_fragment(),3);

    for (int f=0; f<fragments.size(); ++f) {
      g_frag = fragments[f]->g_geom();
      for (int i=0; i<fragments[f]->g_natom(); ++i)
        for (int xyz=0; xyz<3; ++xyz)
          g[g_atom_offset(f)+i][xyz] = g_frag[i][xyz];
      free_matrix(g_frag);
    }
//****AVC****//
double **g_efp_frag;
for(int f=0; f<efp_fragments.size(); f++)
{
  g_efp_frag = efp_fragments[f]->get_xyz_pointer();
  for(int i=0; i<3; i++)
  {
    for(int xyz=0; xyz<3; xyz++)
      g[g_atom_offset(fragments.size()-1)+3*f+i][xyz] = g_efp_frag[i][xyz];
printf("\n%d %d  %15.7f %15.7f %15.7f\n", f, i, g[3*f+i][0], g[3*f+i][1], g[3*f+i][2]);
  }
  free_matrix(g_efp_frag);
}


//****AVC****//

    return g;
  }

  double ** g_grad_2D(void) {
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
  void irc_step(void);
  void nr_step(void);
  void rfo_step(void);
  void prfo_step(void);
  void backstep(void);
  void sd_step(void);
  void linesearch_step(void);

  void apply_intrafragment_step_limit(double * & dq);
  void check_intrafragment_zero_angles(double const * const dq);

  void set_geom_array(double * array_in) {
printf("\nentering set_geom_array() -- fragments.size() is %d and efp_fragments.size() is %d\n", fragments.size(), efp_fragments.size());
printf("\n efp_fragments[%d]->get_libmints_geom_index() = %d \n", 0, efp_fragments[0]->get_libmints_geom_index());
printf("\n efp_fragments[%d]->get_libmints_geom_index() = %d \n", 1, efp_fragments[1]->get_libmints_geom_index());


for(int f=0; f<efp_fragments.size(); f++)
{
  double * xyz_array = &( array_in[3*efp_fragments[f]->get_libmints_geom_index()] );
  printf("\nfragment %d:\n", f);
  for(int i=0; i<3*3; i++)
    printf("\n %15.8f", xyz_array[i]);
}

    for (int f=0; f<fragments.size(); ++f)
      fragments[f]->set_geom_array( &(array_in[3*g_atom_offset(f)]) );
    for (int f=0; f<efp_fragments.size(); f++)
    {
      efp_fragments[f]->set_geom_array( &(array_in[3*efp_fragments[f]->get_libmints_geom_index()]) );
    }
printf("\nleaving set_geom_array()\n");
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

  bool read_intcos(std::ifstream & fin);
  // function to obtain geometry and gradient
  void read_geom_grad(void);

  // Compute constraint matrix
  double ** compute_constraints(void);

  void add_efp_fragments(void);

  // freeze interfragment modes that break symmetry
  void freeze_interfragment_asymm(void);

  // determine whether a linear combination of intcos breaks symmetry
  bool intco_combo_is_symmetric(double *intco_combo, int dim);

  // Apply string list of user-specified internals to be frozen.
  bool apply_input_constraints();

};

}

#endif

