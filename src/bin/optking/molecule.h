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
#if defined(OPTKING_PACKAGE_QCHEM)
  vector<EFP_FRAG *> efp_fragments;   // EFP fixed-body fragment - does not actually
                                      // store atoms (for now at least)
#endif
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
#if defined(OPTKING_PACKAGE_QCHEM)
    for (int i=0; i<efp_fragments.size(); ++i)
      delete efp_fragments[i];
    efp_fragments.clear();
#endif
  }

  // if you have one fragment - make sure all atoms are bonded
  // if not, break up into multiple fragments
  void fragmentize(void);

  void add_interfragment(void);

  int g_nfragment(void) const { return fragments.size(); };

#if defined(OPTKING_PACKAGE_QCHEM)
  int g_nefp_fragment(void) const { return efp_fragments.size(); };
#endif

  int g_natom(void) const { // excludes atoms in efp fragments
    int n = 0;
    for (int f=0; f<fragments.size(); ++f)
      n += fragments[f]->g_natom();
    return n;
  }

  int g_nintco(void) const {
    int n=0;
    for (int f=0; f<fragments.size(); ++f)
      n += fragments[f]->g_nintco();
    for (int i=0; i<interfragments.size(); ++i)
      n += interfragments[i]->g_nintco();
#if defined (OPTKING_PACKAGE_QCHEM)
    for (int e=0; e<efp_fragments.size(); ++e)
      n += efp_fragments[e]->g_nintco();
#endif
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

#if defined (OPTKING_PACKAGE_QCHEM)
  int g_nintco_efp_fragment(void) const {
    int n=0;
    for (int f=0; f<efp_fragments.size(); ++f)
      n += efp_fragments[f]->g_nintco();
    return n;
  }
#endif

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
#if defined (OPTKING_PACKAGE_QCHEM)
  int g_efp_fragment_intco_offset(int index) const {
    int n = g_nintco_intrafragment() + g_nintco_interfragment();
    for (int f=0; f<index; ++f)
      n += efp_fragments[f]->g_nintco();
    return n;
  }
#endif

  double g_energy(void) const { return energy; }

  void update_connectivity_by_distances(void) {
    for (int i=0; i<fragments.size(); ++i)
      fragments[i]->update_connectivity_by_distances();
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

#if defined (OPTKING_PACKAGE_QCHEM)
    for (int i=0; i<efp_fragments.size(); ++i) {
      fprintf(fout,"\t---Fragment %d EFP fragment Coordinates---\n", i+1);
      efp_fragments[i]->print_intcos(fout);
    }
#endif
  }

  // print definition of an internal coordinate from global index 
  std::string get_intco_definition_from_global_index(int coord_index) const;

#if defined (OPTKING_PACKAGE_QCHEM)
  void update_efp_values(void);
#endif

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
printf("adding %d auxiliary bonds\n", n);
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

  void forces(void);
  void project_f_and_H(void);
  void irc_step(void);
  void nr_step(void);
  void rfo_step(void);
  void prfo_step(void);
  void backstep(void);
  void sd_step(void);

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

  bool cartesian_H_to_internals(void) const;

  void set_masses(void) {
    for (int f=0; f<fragments.size(); ++f)
      fragments[f]->set_masses();
  }

  bool read_intcos(std::ifstream & fin);
  // function to obtain geometry and gradient
  void read_geom_grad(void);

  // tell whether internal coordinate is frozen
  double ** compute_constraints(void);

  void add_efp_fragments(void);

  // freeze interfragment modes that break symmetry
  void freeze_interfragment_asymm(void);

};

}

#endif

