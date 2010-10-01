#ifndef _opt_molecule_h_
#define _opt_molecule_h_

#include "frag.h"

#include <fstream>

#define FILENAME_GEOM_GRAD_IN "psi.file11.dat"
#define FILENAME_GEOM_OUT "psi.geom.dat"

namespace opt {

class MOLECULE {

  vector<FRAG *> fragments ;
  double energy;

 public:

  MOLECULE(std::ifstream & fp);   // allocate and read in geometry

  ~MOLECULE() {
    //printf("Destructing molecule\n");
    for (int i=0; i<fragments.size(); ++i)
      delete fragments[i];
    fragments.clear();
  }

  int g_nfragment(void) const { return fragments.size(); };

  int g_natom(void) const {
    int n = 0;
    for (int f=0; f<fragments.size(); ++f)
      n += fragments[f]->g_natom();
    return n;
  }

  int g_nintco(void) const {
    int n=0;
    for (int f=0; f<fragments.size(); ++f)
      n += fragments[f]->g_nintco();
    return n;
  }

  int g_atom_offset(int i) const {
    int n = 0;
    for (int f=1; f<=i; ++f)
      n += fragments[f-1]->g_natom();
    return n;
  }

  int g_intco_offset(int i) const {
    int n = 0;
    for (int f=1; f<=i; ++i)
      n += fragments[i-1]->g_nintco();
    return n;
  }

  double g_energy(void) const { return energy; }

  void update_connectivity_by_distances(double scale_radii = -1) {
    for (int i=0; i<fragments.size(); ++i)
      fragments[i]->update_connectivity_by_distances(scale_radii);
  }

  void print_connectivity(FILE *fout) const {
    for (int i=0; i<fragments.size(); ++i)
      fragments[i]->print_connectivity(fout, i);
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
    for (int i=0; i<fragments.size(); ++i)
      fragments[i]->print_intcos(fout, i);
  }

  void print_intco_dat(FILE *fout) {
    for (int i=0; i<fragments.size(); ++i)
      fragments[i]->print_intco_dat(fout, i);
  }

  int add_simples_by_connectivity(void) {
    int n=0;
    for (int i=0; i<fragments.size(); ++i)
      n += fragments[i]->add_simples_by_connectivity();
    return n;
  }

  // compute intco values from frag member geometries
  double * intco_values(void) const {
    double *q, *q_frag;
    q = init_array(g_nintco());

    for (int f=0; f<fragments.size(); ++f) {
      q_frag = fragments[f]->intco_values();
      for (int i=0; i<fragments[f]->g_nintco(); ++i)
        q[g_intco_offset(f)+i]  = q_frag[i];
      free_array(q_frag);
    }
    return q;
  }

  // compute intco values from given geometry
  double * intco_values(GeomType new_geom) const {
    double *q, *q_frag;
    q = init_array(g_nintco());

    for (int f=0; f<fragments.size(); ++f) {
      q_frag = fragments[f]->intco_values( &(new_geom[g_atom_offset(f)]) );

      for (int i=0; i<fragments[f]->g_nintco(); ++i)
        q[g_intco_offset(f)+i]  = q_frag[i];

      free_array(q_frag);
    }
    return q;
  }

  double * g_grad_array(void) {
    int f, i;
    double *g, *g_frag;

    g = init_array(3*g_natom());
    for (f=0; f<fragments.size(); ++f) {
      g_frag = fragments[f]->g_grad_array();
      for (i=0; i<3*fragments[f]->g_natom(); ++i)
        g[g_atom_offset(f)+i] = g_frag[i];
      free_array(g_frag);
    }
    return g;
  }

  double * g_geom_array(void) {
    double *g, *g_frag;

    g = init_array(3*g_natom());
    for (int f=0; f<fragments.size(); ++f) {
      g_frag = fragments[f]->g_geom_array();
      for (int i=0; i<3*fragments[f]->g_natom(); ++i)
        g[g_atom_offset(f)+i] = g_frag[i];
      free_array(g_frag);
    }
    return g;
  }

  double ** g_geom_2D(void) {
    double **g, *g_frag;

    g = init_matrix(g_natom(),3);
    for (int f=0; f<fragments.size(); ++f) {
      g_frag = fragments[f]->g_geom_array();
      int cnt=0;
      for (int i=0; i<fragments[f]->g_natom(); ++i)
        for (int xyz=0; xyz<3; ++xyz)
          g[g_atom_offset(f)+i][xyz] = g_frag[cnt++];
      free_array(g_frag);
    }
    return g;
  }

  void write_geom_chkpt(void);
  void write_geom_to_active_molecule();

  double ** compute_B(void) const {
    double **B, **B_frag;
    int f, i, j;

    B = init_matrix(g_nintco(), 3*g_natom());
    for (f=0; f<fragments.size(); ++f) {
      B_frag = fragments[f]->compute_B();
      for (i=0; i<fragments[f]->g_nintco(); ++i)
        for (j=0; j<3*fragments[f]->g_natom(); ++j)
          B[i+g_intco_offset(f)][j+g_atom_offset(f)] = B_frag[i][j];
      free_matrix(B_frag);
    }
    return B;
  }

  void H_guess(void) const;
  void forces(void);
  void project_f_and_H(void);
  void nr_step(void);
  void rfo_step(void);

  void apply_intrafragment_step_limit(double * & dq);
  void check_intrafragment_zero_angles(double const * const dq);

  void write_geom(void);

  void set_geom_array(double * array_in) {
    for (int f=0; f<fragments.size(); ++f)
      fragments[f]->set_geom_array( &(array_in[3*g_atom_offset(f)]) );
  }

  void fix_tors_near_180(void) {
    for (int f=0; f<fragments.size(); ++f)
      fragments[f]->fix_tors_near_180();
  }

  void test_B(void);
  void test_derivative_B(void);

/*
  // compute and print B matrix (for debugging)
  void print_B(FILE *fp) const ;

  // check nearness to 180 and save value
  //print s vectors to output file
  void print_s(const FILE *fp, double **geom) const {
    fprintf(const_cast<FILE *>(fp),"\t---S vectors for internals---\n");
    for(int i=0; i<intcos.size(); ++i)
      intcos.at(i)->print_s(fp, geom);
  }
  void displace(double *dq, bool print_disp = false);
  double ** H_guess(void);
*/

};

}

#endif
