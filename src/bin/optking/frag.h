#ifndef _opt_frag_h_
#define _opt_frag_h_

#include <cstdio>
#include <vector>
#include "intcos.h"

namespace opt {

using std::vector;

class FRAG {
  int natom;      // number of atoms in fragment
  double *Z;      // atomic numbers
  double **geom;   // cartesian coordinates
  double **grad;   // cartesian coordinates
  double *mass;   // nuclear masses
  bool   **connectivity; // connectivity matrix
  vector<SIMPLE *> intcos;

 public:

  FRAG(int natom_in, double *Z_in, double **geom_in); // use if Z and geom are known

  ~FRAG(); // free memory

  int g_natom(void) const { return natom; };

  void set_default_masses(void);

  void print_geom(FILE *fp, const int id, bool print_mass = false);
  void print_geom_grad(FILE *fp, const int id, bool print_mass = false);
  void print_intcos(FILE *fp, const int id);
  void print_intco_dat(FILE *fp, const int id);

  void update_connectivity_by_distances(double scale_radii = -1);

  void print_connectivity(FILE *fout, const int id) const ;

  // add simple internals based on connectivity; return number added
  int add_stre_by_connectivity(void);
  int add_bend_by_connectivity(void);
  int add_tors_by_connectivity(void);

  int add_simples_by_connectivity(void) {
    int n;
    n  = add_stre_by_connectivity();
    n += add_bend_by_connectivity();
    n += add_tors_by_connectivity();
    return n;
  }

  // compute B matrix (intcos.size x 3*natom)
  double ** compute_B(void) const ;
  void compute_B(double **) const ; // use prevously allocated memory

  // compute and print B matrix (for debugging)
  void print_B(FILE *fp) const ;

  // check nearness to 180 and save value
  void fix_tors_near_180(void);

  // return number of intrafragment coordinates
  int g_nintco(void) const { return intcos.size(); };

  // don't let angles pass through 0
  void check_zero_angles(double const * const dq);

  //print s vectors to output file
  void print_s(FILE *fp, const double ** const geom) const {
    fprintf(fp,"\t---S vectors for internals---\n");
    for(int i=0; i<intcos.size(); ++i)
      intcos.at(i)->print_s(fp, geom);
  }

  // compute and return values of internal coordinates from member geometry
  double * intco_values(void) const ;

  // compute and return values of internal coordinates from given geometry
  double * intco_values(GeomType new_geom) const ;

  // is simple one already present?
  bool present(const SIMPLE *one) const;

  int find(const SIMPLE *one) const;

  void displace(double *dq, bool print_disp = false);

  double * g_geom_array(void);
  double * g_grad_array(void);
  void set_geom_array(double * geom_array_in);
  void set_grad(double **grad_in);

  void write_geom(FILE *fp_geom); // write cartesian geometry out for next step

  double ** H_guess(void); 

};

}

#endif
