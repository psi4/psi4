/*!
   \ingroup optking
   \file frag.h : header for fragment class
*/
#ifndef _opt_frag_h_
#define _opt_frag_h_

#include <cstdio>
#include <vector>
#include <string>

#include "intcos.h"

namespace opt {

class INTERFRAG;

using std::vector;

/*!
  \ingroup optking
  \class   FRAG
  \brief   A group of atoms, its geometry, and its internal coordinates.
  */
class FRAG {

 protected: // private, except the efp derived class can access
  int natom;      //< number of atoms in fragment
  double *Z;      //< atomic numbers
  double **geom;  //< cartesian coordinates
  double **grad;  //< cartesian coordinates
  double *mass;   //< nuclear masses
  bool   **connectivity; //< connectivity matrix
  vector<SIMPLE *> intcos;
  bool frozen; //< whether to optimize

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

  void print_geom(FILE *fp, const int id, bool print_mass = false);
  void print_geom_grad(FILE *fp, const int id, bool print_mass = false);
  void print_intcos(FILE *fp, int atom_offset=0);
  void print_intco_dat(FILE *fp, int atom_offset=0);

  void update_connectivity_by_distances(void);
  void update_connectivity_by_bonds(void);

  void print_connectivity(FILE *fout, const int id, const int offset = 0) const ;

  // add simple internals based on connectivity; return number added
  int add_stre_by_connectivity(void);
  int add_bend_by_connectivity(void);
  int add_tors_by_connectivity(void);

  int add_simples_by_connectivity(void) {
    int n;
    n  = add_stre_by_connectivity();
    n += add_bend_by_connectivity();
    n += add_tors_by_connectivity(); // but check bond angles
    return n;
  }

  // add connectivity between two atoms; atom numbers are for fragment
  void connect(int i, int j) {
    connectivity[i][j] = connectivity[j][i] = true;
  }

  // compute B matrix (intcos.size x 3*natom)
  double ** compute_B(void) const ;
  void compute_B(double **) const ; // use prevously allocated memory
  void compute_G(double **, bool use_masses=false) const;

  // compute B' matrix for one internal coordinate computed with member geometry
  double ** compute_derivative_B(int intco_index) const;

  // compute B' matrix for one internal coordinate computed with given geometry
  double ** compute_derivative_B(int intco_index, GeomType new_geom) const;

  // compute and print B matrix (for debugging)
  void print_B(FILE *fp) const ;

  // check nearness to 180 and save value
  void fix_tors_near_180(void);

  // check if interior angles of torsion are near 0 or linear
  //bool check_tors_for_bad_angles(void) const;

  // return number of intrafragment coordinates
  int g_nintco(void) const { return intcos.size(); };

  // return natom in definition of intco # intco_index
  int g_intco_natom(const int intco_index) const {
    return intcos.at(intco_index)->g_natom();
  }

  // return atom i in definition of intco # intco_index
  int g_intco_atom(const int intco_index, const int atom) const {
    return intcos.at(intco_index)->g_atom(atom);
  }

  // return array of atomic numbers
  double *g_Z(void) const;
  double *g_Z_pointer(void) { return Z; }

  // don't let angles pass through 0
  void check_zero_angles(double const * const dq);

  //print s vectors to output file
  void print_s(FILE *fp, const double ** const geom) const {
    fprintf(fp,"\t---S vectors for internals---\n");
    for(int i=0; i<intcos.size(); ++i)
      intcos.at(i)->print_s(fp, geom);
    fflush(fp);
  }

  // compute and return values of internal coordinates from member geometry
  double * intco_values(void) const ;

  // compute and return values of internal coordinates from given geometry
  double * intco_values(GeomType new_geom) const ;

  // is simple one already present?
  bool present(const SIMPLE *one) const;

  int find(const SIMPLE *one) const;

  void displace(double *dq, bool print_disp = false, int atom_offset=0);

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

  void print_geom(FILE *fp_geom); // write cartesian geometry out for next step

  double ** H_guess(void);
  bool **g_connectivity(void) const;
  const bool * const * g_connectivity_pointer(void) const;

  bool read_intco(std::vector<std::string> & tokens, int first_atom_offset);

  // return matrix of constraints on coordinates
  double ** compute_constraints(void) const;

  // add any missing hydrogen bonds within the fragment
  // return number added
  int add_hbonds(void);

  void intco_add(SIMPLE * i) {
    intcos.push_back(i);
  }

  double g_mass(int i) { return mass[i]; }

  // relating to frozen fragments
  bool is_frozen(void) const { return frozen; }
  void freeze(void)   { frozen = true; }
  void unfreeze(void) { frozen = false; }

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

};

}

#endif

