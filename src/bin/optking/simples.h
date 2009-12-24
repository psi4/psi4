/*! \file
    \ingroup OPTKING
    \brief Class declaration for simple internal coordinate set
*/

#ifndef _psi3_bin_optking_simples_h_
#define _psi3_bin_optking_simples_h_

#include <libipv1/ip_lib.h>
#include <cov_radii.h>
#include <exception>

#include "stre.h"
#include "bend.h"
#include "tors.h"
#include "out.h"
#include "linb.h"
#include "frag.h"

#include <cctype>
#include <vector>

namespace psi { namespace optking {

class salc_set; // forward declaration for friend functions

using std::vector;

class simples_class {

  vector<stre_class> stre;
  vector<bend_class> bend;
  vector<tors_class> tors;
  vector<out_class> out;
  vector<linb_class> linb;
  vector<frag_class> frag;

  public:

   friend int opt_step(cartesians &, simples_class &, const salc_set &);

   // constructor in frag.cc
   // user_intcos = 1; read in simple coordinates from intco.dat
   //               0; generate simples from geometry
   simples_class(cartesians& carts, int user_intcos);

   simples_class() : stre(), bend(), tors(), out(), linb(), frag() { };

   ~simples_class() {
     stre.clear();
     bend.clear();
     tors.clear();
     out.clear();
     linb.clear();
     frag.clear();
   }

   // print_vals: false => print definitions to a file in intco.dat format
   //             true  => print intcos and their values
   // print_frag_weights: whether to print weights for fragments
   void print(FILE *fp_out, bool print_vals, bool print_frag_weights = false) const;

   // print s vectors
   void print_s(FILE *fp_out) const {
     int i;
     fprintf(fp_out,"S vectors for simple internal coordinates\n");
     for (i=0; i<stre.size(); ++i) stre.at(i).print_s(fp_out);
     for (i=0; i<bend.size(); ++i) bend.at(i).print_s(fp_out);
     for (i=0; i<tors.size(); ++i) tors.at(i).print_s(fp_out);
     for (i=0; i<out.size(); ++i)  out.at(i).print_s(fp_out);
     for (i=0; i<linb.size(); ++i) linb.at(i).print_s(fp_out);
     for (i=0; i<frag.size(); ++i) frag.at(i).print_s(fp_out);
     fprintf(fp_out,"\n");
   }

   // compute values of simple internal coordinates
   void compute(double *geom) {
     int i;
     for (i=0; i<stre.size(); ++i) stre.at(i).compute(geom);
     for (i=0; i<bend.size(); ++i) bend.at(i).compute(geom);
     for (i=0; i<tors.size(); ++i) tors.at(i).compute(geom);
     for (i=0; i<out.size(); ++i)  out.at(i).compute(geom);
     for (i=0; i<linb.size(); ++i) linb.at(i).compute(geom);
     for (i=0; i<frag.size(); ++i) frag.at(i).compute(geom);
     return;
   }

   // compute and store s sectors 
   void compute_s(double *geom) {
     int i;
     for (i=0; i<stre.size(); ++i) stre.at(i).compute_s(geom);
     for (i=0; i<bend.size(); ++i) bend.at(i).compute_s(geom);
     for (i=0; i<tors.size(); ++i) tors.at(i).compute_s(geom);
     for (i=0; i<out.size(); ++i)  out.at(i).compute_s(geom);
     for (i=0; i<linb.size(); ++i) linb.at(i).compute_s(geom);
     for (i=0; i<frag.size(); ++i) frag.at(i).compute_s(geom);
     return;
   }

   // get number of simple internal coordinates
   int get_num(void) const {
     int i, n;
     n = stre.size() + bend.size() + tors.size() + out.size() + linb.size();
     for (i=0; i<frag.size(); ++i)
       n += frag[i].get_dim();
     return n;
   }

   // get number of simple internal coordinates
   int get_num(Intco_type itype, int cnt_sets_as_one = 0) const {
     int i, tot=0;
     if (itype == STRE)      return stre.size();
     else if (itype == BEND) return bend.size();
     else if (itype == TORS) return tors.size();
     else if (itype == OUT)  return out.size();
     else if (itype == LINB) return linb.size();
     else if (itype == FRAG) {
       if (cnt_sets_as_one) {
         return frag.size();
       }
       else {
         for (i=0; i<frag.size(); ++i)
           tot += frag[i].get_dim();
         return tot;
       }
     }
   }

   void fix_near_180(void) {
     int i;
     for (i=0; i<tors.size(); ++i)
       tors[i].fix_near_180();
     for (i=0; i<frag.size(); ++i)
       frag[i].fix_near_180();
   }

   // ** following functions in frag.cc **

   // given id number, return type of coordinate, index within the type, and the
   // sub_index2 (I=1-6) for interfragment coordinates
   void locate_id(int id, Intco_type *itype, int *sub_index, int *sub_index2) const;

   // returns double ** bond connectivity matrix
   double **bond_connectivity_matrix(int natoms) const;

   // given the overall optking index (0-N), returns the user-assigned id number
   // optking index skips over interfragment coordinates turned "off"
   int index_to_id(int index) const;

   // given the id number, returns the overall optking index (0-N)
   int id_to_index(int id) const;

   // return id number that cooresponds to given simple type and atom numbers
   int get_id_from_atoms_stre(int a, int b) const;
   int get_id_from_atoms_bend(int a, int b, int c) const;
   int get_id_from_atoms_tors(int a, int b, int c, int d) const;
   int get_id_from_atoms_out(int a, int b, int c, int d, int *sign) const;
   int get_id_from_atoms_linb(int a, int b, int c, int linval) const;
   int get_id_from_atoms_frag(int a_natom, int b_natom, int *a_atom, int *b_atom) const;

   int get_id(Intco_type itype, int sub_index, int sub_index2 = 0) const;

   //int *atom2fragment(int natom) { return frag.atom2fragment(natom); } 

   bool is_unique (tors_class & t1) const {
     int i;
     bool unique = true;

     for (i=0; i<tors.size(); ++i)
       if (t1 == tors[i])
         unique = false;

     return unique;
   }

   // returns value in angstroms or degrees
   double get_val(Intco_type itype, int sub_index, int sub_index2=0) const;

   // returns value in angstroms or radians
   double get_val_A_or_rad(Intco_type itype, int sub_index, int sub_index2=0) const;

  // itype == STRE, BEND, etc = type of intco
  // subindex == place of coordinate within the set of coordinates of itype type
  // int atom == place of atom within definition (i.e., for STRE 0 or 1; for BEND 0, 1, or 2)
  // X (only applies to interfragment coodinates) == which fragment for which
  //  atom is requested (FRAG_A or FRAG_B); if FRAG_A, then atom belongs to {0, A_natom}
  // returns a = atom included in definition of internal (for which s-vector is non-zero)

   int get_atom(Intco_type itype, int sub_index, int atom, Frag_switch X = FRAG_A) const;

   // get s vector component; last 2 arguments only required for interfragment coordinates
   double get_s(Intco_type itype, int sub_index, int atom, int xyz, int sub_index2 = 0, Frag_switch X = FRAG_A) const;

   // get number of atoms in coordinate (e.g. 3 for bend), i.e. the number with nonzero s vectors
   int get_natom(Intco_type itype, int sub_index, Frag_switch X = FRAG_A) const;

   int linb_get_linval(int sub_index) const {
     return linb[sub_index].get_linval();
   }

   int frag_get_coord_on(int sub_index, int sub_index2) const {
     return frag[sub_index].get_coord_on(sub_index2);
   }

};

}} /* namespace psi::optking */

#endif
