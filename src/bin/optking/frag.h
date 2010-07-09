/*! \file
    \ingroup OPTKING
    \brief frag_class describes a set of interfragment coordinates connecting two fragments
*/

#ifndef _psi3_bin_optking_fragment_h_
#define _psi3_bin_optking_fragment_h_

#include <libqt/qt.h>
#include <libciomr/libciomr.h>

namespace psi { //namespace optking {

/*
Fragment_class describes a SET of interfragment coordinates connecting two fragments
(A and B).  Each member of the set shares the definition of the fragments:
(A_natom, *A_atom, *A_weight, B_natom, *B_atom, *B_weight) and id number.
Each member gets its own s-vectors.
The following index order of vals and s-vectors (A_s and B_s) is assumed:
RAB, theta_A, theta_B, tau, phi_A, phi_B.
*/

class frag_class {

  private:

    bool coord_on[6]; /* indicates on/off for which coordinates [1-6] are present */
    int id;        /* unique id number for the set */
    int A_natom;   /* number of atoms in fragment A */
    int B_natom;   /* number of atoms in fragment B */
    int *A_atom;   /* list of atoms in fragment A */
    int *B_atom;   /* list of atoms in fragment B */
    /* P=1 (for an atom); P=2 (for a linear fragment); P=3 (for non-linear fragments) */
    int A_P;       /* the # of reference points for fragment A to worry about */
    int B_P;       /* the # of reference points for fragment B to worry about; A_P >= B_P */
    double **A_weight; /* weights fix reference points via a linear combination [3][A_natom] */
    double **B_weight; /* weights fix reference points via a linear combination [3][B_natom] */
    double val[6];  /* D = val of coordinate [6] */
    double **A_s; /* S vectors fragment A [6][A_P*3] */
    double **B_s; /* S vectors fragment B [6][B_P*3] */
    int near_180[6]; /* [6] ; +1 => val>FIX_NEAR_180, -1 => val<FIX_NEAR_180, 0 otherwise */

  public:

    friend class simples_class;

    frag_class(int id_in, int A_natom_in, int B_natom_in, int A_P_in, int B_P_in) { 
      id = id_in;
      A_natom = A_natom_in;
      B_natom = B_natom_in;
      A_P = A_P_in;
      B_P = B_P_in;
      A_atom   = new int[A_natom];
      B_atom   = new int[B_natom];
      A_weight = block_matrix(3,A_natom);
      B_weight = block_matrix(3,B_natom);
      A_s = block_matrix(6,A_P*3);
      B_s = block_matrix(6,B_P*3);
    }

    ~frag_class() {
      delete [] A_atom;
      delete [] B_atom;
      free_block(A_weight);
      free_block(B_weight);
      free_block(A_s);
      free_block(B_s);
    }

    /* functions in frag.cc */
    void print(FILE *fp_out, bool print_vals, bool print_weights) const;
    void print_s(FILE *fp_out) const;
    void compute(double *geom);
    void compute_s(double *geom);

    // functions to set member variables
    void set_id(int i){ id = i;}
    void set_A_natom(int i) { A_natom = i;}
    void set_B_natom(int i) { B_natom = i;}
    void set_A_P(int i) { A_P = i;}
    void set_B_P(int i) { B_P = i;}

    void set_A_atom(int i, int j) {
      if (i >= A_natom) throw("fragment.set_A_atom() : index is too large");
      A_atom[i] = j;
    }
    void set_B_atom(int i, int j) {
      if (i >= B_natom) throw("fragment.set_B_atom() : index is too large");
      B_atom[i] = j;
    }

    void set_A_weight(int ref_atom, int frag_index, double new_weight) {
      if (ref_atom >= A_P)
        throw("fragment.get_A_weight() : ref_atom is greater than A_P\n");
      else if (frag_index >= A_natom)
        throw("fragment.get_A_weight() : frag_index is greater than A_natom\n");
      A_weight[ref_atom][frag_index] = new_weight;
    }
    void set_B_weight(int ref_atom, int frag_index, double new_weight) {
      if (ref_atom >= B_P)
        throw("fragment.get_B_weight() : ref_atom is greater than B_P\n");
      else if (frag_index >= B_natom)
        throw("fragment.get_B_weight() : frag_index is greater than B_natom\n");
      B_weight[ref_atom][frag_index] = new_weight;
    }

    void set_near_180(int I, int new_val) {
      if (I<0 || I>5) throw("fragment.set_near_180() : expecting id between 0 and 5");
      if ((new_val != -1) && (new_val != 1) && (new_val != 0))
        throw("fragment.set_near_180() : new val not understood");
      near_180[I] = new_val;
    }

    void set_coord_on(int I, bool on_or_off) {
      if (I<0 || I>5) throw("fragment.set_coord_on() : expecting id between 0 and 5");
      coord_on[I] = on_or_off;
    }

    void set_val(int I, double new_val) {
      if (I<0 || I>5) throw("fragment.set_val() : expecting id between 0 and 5");
      val[I] = new_val;
    }

    // functions to get member variables
    int get_id(void) const { return id;}
    int get_A_natom(void) const { return A_natom;}
    int get_B_natom(void) const { return B_natom;}
    int get_A_P(void) const { return A_P;}
    int get_B_P(void) const { return B_P;}

    int get_A_atom(int i) const {
      if (i < 0 || i >= A_natom) throw("fragment.get_A_atom() : bad atom index");
      return A_atom[i];
    }
    int get_B_atom(int i) const {
      if (i < 0 || i >= B_natom) throw("fragment.get_B_atom() : bad atom index");
      return B_atom[i];
    }

    int get_atom(Frag_switch X, int a) const  {
      int atom;
      if (X == FRAG_A)
        atom = get_A_atom(a);
      else if (X == FRAG_B)
        atom = get_B_atom(a);
      return atom;
    }

    double get_s(Frag_switch X, int sub_index2, int atom, int xyz) const  {
      int K;
      double tval = 0.0;
      if ( xyz < 0 || xyz > 2) throw ("frag_class::get_s() : xyz must be 0, 1 or 2");
      if ( sub_index2 < 0 || sub_index2 > 5) throw ("frag_class::get_s() : bad sub_index2 value, must be {0,5}");
 
      // A_s and B_s contain s vectors on reference atoms
      // compute correlating s vectors on fragment atoms

      if (X == FRAG_A) {
        if ( atom < 0 || atom > A_natom) throw ("frag_class::get_s() : bad atom number");
        for (K=0; K<A_P; ++K)             // loop over reference atoms of A
          tval += A_weight[K][atom] * A_s[sub_index2][3*K+xyz];
      }
      else if (X == FRAG_B) {
        if ( atom < 0 || atom > B_natom) throw ("frag_class::get_s() : bad atom number");
        for (K=0; K<B_P; ++K)             // loop over reference atoms of A
          tval += B_weight[K][atom] * B_s[sub_index2][3*K+xyz];
      }
      return tval;
    }

    /* ref_atom = 0-2, frag_index 0-3*A_natom */
    double get_A_weight(int ref_atom, int frag_index) const {
      if (ref_atom >= A_P)
        throw("fragment.get_A_weight() : ref_atom is greater than A_P\n");
      else if (frag_index >= A_natom)
        throw("fragment.get_A_weight() : frag_index is greater than A_natom\n");
      return A_weight[ref_atom][frag_index];
    }
    double get_B_weight(int ref_atom, int frag_index) const {
      if (ref_atom >= B_P)
        throw("fragment.get_B_weight() : ref_atom is greater than B_P\n");
      else if (frag_index >= B_natom)
        throw("fragment.get_B_weight() : frag_index is greater than B_natom\n");
      return B_weight[ref_atom][frag_index];
    }

    int get_near_180(int I) const {
      if (I<0 || I>5) throw("fragment.get_near_180() : expecting id between 0 and 5");
      return near_180[I];
    }

    bool get_coord_on(int I) const {
      if (I<0 || I>5) throw("fragment.get_coord_on() : expecting id between 0 and 5");
      return coord_on[I];
    }

    double get_val(int I) const {
      if (I<0 || I>5) throw("fragment.get_val() : expecting id between 0 and 5");
      return val[I];
    }
    // returns bond lengths in normal angstroms but angles in degrees
    double get_val_A_or_rad(int I) const {
      if (I==0)
        return get_val(0);
      else 
        return (get_val(I) * _pi / 180.0);
    }

    double get_A_s(int I, int ref_atom_xyz) const {
      if (I<0 || I>5)
        throw("fragment.get_A_s() : expecting id between 0 and 5");
      else if (ref_atom_xyz >= 3*A_P)
        throw("fragment.get_A_s() : ref_atom_xyz >= 3*A_natom\n");
      return A_s[I][ref_atom_xyz];
    }
    double get_B_s(int I, int ref_atom_xyz) const {
      if (I<0 || I>5)
        throw("fragment.get_B_s() : expecting id between 0 and 5");
      else if (ref_atom_xyz >= 3*B_P)
        throw("fragment.get_B_s() : ref_atom_xyz >= 3*B_natom\n");
      return B_s[I][ref_atom_xyz];
    }

    int  get_dim(void) const {
      int I, dim=0;
      for (I=0; I<6; ++I) 
        if (coord_on[I]) ++dim;
      return dim;
    }

    // fix torsional angles by recording nearness (or not) to +180 or -180
    void fix_near_180(void) {
      int I;
      for (I=3; I<6; ++I) {
        if (val[I] > FIX_NEAR_180)
          near_180[I] = +1;
        else if (val[I] < -1*FIX_NEAR_180)
          near_180[I] = -1;
        else 
          near_180[I] = 0;
      }
    }

    bool operator==(const frag_class & s2) const {
      int i, k;

      if ( (this->A_natom != s2.A_natom) || (this->B_natom != s2.B_natom) )
        return false;

      for (i=0; i<A_natom; ++i) {
        if ( this->A_atom[i] != s2.A_atom[i] ) return false;
      }

      for (i=0; i<B_natom; ++i) {
        if ( this->B_atom[i] != s2.B_atom[i] ) return false;
      }

      for (k=0; k<A_P; ++k) {
        for (i=0; i<A_natom; ++i)
          if ( this->A_weight[k][i] != s2.A_weight[k][i] ) return false;
      }

      for (k=0; k<B_P; ++k) {
        for (i=0; i<B_natom; ++i)
          if ( this->B_weight[k][i] != s2.B_weight[k][i] ) return false;
      }

      return true;
    };
};

}//} /* namespace psi::optking */

#endif
