/*! \file stretch.h
    \ingroup OPT09
    \brief STRETCH class declaration
*/

#ifndef _psi4_src_bin_opt09_stretch_h_
#define _psi4_src_bin_opt09_stretch_h_

#include <cstdio>
#include <psi4-dec.h>

#include "simple.h"
#include "v3d.h"

using namespace psi::v3d;

namespace psi {

class STRETCH : public SIMPLE {
    int    A;    // atom A index
    int    B;    // atom B index
    // The "s vector for atom A" is defined as d(val)/dR_A.
    // It points in the direction of maximum increase in the coordinate.
    double sA[3]; // s vector for A
    double sB[3]; // s vector for B

  public:

    STRETCH(int id_in, int A_in, int B_in);

    // ~STRETCH() { fprintf(outfile,"destructing stretch\n"); }

    double get_val(void) const {return val;}
    int get_id(void) const {return id;}

    int  get_A(void) const { return A;}
    int  get_B(void) const { return B;}
    int get_atom(int a) const {
      if (a == 0) return A;
      else if (a == 1) return B;
      else throw("STRETCH::get_atom() - only 2 atoms in stretch");
    }
    //void set(int id, int A,int B); // later make id automatic

    void compute(double **geom);
    void compute_s(double **geom);
    void print(int print_flag) const;
    void print_s(void) const;
    double get_s(int atom_index, int xyz) const;

    bool operator==(const STRETCH & s2) const;
};

}

#endif
