/*! \file bend.h
    \ingroup OPT09
    \brief STRETCH class declaration
*/

#ifndef _psi3_src_bin_opt09_bend_h_
#define _psi3_src_bin_opt09_bend_h_

#include <cstdio>
extern "C" { extern FILE *outfile;}

#include "simple.h"
#include "v3d.h"

namespace psi { namespace opt09 {

using namespace psi::v3d;

class BEND : public SIMPLE {
    int    A;    // atom A index
    int    B;    // atom B index
    int    C;    // atom C index
    // The "s vector for atom A" is defined as d(val)/dR_A.
    // It points in the direction of maximum increase in the coordinate.
    double sA[3]; // s vector for A
    double sB[3]; // s vector for B
    double sC[3]; // s vector for B

  public:

    BEND(int id_in, int A_in, int B_in, int C_in);

    // ~BEND() { fprintf(outfile,"destructing bend\n"); }

    double get_val(void) const {return val;}
    int get_id(void) const {return id;}

    //double get_s(int atom_index, int xyz) const;
    //int get_atom(int atom_index) const;
    int  get_A(void) const { return A;}
    int  get_B(void) const { return B;}
    int  get_C(void) const { return C;}
    int get_atom(int a) const {
      if (a == 0) return A;
      else if (a == 1) return B;
      else if (a == 2) return C;
      else throw("BEND::get_atom() - only 3 atoms in stretch");
    }

    void compute(double **geom);
    void compute_s(double **geom);

    void print(int print_flag) const;
    void print_s(void) const;
    double get_s(int atom_index, int xyz) const;
  
    //void set(int id, int A,int B); // later make id automatic
    bool operator==(const BEND & s2) const;
};

}}
#endif

