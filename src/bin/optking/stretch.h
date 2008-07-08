/*! \file stretch.h
    \ingroup OPT08
    \brief Stretch class declaration
*/

#ifndef _psi3_src_bin_optking_stretch_h_
#define _psi3_src_bin_optking_stretch_h_

#include <cstdio>
extern "C" { extern FILE *outfile;}

#include "v3d.h"
using namespace psi::v3d;

namespace psi { namespace optking {

class Stretch {
    int   id;    // identifier
    int    A;    // atom A index
    int    B;    // atom B index
    double val;  // value of distance from A to B
    // The "s vector for atom A" is defined as d(val)/dR_A.
    // It points in the direction of maximum increase in the coordinate.
    double sA[3]; // s vector for A
    double sB[3]; // s vector for B
    static const int natom_coord = 2; // number of atoms involved in coordinate

  public:
    Stretch(int id_in = -1, int A_in = -1, int B_in = -1);
    ~Stretch() { printf("destructing stretch\n"); }

    int get_natom(void) const { return natom_coord; }
    double get_s(int atom_index, int xyz) const;
    int get_atom(int atom_index) const;
    int  get_A(void) const { return A;}
    int  get_B(void) const { return B;}

    // if print_flag == 0, print intco formatted definition of coordinate
    // if print_flag == 1, print format for output with value
    void print(int print_flag) const;
    void print_s(void) const;
  
    void set(int id, int A,int B); // later make id automatic
    void compute(double **geom);
    void compute_s(double **geom);

    bool operator==(const Stretch & s2) const;

};


}}

#endif

