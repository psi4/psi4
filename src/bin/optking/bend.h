/*! \file stretch.h
    \ingroup OPT08
    \brief Bend class declaration
*/

#ifndef _psi3_src_bin_optking_bend_h_
#define _psi3_src_bin_optking_bend_h_

#include <cstdio>
extern "C" { extern FILE *outfile;}
#include <cmath>
#include "physconst.h"

#include "v3d.h"
using namespace psi::v3d;

namespace psi { namespace optking {

class Bend {
    int   id;   
    int    A;    // atom A
    int    B;    // atom B
    int    C;    // atom C
    double val;  // value of angle A-B-C
    // The "s vector for atom A" is defined as d(val)/dR_A.
    // It points in the direction of maximum increase in the coordinate.
    double sA[3]; // s vector for A
    double sB[3]; // s vector for B
    double sC[3]; // s vector for C
    static const int natom_coord = 3; // number of atoms involved in coordinate

  public:
    Bend(int id_in = -1, int A_in = -1, int B_in = -1, int C_in = -1);
    ~Bend() { printf("destructing bend\n"); }

    int get_natom(void) const {return natom_coord;}
    double get_s(int atom_index, int xyz) const;
    int get_atom(int atom_index) const;

    // if print_flag == 0, print intco formatted definition of coordinate
    // if print_flag == 1, print format for output with value
    void print(int print_flag) const; 
    void print_s(void) const;
  
    void set(int id, int A, int B, int C); // later make id automatic?
    void compute(double **geom);
    void compute_s(double **geom);

    bool operator==(const Bend & s2) const;
};

}}

#endif

