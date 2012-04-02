/*! \file bend.h
    \ingroup OPTKING
    \brief BEND class declaration
*/

#ifndef _opt_bend_h_
#define _opt_bend_h_

#include "simple.h"

namespace opt {

class BEND : public SIMPLE {

    // if true, this bend is a secondary complement to another linear bend
    bool linear_bend;

  public:

    BEND(int A_in, int B_in, int C_in, bool freeze_in=false);

    ~BEND() { } // also calls ~SIMPLE()

    double value(GeomType geom) const;

    // compute and return array of first derivative (B marix elements)
    double **DqDx(GeomType geom) const;

    // compute and return array of second derivative (B' matrix elements)
    double **Dq2Dx2(GeomType geom) const;

    void print(FILE *fp, GeomType geom, int atom_offset=0) const;
    void print_intco_dat(FILE *fp, int atom_offset=0) const;
    void print_s(FILE *fp, GeomType geom) const;
    void print_disp(FILE *fp, const double old_q, const double f_q,
      const double dq, const double new_q, int atom_offset=0) const;
    bool operator==(const SIMPLE & s2) const;
    std::string get_definition_string(int atom_offset=0) const;

    void make_linear_bend(void) { linear_bend = true; }
    bool is_linear_bend(void) const { return linear_bend; }

};

}

#endif

