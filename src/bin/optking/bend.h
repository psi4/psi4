/*! \file bend.h
    \ingroup OPT10
    \brief BEND class declaration
*/

#ifndef _opt_bend_h_
#define _opt_bend_h_

#include "simple.h"

namespace opt {

class BEND : public SIMPLE {

  public:

    BEND(int A_in, int B_in, int C_in);

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
      const double dq, const double new_q) const;
    bool operator==(const SIMPLE & s2) const;

};

}

#endif

