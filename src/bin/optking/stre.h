/*! \file stre.h
    \ingroup OPT10
    \brief STRE class declaration
*/

#ifndef _opt_stretch_h_
#define _opt_stretch_h_

#include "simple.h"

namespace opt {

class STRE : public SIMPLE {

  public:

    STRE(int A_in, int B_in);

    ~STRE() { } // also calls ~SIMPLE()

    double value(GeomType geom) const;

    // compute and return array of first derivative (B matrix elements)
    // returned matrix is [atom][x,y,z]
    double **DqDx(GeomType geom) const;

    // compute and return array of second derivative (B' matrix elements)
    // returned matrix is order 3N cart by 3N cart
    double **Dq2Dx2(GeomType geom) const;

    void print(FILE *fp, GeomType geom) const;
    void print_s(FILE *fp, GeomType geom) const;
    void print_intco_dat(FILE *fp) const;
    void print_disp(FILE *fp, const double old_q, const double f_q,
      const double dq, const double new_q) const;
    bool operator==(const SIMPLE & s2) const;

};

}

#endif

