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

    void compute_val(double **geom);
    double * g_s(double **geom, const int iatom) const;
    void print(const FILE *fp) const;
    void print_s(const FILE *fp, double **geom) const;
    void print_disp(const FILE *fp, const double old_q, const double f_q,
      const double dq, const double new_q) const;
    bool operator==(const SIMPLE & s2) const;

};

}

#endif

