/*! \file tors.h
    \ingroup OPT10
    \brief TORS class declaration
*/

#ifndef _opt_tors_h_
#define _opt_tors_h_

#include "simple.h"

namespace opt {

class TORS : public SIMPLE {

  private:
    int near_180;
        // +1 if positive and approaching 180
        // -1 if negative and approaching -180
        //  0 otherwise

  public:

    TORS(int A_in, int B_in, int C_in, int D_in);

    ~TORS() { } // also calls ~SIMPLE

    void compute_val(double **geom);
    double * g_s(double **geom, const int iatom) const;
    void print(const FILE *fp) const;
    void print_s(const FILE *fp, double **geom) const;
    void print_disp(const FILE *fp, const double old_q, const double f_q, 
      const double dq, const double new_q) const;
    bool operator==(const SIMPLE & s2) const;

    void fix_near_180(void);
};

}

#endif

