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

    //STRE(const STRE & from) : SIMPLE (from)  { fprintf(stdout,"stre copy constructor\n"); } 

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

