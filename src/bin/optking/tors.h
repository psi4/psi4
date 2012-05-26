/*! \file tors.h
    \ingroup optking
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

    TORS(int A_in, int B_in, int C_in, int D_in, bool freeze_in=false);

    ~TORS() { } // also calls ~SIMPLE

    double value(GeomType geom) const;
    //bool fix_tors_value_corrected(GeomType geom) const;

    // compute and return array of first derivative (B marix elements)
    double ** DqDx(GeomType geom) const;

    // compute and return array of second derivative (B' matrix elements)
    double ** Dq2Dx2(GeomType geom) const;

    void print(FILE *fp, GeomType geom, int atom_offset=0) const;
    void print_intco_dat(FILE *fp, int atom_offset=0) const;
    void print_s(FILE *fp, GeomType geom) const;
    void print_disp(FILE *fp, const double old_q, const double f_q, 
      const double dq, const double new_q, int atom_offset=0) const;
    bool operator==(const SIMPLE & s2) const;
    std::string get_definition_string(int atom_offset=0) const;

    void fix_tors_near_180(GeomType geom);

    //bool check_tors_for_bad_angles(GeomType geom) const;

};

}

#endif

