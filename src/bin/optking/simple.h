/*! \file simple.h
    \ingroup opt
    \brief simple internal coordinate base class
*/

#ifndef _opt_simple_h_
#define _opt_simple_h_

#include <cstdio>
#include "mem.h"

namespace opt {

enum INTCO_TYPE {min_type, stre_type, bend_type, tors_type, max_type};

class SIMPLE {

  protected:

    INTCO_TYPE s_type; // type of simple
    int s_natom;       // # of atoms in internal definition, i.e., # of
                      // atoms with non-zero s-vectors
    int *s_atom;      // atoms indices in internal definition
    double s_val;     // value of coordinate

  public:

    SIMPLE (INTCO_TYPE s_type_in, int s_natom_in) {
      s_type = s_type_in;
      s_natom = s_natom_in;
      s_atom = init_int_array(s_natom);
    };

    virtual ~SIMPLE() { // derived class destructors called first; then this one
      free_int_array(s_atom);
    }

    INTCO_TYPE g_type(void) const { return s_type; }
    int     g_natom(void) const { return s_natom; }
    int     g_atom(int a) const { return s_atom[a]; }
    double  g_val(void)   const { return s_val; }

    // do-nothing function overridden only by torsion class
    virtual void fix_near_180(void) { return; }

    // each internal coordinate type must provide the following virtual functions:

    // compute value of internal coordinate and return
    virtual void compute_val(double **geom) = 0;

    // compute s vector (xyz) for atom a and return
    virtual double * g_s(double **geom, const int iatom) const = 0;

    // print coordinates and value to output file
    virtual void print(const FILE *fp) const = 0;

    // print coordinates and displacements 
    virtual void print_disp(const FILE *fp, const double old_q, const double f_q,
      const double dq, const double new_q) const = 0;

    // for debugging, print s vectors to output file
    virtual void print_s(const FILE *fp, double **geom) const = 0;

    // each derived class should have an equality operator of this type that
    // first checks types with g_type() and then, if true, compares
    // the internal coordinates
    virtual bool operator==(const SIMPLE & s2) const  = 0;


};

}

#endif
