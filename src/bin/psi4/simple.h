/*! \file simple.h
    \ingroup OPT09
    \brief simple internal coordinate base class
*/

#ifndef _psi3_src_bin_opt09_simple_h_
#define _psi3_src_bin_opt09_simple_h_

#include <cstdio>
#include <typeinfo>
extern "C" { extern FILE *outfile; }

namespace psi { namespace opt09 {

enum intco_type {min, stre, bend, max};

class SIMPLE {
  protected:
    intco_type itype; // type of simple
    int id;          // identifier; a natural number
    double val;      // value
    int na; // # of atoms in internal definition

    SIMPLE (int id_in) {
      // fprintf(outfile,"constructing SIMPLE\n");
      id = id_in;
      val = 0.0; // value in angstroms or degrees
    };

  public:
    //virtual ~SIMPLE() { fprintf(outfile,"destructing SIMPLE\n"); };

    intco_type get_itype(void) const {return itype; }
    int get_na(void) const {return na; }

    virtual double get_val(void) const = 0;
    virtual int get_id(void) const = 0 ;

    // return atom number for A,B,C,... in internal
    virtual int get_atom(int a) const = 0;

    // compute value of simple
    virtual void compute(double **geom) = 0;

    // compute s vectors of simple
    virtual void compute_s(double **geom) = 0;

    // return s-vector value for atom 0,1,2,... and xyz
    virtual double get_s(int atom_index, int xyz) const = 0;

    // if print_flag == 0, print intco formatted definition of coordinate
    // if print_flag == 1, print format for output with value
    virtual void print(int print_flag) const = 0;

    // print s vectors
    virtual void print_s(void) const = 0;

    // also define an equality operator for derived classes like
    // virtual bool operator==(const STRETCH & s2) const;
};

}}

#endif
