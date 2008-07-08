/*! \file
    \ingroup OPT08
    \brief Class for Salcs
*/

#ifndef _psi3_bin_optking_salc_h_
#define _psi3_bin_optking_salc_h_

#include <vector>
#include <iostream>
#include <string>
#include <cmath>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include "except.h"

extern "C" { extern FILE *outfile; }

using std::string;
using std::vector;

namespace psi { namespace optking {

const int MAX_LINELENGTH = 133;

inline double sqr(double val) { return (val*val); }
inline char * c_string(string);

class One_salc {
  string lbl;
  double prefactor;
  int length;
  vector<int> simple;
  vector<double> coeff;

  void set_prefactor(void);

  public:
    friend class Salc_set;
    One_salc(string lbl_in, int length_in, vector<int> simple_in,
      vector<double> coeff_in)
    {
      lbl = lbl_in;
      length = length_in;
      simple.resize(length);
      coeff.resize(length);
      simple = simple_in;
      coeff = coeff_in;
      if (length>0) set_prefactor();
    }
    ~One_salc() {
      simple.clear();
      coeff.clear();
    }
    One_salc(const One_salc & os)
    {
      lbl = os.lbl;
      length = os.length;
      simple.resize(length);
      coeff.resize(length);
      simple = os.simple;
      coeff  = os.coeff;
      prefactor = os.prefactor;
    }
    void print(int print_flag=0) const;

};

class Salc_set {
 string keyword;
 int irrep; // set to -1 if not symmetry adapted?
 vector<One_salc> salc;
 int natom;
 int nsimples;
 int matrix_present;
 double **matrix; // num of salcs by num of cartesians

 public:

   Salc_set(string key_in = "SYMM", int irrep_in = 0,
     int natom=0, int nsimples=0) throw(bad_intco_io);
   ~Salc_set();
   void print(int print_flag = 0) const;
   Salc_set(const Salc_set &);
   int get_natom(void) const {return natom;}
   int size(void) const { return salc.size();}
   void build_matrix(void);
   void print_matrix(void) const;
   double **get_matrix(void) { return matrix; }
};

}} /* namespace psi::optking */

#endif
