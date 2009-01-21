#ifndef _psi_include_molecular_system_h_
#define _psi_include_molecular_system_h_

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <fragment.h>
#include <psi4-dec.h> // for outfile

namespace psi {

using std::vector;
using std::string;
using std::ostringstream;

class Molecular_system {
  vector<Fragment> fragment;
  int charge;

 public:
  Molecular_system(double conv_factor=1.0);
  ~Molecular_system() { fragment.clear(); }
  Molecular_system(const Molecular_system & sys);

  int get_natoms(void) const;
  int get_nfragments(void) const { return fragment.size(); }
  double **get_geom(void) const;
  double *get_Z(void) const;
  string *get_atom_label(void) const;
  void print(void) const;
};

}

#endif
