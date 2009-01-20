#ifndef _psi_include_molecular_system_h_
#define _psi_include_molecular_system_h_

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "fragment.h"

namespace psi { namespace opt09 {

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

  int get_natoms(void) {
    int i, num=0;
    for (i=0; i<fragment.size(); ++i)
      num += fragment.at(i).get_natoms();
    return num;
  };

  int get_nfragments(void) { return fragment.size(); }
  double **get_geom(void);
  double *get_Z(void);
  string *get_atom_label(void);
  void print(void) const;
};

}}

#endif
