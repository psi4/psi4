#ifndef _psi_include_molecular_system_h_
#define _psi_include_molecular_system_h_

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "Fragment.h"

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
  Molecular_system(const Molecular_system & sys); // deep copy constructor

  int get_num_atoms(void) {
    int num=0;
    for (int i=0; i<fragment.size(); ++i)
      num += fragment.at(i).num_atoms;
    return num;
  };

  int get_num_fragments(void) { return fragment.size(); }
  double **get_geom(void);
  double *get_Z(void);
  string *get_atom_label(void);
  void print(void) const;

};

} // psi

#endif

