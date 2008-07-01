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
  int num_atoms;
  int num_fragments;
  vector<Fragment> fragment;
  int charge;

 public:
  Molecular_system();
  ~Molecular_system() { fragment.clear(); }
  Molecular_system(const Molecular_system & sys); // deep copy constructor

  int get_num_atoms(void) { return num_atoms; }
  int get_num_fragments(void) { return num_fragments; }
  int get_num_fragments2(void) { return fragment.size(); }
  double **get_geom(void);
  double *get_Z(void);
  string *get_atom_label(void);
  void print(void) const;

};

} // psi

#endif

