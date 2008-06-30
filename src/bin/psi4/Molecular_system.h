#ifndef _psi4_src_lib_molecule_h_
#define _psi4_src_lib_molecule_h_

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
  int num_frag;
  vector<Fragment> frag;
  int charge;

 public:
  Molecular_system();
  ~Molecular_system() { frag.clear(); }
};

// default constructor: get molecule from input file
Molecular_system::Molecular_system() {
  int i, errcod = 0;
  string geom_label;
  bool cart;

  // read how many fragments there are
  num_frag = 1;
  errcod = ip_data("NUM_FRAGMENTS","%d",&i,0);
  if (errcod == IPE_OK) num_frag = i;

  // Cartesian or z-matrix input?
  cart = false;
  if (ip_exist("GEOMETRY",0))
    cart = true;
  else if (ip_exist("ZMAT",0))
    throw("zmatrix not yet implemented");
  else
    throw("could not find GEOMETRY or ZMAT in input!");

  if (cart) {
    Fragment lfrag;
    ostringstream outstr;
    for (i=0; i<num_frag; ++i) {
      outstr << "GEOMETRY" << i;
      geom_label = outstr.str();
      // outstr.clear()?
      lfrag.read_cartesian_from_input(geom_label);
      frag.push_back(lfrag);
    }
  }

}

}

#endif
