#ifndef _psi4_src_lib_molecule_h_
#define _psi4_src_lib_molecule_h_

#include <string>

namespace psi {

using std::vector;
using std::cout;
using std::string;
using std::endl;

class Molecular_System {
  int num_frag;
  vector<Fragment> frag
  int charge;

 public:
  Molecule();
  ~Molecule();
}

// default constructor: get molecule from input file
Molecule::Molecule(void) {
  int i, errcod = 0, char *label;

  label = new char [20];

  // read how many fragments there are
  num_frag = 1;
  errcod = ip_data("NUM_FRAGMENTS","%d",&i,0)
  if (errcod == IPE_OK) num_frag = i;

  // Cartesian or z-matrix input?
  if (ip_exist("GEOMETRY",0)) {
    Fragment lfrag;
    for (i=0; i<nfragments; ++i) {
      sprintf(label, "GEOMETRY%d", i);
      lfrag.read_cartesian_from_input(label);
    }

  }
  else if (ip_exist("ZMAT",0))
  else throw("could not find GEOMETRY or ZMAT in input!");

  frag.push_back(lfrag);

  delete [] label;

}

}

#endif
