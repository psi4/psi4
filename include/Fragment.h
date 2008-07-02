#ifndef _psi_include_fragment_h_
#define _psi_include_fragment_h_

#include <string>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <Element_to_Z.h>

namespace psi {

using std::string;

class Fragment {
  friend class Molecular_system;
  int num_atoms;  // number of atoms in fragment
  double *Z;         // atomic numbers
  double **geom;   // cartesian coordinates
  string *atom_label;

 public:
  Fragment() {};  // default constructor is empty
  ~Fragment() {
    delete [] Z;
    delete [] atom_label;
    if (geom != NULL) free_block(geom);
  }
  Fragment(const Fragment &);             // deep copy constructor
  Fragment & operator=(const Fragment &); // deep assignment

  // read cart geom from input file and geom_label
  void read_cartesian_from_input(string geom_label,double conv_factor=1.0);
  int get_num_atoms(void) {return num_atoms;};
};

} // psi

#endif

