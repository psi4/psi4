#ifndef _psi4_src_lib_fragment_h_
#define _psi4_src_lib_fragment_h_

#include <string>
#include <libipv1/ip_lib.h>
#include <element_to_Z.h>

namespace psi {

using std::string;

class Fragment {
  int num_atoms;  // number of atoms in fragment
  double *Z;         // atomic numbers
  double **geom;   // cartesian coordinates
 public:
  Fragment() {};  // default constructor is empty
  ~Fragment() {
    delete [] Z;
    if (geom != NULL) free_block(geom);
  }

  // read cart geom from input file and geom_label
  void read_cartesian_from_input(string geom_label);
};

// lets obselete non-"simple" geometry format:
// geometry = ( (atom 1) (atom 2) ... )
void Fragment::read_cartesian_from_input(string geom_label) {
  int i, j, errcod, num_elem;
  double tmp;
  char *c_geom_label, *c_atom_label;
  string error_message, string_atom_label;

  c_geom_label = const_cast<char *>(geom_label.c_str());

  // Determine number of atoms in fragment
  num_elem = 0;
  ip_count(c_geom_label, &num_elem, 0);
  if (!num_atoms) {
    error_message = "Fragment " + geom_label + " has no atoms.";
    throw(error_message);
  }
  if (num_elem % 4) {
    error_message = "Problem with number of entries in " + geom_label;
    throw(error_message);
  }
  num_atoms = num_elem/4;

  // allocate memory for Z, geom
  Z = init_array(num_atoms);
  geom = block_matrix(num_atoms,3);

  for (i=0; i<num_atoms; ++i) {
    // read symbol and determine atomic number
    errcod = ip_string(c_geom_label,&c_atom_label,1,4*i);
    if (errcod != IPE_OK) {
      error_message = "Problem reading the " + geom_label + " array.";
      throw(error_message);
    }
    Element_to_Z elem_map;
    Z[i] = elem_map[string_atom_label = c_atom_label];
    free(c_atom_label);

    // read geometry
    tmp = 0.0;
    errcod = ip_data(c_geom_label,"%lf", &tmp,1,4*i+j+1);
    if (errcod != IPE_OK) {
      error_message = "Problem reading the " + geom_label + " array.";
      throw(error_message);
    }
    geom[i][j] = tmp;
  }

  return;
}

} // psi

#endif
