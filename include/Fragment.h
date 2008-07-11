#ifndef _psi_include_fragment_h_
#define _psi_include_fragment_h_

#include <string>
#include <iostream>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <Element_to_Z.h>

namespace psi {

class bad_masses_io
{
  std::string keyword;
  int row;
public:
  bad_masses_io(std::string key_in=0, int row_in=0) : keyword(key_in), row(row_in) {}
  void mesg(void) {
    std::cout << "bad_masses_io exception: cannot read " << keyword <<
     " entry row " << row << std::endl;
  }
};

using std::string;

class Fragment {
  friend class Molecular_system;
  int num_atoms;   // number of atoms in fragment
  double *Z;       // atomic numbers
  double **geom;   // cartesian coordinates
  double *masses;  // isotopic masses
  string *atom_label;

 public:
  Fragment() {};  // default constructor is empty
  ~Fragment() {
    delete [] Z;
    delete [] atom_label;
    if (geom != NULL) free_block(geom);
    delete [] masses;
  }
  Fragment(const Fragment &);             // deep copy constructor
  Fragment & operator=(const Fragment &); // deep assignment

  int get_num_atoms(void) {return num_atoms;};

  // read cart geom from input file and geom_label
  void read_cartesian_from_input(string geom_label,double conv_factor=1.0);
  // read masses from input file - return 1 if successful
  int read_masses_from_input(string geom_label) throw(bad_masses_io);
  // use atomic numbers to set mass values
  void load_default_masses(string geom_label);
};

} // psi

#endif

