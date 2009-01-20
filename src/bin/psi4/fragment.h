#ifndef _psi_include_fragment_h_
#define _psi_include_fragment_h_

#include <string>
#include <iostream>
#include <sstream>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include "element_to_Z.h"
#include <masses.h>

extern "C" { extern FILE *outfile; }

namespace psi { namespace opt09 {

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
using std::ostringstream;

class Fragment {
  friend class Molecular_system;
  int natoms;   // number of atoms in fragment
  double *Z;       // atomic numbers
  double *masses;  // isotopic masses
  string *atom_label;
  double **geom;   // cartesian coordinates

 public:
  Fragment() {};  // default constructor is empty
  ~Fragment() {
    delete [] Z;
    delete [] masses;
    delete [] atom_label;
    if (geom != NULL) { free_block(geom); geom = NULL; }
  }
  Fragment(const Fragment &);             // deep copy constructor
  Fragment & operator=(const Fragment &); // deep assignment

  int get_natoms(void) const {return natoms;};

  // read cart geom from input file and geom_label
  void read_cartesian_from_input(string geom_label,double conv_factor=1.0);

  // read masses from input file - return 1 if successful
  int read_masses_from_input(int fragment_id) throw(bad_masses_io);

  // use atomic numbers to set mass values
  void read_default_masses(void);

  // use atomic numbers to set mass values
  void print(void) const {
    for (int i=0; i<natoms; ++i) {
      fprintf(outfile,"%5.2lf %12.8lf %s %13.8lf %13.8lf %13.8lf \n", Z[i], masses[i], atom_label[i].c_str(), \
        geom[i][0], geom[i][1], geom[i][2]);
    }
  };
};

}}

#endif

