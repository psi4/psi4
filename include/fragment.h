/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#ifndef _psi_include_fragment_h_
#define _psi_include_fragment_h_

#include <string>
#include <iostream>
#include <sstream>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <element_to_Z.h>
#include <masses.h>

#include <psi4-dec.h>

using std::ostringstream;
using std::string;

namespace psi {

/* Throw an InputException instead of this.
class bad_masses_io {
  std::string keyword;
  int row;
public:
  bad_masses_io(std::string key_in=0, int row_in=0) : keyword(key_in), row(row_in) {}
  void mesg(void) {
    std::cout << "bad_masses_io exception: cannot read " << keyword <<
     " entry row " << row << std::endl;
  }
};
*/

/* Use the Molecule class in libmints instead
class Fragment {
  friend class Molecular_system;
  int natoms;         // number of atoms in fragment
  double *Z;          // atomic numbers
  double *masses;     // isotopic masses
  string *atom_label; // atom labels
  double **geom;      // cartesian coordinates

 public:
  Fragment(void) {} ; //default constructor is empty
  ~Fragment() {
    delete [] Z;
    delete [] masses;
    delete [] atom_label;
    if (geom != NULL) { free_block(geom); geom = NULL; }
  }

  Fragment(const Fragment &);             // deep copy constructor
  Fragment & operator=(const Fragment &); // deep assignment

  int get_natoms(void) const { return natoms; }

  // read cart geom from input file and geom_label
  void read_cartesian_from_input(string geom_label,double conv_factor=1.0);

  // read masses from input file - return 1 if successful
  int read_masses_from_input(int fragment_id) throw(bad_masses_io);

  // use atomic numbers to set mass values
  void read_default_masses(void);

  void print(void) const ;
  };
}
*/

#endif
