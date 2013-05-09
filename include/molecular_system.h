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

#ifndef _psi_include_molecular_system_h_
#define _psi_include_molecular_system_h_

#include <string>
#include <vector>
#include <libmints/molecule.h>

namespace psi {

using std::vector;
using std::string;
using std::ostringstream;

class MolecularSystem {
  vector<Molecule> fragment_;
  int charge_;

 public:
  MolecularSystem(double conv_factor=1.0);
  virtual ~MolecularSystem() { }
  MolecularSystem(const MolecularSystem & sys);

  int natom() const;
  int nfragment() const { return fragment_.size(); }
  double **get_geom() const;
  double *get_Z() const;
  string* get_atom_label() const;
  void print() const;
};

}

#endif
