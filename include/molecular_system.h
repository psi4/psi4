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
