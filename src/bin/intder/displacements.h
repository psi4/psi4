/*! \file
    \ingroup INTDER
    \brief Enter brief description of file here 
*/
#ifndef _psi_bin_intder_displacements_h_
#define _psi_bin_intder_displacements_h_

#include "molecule.h"
#include <vector>

namespace psi { namespace intder {

class Displacements
{
  std::vector<Molecule> vectorMolecules;

public:
  bool loadFromOptKing();
  bool loadFromCheckPoint();
  bool loadFromInput();

  void printGeometries();
  void moveToCenterOfMass();

  void atom_num(char*, double*);
  void addDisplacement(Molecule& mol);
  void useMasses(double *mass);
  Molecule* displacement(int disp);
};

}} // namespace psi::intder

#endif // header guard

