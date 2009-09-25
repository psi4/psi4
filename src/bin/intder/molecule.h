/*! \file
    \ingroup INTDER
    \brief Enter brief description of file here 
*/
#ifndef _psi_bin_intder_molecule_h_
#define _psi_bin_intder_molecule_h_

#include "atom.h"
#include <vector>

namespace psi { namespace intder {

class Molecule
{
  std::vector<Atom> vectorAtoms;

public:
  Molecule();
  ~Molecule();

  void addAtom(int an, double aw, double ax, double ay, double az);
  void useMasses(double *mass);
  void moveToCenterOfMass();

  void printGeometry();
  Atom* atom(int an);
};

}} // namespace psi::intder

#endif // header guard

