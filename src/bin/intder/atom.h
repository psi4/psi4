/*! \file
    \ingroup INTDER
    \brief Enter brief description of file here 
*/
#ifndef _psi_bin_intder_atom_h_
#define _psi_bin_intder_atom_h_

#include "cartesian.h"

namespace psi { namespace intder {

class Atom : public Cartesian
{
public:
  int atomicNumber;
  double atomicWeight;
  Atom(int an, double aw, double ax, double ay, double az);
  Atom(int an, double aw, Cartesian& A);
  Atom(int an, double aw);
  Atom();

  int getAtomicNumber()
    { return atomicNumber; }
  void setAtomicNumber(int an)
    { atomicNumber = an; }
  double getAtomicWeight()
    { return atomicWeight; }
  void setAtomicWeight(double aw)
    { atomicWeight = aw; }
};

}} // namespace psi::intder

#endif // header guard

