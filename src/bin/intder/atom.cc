/*! \file
    \ingroup INTDER
    \brief Enter brief description of file here 
*/
#include "atom.h"

using namespace psi::intder;

Atom::Atom(int an, double aw, double ax, double ay, double az) 
    : Cartesian(ax, ay, az)
{
  atomicNumber = an;
  atomicWeight = aw;
}

Atom::Atom(int an, double aw, Cartesian& A)
    : Cartesian(A)
{
  atomicNumber = an;
  atomicWeight = aw;
}

Atom::Atom(int an, double aw)
  : Cartesian()
{
  atomicNumber = an;
  atomicWeight = aw;
}

Atom::Atom()
  : Cartesian()
{
  atomicNumber = 0;
  atomicWeight = 0.0;
}

