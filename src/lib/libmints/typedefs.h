#ifndef libmints_typedefs_h
#define libmints_typedefs_h

#include <compiler.h>

// Handy mints timer macros, requires libqt to be included
#ifdef MINTS_TIMER
#   define mints_timer_on(a) timer_on((a));
#   define mints_timer_off(a) timer_off((a));
#else
#   define mints_timer_on(a)
#   define mints_timer_off(a)
#endif


// Forward declare boost
namespace boost {
template <class T>
class shared_ptr;
}

// Dimension is lightweight
#include "dimension.h"

// Forward declare psi
namespace psi {

// libmints objects
class BasisSet;
class CdSalcList;
class CoordEntry;
class CoordValue;
class GaussianShell;
class IntegralFactory;
class Matrix;
class Molecule;
class ObaraSaikaTwoCenterRecursion;
class ObaraSaikaTwoCenterVIDeriv2Recursion;
class OneBodyAOInt;
class TwoBodyAOInt;
class PointGroup;
class SimpleVector;
class SOBasisSet;
class SphericalTransform;
class Vector;
class Vector3;

// objects from other libraries
class Chkpt;
class PSIO;

typedef boost::shared_ptr<psi::Matrix> SharedMatrix;
typedef boost::shared_ptr<psi::Vector> SharedVector;

// Useful when working with SO-TEIs
template<typename T>
void swap_index(T& a, T& b) {
    T temp;
    temp = b;
    b = a;
    a = temp;
}

#define SWAP_INDEX(a, b) swap_index(a ## abs, b ## abs); swap_index(a ## rel, b ## rel); swap_index(a ## irrep, b ## irrep);

}


#endif // libmints_typedefs_h
