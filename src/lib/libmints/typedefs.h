#ifndef libmints_typedefs_h
#define libmints_typedefs_h

// Forward declare boost
namespace boost {
template <class T>
class shared_ptr;
}

// Forward declare psi
namespace psi {

// libmints objects
class BasisSet;
class CoordEntry;
class CoordValue;
class GaussianShell;
class IntegralFactory;
class Matrix;
class Molecule;
class ObaraSaikaTwoCenterRecursion;
class ObaraSaikaTwoCenterVIDeriv2Recursion;
class OneBodyAOInt;
class PointGroup;
class SimpleVector;
class SphericalTransform;
class Vector;
class Vector3;

// objects from other libraries
class Chkpt;
class PSIO;
}

typedef boost::shared_ptr<psi::Matrix> SharedMatrix;
typedef boost::shared_ptr<psi::Vector> SharedVector;

#endif // libmints_typedefs_h
