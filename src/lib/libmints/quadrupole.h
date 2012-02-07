#ifndef _psi_src_lib_libmints_quadrupole_h_
#define _psi_src_lib_libmints_quadrupole_h_

#include <vector>
#include "typedefs.h"

namespace psi {

class OneBodyAOInt;
class ObaraSaikaTwoCenterRecursion;
class GaussianShell;
class SphericalTransform;
class BasisSet;
class Matrix;

/*! \ingroup MINTS
 *  \class QuadrupoleInt
 *  \brief Computes quadrupole integrals. At last check this may not be working.
 *  Use an IntegralFactory to create this object.
 */
class QuadrupoleInt : public OneBodyAOInt
{
    ObaraSaikaTwoCenterRecursion overlap_recur_;

    // This the work horse function.
    void compute_pair(const GaussianShell&, const GaussianShell&);
public:
    QuadrupoleInt(std::vector<SphericalTransform>&, boost::shared_ptr<BasisSet>, boost::shared_ptr<BasisSet>);
    virtual ~QuadrupoleInt();

    static SharedVector nuclear_contribution(boost::shared_ptr<Molecule> mol);

};

}

#endif
