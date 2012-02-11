#ifndef _psi_src_lib_libmints_tracelessquadrupole_h_
#define _psi_src_lib_libmints_tracelessquadrupole_h_

#include <vector>
#include <boost/shared_ptr.hpp>

namespace psi {

class OneBodyAOInt;
class ObaraSaikaTwoCenterRecursion;
class GaussianShell;
class SphericalTransform;
class BasisSet;
class Matrix;

/*! \ingroup MINTS
 *  \class TracelessQuadrupoleInt
 *  \brief Computes quadrupole integrals. At last check this may not be working.
 *  Use an IntegralFactory to create this object.
 */
class TracelessQuadrupoleInt : public OneBodyAOInt
{
    ObaraSaikaTwoCenterRecursion overlap_recur_;

    // This the work horse function.
    void compute_pair(const GaussianShell&, const GaussianShell&);
public:
    TracelessQuadrupoleInt(std::vector<SphericalTransform>&, boost::shared_ptr<BasisSet>, boost::shared_ptr<BasisSet>);
    virtual ~TracelessQuadrupoleInt();
};

}

#endif
