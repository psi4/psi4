#ifndef _psi_src_lib_libmints_overlap_h_
#define _psi_src_lib_libmints_overlap_h_

#include <vector>
#include <libmints/onebody.h> // needed because we derive from OneBodyAOInt

namespace psi {

    class BasisSet;
    class GaussianShell;
    class ObaraSaikaTwoCenterRecursion;
    class OneBodyAOInt;
    class IntegralFactory;
    class SphericalTransform;
    class Matrix;

/*! \ingroup MINTS
 *  \class OverlapInt
 *  \brief This class computes overlap integrals and soon overlap integral derivatives.
 *  Use an IntegralFactory to create this object.
 */
class OverlapInt : public OneBodyAOInt
{
    /// Generic Obara Saika recursion object.
    ObaraSaikaTwoCenterRecursion overlap_recur_;

    /// Computes the overlap between a given shell pair.
    void compute_pair(const GaussianShell& , const GaussianShell&);
    void compute_pair_deriv1(const GaussianShell& s1, const GaussianShell& s2);
    void compute_pair_deriv2(const GaussianShell&, const GaussianShell&);

public:
    /// Constructor, it assumes you are not computing derivatives by default
    OverlapInt(std::vector<SphericalTransform>&, boost::shared_ptr<BasisSet>, boost::shared_ptr<BasisSet>, int deriv=0);
    virtual ~OverlapInt();

    /// Does the method provide first derivatives?
    bool has_deriv1() { return true; }
    /// Does the method provide second derivatives?
    bool has_deriv2() { return true; }
};

}

#endif
