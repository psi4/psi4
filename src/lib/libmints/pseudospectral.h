#ifndef _psi_src_lib_libmints_pseudospectral_h_
#define _psi_src_lib_libmints_pseudospectral_h_

#include <vector>
#include <boost/shared_ptr.hpp>

namespace psi {

    class BasisSet;
    class GaussianShell;
    class ObaraSaikaTwoCenterVIRecursion;
    class ObaraSaikaTwoCenterVIDerivRecursion;
    class OneBodyAOInt;
    class IntegralFactory;
    class SphericalTransform;

/*! \ingroup MINTS
 *  \class PotentialInt
 *  \brief Computes pseudospectral integrals.
 * Use an IntegralFactory to create this object.
 */
class PseudospectralInt : public OneBodyAOInt
{

    /// Computes integrals between two shell objects.
    void compute_pair(const boost::shared_ptr<GaussianShell>&, const boost::shared_ptr<GaussianShell>&);
    /// Computes integrals between two shell objects.
    void compute_pair_deriv1(boost::shared_ptr<GaussianShell>, boost::shared_ptr<GaussianShell>);

protected:

    /// The integration point 
    double C_[3];
    /// Recursion object that does the heavy lifting.
    ObaraSaikaTwoCenterVIRecursion potential_recur_;
    /// Recursion object that does the heavy lifting.
    ObaraSaikaTwoCenterVIDerivRecursion potential_deriv_recur_;

public:
    /// Constructor
    PseudospectralInt(std::vector<SphericalTransform>&, boost::shared_ptr<BasisSet>, boost::shared_ptr<BasisSet>, int deriv=0);
    ~PseudospectralInt();

    /// Computes integrals between two shells.
    void compute_shell_deriv1(int, int);

    /// Set integration point
    void set_point(double x, double y, double z) { C_[0] = x; C_[1] = y; C_[2] = z; }

    /// Does the method provide first derivatives?
    bool has_deriv1() { return false; }
};

}

#endif
