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
    void compute_pair(const GaussianShell&, const GaussianShell&);
    /// Computes integrals between two shell objects.
    void compute_pair_deriv1(const GaussianShell&, const GaussianShell&);

protected:

    /// Use range-separation or not? Defaults to false. If so, produce <m|erf(\omega r) / r|n> integrals
    bool use_omega_;

    /// The range-separation parameter. Defaults to 0.0
    double omega_;

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

    /// Set omega value, turns use_omega_ to true
    void set_omega(double omega) { use_omega_ = (omega != 0.0); omega_ = omega; }

    /// Set the value of the use_omega_ flag
    void use_omega(bool yes) { use_omega_ = yes; }

    /// Does the method provide first derivatives?
    bool has_deriv1() { return false; }
};

}

#endif
