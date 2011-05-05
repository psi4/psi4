#ifndef _psi_src_lib_libmints_potential_h_
#define _psi_src_lib_libmints_potential_h_

#include <vector>
#include <boost/shared_ptr.hpp>

namespace psi {

    class Matrix;
    class BasisSet;
    class GaussianShell;
    class ObaraSaikaTwoCenterVIRecursion;
    class ObaraSaikaTwoCenterVIDerivRecursion;
    class OneBodyAOInt;
    class IntegralFactory;
    class SphericalTransform;

/*! \ingroup MINTS
 *  \class PotentialInt
 *  \brief Computes potential integrals.
 * Use an IntegralFactory to create this object.
 */
class PotentialInt : public OneBodyAOInt
{
    /// Computes integrals between two shell objects.
    void compute_pair(const boost::shared_ptr<GaussianShell>&, const boost::shared_ptr<GaussianShell>&);
    /// Computes integrals between two shell objects.
    void compute_pair_deriv1(const boost::shared_ptr<GaussianShell>&, const boost::shared_ptr<GaussianShell>& );

protected:
    /// Recursion object that does the heavy lifting.
    ObaraSaikaTwoCenterVIRecursion potential_recur_;
    /// Recursion object that does the heavy lifting.
    ObaraSaikaTwoCenterVIDerivRecursion potential_deriv_recur_;

    /// Matrix of coordinates/charges of partial charges
    boost::shared_ptr<Matrix> xyzZ_;

public:
    /// Constructor. Assumes nuclear centers/charges as the potential
    PotentialInt(std::vector<SphericalTransform>&, boost::shared_ptr<BasisSet>, boost::shared_ptr<BasisSet>, int deriv=0);
    ~PotentialInt();

    /// Computes the first derivatives and stores them in result
    virtual void compute_deriv1(std::vector<boost::shared_ptr<SimpleMatrix> > &result);

    /// Set the field of charges
    void setChargeField(boost::shared_ptr<Matrix> xyzZ) { xyzZ_ = xyzZ; }

    /// Get the field of charges
    boost::shared_ptr<Matrix> chargeField() const { return xyzZ_; }

    /// Does the method provide first derivatives?
    bool has_deriv1() { return true; }
};

}

#endif
