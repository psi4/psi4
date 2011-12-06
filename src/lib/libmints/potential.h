#ifndef _psi_src_lib_libmints_potential_h_
#define _psi_src_lib_libmints_potential_h_

#include <vector>
#include "typedefs.h"

namespace psi {

    class Matrix;
    class BasisSet;
    class GaussianShell;
    class ObaraSaikaTwoCenterVIRecursion;
    class ObaraSaikaTwoCenterVIDerivRecursion;
    class ObaraSaikaTwoCenterVIDeriv2Recursion;
    class OneBodyAOInt;
    class IntegralFactory;
    class SphericalTransform;
    class OneBodySOInt;
    class CdSalcList;

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
    void compute_pair_deriv2(const boost::shared_ptr<GaussianShell>&, const boost::shared_ptr<GaussianShell>& );

protected:
    /// Recursion object that does the heavy lifting.
    ObaraSaikaTwoCenterVIRecursion* potential_recur_;

    /// Matrix of coordinates/charges of partial charges
    SharedMatrix Zxyz_;

public:
    /// Constructor. Assumes nuclear centers/charges as the potential
    PotentialInt(std::vector<SphericalTransform>&, boost::shared_ptr<BasisSet>, boost::shared_ptr<BasisSet>, int deriv=0);
    virtual ~PotentialInt();

    /// Computes the first derivatives and stores them in result
    virtual void compute_deriv1(std::vector<SharedMatrix > &result);

    /// Set the field of charges
    void set_charge_field(SharedMatrix Zxyz) { Zxyz_ = Zxyz; }

    /// Get the field of charges
    SharedMatrix charge_field() const { return Zxyz_; }

    /// Does the method provide first derivatives?
    bool has_deriv1() { return true; }
};

class PotentialSOInt : public OneBodySOInt
{
    int natom_;
public:
    PotentialSOInt(const boost::shared_ptr<OneBodyAOInt>& , const boost::shared_ptr<IntegralFactory> &);
    PotentialSOInt(const boost::shared_ptr<OneBodyAOInt>& , const IntegralFactory*);

    /**
     * Computes one-electron integral derivative matrices.
     * Specifically handles CdSalc SO potential integral derivatives.
     *
     * \param result Where the integral derivatives are going.
     * \param cdsalcs The Cartesian displacement SALCs that you are interested in.
     */
    void compute_deriv1(std::vector<SharedMatrix > result,
                        const CdSalcList& cdsalcs);
};

}

#endif
