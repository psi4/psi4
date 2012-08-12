#ifndef _psi_src_lib_libmints_multipoles_h_
#define _psi_src_lib_libmints_multipoles_h_

#include <vector>
#include "typedefs.h"

namespace psi {

class ObaraSaikaTwoCenterMIRecursion;

/*! \ingroup MINTS
 *  \class MultipoleInt
 *  \brief Computes arbitrary-order multipole integrals.
 *
 * Use an IntegralFactory to create this object. */
class MultipoleInt : public OneBodyAOInt
{
    //! Obara and Saika recursion object to be used.
    ObaraSaikaTwoCenterMIRecursion mi_recur_;

    //! Computes the multipole integrals between two gaussian shells.
    void compute_pair(const GaussianShell&, const GaussianShell&);
    
    //! Computes the multipole derivative between two gaussian shells.
    void compute_pair_deriv1(const GaussianShell&, const GaussianShell&){ throw PSIEXCEPTION("NYI"); }

    //! The order of multipole moment to compute
    int order_;
public:
    //! Constructor. Do not call directly. Use an IntegralFactory.
    MultipoleInt(std::vector<SphericalTransform>&, boost::shared_ptr<BasisSet>,
                 boost::shared_ptr<BasisSet>, int order, int deriv=0);
    //! Virtual destructor
    virtual ~MultipoleInt();

    //! Does the method provide first derivatives?
    bool has_deriv1() { return false; }

    /// Returns the nuclear contribution to the multipole moments, with angular momentum up to order
    static SharedVector nuclear_contribution(boost::shared_ptr<Molecule> mol, int order);
};

}
#endif
