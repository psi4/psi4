#ifndef _psi_src_lib_libmints_dipole_h_
#define _psi_src_lib_libmints_dipole_h_

#include <vector>
#include "typedefs.h"

namespace psi {

/*! \ingroup MINTS
 *  \class DipoleInt
 *  \brief Computes dipole integrals.
 *
 * Use an IntegralFactory to create this object. */
class DipoleInt : public OneBodyAOInt
{
    //! Obara and Saika recursion object to be used.
    ObaraSaikaTwoCenterRecursion overlap_recur_;

    //! Computes the dipole between two gaussian shells.
    void compute_pair(const boost::shared_ptr<GaussianShell>&, const boost::shared_ptr<GaussianShell>&);
    //! Computes the dipole derivative between two gaussian shells.
    void compute_pair_deriv1(const boost::shared_ptr<GaussianShell>&, const boost::shared_ptr<GaussianShell>&);

public:
    //! Constructor. Do not call directly use an IntegralFactory.
    DipoleInt(std::vector<SphericalTransform>&, boost::shared_ptr<BasisSet>, boost::shared_ptr<BasisSet>, int deriv=0);
    //! Virtual destructor
    virtual ~DipoleInt();

    //! Does the method provide first derivatives?
    bool has_deriv1() { return true; }

    /// Returns the nuclear contribution to the dipole moment
    static SharedVector nuclear_contribution(boost::shared_ptr<Molecule> mol);
};

}

#endif
