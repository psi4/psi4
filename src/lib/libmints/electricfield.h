#ifndef _psi_src_lib_libmints_electricfield_h_
#define _psi_src_lib_libmints_electricfield_h_

#include <vector>
#include "typedefs.h"

namespace psi {

/*! \ingroup MINTS
 *  \class ElectricFieldInt
 *  \brief Computes electric field integrals.
 *
 *  Use an IntegralFactory to create this object.
 */
class ElectricFieldInt : public OneBodyAOInt
{
    //! Obara and Saika recursion object to be used.
    ObaraSaikaTwoCenterVIDeriv2Recursion efield_recur_;

    //! Number of atoms.
    int natom_;

    //! Computes the electric field between two gaussian shells.
    void compute_pair(const boost::shared_ptr<GaussianShell>&, const boost::shared_ptr<GaussianShell>&);

    //! Computes the electric field gradient between two gaussian shells.
    void compute_pair_deriv1(const boost::shared_ptr<GaussianShell>&, const boost::shared_ptr<GaussianShell>&);

public:
    //! Constructor. Do not call directly use an IntegralFactory.
    ElectricFieldInt(std::vector<SphericalTransform>&, boost::shared_ptr<BasisSet>, boost::shared_ptr<BasisSet>, int deriv=0);
    //! Virtual destructor
    virtual ~ElectricFieldInt();

    //! Does the method provide first derivatives?
    bool has_deriv1() { return true; }

    static SharedMatrix nuclear_contribution(boost::shared_ptr<Molecule> mol);
};

}

#endif

