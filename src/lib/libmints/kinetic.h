#ifndef _psi_src_lib_libmints_kinetic_h_
#define _psi_src_lib_libmints_kinetic_h_

#include <boost/shared_ptr.hpp>
#include <vector>

namespace psi {

    class BasisSet;
    class GaussianShell;
    class ObaraSaikaTwoCenterRecursion;
    class OneBodyInt;
    class IntegralFactory;
    class SphericalTransform;
    class SimpleMatrix;

/*! \ingroup MINTS
 *  \class KineticInt
 *  \brief Computes kinetic integrals.
 *
 * Use an IntegralFactory to create this object.
 */
class KineticInt : public OneBodyInt
{
    //! Obara and Saika recursion object to be used.
    ObaraSaikaTwoCenterRecursion overlap_recur_;

    //! Computes the kinetic integral between two gaussian shells.
    void compute_pair(const boost::shared_ptr<GaussianShell>&, const boost::shared_ptr<GaussianShell>&);
    //! Computes the kinetic derivatve between two gaussian shells.
    void compute_pair_deriv1(boost::shared_ptr<GaussianShell>, boost::shared_ptr<GaussianShell>);

public:
    //! Constructor. Do not call directly, use an IntegralFactory.
    KineticInt(std::vector<SphericalTransform>&, boost::shared_ptr<BasisSet>, boost::shared_ptr<BasisSet>, int deriv=0);
    //! Virtual destructor.
    virtual ~KineticInt();

    //! Compute kinetic derivative between two shells, result stored in buffer_.
    void compute_shell_deriv1(int, int);

    /// Does the method provide first derivatives?
    bool has_deriv1() { return true; }
};

}

#endif
