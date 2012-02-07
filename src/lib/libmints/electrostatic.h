#ifndef _psi_src_lib_libmints_electrostatic_h_
#define _psi_src_lib_libmints_electrostatic_h_

namespace boost {
template<class T> class shared_ptr;
}
namespace psi {

class BasisSet;
class GaussianShell;
class ObaraSaikaTwoCenterRecursion;
class OneBodyAOInt;
class PotentialInt;
class IntegralFactory;
class SphericalTransform;
class Vector3;

/*! \ingroup MINTS
 *  \class PotentialInt
 *  \brief Computes potential integrals.
 * Use an IntegralFactory to create this object.
 */
class ElectrostaticInt : public PotentialInt
{
    void compute_pair(const GaussianShell&, const GaussianShell&)
    {}

public:
    /// Constructor
    ElectrostaticInt(std::vector<SphericalTransform>&, boost::shared_ptr<BasisSet>, boost::shared_ptr<BasisSet>, int deriv=0);
    ~ElectrostaticInt();

    // Intel C++ 12 thinks we're trying to overload the "void compute_shell(int, int)" and warns us about it.
    // The following line is to shut it up.
    #pragma warning disable 1125
    /// Computes integrals between two shells.
    void compute_shell(int, int, Vector3&);
    /// Computes integrals between two shells.
    void compute_pair(const GaussianShell&, const GaussianShell&, Vector3&);

    /// Does the method provide first derivatives?
    bool has_deriv1() { return false; }
};

}

#endif
