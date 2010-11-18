#ifndef _psi_src_lib_libmints_electrostatic_h_
#define _psi_src_lib_libmints_electrostatic_h_

namespace psi {

    class BasisSet;
    class GaussianShell;
    class ObaraSaikaTwoCenterRecursion;
    class OneBodyInt;
    class PotentialInt;
    class IntegralFactory;
    class SphericalTransform;
    class SimpleMatrix;
    class Vector3;

/*! \ingroup MINTS
 *  \class PotentialInt
 *  \brief Computes potential integrals.
 * Use an IntegralFactory to create this object.
 */
class ElectrostaticInt : public PotentialInt
{
public:
    /// Constructor
    ElectrostaticInt(std::vector<SphericalTransform>&, shared_ptr<BasisSet>, shared_ptr<BasisSet>, int deriv=0);
    ~ElectrostaticInt();

    // Intel C++ 12 thinks we're trying to overload the "void compute_shell(int, int)" and warns us about it.
    // The following line is to shut it up.
    #pragma warning disable 1125
    /// Computes integrals between two shells.
    void compute_shell(int, int, Vector3&);
    /// Computes integrals between two shells.
    void compute_pair(shared_ptr<GaussianShell>, shared_ptr<GaussianShell>, Vector3&);

    /// Does the method provide first derivatives?
    bool has_deriv1() { return false; }
};

}

#endif
