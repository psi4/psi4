#ifndef _psi_src_lib_libmints_electrostatic_h_
#define _psi_src_lib_libmints_electrostatic_h_

#include <libmints/basisset.h>
#include <libmints/gshell.h>
#include <libmints/osrecur.h>
#include <libmints/onebody.h>
#include <libmints/integral.h>
#include <libmints/potential.h>
#include <libmints/vector3.h>

namespace psi {

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
    
    /// Computes integrals between two shells.
    void compute_shell(int, int, Vector3&);
    /// Computes integrals between two shells.
    void compute_pair(shared_ptr<GaussianShell>, shared_ptr<GaussianShell>, Vector3&);
    
    /// Does the method provide first derivatives?
    bool has_deriv1() { return false; }
};

}

#endif
