#ifndef _psi_src_lib_libmints_potential_h_
#define _psi_src_lib_libmints_potential_h_

/*!
    \file libmints/potential.h
    \ingroup MINTS
*/

#include <libutil/ref.h>

#include <libmints/basisset.h>
#include <libmints/gshell.h>
#include <libmints/osrecur.h>
#include <libmints/onebody.h>
#include <libmints/integral.h>

namespace psi {
    
/// Computes potential integrals.
/// Use an IntegralFactory to create this object.
class PotentialInt : public OneBodyInt
{
    /// Recursion object that does the heavy lifting.
    ObaraSaikaTwoCenterVIRecursion potential_recur_;
    /// Recursion object that does the heavy lifting.
    ObaraSaikaTwoCenterVIDerivRecursion potential_deriv_recur_;
    
    /// Computes integrals between two shell objects.
    void compute_pair(GaussianShell*, GaussianShell*);
    /// Computes integrals between two shell objects.
    void compute_pair_deriv1(GaussianShell*, GaussianShell*);
    
public:
    /// Constructor
    PotentialInt(IntegralFactory*, BasisSet*, BasisSet*, int deriv=0);
    ~PotentialInt();
    
    /// Computes integrals between two shells.
    void compute_shell(int, int);
    /// Computes integrals between two shells.
    void compute_shell_deriv1(int, int);
    
    /// Does the method provide first derivatives?
    bool has_deriv1() { return true; }
};

}

#endif
