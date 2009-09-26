#ifndef _psi_src_lib_libmints_overlap_h_
#define _psi_src_lib_libmints_overlap_h_

/*!
    \file libmints/overlap.h
    \ingroup MINTS
*/

#include <libutil/ref.h>

#include <libmints/basisset.h>
#include <libmints/gshell.h>
#include <libmints/osrecur.h>
#include <libmints/onebody.h>
#include <libmints/integral.h>

namespace psi {
    
/// This class computes overlap integrals and soon overlap integral derivatives.
/// Use an IntegralFactory to create this object.
class OverlapInt : public OneBodyInt
{
    /// Generic Obara Saika recursion object.
    ObaraSaikaTwoCenterRecursion overlap_recur_;
    
    /// Computes the overlap between a given shell pair.
    void compute_pair(GaussianShell* , GaussianShell*);
    void compute_pair_deriv1(GaussianShell*, GaussianShell*);
    void compute_pair_deriv2(GaussianShell*, GaussianShell*);
    
public:
    /// Constructor, it assumes you are not computing derivatives by default
    OverlapInt(IntegralFactory*, BasisSet*, BasisSet*, int deriv=0);
    ~OverlapInt();
    
    /// Compute overlap between 2 shells. Result is stored in buffer.
    void compute_shell(int, int);
    void compute_shell_deriv1(int, int);
    // void compute_shell_deriv2(int, int);
    
    /// Does the method provide first derivatives?
    bool has_deriv1() { return true; }
    /// Does the method provide second derivatives?
    bool has_deriv2() { return true; }
};

}

#endif
