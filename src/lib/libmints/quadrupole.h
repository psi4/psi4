#ifndef _psi_src_lib_libmints_quadrupole_h_
#define _psi_src_lib_libmints_quadrupole_h_

#include <libmints/basisset.h>
#include <libmints/gshell.h>
#include <libmints/osrecur.h>
#include <libmints/onebody.h>
#include <libmints/integral.h>

namespace psi {

/*! \ingroup MINTS
 *  \class QuadrupoleInt
 *  \brief Computes quadrupole integrals. At last check this may not be working.
 *  Use an IntegralFactory to create this object.
 */
class QuadrupoleInt : public OneBodyInt
{
    ObaraSaikaTwoCenterRecursion overlap_recur_;
    
    void compute_pair(shared_ptr<GaussianShell>, shared_ptr<GaussianShell>);
    
public:
    QuadrupoleInt(std::vector<SphericalTransform>&, shared_ptr<BasisSet>, shared_ptr<BasisSet>);
    virtual ~QuadrupoleInt();
    
    void compute_shell(int, int);
    
    /// Computes all quadrupole integrals (Qxx, Qxy, Qxz, Qyy, Qyz, Qzz) result must be an array of enough
    /// size to contain it.
    void compute(std::vector<shared_ptr<Matrix> > &result);
    void compute(std::vector<shared_ptr<SimpleMatrix> > &result);
    
    virtual void spherical_transform(shared_ptr<GaussianShell> , shared_ptr<GaussianShell>);
};

}

#endif
