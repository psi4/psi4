#ifndef _psi_src_lib_libmints_eri_h
#define _psi_src_lib_libmints_eri_h

/*!
    \file libmints/eri.h
    \ingroup MINTS
*/

#include <libutil/ref.h>

#include <libmints/basisset.h>
#include <libmints/gshell.h>
#include <libmints/osrecur.h>
#include <libmints/onebody.h>
#include <libmints/twobody.h>
#include <libmints/integral.h>

#include <libint/libint.h>
#include <libderiv/libderiv.h>

namespace psi {
    
class ERI : public TwoBodyInt
{
    //! Libint object.
    Libint_t libint_;
    //! Libderiv object
    Libderiv_t libderiv_;

    //! Maximum cartesian class size.
    int max_cart_;
    double **d_;
    double *denom_;
    double wval_infinity_;
    int itable_infinity_;
    
    void init_fjt(int);
    void int_fjt(double *, int, double);

    //! Computes the ERIs between four shells.
    void compute_quartet(int, int, int, int);

    //! Computes the ERI derivatives between four shells.
    void compute_quartet_deriv1(int, int, int, int);
public:
    //! Constructor. Use an IntegralFactory to create this object.
    ERI(shared_ptr<BasisSet>, shared_ptr<BasisSet>, shared_ptr<BasisSet>, shared_ptr<BasisSet>, int deriv=0);
    ~ERI();
    
    /// Compute ERIs between 4 shells. Result is stored in buffer.
    void compute_shell(int, int, int, int);

    /// Compute ERI derivatives between 4 shells. Result is stored in buffer.
    void compute_shell_deriv1(int, int, int, int);
};

}

#endif
