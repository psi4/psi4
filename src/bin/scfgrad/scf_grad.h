#ifndef SCF_GRAD_H
#define SCF_GRAD_H

#include <libmints/wavefunction.h>
#include <libmints/typedefs.h>

namespace psi {

namespace scfgrad {

class SCFGrad : public Wavefunction {

protected:

    /// Common initialization
    void common_init();
    
public:
    SCFGrad();
    virtual ~SCFGrad();
    
    double compute_energy() { throw PSIEXCEPTION("SCFGrad needs a rehash, call Rob."); }
   
    SharedMatrix compute_gradient(); 

    SharedMatrix compute_hessian();
};

}} // Namespaces

#endif
