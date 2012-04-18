#ifndef SCF_GRAD_H
#define SCF_GRAD_H

#include <libmints/wavefunction.h>
#include <libmints/typedefs.h>

namespace psi {

namespace scfgrad {

class SCFGrad : public Wavefunction {

protected:

    /// Gradient components
    std::map<std::string, SharedMatrix> gradients_;

    /// Common initialization
    void common_init();
    /// Print the specialization 
    void print_header() const;
    
public:
    SCFGrad();
    virtual ~SCFGrad();
    
    double compute_energy() { return 0.0; }
   
    SharedMatrix compute_gradient(); 
};

}} // Namespaces

#endif
