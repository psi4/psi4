#ifndef RKS_H
#define RKS_H
/*
 *  rks.h
 *
 *  Created by Rob Parrish on 01/20/2010
 *
 */


#include <libpsio/psio.hpp>
#include "hf.h"
#include "rhf.h"
#include <libfunctional/superfunctional.h>
#include "integrator.h"
#include <libmints/properties.h>
using namespace psi;
using namespace psi::functional;

namespace psi { namespace scf {
     
/*! \ingroup SCF */
/* Class RKS contains the member data and methods needed for
* Restricted Kohn-Sham DFT computations
*/
class RKS : public RHF {
protected:
    /// Vxc matrix
    SharedMatrix V_;
    /// Exchange-Correlation Functional
    shared_ptr<SuperFunctional> functional_;
    /// Integrator
    SharedIntegrator integrator_;
    /// Properties evaluator
    SharedProperties properties_;
    /// Density check values
    double densityCheck_;   
    double dipoleCheckX_;   
    double dipoleCheckY_;   
    double dipoleCheckZ_;  

    /// Correlation Functional Energy
    double functional_energy_;
 
public:
    /// Constructor, same as RHF, which it derives from
    RKS(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    /// Destructor, frees member data
    virtual ~RKS();
    /// Computes KS-DFT energy and C/D matrices
    double compute_energy();
    /// Computes the Vxc matrix
    void form_V();
    void form_F();
    void form_J();
    void form_K();
    
    double compute_E();
    /// Save DFT grid points and weights
    void save_DFT_grid();
};

}}

#endif
