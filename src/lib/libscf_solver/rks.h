/*
 *  rks.h
 *
 *  Created by Rob Parrish on 01/20/2010
 *
 */

#ifndef RKS_H
#define RKS_H

#include <libpsio/psio.hpp>
#include "hf.h"
#include "rhf.h"

using namespace psi;

namespace psi { namespace scf {
     
class RKS : public RHF {
protected:
    
public:
    RKS(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    virtual ~RKS();
    
    double compute_energy();
};

}}

#endif
