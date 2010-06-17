/*
 *  dfmp2_energy.h
 *  
 *
 *
 */

#ifndef DFMP2ENERGY_H
#define DFMP2ENERGY_H

#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <psi4-dec.h>

namespace psi { namespace dfmp2 {
 
class DFMP2Energy : public Wavefunction {
public:
    DFMP2Energy(Options & options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    virtual ~DFMP2Energy() {}
    
    double compute_E();
};

}}

#endif
