/*
 *  hfenergy.h
 *  matrix
 *
 *  Created by Justin Turney on 4/9/08.
 *  Copyright 2008 by Justin M. Turney, Ph.D.. All rights reserved.
 *
 */

#ifndef HFENERGY_H
#define HFENERGY_H

#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>

 using namespace psi;
 
class HFEnergy : public Wavefunction {
public:
    HFEnergy(PSIO *psio, Chkpt *chkpt = 0);
    HFEnergy(PSIO &psio, Chkpt &chkpt);
    virtual ~HFEnergy() {}
    
    double compute_energy();
};

#endif
