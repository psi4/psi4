/*
 *  hfenergy.cpp
 *  matrix
 *
 *  Created by Justin Turney on 4/9/08.
 *  Copyright 2008 by Justin M. Turney, Ph.D.. All rights reserved.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <string>

#include <libipv1/ip_lib.h>

#include "hfenergy.h"
#include "libscf_solver/rhf.cc"
#include "libscf_solver/rohf.h" 
#include "libscf_solver/uhf.h"
#include "libscf_solver/rks.h"
#include "libscf_solver/uks.h"

using namespace std;
using namespace psi;

namespace psi { namespace scf {
     
HFEnergy::HFEnergy(Options & options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt) 
    : Wavefunction(options, psio, chkpt)
{
    
}

double HFEnergy::compute_energy()
{
    // Check the requested reference in the input file
    string reference;
    double energy;
    
    reference = options_.get_str("REFERENCE");
    
    if (reference == "RHF") {
        RHF rhf_energy(options_, psio_, chkpt_);
        energy = rhf_energy.compute_energy();
    }
    else if (reference == "ROHF") {
        ROHF rohf_energy(options_, psio_, chkpt_);
        energy = rohf_energy.compute_energy();
    }
    else if (reference == "UHF") {
        UHF uhf_energy(options_, psio_, chkpt_);
        energy = uhf_energy.compute_energy();
    }
    else if (reference == "RKS") {
        RKS rks_energy(options_, psio_, chkpt_);
        energy = rks_energy.compute_energy();
    }
    else if (reference == "UKS") {
        UKS uks_energy(options_, psio_, chkpt_);
        energy = uks_energy.compute_energy();
    }
    else {
        throw InputException("Unknown reference " + reference, "REFERENCE", __FILE__, __LINE__);
    	energy = 0.0;
    }
    
    return energy;
}

}}
