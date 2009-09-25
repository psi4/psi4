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
#include <cstring>

#include <libipv1/ip_lib.h>

#include "hfenergy.h"
#include "rhf.h"
// #include "rohf.h"
// #include "uhf.h"

using namespace psi;
 
extern FILE *outfile;

HFEnergy::HFEnergy(PSIO *psio, Chkpt *chkpt) : Wavefunction(psio, chkpt)
{
    
}

HFEnergy::HFEnergy(PSIO &psio, Chkpt &chkpt) : Wavefunction(psio, chkpt)
{
    
}

double HFEnergy::compute_energy()
{
    // Check the requested reference in the input file
    char *reference;
    double energy;
    
    ip_string(const_cast<char*>("REFERENCE"), &reference, 0);
    
    if (strcmp(reference, "RHF") == 0) {
        RHF rhf_energy(psio_, chkpt_);
        energy = rhf_energy.compute_energy();
    }
    // else if (strcmp(reference, "ROHF") == 0) {
    //     ROHF rohf_energy(psio_, chkpt_);
    //     energy = rohf_energy.compute_energy();
    // }
    // else if (strcmp(reference, "UHF") == 0) {
    //  UHF uhf_energy(psio_, chkpt_);
    //  energy = uhf_energy.compute_energy();
    // }
    else {
    	fprintf(outfile, "ERROR: Unrecognized reference wavefunction.\n");
    	energy = 0.0;
    }
    // Free the memory allocated by ip_string
    free(reference);
    
    return energy;
}
