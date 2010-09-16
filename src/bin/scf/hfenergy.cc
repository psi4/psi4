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
#include <libparallel/parallel.h>
#include "libscf_solver/rhf.cc"
#include "libscf_solver/rohf.h" 
#include "libscf_solver/uhf.h"
#include "libscf_solver/rks.h"
#include "libscf_solver/uks.h"

using namespace std;
using namespace psi;

namespace psi { namespace scf {
     
HFEnergy::HFEnergy(Options & options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt) 
    : options_(options), psio_(psio), chkpt_(chkpt)
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

#if HAVE_MADNESS == 1
    double HFEnergy::compute_parallel_energy()
    {
        // Check the requested reference in the input file
        string reference, type;
        double energy;

        type = options_.get_str("SCF_TYPE");
        reference = options_.get_str("REFERENCE");

        if(Communicator::world->me() == 0)
            fprintf(outfile, "\n Running in parallel with the MADNESS runtime library\n\n");
        if(type == "DIRECT") {
            if (reference == "RHF") {
                RHF rhf_energy(options_, psio_, chkpt_);
                energy = rhf_energy.compute_energy();
            }
            else {
                throw InputException("Parallel SCF only works for RHF reference" , "REFERENCE", __FILE__, __LINE__);
                energy = 0.0;
            }
        }
        else {
            throw InputException("Parallel SCF is direct only. Please set SCF_TYPE=direct in input", "SCF_TYPE", __FILE__, __LINE__);
            energy = 0.0;
        }


        return energy;
    }
#endif

}}
