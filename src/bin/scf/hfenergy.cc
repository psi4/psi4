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

    shared_ptr<Wavefunction> scf;

    if (reference == "RHF") {
        scf = shared_ptr<Wavefunction>(new RHF(options_, psio_, chkpt_));
    }
    else if (reference == "ROHF") {
        scf = shared_ptr<Wavefunction>(new ROHF(options_, psio_, chkpt_));
    }
    else if (reference == "UHF") {
        scf = shared_ptr<Wavefunction>(new UHF(options_, psio_, chkpt_));
    }
    else if (reference == "RKS") {
        scf = shared_ptr<Wavefunction>(new RKS(options_, psio_, chkpt_));
    }
    else if (reference == "UKS") {
        scf = shared_ptr<Wavefunction>(new UKS(options_, psio_, chkpt_));
    }
    else {
        throw InputException("Unknown reference " + reference, "REFERENCE", __FILE__, __LINE__);
        energy = 0.0;
    }
    energy = scf->compute_energy();
    
    Process::environment.set_reference_wavefunction(scf);

    return energy;
}

#if HAVE_MPI == 1
    double HFEnergy::compute_parallel_energy() {

        // Check the requested reference in the input file
        string reference, type;
        double energy;

        type = options_.get_str("SCF_TYPE");
        reference = options_.get_str("REFERENCE");

        Communicator::world->print(outfile);

        if(type == "DIRECT") {
            if (reference == "RHF") {
                RHF rhf_energy(options_, psio_, chkpt_);
                energy = rhf_energy.compute_energy_parallel();
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
