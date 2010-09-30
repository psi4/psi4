#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>

#include <libsapt_solver/sapt0.h>
#include <libsapt_solver/sapt_dft.h>
#include <libsapt_solver/sapt2.h>
#include <libsapt_solver/sapt2p.h>
#include <libsapt_solver/sapt2p3.h>
#include <libmints/mints.h>
#include "wrapper.h"

#include <psi4-dec.h>

namespace psi { namespace sapt {

std::string to_string(const int val);   // In matrix.cpp

PsiReturnType sapt(Options & options)
{
    tstart();
    
    Wavefunction::initialize_singletons();
    
    shared_ptr<PSIO> psio(new PSIO);
    psiopp_ipv1_config(psio);
    
    shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));
    
    // Initialize the psi3 timer library.
    timer_init();

    if (options.get_str("SAPT_LEVEL") == "SAPT0") {
        SAPT0 sapt(options, psio, chkpt);
        sapt.compute_energy();
    }

    if (options.get_str("SAPT_LEVEL") == "SAPT_DFT") {
        fprintf(outfile,"  SAPT(DFT) is not currently implemented\n");
//      SAPT_DFT sapt(options, psio, chkpt);
//      sapt.compute_energy();
    }

    if (options.get_str("SAPT_LEVEL") == "SAPT2") {
        SAPT2 sapt(options, psio, chkpt);
        sapt.compute_energy();
    }
    
    if (options.get_str("SAPT_LEVEL") == "SAPT2+") {
        SAPT2p sapt(options, psio, chkpt);
        sapt.compute_energy();
    }

    if (options.get_str("SAPT_LEVEL") == "SAPT2+3") {
        SAPT2p3 sapt(options, psio, chkpt);
        sapt.compute_energy();
    }
    
    // Shut down psi. 
    timer_done();

    tstop();
    
    return Success;
}

}}
