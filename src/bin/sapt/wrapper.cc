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
    
    // Shut down psi. 
    timer_done();

    tstop();
    
    return Success;
}

}}
