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

#include <libmints/mints.h>
#include "hfenergy.h"

#include <psi4-dec.h>

namespace psi { namespace scf {

PsiReturnType scf(Options & options)
{
    tstart();

    Wavefunction::initialize_singletons();

    shared_ptr<PSIO> psio(new PSIO);
    psiopp_ipv1_config(psio);

    shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));

    // Initialize the psi3 timer library.
    timer_init();

    // Compute the Hartree-Fock energy
    HFEnergy hf(options, psio, chkpt);
    double energy = hf.compute_energy();

    Process::environment.globals["SCF ENERGY"] = energy;
    Process::environment.globals["CURRENT ENERGY"] = energy;

    // Shut down psi.
    timer_done();

    tstop();

    return Success;
}

}}
