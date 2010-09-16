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
#include <libparallel/parallel.h>
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

    shared_ptr<PSIO> psio(new PSIO);
    psiopp_ipv1_config(psio);

    shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));

    // Initialize the psi3 timer library.
    timer_init();

    double energy;
    bool parallel;

#if HAVE_MADNESS == 1
    int use_madness;
    std::string mad = Process::environment("MADNESS");
    if ((mad == "true") || (mad == "True") || (mad == "TRUE")) {
        use_madness = 1;
    }
    else {
        use_madness = 0;
    }

    parallel = options.get_bool("PARALLEL");


    if (parallel) {
        if(use_madness) {
            // Compute the Hartree-Fock energy
            HFEnergy hf(options, psio, chkpt);
            energy = hf.compute_parallel_energy();
        }
        else {
            throw InputException("If you want to run SCF in parallel please set MADNESS=true in environment." , "PARALLEL", __FILE__, __LINE__);
        }
    }
    else {
#endif
#if HAVE_MPI == 1
        if(Communicator::world->me() == 0) {
#endif

            // Compute the Hartree-Fock energy
            HFEnergy hf(options, psio, chkpt);
            energy = hf.compute_energy();

#if HAVE_MPI == 1
        }
#endif
#if HAVE_MADNESS == 1
    }
#endif

    
    Process::environment.globals["SCF ENERGY"] = energy;
    Process::environment.globals["CURRENT ENERGY"] = energy;

    // Shut down psi.
    timer_done();

    tstop();

    return Success;
}

}}
