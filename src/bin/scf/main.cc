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
    bool parallel = options.get_bool("PARALLEL");

    if(parallel) {
#if HAVE_MPI == 1
        HFEnergy hf(options, psio, chkpt);
        // Compute the Hartree-Fock energy in parallel
        energy = hf.compute_parallel_energy();
#else
        throw InputException("Parallel SCF requires MPI or MADNESS. If you want to run in parallel, please install one or both." , "PARALLEL", __FILE__, __LINE__);
#endif
    }
    else {
            HFEnergy hf(options, psio, chkpt); 
            // Compute the Hartree-Fock energy
            energy = hf.compute_energy();
    }
    
    Process::environment.globals["SCF ENERGY"] = energy;
    Process::environment.globals["CURRENT ENERGY"] = energy;

    // Shut down psi.
    timer_done();

    tstop();

    return Success;
}

}}
