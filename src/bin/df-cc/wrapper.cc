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
#include <libiwl/iwl.h>
#include <libqt/qt.h>

#include <libdfcc_solver/ccd.h>
#include <libmints/mints.h>
#include "wrapper.h"

#include <psi4-dec.h>

namespace psi { namespace dfcc {

std::string to_string(const int val); 

PsiReturnType dfcc(Options & options)
{
    tstart();

    Wavefunction::initialize_singletons();

    shared_ptr<PSIO> psio(new PSIO);
    shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));

    timer_init();

    CCD ccd(options, psio, chkpt);
    ccd.compute_energy();

    timer_done();

    tstop();

    return Success;
}

}}
