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

#include <libmints/mints.h>
#include "dfmp2.h"

#include <psi4-dec.h>

using namespace boost;

namespace psi { namespace dfmp2 {

std::string to_string(const int val);   // In matrix.cpp

PsiReturnType dfmp2(Options & options)
{
    tstart();

    shared_ptr<PSIO> psio(new PSIO);
//    psiopp_ipv1_config(psio);

    shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));

    // Initialize the psi3 timer library.
    timer_init();

    DFMP2 df(options, psio, chkpt);
    df.compute_energy();

    // Shut down psi.
    timer_done();

    tstop();

    return Success;
}

}}
