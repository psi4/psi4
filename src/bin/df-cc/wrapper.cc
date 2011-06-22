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
#include <libdfcc_solver/mp2.h>
#include <libdfcc_solver/mp3.h>
#include <libdfcc_solver/drpa.h>
#include <libmints/mints.h>
#include "wrapper.h"

#include <psi4-dec.h>

using namespace boost;

namespace psi { namespace dfcc {

std::string to_string(const int val);

PsiReturnType dfcc(Options & options)
{
    tstart();

    Wavefunction::initialize_singletons();

    boost::shared_ptr<PSIO> psio(new PSIO);
    boost::shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));

    if (options.get_str("WAVEFUNCTION") == "CCD") {

    CCD ccd(options, psio, chkpt);
    ccd.compute_energy();

    } else if (options.get_str("WAVEFUNCTION") == "MP2") {

    MP2 mp2(options, psio, chkpt);
    mp2.compute_energy();

    } else if (options.get_str("WAVEFUNCTION") == "MP3") {

    MP3 mp3(options, psio, chkpt);
    mp3.compute_energy();

    } else if (options.get_str("WAVEFUNCTION") == "DRPA") {

    dRPA drpa(options, psio, chkpt);
    drpa.compute_energy();

    }

    tstop();

    return Success;
}

}}
