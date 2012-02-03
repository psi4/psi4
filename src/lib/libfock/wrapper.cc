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
#include <psi4-dec.h>

#include "apps.h"

using namespace boost;

namespace psi { 
namespace libfock {

PsiReturnType libfock(Options & options)
{
    tstart();

    if (options.get_str("MODULE") == "RCPHF") {
        boost::shared_ptr<RCPHF> rcphf(new RCPHF());
        for (int i = 0; i < options["CPHF_TASKS"].size(); i++) {
            rcphf->add_task(options["CPHF_TASKS"][i].to_string());
        }
        rcphf->compute_energy();
    }

    tstop();

    return Success;
}

}}

