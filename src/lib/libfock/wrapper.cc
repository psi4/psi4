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

    boost::shared_ptr<RBase> wfn;

    if (options.get_str("MODULE") == "RCPHF") {
        RCPHF* cphf(new RCPHF());
        for (int i = 0; i < options["CPHF_TASKS"].size(); i++) {
            cphf->add_task(options["CPHF_TASKS"][i].to_string());
        }
        wfn = boost::shared_ptr<RBase>(cphf);
    } else if (options.get_str("MODULE") == "RCIS") {
        wfn = boost::shared_ptr<RBase>(new RCIS());
    } else if (options.get_str("MODULE") == "RTDHF") {
        wfn = boost::shared_ptr<RBase>(new RTDHF());
    } else if (options.get_str("MODULE") == "RCPKS") {
        RCPKS* cphf(new RCPKS());
        for (int i = 0; i < options["CPHF_TASKS"].size(); i++) {
            cphf->add_task(options["CPHF_TASKS"][i].to_string());
        }
        wfn = boost::shared_ptr<RBase>(cphf);
    } else if (options.get_str("MODULE") == "RTDA") {
        wfn = boost::shared_ptr<RBase>(new RTDA());
    } else if (options.get_str("MODULE") == "RTDDFT") {
        wfn = boost::shared_ptr<RBase>(new RTDDFT());
    } else {
        throw PSIEXCEPTION("Libfock: Applications module not recognized");
    }

    wfn->compute_energy();

    tstop();

    return Success;
}

}}

