/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>

#include "psi4/psifiles.h"
#include "psi4/src/lib/libciomr/libciomr.h"
#include "psi4/src/lib/libpsio/psio.h"
#include "psi4/src/lib/libpsio/psio.hpp"
#include "psi4/src/lib/libiwl/iwl.h"
#include "psi4/src/lib/libqt/qt.h"

#include "psi4/psi4-dec.h"

#include "apps.h"

using namespace boost;

namespace psi {
namespace libfock {

SharedWavefunction libfock(SharedWavefunction ref_wfn, Options & options)
{
    tstart();

    boost::shared_ptr<RBase> wfn;

    if (options.get_str("MODULE") == "RCPHF") {
        RCPHF* cphf(new RCPHF(ref_wfn, options));
        for (unsigned int i = 0; i < options["CPHF_TASKS"].size(); i++) {
            cphf->add_task(options["CPHF_TASKS"][i].to_string());
        }
        wfn = boost::shared_ptr<RBase>(cphf);
    } else if (options.get_str("MODULE") == "RCIS") {
        wfn = boost::shared_ptr<RBase>(new RCIS(ref_wfn, options));
    } else if (options.get_str("MODULE") == "RTDHF") {
        wfn = boost::shared_ptr<RBase>(new RTDHF(ref_wfn, options));
    } else if (options.get_str("MODULE") == "RCPKS") {
        RCPKS* cphf(new RCPKS(ref_wfn, options));
        for (unsigned int i = 0; i < options["CPHF_TASKS"].size(); i++) {
            cphf->add_task(options["CPHF_TASKS"][i].to_string());
        }
        wfn = boost::shared_ptr<RBase>(cphf);
    } else if (options.get_str("MODULE") == "RTDA") {
        wfn = boost::shared_ptr<RBase>(new RTDA(ref_wfn, options));
    } else if (options.get_str("MODULE") == "RTDDFT") {
        wfn = boost::shared_ptr<RBase>(new RTDDFT(ref_wfn, options));
    } else {
        throw PSIEXCEPTION("Libfock: Applications module not recognized");
    }

//    wfn->compute_energy();

    tstop();

    return wfn;
//    return Success;
}

}}
