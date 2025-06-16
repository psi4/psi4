/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "apps.h"

#include "psi4/psifiles.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.h"
#include "psi4/libqt/qt.h"

#include "psi4/psi4-dec.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>

namespace psi {
namespace libfock {

SharedWavefunction libfock(SharedWavefunction ref_wfn, Options& options) {
    tstart();

    std::shared_ptr<RBase> wfn;

    if (options.get_str("MODULE") == "RCPHF") {
        RCPHF* cphf(new RCPHF(ref_wfn, options));
        for (size_t i = 0; i < options["CPHF_TASKS"].size(); i++) {
            cphf->add_task(options["CPHF_TASKS"][i].to_string());
        }
        wfn = std::shared_ptr<RBase>(cphf);
    } else {
        throw PSIEXCEPTION("Libfock: Applications module not recognized");
    }

    //    wfn->compute_energy();

    tstop();

    return wfn;
    //    return Success;
}
}  // namespace libfock
}  // namespace psi
