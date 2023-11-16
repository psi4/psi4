/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>

#include "psi4/psifiles.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/exception.h"

#include "psi4/psi4-dec.h"

#include "mp2f12.h"

namespace psi { namespace mp2f12 {

SharedWavefunction mp2f12(SharedWavefunction ref_wfn, Options& options)
{

    std::shared_ptr<Wavefunction> mp2f12;
    if (options.get_str("F12_TYPE").find("DISK") != std::string::npos) {
        mp2f12 = std::make_shared<DiskMP2F12>(ref_wfn, options, psio);
    } else if {
        mp2f12 = std::make_shared<MP2F12>(ref_wfn, options, psio);
    } else {
        throw PSIEXCEPTION("MP2F12: Unrecognized F12_TYPE.");
    }

    return mp2f12;
}

}} // End namespaces
