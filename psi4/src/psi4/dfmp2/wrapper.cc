/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
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
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.h"
#include "psi4/libqt/qt.h"

#include "psi4/psi4-dec.h"

#include "mp2.h"



namespace psi { namespace dfmp2 {

SharedWavefunction dfmp2(SharedWavefunction ref_wfn, Options & options)
{

    std::shared_ptr<PSIO> psio(new PSIO);

    std::shared_ptr<Wavefunction> dfmp2;
    if (options.get_str("REFERENCE") == "RHF" || options.get_str("REFERENCE") == "RKS") {
        dfmp2 = std::shared_ptr<Wavefunction>(new RDFMP2(ref_wfn, options, psio));
    } else if (options.get_str("REFERENCE") == "UHF" || options.get_str("REFERENCE") == "UKS") {
        dfmp2 = std::shared_ptr<Wavefunction>(new UDFMP2(ref_wfn, options, psio));
    } else if (options.get_str("REFERENCE") == "ROHF") {
        dfmp2 = std::shared_ptr<Wavefunction>(new RODFMP2(ref_wfn, options, psio));
    } else {
        throw PSIEXCEPTION("DFMP2: Unrecognized reference");
    }

    return dfmp2;
}


}}
