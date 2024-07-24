/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#include "dlpno.h"

#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/exception.h"

namespace psi {
namespace dlpno {

SharedWavefunction dlpno(SharedWavefunction ref_wfn, Options& options) {

    std::shared_ptr<Wavefunction> dlpno;
    if (options.get_str("REFERENCE") == "RHF") {
        if (options.get_str("DLPNO_ALGORITHM") == "MP2") {
            dlpno = std::make_shared<DLPNOMP2>(ref_wfn, options);
        } else {
            throw PSIEXCEPTION("Requested DLPNO method is not yet available!");
        }
    } else {
        throw PSIEXCEPTION("DLPNO requires closed-shell reference"); 
    }

    return dlpno;
}
}  // namespace dlpno
}  // namespace psi
