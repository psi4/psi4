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

/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/

#include "MOInfo.h"
#include "psi4/cc/ccwave.h"

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <cstdio>

namespace psi {
namespace ccenergy {

void CCEnergyWavefunction::update() {
    outfile->Printf("  %4d      %20.15f    %4.3e    %7.6f    %7.6f    %7.6f    %7.6f\n", moinfo_.iter, moinfo_.ecc,
                    moinfo_.conv, moinfo_.t1diag, moinfo_.d1diag, moinfo_.new_d1diag, moinfo_.d2diag);
}
}  // namespace ccenergy
}  // namespace psi
