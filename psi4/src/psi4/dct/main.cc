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

#include "dct.h"

#include <cstdio>
#include <cstdlib>

#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"

using namespace psi;

namespace psi {
namespace dct {

SharedWavefunction dct(SharedWavefunction ref_wfn, Options& options) {
    // Start the timers
    tstart();

    outfile->Printf("\n\n\t***********************************************************************************\n");
    outfile->Printf("\t*                             Density Cumulant Theory                             *\n");
    outfile->Printf("\t*                by Alexander Sokolov, Andy Simmonett, and Xiao Wang              *\n");
    outfile->Printf("\t***********************************************************************************\n");

    auto dct = std::make_shared<DCTSolver>(ref_wfn, options);
    dct->compute_energy();

    // Shut down the timers
    tstop();
    return dct;
}

}  // namespace dct
}  // namespace psi
