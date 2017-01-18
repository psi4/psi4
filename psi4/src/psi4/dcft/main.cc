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

#include <stdio.h>
#include <stdlib.h>
#include "psi4/psi4-dec.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.h"
#include "defines.h"
#include "dcft.h"

using namespace psi;


namespace psi{ namespace dcft{

SharedWavefunction dcft(SharedWavefunction ref_wfn, Options& options)
{
    // Start the timers
    tstart();

    outfile->Printf("\n\n\t***********************************************************************************\n");
    outfile->Printf(    "\t*                        Density Cumulant Functional Theory                       *\n");
    outfile->Printf(    "\t*                by Alexander Sokolov, Andy Simmonett, and Xiao Wang              *\n");
    outfile->Printf(    "\t***********************************************************************************\n");

    SharedWavefunction dcft = SharedWavefunction(new DCFTSolver(ref_wfn, options));
    dcft->compute_energy();

    // Shut down the timers
    tstop();
    return dcft;
}

}} // End Namespaces
