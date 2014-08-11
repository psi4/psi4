/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#include <stdio.h>
#include <stdlib.h>
#include "psi4-dec.h"
#include "libqt/qt.h"
#include "libchkpt/chkpt.hpp"
#include "libpsio/psio.h"
#include "libchkpt/chkpt.h"
#include "defines.h"
#include "dcft.h"

using namespace psi;
using namespace boost;

namespace psi{ namespace dcft{

PsiReturnType
dcft(Options &options)
{

    // Start the timers
    tstart();

    outfile->Printf("\n\n\t***********************************************************************************\n");
    outfile->Printf(    "\t*                        Density Cumulant Functional Theory                       *\n");
    outfile->Printf(    "\t*                     by Alexander Sokolov and Andy Simmonett                     *\n");
    outfile->Printf(    "\t***********************************************************************************");

    boost::shared_ptr<Wavefunction> dcft = boost::shared_ptr<Wavefunction>(new DCFTSolver(Process::environment.wavefunction(), options));
    Process::environment.set_wavefunction(dcft);

    dcft->compute_energy();

    // Shut down the timers
    tstop();

    return Success;
}

}} // End Namespaces
