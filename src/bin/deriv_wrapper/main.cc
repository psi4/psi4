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

#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libmints/mints.h>
#include <libmints/cdsalclist.h>
#include <libmints/deriv.h>

#include <psi4-dec.h>

using namespace boost;

namespace psi { namespace deriv {

PsiReturnType deriv(Options &)
{
    tstart();

    outfile->Printf( " DERIV: Derivative code.\n   by Justin Turney\n\n");

    Deriv test(Process::environment.wavefunction(),
               0x1,
               false,
               false);
    test.compute();

    // Shut down psi
    tstop();

    return Success;
}

}} // namespaces
