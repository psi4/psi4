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

/*!
** \file
** \brief Get the amount of core memory available from input
** \ingroup CIOMR
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libparallel/parallel.h>
#include <psi4-dec.h>

namespace psi {

#define DEF_MAXCRR (256000000)  // default maxcor 256 M bytes

void set_memory(std::string OutFileRMR)
{
    long int maxcrr;

    maxcrr = DEF_MAXCRR; // set to default

    if (maxcrr < 1e9) {
       outfile->Printf("    Memory level set to %.3lf MB\n", maxcrr / 1e6 );
    }
    else {
        outfile->Printf("    Memory level set to %.3lf GB\n", maxcrr / 1e9 );
    }

    Process::environment.set_memory(maxcrr);

    return;
}

}
