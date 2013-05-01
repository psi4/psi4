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

/*
 *
 *
 */

#include <cstdio>
#include <cstdlib>
#include <string>

#include "dfmp2_energy.h"
#include "dfmp2.h"

using namespace boost;
using namespace psi;

namespace psi { namespace dfmp2 {

DFMP2Energy::DFMP2Energy(Options & options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt)
    : Wavefunction(options, psio, chkpt)
{

}

double DFMP2Energy::compute_E()
{
    double energy;

    psi::dfmp2::DFMP2 dfmp2_energy(options_, psio_, chkpt_);
    energy = dfmp2_energy.compute_energy();

    return energy;
}

}}
