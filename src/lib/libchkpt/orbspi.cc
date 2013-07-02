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
  \file
*/

#include <cstdio>
#include <cstdlib>
#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

int *Chkpt::rd_orbspi(void)
{
        int nirreps;
        int *orbspi;
        char *keyword;
        keyword = build_keyword("MO's per irrep");

        nirreps = rd_nirreps();
        orbspi  = array<int>(nirreps);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) orbspi, nirreps*sizeof(int));

        free(keyword);
        return orbspi;
}

void Chkpt::wt_orbspi(int *orbspi)
{
        int nirreps;
        char *keyword;
        keyword = build_keyword("MO's per irrep");

        nirreps = rd_nirreps();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) orbspi, nirreps*sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_orbspi():  Reads in the number of molecular orbitals in each irrep.
**
**   takes no arguments.
**
**   returns:
**     int *orbspi  an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each
**                 element contains the number of molecular orbitals for
**                 that irrep. Also, see chkpt_rd_sopi().
*/
        int *chkpt_rd_orbspi(void)
        {
                return _default_chkpt_lib_->rd_orbspi();
        }

/*!
** chkpt_wt_orbspi():  Writes the number of molecular orbitals in each irrep.
**
** \param orbspi = an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each
**                 element contains the number of molecular orbitals for
**                 that irrep. Also, see chkpt_rd_sopi().
**
** returns: none
*/
        void chkpt_wt_orbspi(int *orbspi)
        {
                _default_chkpt_lib_->wt_orbspi(orbspi);
        }
}
