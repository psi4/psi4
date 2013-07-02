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
  \ingroup CHKPT
*/

#include <cstdio>
#include <cstdlib>
#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

int *Chkpt::rd_statespi(void)
{
        int nirreps;
        int *statespi;
        char *keyword;
        keyword = build_keyword("States per irrep");

        nirreps = rd_nirreps();
        statespi = array<int>(nirreps);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) statespi,
                nirreps*sizeof(int));

        free(keyword);
        return statespi;
}

void Chkpt::wt_statespi(int *statespi)
{
        int nirreps;
        char *keyword;
        keyword = build_keyword("States per irrep");

        nirreps = rd_nirreps();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) statespi,
                nirreps*sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_statespi():  Reads in the number of excited-states for each
**   irrep.
**
**   takes no arguments.
**
**   returns:
**     int *statespi  an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each
**                 element contains the number of excited states of that
**                 irrep to be studied.
** \ingroup CHKPT
*/
        int *chkpt_rd_statespi(void)
        {
                return _default_chkpt_lib_->rd_statespi();
        }

/*!
** chkpt_wt_statespi():  Writes the number of excited states in each irrep.
**
** \param statespi = an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each
**                 element contains the number of excited states of that
**                 irrep to be studied.
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_statespi(int *statespi)
        {
                _default_chkpt_lib_->wt_statespi(statespi);
        }
}
