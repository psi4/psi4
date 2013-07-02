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

#include <cstdlib>
#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

int Chkpt::rd_phase_check(void)
{
        int pcheck;
        char *keyword;
        keyword = build_keyword("Phase check");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &pcheck, sizeof(int));

        free(keyword);
        return pcheck;
}

void Chkpt::wt_phase_check(int pcheck)
{
        char *keyword;
        keyword = build_keyword("Phase check");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &pcheck, sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_phase_check()
**
** Reads a boolean flag indicating whether the SCF code was able to correct
** the phases of the molecular orbitals relative to the guess orbitals.  This
** is important for restarting correlated wfn calculations from earlier vectors.
**
** arguments: none
**
** returns: pcheck = Phase check flag (1 if phase has been checked, else 0)
**
** \ingroup CHKPT
*/
        int chkpt_rd_phase_check(void)
        {
                return _default_chkpt_lib_->rd_phase_check();
        }

/*!
** void chkpt_wt_phase_check(int)
**
** Reads a boolean flag indicating whether the SCF code was able to correct
** the phases of the molecular orbitals relative to the guess orbitals.  This
** is important for restarting correlated wfn calculations from earlier vectors.
**
** \param pcheck = Phase check flag (1 if phase has been checked, else 0)
**
** returns: none
**
** \ingroup CHKPT
*/
        void chkpt_wt_phase_check(int pcheck)
        {
                _default_chkpt_lib_->wt_phase_check(pcheck);
        }
}
