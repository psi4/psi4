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

int Chkpt::rd_disp(void)
{
        int disp;
        char *keyword;
        keyword = build_keyword("Current displacement");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &disp,
                sizeof(int));

        free(keyword);
        return disp;
}

void Chkpt::wt_disp(int disp)
{
        char *keyword;
        keyword = build_keyword("Current displacement");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &disp,
                sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_disp():  Reads in the current geometry displacement number.
**
** Parameters: none
**
** not used by OPTKING; used by anybody else ???
**
** Returns: int disp, the current geometry displacement number
** \ingroup CHKPT
*/
        int chkpt_rd_disp(void)
        {
                int disp;
                disp = _default_chkpt_lib_->rd_disp();
                return disp;
        }


/*!
** chkpt_wt_disp():  Writes out the current geometry displacement number.
**
** \param int disp   the current geometry displacement number
**
** Returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_disp(int disp)
        {
                _default_chkpt_lib_->wt_disp(disp);
        }
}

