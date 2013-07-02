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

int Chkpt::rd_disp_irrep()
{
        int h;
        char *keyword;
        keyword = build_keyword("Current Displacement Irrep");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &h, sizeof(int));

        free(keyword);
        return h;
}

void Chkpt::wt_disp_irrep(int disp_irrep)
{
        char *keyword;
        keyword = build_keyword("Current Displacement Irrep");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &disp_irrep, sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_disp_irrep()
**
** Reads in the irrep of the current displaced geometry assuming
** Cotton ordering of irreps - to be used by input to determine
** docc and socc
**
** Returns:
**   disp_irrep = irrep of current displaced geometry
** \ingroup CHKPT
*/
        int chkpt_rd_disp_irrep(void)
        {
                int h;
                h = _default_chkpt_lib_->rd_disp_irrep();
                return h;
        }

/*!
** void chkpt_wt_disp_irrep()
**
** Writes the irrep of the current displaced geometry assuming
** Cotton ordering of irreps - to be used by input to determine
** docc and socc
**
** \param disp_irrep = irrep of current displaced geometry
**
** Returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_disp_irrep(int disp_irrep)
        {
                _default_chkpt_lib_->wt_disp_irrep(disp_irrep);
        }
}

