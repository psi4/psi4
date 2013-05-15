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

int Chkpt::rd_nallatom(void)
{
        int num_allatoms;
        char *keyword;
        keyword = build_keyword("Num. all atoms");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &num_allatoms,
          sizeof(int));

        free(keyword);
        return num_allatoms;
}

void Chkpt::wt_nallatom(int num_allatoms)
{
        char *keyword;
        keyword = build_keyword("Num. all atoms");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &num_allatoms,
          sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_nallatom()
**
** Reads number of all atoms (including dummy atoms)
**
** Parameters: none
**
** Returns:
**   nallatom = number of all atoms (including dummies).
** \ingroup CHKPT
*/
        int chkpt_rd_nallatom(void)
        {
                return _default_chkpt_lib_->rd_nallatom();
        }


/*!
** chkpt_wt_nallatom()
**
** Writes the number of all atoms (including dummy atoms)
**
** Parameters:
**   \param nallatom  = number of all atoms (including dummies).
**
** Returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_nallatom(int num_allatoms)
        {
                _default_chkpt_lib_->wt_nallatom(num_allatoms);
        }
}

