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

int Chkpt::rd_max_am(const char *key2)
{
        int max_am;
        char *keyword;
        keyword = build_keyword("Max. AM", key2);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &max_am, sizeof(int));

        free(keyword);
        return max_am;
}

void Chkpt::wt_max_am(int max_am, const char *key2)
{
        char *keyword;
        keyword = build_keyword("Max. AM", key2);

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &max_am, sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_max_am()
** Reads in the maximum orbital quantum number
** of AOs in the basis.
**
** Returns: max_am = (0 corresponds to s-functions,
**                    1 - to up to p-functions, etc.)
** \ingroup CHKPT
*/
        int chkpt_rd_max_am(void)
        {
                int max_am;
                max_am = _default_chkpt_lib_->rd_max_am();
                return max_am;
        }

/*!
** void chkpt_wt_max_am()
** Writes out the maximum orbital quantum number
** of AOs in the basis.
**
** \param max_am = (0 corresponds to s-functions,
**                   1 - to up to p-functions, etc.)
**
** Returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_max_am(int max_am, const char *key2)
        {
                _default_chkpt_lib_->wt_max_am(max_am, key2);
        }
}

