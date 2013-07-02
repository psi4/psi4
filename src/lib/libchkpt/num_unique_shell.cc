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

int Chkpt::rd_num_unique_shell(const char *key2)
{
        int nunique;
        char *keyword;
        keyword = build_keyword("Num. unique shells", key2);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &nunique, sizeof(int));

        free(keyword);
        return nunique;
}

void Chkpt::wt_num_unique_shell(int nunique, const char *key2)
{
        char *keyword;
        keyword = build_keyword("Num. unique shells", key2);

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &nunique, sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_num_unique_shell()
** Reads in the number of symmetry unique shells.
**
** returns: nunique = number of symmetry unique shells.
** \ingroup CHKPT
*/
        int chkpt_rd_num_unique_shell(void)
        {
                return _default_chkpt_lib_->rd_num_unique_shell();
        }

/*!
** void chkpt_wt_num_unique_shell(int)
** Writes out the number of symmetry unique shells.
**
** \param nunique = number of symmetry unique shells.
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_num_unique_shell(int nunique, const char *key2)
        {
                _default_chkpt_lib_->wt_num_unique_shell(nunique, key2);
        }
}
