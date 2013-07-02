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

int *Chkpt::rd_us2s(const char *key2)
{
        int *us2s;
        int num_unique_shells;
        char *keyword;
        keyword = build_keyword("Unique shell -> full shell map", key2);

        num_unique_shells = rd_num_unique_shell(key2);
        us2s = array<int>(num_unique_shells);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) us2s, num_unique_shells*sizeof(int));

        free(keyword);
        return us2s;
}

void Chkpt::wt_us2s(int *us2s, const char *key2)
{
        int num_unique_shells;
        char *keyword;
        keyword = build_keyword("Unique shell -> full shell map", key2);

        num_unique_shells = rd_num_unique_shell(key2);

        psio->write_entry(PSIF_CHKPT, keyword, (char *) us2s, num_unique_shells*sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int *chkpt_rd_us2s()
** Read in a mapping array betwen unique shell and
** full shell lists
**
** returns: us2s = Read in an array num_unique_shell
**
** \ingroup CHKPT
*/
        int *chkpt_rd_us2s(void)
        {
                return _default_chkpt_lib_->rd_us2s();
        }

/*!
** void chkpt_wt_us2s(int *)
** Writes out a mapping array betwen unique shell and
** full shell lists.
**
** \param us2s = An array num_unique_shell
**
** returns: none
**
** \ingroup CHKPT
*/
        void chkpt_wt_us2s(int *us2s, const char *key2)
        {
                _default_chkpt_lib_->wt_us2s(us2s, key2);
        }
}
