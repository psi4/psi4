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

int *Chkpt::rd_shells_per_am(const char *key2)
{
        int *shells_per_am;
        int max_am;
        char *keyword;
        keyword = build_keyword("Shells per am", key2);

        max_am = rd_max_am(key2);
        shells_per_am = array<int>(max_am+1);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) shells_per_am,
                (max_am+1)*sizeof(int));

        free(keyword);
        return shells_per_am;
}

void Chkpt::wt_shells_per_am(int *shells_per_am, const char *key2)
{
        int max_am;
        char *keyword;
        keyword = build_keyword("Shells per am", key2);

        max_am = rd_max_am(key2);

        psio->write_entry(PSIF_CHKPT, keyword, (char *) shells_per_am,
                (max_am+1)*sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int *chkpt_rd_shells_per_am()
** Reads in the numbers of shells of each angular momentum.
**
** returns: shells_per_am = array of shells per angular momentum
**
** \ingroup CHKPT
*/
        int *chkpt_rd_shells_per_am(void)
        {
                return _default_chkpt_lib_->rd_shells_per_am();
        }

/*!
** void chkpt_wt_shells_per_am(int *)
** Writes out the numbers of shells of each angular momentum.
**
** \param shells_per_am = array of shells per angular momentum
**
** returns: none
**
** \ingroup CHKPT
*/
        void chkpt_wt_shells_per_am(int *shells_per_am, const char *key2)
        {
                _default_chkpt_lib_->wt_shells_per_am(shells_per_am,key2);
        }
}
