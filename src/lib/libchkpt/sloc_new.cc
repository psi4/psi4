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

int *Chkpt::rd_sloc_new(const char *key2)
{
        int *sloc_new;
        int nshell;
        char *keyword;
        keyword = build_keyword("First BF per shell", key2);

        nshell = rd_nshell(key2);
        sloc_new = array<int>(nshell);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) sloc_new, nshell*sizeof(int));

        free(keyword);
        return sloc_new;
}

void Chkpt::wt_sloc_new(int *sloc_new, const char *key2)
{
        int nshell;
        char *keyword;
        keyword = build_keyword("First BF per shell", key2);

        nshell = rd_nshell(key2);

        psio->write_entry(PSIF_CHKPT, keyword, (char *) sloc_new, nshell*sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int *chkpt_rd_sloc_new()
** Read in an array of the numbers of the first basis
** functions (not AOs as rd_sloc does)  from the shells.
**
** returns:
**   sloc = Read in an array nshell long of the numbers of
**          the first basis functions from the shells.
*/
        int *chkpt_rd_sloc_new(void)
        {
                return _default_chkpt_lib_->rd_sloc_new();
        }

/*!
** void chkpt_wt_sloc_new(int *)
** Writes out an array of the numbers of the first basis
** functions (not AOs as rd_sloc does)  from the shells.
**
** \param sloc = An array nshell long of the numbers of
**               the first basis functions from the shells.
**
** returns: none
*/
        void chkpt_wt_sloc_new(int *sloc_new, const char *key2)
        {
                _default_chkpt_lib_->wt_sloc_new(sloc_new, key2);
        }
}
