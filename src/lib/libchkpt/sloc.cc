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

int *Chkpt::rd_sloc(const char *key2)
{
        int *sloc;
        int nshell;
        char *keyword;
        keyword = build_keyword("First AO per shell", key2);

        nshell = rd_nshell(key2);

        sloc = array<int>(nshell);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) sloc, nshell*sizeof(int));

        free(keyword);
        return sloc;
}

void Chkpt::wt_sloc(int *sloc, const char *key2)
{
        int nshell;
        char *keyword;
        keyword = build_keyword("First AO per shell", key2);

        nshell = rd_nshell(key2);

        psio->write_entry(PSIF_CHKPT, keyword, (char *) sloc, nshell*sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_sloc():	Read in an array of the numbers of the first AO
**			from the shells.
**
**  takes no arguments.
**
**  returns: sloc = An array nshell long of the numbers of
**                  the first AOs from the shells.
**
** \ingroup CHKPT
*/
        int *chkpt_rd_sloc(void)
        {
                return _default_chkpt_lib_->rd_sloc();
        }

/*!
** chkpt_wt_sloc():
**  Writes out an array of the numbers of the first AO from the shells.
**
**  \param sloc = An array nshell long of the numbers of the first AOs
**                from the shells.
**  returns: none
**
** \ingroup CHKPT
*/
        void chkpt_wt_sloc(int *sloc, const char *key2)
        {
                _default_chkpt_lib_->wt_sloc(sloc, key2);
        }
}
