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

int *Chkpt::rd_sprim(const char *key2)
{
        int *sprim;
        int nshell;
        char *keyword;
        keyword = build_keyword("First primitive per shell", key2);

        nshell = rd_nshell(key2);

        sprim = array<int>(nshell);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) sprim, nshell*sizeof(int));

        free(keyword);
        return sprim;
}

void Chkpt::wt_sprim(int *sprim, const char *key2)
{
        int nshell;
        char *keyword;
        keyword = build_keyword("First primitive per shell", key2);

        nshell = rd_nshell(key2);

        psio->write_entry(PSIF_CHKPT, keyword, (char *) sprim, nshell*sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_sprim(): Reads in array of the numbers of first primitives
**                   from the shells.
**
**  takes no arguments.
**
**  returns: sprim = an array of the numbers of first primitives
**           from the shells.
**
** \ingroup CHKPT
*/
        int *chkpt_rd_sprim(void)
        {
                return _default_chkpt_lib_->rd_sprim();
        }

/*!
** chkpt_wt_sprim():	Writes out an array of the numbers of first primitives
**			from the shells.
**
**  \param sprim = an array of the numbers of first primitives
**                 from the shells.
**
**  returns: none
**
** \ingroup CHKPT
*/
        void chkpt_wt_sprim(int *sprim, const char *key2)
        {
                _default_chkpt_lib_->wt_sprim(sprim, key2);
        }
}
