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

int *Chkpt::rd_snuc(const char *key2)
{
        int *snuc;
        int nshell;
        char *keyword;
        keyword = build_keyword("Shell nucleus", key2);

        nshell = rd_nshell(key2);

        snuc = array<int>(nshell);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) snuc, nshell*sizeof(int));

        free(keyword);
        return snuc;
}

void Chkpt::wt_snuc(int *snuc, const char *key2)
{
        int nshell;
        char *keyword;
        keyword = build_keyword("Shell nucleus", key2);

        nshell = rd_nshell(key2);

        psio->write_entry(PSIF_CHKPT, keyword, (char *) snuc, nshell*sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_snuc(): Reads in array of the nuclei numbers shells belong to.
**
**  takes no arguments.
**
**  returns: snuc = an array of the nuclei numbers to which shells
**                  belong to.
**
** \ingroup CHKPT
*/
        int *chkpt_rd_snuc(void)
        {
                return _default_chkpt_lib_->rd_snuc();
        }

/*!
** chkpt_wt_snuc(): Writes out array of the nuclei numbers shells belong to.
**
**  \param snuc = an array of the nuclei numbers to which shells belong to
**
**  returns: none
**
** \ingroup CHKPT
*/
        void chkpt_wt_snuc(int *snuc, const char *key2)
        {
                _default_chkpt_lib_->wt_snuc(snuc, key2);
        }
}
