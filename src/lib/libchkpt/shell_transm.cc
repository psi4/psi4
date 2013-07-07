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

int **Chkpt::rd_shell_transm(const char *key2)
{
        int i, nshell, nirreps;
        int **shell_transm;
        psio_address ptr;
        char *keyword;
        keyword = build_keyword("Shell transmat", key2);

        nshell = rd_nshell(key2);
        nirreps = rd_nirreps();

        shell_transm = matrix<int>(nshell,nirreps);
        ptr = PSIO_ZERO;
        for(i=0; i < nshell; i++)
                psio->read(PSIF_CHKPT, keyword, (char *) shell_transm[i],
                nirreps*sizeof(int), ptr, &ptr);

        free(keyword);
        return shell_transm;
}

void Chkpt::wt_shell_transm(int **shell_transm, const char *key2)
{
        int i, nshell, nirreps;
        psio_address ptr;
        char *keyword;
        keyword = build_keyword("Shell transmat", key2);

        nshell = rd_nshell(key2);
        nirreps = rd_nirreps();

        ptr = PSIO_ZERO;
        for(i=0; i < nshell; i++) {
                psio->write(PSIF_CHKPT, keyword, (char *) shell_transm[i],
                        nirreps*sizeof(int), ptr, &ptr);
        }

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_shell_transm():	Read in a matrix of nshell*nirreps integers
**			        that contains symmetry information.
**
**  takes no arguments.
**
**  returns:
**    shell_transm = matrix of nshell*nirrpes ints w/ symmetry info
**
** \ingroup CHKPT
*/
        int **chkpt_rd_shell_transm(void)
        {
                return _default_chkpt_lib_->rd_shell_transm();
        }

/*!
** chkpt_wt_shell_transm():	Write out a matrix of nshell*nirreps integers
**			        that contains symmetry information.
**
** \param shell_transm = matrix of nshell*nirreps ints w/ symmetry info
**
** returns: none
**
** \ingroup CHKPT
*/
        void chkpt_wt_shell_transm(int **shell_transm, const char *key2)
        {
                _default_chkpt_lib_->wt_shell_transm(shell_transm, key2);
        }
}
