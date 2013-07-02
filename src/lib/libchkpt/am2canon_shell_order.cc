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

/*! \defgroup CHKPT libchkpt: The Checkpoint Interface library */

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

/*!
** int *Chkpt::rd_am2canon_shell_order()
**
** Reads in the mapping array from the am-ordered
** to the canonical (in the order of appearance)
** list of shells.
**
** Returns: int *am2can_shell_order
** \ingroup CHKPT
*/
int *Chkpt::rd_am2canon_shell_order(const char *key2)
{
        int *am2can_sh_ord, nshell;
        char *keyword;
        keyword = build_keyword("AM -> canonical shell map", key2);

        nshell = rd_nshell(key2);
        am2can_sh_ord = array<int>(nshell);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) am2can_sh_ord,
                nshell*sizeof(int));

        free(keyword);
        return am2can_sh_ord;
}

/*!
** void Chkpt::wt_am2canon_shell_order()
**
** Writes out the mapping array from the am-ordered
** to the canonical (in the order of appearance)
** list of shells.
**
** \param am2can_shell_order = array to store the mapping array
**
** Returns: none
** \ingroup CHKPT
*/
void Chkpt::wt_am2canon_shell_order(int *am2can_sh_ord, const char *key2)
{
        int nshell;
        char *keyword;
        keyword = build_keyword("AM -> canonical shell map", key2);

        nshell = rd_nshell(key2);

        psio->write_entry(PSIF_CHKPT, keyword, (char *) am2can_sh_ord,
                nshell*sizeof(int));

  free(keyword);
}

extern "C" {
/*!
** int *chkpt_rd_am2canon_shell_order()
**
** Reads in the mapping array from the am-ordered
** to the canonical (in the order of appearance)
** list of shells.
**
** Returns: int *am2can_shell_order
** \ingroup CHKPT
*/
        int *chkpt_rd_am2canon_shell_order(void)
        {
                return _default_chkpt_lib_->rd_am2canon_shell_order();
        }

/*!
** void chkpt_wt_am2canon_shell_order()
**
** Writes out the mapping array from the am-ordered
** to the canonical (in the order of appearance)
** list of shells.
**
** \param am2can_shell_order = array to store the mapping array
**
** Returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_am2canon_shell_order(int *am2can_sh_ord, const char *key2)
        {
                _default_chkpt_lib_->wt_am2canon_shell_order(am2can_sh_ord,
                  key2);
        }
}

