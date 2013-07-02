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

int Chkpt::rd_nfzv(void)
{
        int nfzv;
        char *keyword;
        keyword = build_keyword("Num. Frozen UOCC");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &nfzv, sizeof(int) );

        free(keyword);
        return nfzv;
}

void Chkpt::wt_nfzv(int nfzv)
{
        char *keyword;
        keyword = build_keyword("Num. Frozen UOCC");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &nfzv, sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_nfzv()
** Reads in the total number of frozen unoccupied molecular orbitals.
**
** returns: nfzv = total number of frozen unoccupied molecular orbitals.
** \ingroup CHKPT
*/
        int chkpt_rd_nfzv(void)
        {
                return _default_chkpt_lib_->rd_nfzv();
        }

/*!
** void chkpt_wt_nfzv(int)
** Writes out the total number of frozen unoccupied molecular orbitals.
**
** \param nfzv = total number of frozen unoccupied molecular orbitals.
**
** \ingroup CHKPT
*/
        void chkpt_wt_nfzv(int nfzv)
        {
                _default_chkpt_lib_->wt_nfzv(nfzv);
        }
}
