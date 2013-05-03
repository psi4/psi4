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

#include <stdlib.h>
#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

int Chkpt::rd_nfragment(void)
{
        int nfragment;
        char *keyword;
        keyword = build_keyword("Num. fragments");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &nfragment, sizeof(int));

        free(keyword);
        return nfragment;
}

void Chkpt::wt_nfragment(int nfragment)
{
        char *keyword;
        keyword = build_keyword("Num. fragments");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &nfragment, sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_nfragment()
** Reads in the total number of irreducible representations
** in the point group in which the molecule is being considered.
**
** returns: nfragment = total number of irreducible representations.
** \ingroup CHKPT
*/
        int chkpt_rd_nfragment(void)
        {
                return _default_chkpt_lib_->rd_nfragment();
        }

/*!
** void chkpt_wt_nfragment(int)
** Writes out the total number of irreducible representations
** in the point group in which the molecule is being considered.
**
** \param nfragment = total number of irreducible representations.
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_nfragment(int nfragment)
        {
                _default_chkpt_lib_->wt_nfragment(nfragment);
        }
}
