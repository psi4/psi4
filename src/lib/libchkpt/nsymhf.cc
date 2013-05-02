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

int Chkpt::rd_nsymhf(void)
{
        int nsymhf;
        char *keyword;
        keyword = build_keyword("Num. HF irreps");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &nsymhf, sizeof(int));

        free(keyword);
        return nsymhf;
}

void Chkpt::wt_nsymhf(int nsymhf)
{
        char *keyword;
        keyword = build_keyword("Num. HF irreps");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &nsymhf, sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_nsymhf()
** Reads in the total number of irreps in the point group
** in which the molecule is being considered which
** have non-zero number of basis functions.
**
** returns: nirreps = total number of irreducible representations
**      with a non-zero number of basis functions. For STO or DZ water, for
**      example, this is three, even though nirreps is 4 (see rd_nirreps()).
** \ingroup CHKPT
*/
        int chkpt_rd_nsymhf(void)
        {
                return _default_chkpt_lib_->rd_nsymhf();
        }

/*!
** void chkpt_wt_nsymhf(int)
** Writes out the total number of irreps in the point group
** in which the molecule is being considered which
** have non-zero number of basis functions.
**
** \param nirreps = total number of irreducible representations
**      with a non-zero number of basis functions. For STO or DZ water, for
**      example, this is three, even though nirreps is 4 (see rd_nirreps()).
** \ingroup CHKPT
*/
        void chkpt_wt_nsymhf(int nsymhf)
        {
                _default_chkpt_lib_->wt_nsymhf(nsymhf);
        }
}
