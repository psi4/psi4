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

int Chkpt::rd_ncalcs(void)
{
        char *keyword_mo, *keyword_alpha_mo;
        keyword_mo = build_keyword("MO coefficients");
        keyword_alpha_mo = build_keyword("Alpha MO coefficients");

        if (psio->tocscan(PSIF_CHKPT, keyword_mo) == NULL &&
                psio->tocscan(PSIF_CHKPT, keyword_alpha_mo) == NULL)
                return 0;
        else
                return 1;

        free(keyword_mo);
        free(keyword_alpha_mo);
}

extern "C" {
/*!
** int chkpt_rd_ncalcs()
** Reads in the total number of HF wave functions.
**
** returns: ncalcs = total number of HF wave functions in checkpoint
** \ingroup CHKPT
*/
        int chkpt_rd_ncalcs(void)
        {
                return _default_chkpt_lib_->rd_ncalcs();
        }
}

