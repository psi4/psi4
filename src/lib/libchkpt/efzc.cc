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

double Chkpt::rd_efzc(void)
{
        double efzc;
        char *keyword;
        keyword = build_keyword("Frozen core energy");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &efzc, sizeof(double));

        free(keyword);
        return efzc;
}

void Chkpt::wt_efzc(double efzc)
{
        char *keyword;
        keyword = build_keyword("Frozen core energy");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &efzc,
                sizeof(double));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_efzc(): Reads in the frozen-core energy.
**
**   takes no arguments.
**
**   returns: double efzc  the frozen-core energy.
** \ingroup CHKPT
*/
        double chkpt_rd_efzc(void)
        {
                double efzc;
                efzc = _default_chkpt_lib_->rd_efzc();
                return efzc;
        }

/*!
** chkpt_wt_efzc(): Writes out the frozen-core energy.
**
** \param efzc = the frozen-core energy.
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_efzc(double efzc)
        {
                _default_chkpt_lib_->wt_efzc(efzc);
        }
}
