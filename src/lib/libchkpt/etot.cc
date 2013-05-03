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

double Chkpt::rd_etot(void)
{
        double etot;
        char *keyword;
        keyword = build_keyword("Total energy");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &etot, sizeof(double));

        free(keyword);
        return etot;
}

void Chkpt::wt_etot(double etot)
{
        char *keyword;
        keyword = build_keyword("Total energy");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &etot, sizeof(double));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_etot(): Reads in the total energy.
**
**  takes no arguments.
**
**  returns: double etot  the total energy.
**  \ingroup CHKPT
*/
        double chkpt_rd_etot(void)
        {
                double etot;
                etot = _default_chkpt_lib_->rd_etot();
                return etot;
        }

/*!
** chkpt_wt_etot(): Writes out the total energy.
**
**  arguments:
**   \param double etot  the total energy.
**
**  returns: none
**  \ingroup CHKPT
*/
        void chkpt_wt_etot(double etot)
        {
                _default_chkpt_lib_->wt_etot(etot);
        }
}
