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

double Chkpt::rd_eref(void)
{
        double eref;
        char *keyword;
        keyword = build_keyword("Reference energy");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &eref, sizeof(double));

        free(keyword);
        return eref;
}

void Chkpt::wt_eref(double eref)
{
        char *keyword;
        keyword = build_keyword("Reference energy");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &eref,
                sizeof(double));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_eref(): Reads in the reference energy.
**
**   takes no arguments.
**
**   returns: double eref  the reference energy.
**
** \ingroup CHKPT
*/
        double chkpt_rd_eref(void)
        {
                double eref;
                eref = _default_chkpt_lib_->rd_eref();
                return eref;
        }


/*!
** chkpt_wt_eref(): Writes out the reference energy.
**
** \param double eref = the reference energy.
**
** returns: none
**
** \ingroup CHKPT
*/
        void chkpt_wt_eref(double eref)
        {
                _default_chkpt_lib_->wt_eref(eref);
        }
}
