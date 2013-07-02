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

double Chkpt::rd_enuc(void)
{
        double enuc;
        char *keyword;
        keyword = build_keyword("Nuclear rep. energy");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &enuc, sizeof(double));

        free(keyword);
        return enuc;
}

void Chkpt::wt_enuc(double enuc)
{
        char *keyword;
        keyword = build_keyword("Nuclear rep. energy");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &enuc, sizeof(double));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_enuc(): Reads in the nuclear repulsion energy
**
**   takes no arguments.
**
**   returns: double enuc  the nuclear repulsion energy.
**
** \ingroup CHKPT
*/
        double chkpt_rd_enuc(void)
        {
                double enuc;
                enuc = _default_chkpt_lib_->rd_enuc();
                return enuc;
        }

/*!
** chkpt_wt_enuc(): Writes out the nuclear repulsion energy
**
** \param enuc = the nuclear repulsion energy.
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_enuc(double enuc)
        {
                _default_chkpt_lib_->wt_enuc(enuc);
        }
}
