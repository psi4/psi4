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

double Chkpt::rd_ecorr(void)
{
        double ecorr;
        char *keyword;
        keyword = build_keyword("Correlation energy");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &ecorr,
                sizeof(double));

        free(keyword);
        return ecorr;
}

void Chkpt::wt_ecorr(double ecorr)
{
        char *keyword;
        keyword = build_keyword("Correlation energy");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &ecorr,
                sizeof(double));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_ecorr():  Reads in the correlated energy.
**
** takes no arguments.
**
** returns: e_corr = the correlated energy.  To get some
**        information (a label) on the type of correlated wavefunction
**        used to get this energy, see rd_corr_lab().
** \ingroup CHKPT
*/
        double chkpt_rd_ecorr(void)
        {
                double energy;
                energy = _default_chkpt_lib_->rd_ecorr();
                return energy;
        }

/*!
** chkpt_wt_ecorr():  Writes out the correlated energy.
**
** \param e_corr = the correlated energy.  To get some
**        information (a label) on the type of correlated wavefunction
**        used to get this energy, see rd_corr_lab().
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_ecorr(double ecorr)
        {
                _default_chkpt_lib_->wt_ecorr(ecorr);
        }
}
