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

double *Chkpt::rd_exps(const char *key2)
{
        double *exps;
        int nprim = 0;
        char *keyword;
        keyword = build_keyword("Exponents",key2);

        nprim = rd_nprim(key2);
        exps = array<double>(nprim);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) exps,
                nprim*sizeof(double));

        free(keyword);
        return exps;
}

void Chkpt::wt_exps(double *exps, const char *key2)
{
        int nprim;
        char *keyword;
        keyword = build_keyword("Exponents",key2);

        nprim = rd_nprim(key2);

        psio->write_entry(PSIF_CHKPT, keyword, (char *) exps,
                nprim*sizeof(double));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_exps():
** Reads in the exponents of the primitive Gaussian functions.
**
** takes no arguments.
**
** returns: double *exps
** The exponents are returned as an array of doubles.
** \ingroup CHKPT
*/
        double *chkpt_rd_exps(void)
        {
                double *exps;
                exps = _default_chkpt_lib_->rd_exps();
                return exps;
        }

/*!
** chkpt_wt_exps():
** Writes out the exponents of the primitive Gaussian functions.
**
** arguments:
**  \param exps = The exponents are returned as an array of doubles.
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_exps(double *exps, const char *key2)
        {
                _default_chkpt_lib_->wt_exps(exps, key2);
        }
}
