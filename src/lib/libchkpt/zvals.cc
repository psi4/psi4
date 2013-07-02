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
*/

#include <cstdio>
#include <cstdlib>
#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

double *Chkpt::rd_zvals(void)
{
        int natom;
        double *zvals;
        char *keyword;
        keyword = build_keyword("Nuclear charges");

        natom = rd_natom();
        zvals = array<double>(natom);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) zvals,
                natom*sizeof(double));

        free(keyword);
        return zvals;
}

void Chkpt::wt_zvals(double *zvals)
{
        int natom;
        char *keyword;
        keyword = build_keyword("Nuclear charges");

        natom = rd_natom();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) zvals,
                natom*sizeof(double));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_zvals()
** Reads the nuclear charges from the checkpoint file.
**
** arguments: none
**
** returns:
**   double *zvals: An array of the charges
*/
        double *chkpt_rd_zvals(void)
        {
                return _default_chkpt_lib_->rd_zvals();
        }

/*!
** chkpt_wt_zvals()
** Writes the nuclear charges to the checkpoint file.
**
** \param double *zvals: An array of the charges
**
** returns: nothing
*/
        void chkpt_wt_zvals(double *zvals)
        {
                _default_chkpt_lib_->wt_zvals(zvals);
        }
}
