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

double **Chkpt::rd_fgeom(void)
{
        int nallatom;
        double **fgeom;
        char *keyword;
        keyword = build_keyword("Full cartesian geometry");

        nallatom = rd_nallatom();
        fgeom = matrix<double>(nallatom,3);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) fgeom[0],
                                        (int) 3*nallatom*sizeof(double));

        free(keyword);
        return  fgeom;
}

void Chkpt::wt_fgeom(double **fgeom)
{
        int nallatom;
        char *keyword;
        keyword = build_keyword("Full cartesian geometry");

        nallatom = rd_nallatom();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) fgeom[0],
                (int) 3*nallatom*sizeof(double));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_fgeom():  Reads in full cartesian geometry including dummy atoms
**
** takes no arguments.
** returns: double **full_geom;
** \ingroup CHKPT
*/
        double **chkpt_rd_fgeom(void)
        {
                return _default_chkpt_lib_->rd_fgeom();
        }

/*!
** chkpt_wt_fgeom():  Writes out full cartesian geometry including dummy atoms
**
** arguments:
**   \param full_geom = Matrix for cartesian coordinates
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_fgeom(double **fgeom)
        {
                _default_chkpt_lib_->wt_fgeom(fgeom);
        }
}
