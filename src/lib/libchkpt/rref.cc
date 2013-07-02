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

double **Chkpt::rd_rref(void)
{
//	char *key;
        double **Rref;
        char *keyword;
        keyword = build_keyword("Transmat to reference frame");

        Rref = matrix<double>(3,3);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) Rref[0], sizeof(double)*9);

        free(keyword);
        return Rref;
}

void Chkpt::wt_rref(double **Rref)
{
        char *keyword;
        keyword = build_keyword("Transmat to reference frame");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) Rref[0], sizeof(double)*9);

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_rref()
** Reads in a 3x3 matrix used to rotate back to the reference frame.
**
**   takes no arguments.
**
**   returns: rref = A 3x3 matrix describing the rotation back to the
**            reference frame.  The reference frame is a coordinate system
**            defined by the "raw" geometry specification (either Z-matrix
**            or geometry array in input.dat or chkpt). Can be used to
**            transform quantities corresponding to different but similar
**            calculations (gradients at displaced geometries) to a
**            common frame.
**
** \ingroup CHKPT
*/
        double **chkpt_rd_rref(void)
        {
                return _default_chkpt_lib_->rd_rref();
        }

/*!
** chkpt_wt_rref()
** Writes out a 3x3 matrix used to rotate back to the reference frame.
**
** \params rref = A 3x3 matrix describing the rotation back to the reference
**                frame.  The reference frame is a coordinate system defined
**                by the "raw" geometry specification (either Z-matrix or
**                geometry array in input.dat or chkpt). Can be used to
**                transform quantities corresponding to different but
**                similar calculations (gradients at displaced geometries)
**                to a common frame.
**
** returns: none
**
** \ingroup CHKPT
*/
        void chkpt_wt_rref(double **Rref)
        {
                _default_chkpt_lib_->wt_rref(Rref);
        }
}
