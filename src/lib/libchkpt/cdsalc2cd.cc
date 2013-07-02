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

double **Chkpt::rd_cdsalc2cd(void)
{
  const int num_cd = 3*rd_natom();
  double **cdsalc2cd = matrix<double>(num_cd,num_cd);
  psio_address ptr = PSIO_ZERO;
  char *keyword = build_keyword("cartdisp->SALC matrix");

  psio->read(PSIF_CHKPT, keyword, (char *) cdsalc2cd[0],
    num_cd*num_cd*sizeof(double), ptr, &ptr);

  free(keyword);
  return cdsalc2cd;
}

void Chkpt::wt_cdsalc2cd(const double **cdsalc2cd)
{
  const int num_cd = 3*rd_natom();
  psio_address ptr = PSIO_ZERO;
  char *keyword = build_keyword("cartdisp->SALC matrix");

  psio->write(PSIF_CHKPT, keyword, (char *) cdsalc2cd[0],
    num_cd*num_cd*sizeof(double), ptr, &ptr);

  free(keyword);
}

extern "C" {
/*!
** chkpt_rd_cdsalc2cd(): Read in (normalized) SALCs of cartesian displacements
**
** Parameters: none
**
** Returns: cdsalc2cd = A natom*3 by natom*3 blocked matrix of doubles.
**   Columns correpond to symmetry-blocked SALCs
**
** \ingroup CHKPT
*/
        double **chkpt_rd_cdsalc2cd(void)
        {
                double **cdsalc2cd;
                cdsalc2cd = _default_chkpt_lib_->rd_cdsalc2cd();
                return cdsalc2cd;
        }


/*!
** chkpt_wt_cdsalc2cd(): Writes (normalized) SALCs of cartesian displacements
**
** \param cdsalc2cd = A natom*3 by natom*3 blocked matrix of doubles.
**   Columns correpond to symmetry-blocked SALCs
**
** Returns: none
**
** \ingroup CHKPT
*/
        void chkpt_wt_cdsalc2cd(const double **cdsalc2cd)
        {
                _default_chkpt_lib_->wt_cdsalc2cd(cdsalc2cd);
        }
}

