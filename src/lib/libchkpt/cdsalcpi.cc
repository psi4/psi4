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

int *Chkpt::rd_cdsalcpi(void)
{
        const int nirreps = rd_nirreps();
        int *cdsalcpi = array<int>(nirreps);
        psio_address ptr = PSIO_ZERO;
        char *keyword = build_keyword("cartdisp SALCs per irrep");

        psio->read(PSIF_CHKPT, keyword, (char *) cdsalcpi, nirreps*sizeof(int), ptr, &ptr);

        free(keyword);
        return cdsalcpi;
}

void Chkpt::wt_cdsalcpi(const int *cdsalcpi)
{
        const int nirreps = rd_nirreps();
        psio_address ptr = PSIO_ZERO;
        char *keyword = build_keyword("cartdisp SALCs per irrep");

        psio->write(PSIF_CHKPT, keyword, (char *) cdsalcpi, nirreps*sizeof(int), ptr, &ptr);

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_cdsalcpi(): Read in number of SALCs per irrep
**
** Parameters: none
**
** Returns: cdsalcpi = An array of nirreps integers.
**
** \ingroup CHKPT
*/
        int *chkpt_rd_cdsalcpi(void)
        {
                return _default_chkpt_lib_->rd_cdsalcpi();
        }

/*!
** chkpt_wt_cdsalcpi(): Writes out number of SALCs per irrep
**
** \param cdsalcpi = An array of nirreps integers
**
** Returns: none
**
** \ingroup CHKPT
*/
        void chkpt_wt_cdsalcpi(const int *cdsalcpi)
        {
                _default_chkpt_lib_->wt_cdsalcpi(cdsalcpi);
        }
}

