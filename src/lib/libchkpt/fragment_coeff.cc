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

#include <stdio.h>
#include <stdlib.h>
#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>
#include <libciomr/libciomr.h>

using namespace psi;

double ***Chkpt::rd_fragment_coeff(void)
{
        int nfragment, *natom_per_fragment, *nref_per_fragment, i, j;
    psio_address ptr;
        double ***fragment_coeff;
        char *keyword;
        keyword = build_keyword("Fragment coeff");

        nfragment = rd_nfragment();
        natom_per_fragment = rd_natom_per_fragment();
        nref_per_fragment = rd_nref_per_fragment();

        fragment_coeff = array<double **>(nfragment);
    for (i=0; i<nfragment; ++i)
          fragment_coeff[i] = matrix<double>(nref_per_fragment[i],natom_per_fragment[i]);

    ptr = PSIO_ZERO;
    for (i=0; i<nfragment; ++i)
      for (j=0; j<nref_per_fragment[i]; ++j)
            psio->read(PSIF_CHKPT, keyword, (char *) fragment_coeff[i][j],
            (int) natom_per_fragment[i]*sizeof(double), ptr, &ptr);

    free(natom_per_fragment);
    free(nref_per_fragment);
        free(keyword);
        return fragment_coeff;
}

void Chkpt::wt_fragment_coeff(double ***fragment_coeff)
{
        int nfragment, *natom_per_fragment, *nref_per_fragment, i, j;
    psio_address ptr;
        char *keyword;
        keyword = build_keyword("Fragment coeff");

        nfragment = rd_nfragment();
    natom_per_fragment = rd_natom_per_fragment();
    nref_per_fragment = rd_nref_per_fragment();

    ptr = PSIO_ZERO;
    for (i=0; i<nfragment; ++i)
      for (j=0; j<nref_per_fragment[i]; ++j)
            psio->write(PSIF_CHKPT, keyword, (char *) fragment_coeff[i][j],
            (int) natom_per_fragment[i]*sizeof(double), ptr, &ptr);

    free(natom_per_fragment);
    free(nref_per_fragment);
        free(keyword);
    return;
}

extern "C" {
/*!
** chkpt_rd_fragment_coeff():  Reads in the coefficients specifying reference points
** for molecular fragments
**
**   takes no arguments.
**
**   returns:
**     double ***fragment_coeff[fragment][reference point][atom in fragment]
** \ingroup CHKPT
*/
        double ***chkpt_rd_fragment_coeff(void)
        {
                double ***fragment_coeff;
                fragment_coeff = _default_chkpt_lib_->rd_fragment_coeff();
                return fragment_coeff;
        }


/*!
** chkpt_wt_fragment_coeff():  Writes out the coefficients specifying the reference
** points for molecular fragments
**
** \param double ***fragment_coeff[fragment][reference point][atom in fragment]
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_fragment_coeff(double ***fragment_coeff)
        {
                _default_chkpt_lib_->wt_fragment_coeff(fragment_coeff);
        }
}
