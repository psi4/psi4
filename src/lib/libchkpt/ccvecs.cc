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
#include <cmath>
#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

double **Chkpt::rd_ccvecs(void)
{
        int nmo, ccvec_length;
        double **ccvecs;
        char *keyword;

        nmo = rd_nmo();
        ccvec_length = abs(rd_iopen());

        if (ccvec_length > 0) {
                ccvecs = matrix<double>(2,ccvec_length);

                keyword = build_keyword("SCF coupling coefficients");

                psio->read_entry(PSIF_CHKPT, keyword, (char *) ccvecs[0],
                        2*ccvec_length*sizeof(double));

                free(keyword);

                return ccvecs;
        }
        else return NULL;
}

void Chkpt::wt_ccvecs(double **ccvecs)
{
        int nmo, ccvec_length;
        char *keyword;

        nmo = rd_nmo();
        ccvec_length = abs(rd_iopen());

        if (ccvec_length > 0) {
                keyword = build_keyword("SCF coupling coefficients");

                psio->write_entry(PSIF_CHKPT, keyword, (char *) ccvecs[0],
                        2*ccvec_length*sizeof(double));

                free(keyword);
        }
}

extern "C" {
/*!
** chkpt_rd_ccvecs()
**
** Reads in a matrix, rows of which are ALPHA (ccvecs[0])
** and BETA (ccvecs[1]) matrices of coupling coefficients
** for open shells stored in lower triangular form.
** Coupling coefficients are defined NOT as in
** C.C.J.Roothaan Rev. Mod. Phys. 32, 179 (1960) as it's
** stated in the manual pages for CSCF, but according to
** Pitzer (...) and are **different** from those in
** Yamaguchi, Osamura, Goddard, and Schaefer's book
** "Analytic Derivative Methods in Ab Initio Molecular
** Electronic Structure Theory".
**
** The relationship between Pitzer's and Yamaguchi's
** conventions are follows :
** ALPHA = 1-2*a , BETA = -1-4*b , where a and b are
** alpha's and beta's for open shells defined on pp. 69-70
** of Dr. Yamaguchi's book.
**
** Parameters: none
**
** Returns:
** double **ccvecs, a matrix 2 by abs(IOPEN) rows of which are coupling
** coefficient matrices for open-shells in packed form.
** \ingroup CHKPT
*/
        double **chkpt_rd_ccvecs(void)
        {
                return _default_chkpt_lib_->rd_ccvecs();
        }


/*!
** chkpt_wt_ccvecs()
**
** Writes a matrix of coupling coefficients.  See the
** comments chkpt_rd_ccvecs() above.
**
** \param ccvecs = a matrix 2 by abs(IOPEN) rows of which are coupling
**   coefficient matrices for open-shells in packed form.
**
** Returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_ccvecs(double **ccvecs)
        {
                _default_chkpt_lib_->wt_ccvecs(ccvecs);
        }
}

