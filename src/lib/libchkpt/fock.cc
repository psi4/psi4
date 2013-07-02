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

double *Chkpt::rd_fock(void)
{
        double *fmat;
        char *keyword;
        keyword = build_keyword("Fock Matrix");

        int nso = rd_nso();

        fmat = array<double>((nso*nso+nso)/2);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) fmat,
                           ((nso*nso+nso)/2)*sizeof(double));
        free(keyword);
        return fmat;
}

void Chkpt::wt_fock(double *fmat)
{
        char *keyword;
        keyword = build_keyword("Fock Matrix");

        int nso = rd_nso();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) fmat,
                ((nso*nso+nso)/2)*sizeof(double));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_fock():  Reads in the Fock Matrix
**
**  takes no arguments.
**
**  returns: double *fmat  an array lower triangle closed shell fock matrix
**      ordered by irrep.
** \ingroup CHKPT
*/
        double *chkpt_rd_fock(void)
        {
                double *fmat;
                fmat = _default_chkpt_lib_->rd_fock();
                return fmat;
        }

/*!
** chkpt_wt_fock():  Writes the Fock Matrix
**
** arguments:
**  \param evals = an array of the lower triangle part of the fock matrix
**      ordered by irrep.
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_fock(double *fmat)
        {
                _default_chkpt_lib_->wt_fock(fmat);
        }

}
