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

double **Chkpt::rd_contr_full(const char */*key2*/)
{
#if 0
        double **contr, *temp_contr;
        int nprim, i, j, ij = 0;
        char *keyword;
        keyword = build_keyword("Contraction coefficients", key2);

        nprim = rd_nprim(key2);

        temp_contr = array<double>(MAXANGMOM*nprim);
        contr = matrix<double>(nprim,MAXANGMOM);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) temp_contr,
                MAXANGMOM*nprim*sizeof(double));

/* Picking non-zero coefficients to the "master" array contr */
        for(i=0,ij=0; i < MAXANGMOM; i++)
        for(j=0; j < nprim; j++, ij++) {
                contr[j][i] = temp_contr[ij];
        }

        free(temp_contr);
        free(keyword);
        return contr;
#endif
    return NULL;
}

extern "C" {
/*!
** chkpt_rd_contr_full(): Reads in the normalized contraction coefficients.
**
** Parameters: none
**
** Returns:
** double **contr, Normalized contraction coefficients are returned
** as a matrix of doubles.
** \ingroup CHKPT
*/

        double **chkpt_rd_contr_full(void)
        {
                return _default_chkpt_lib_->rd_contr_full();
        }
}

