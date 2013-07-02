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

double **Chkpt::rd_usotbf(const char *key2)
{
  if (rd_puream(key2)) {
        double **usotbf;
        int num_so, i;
        psio_address ptr;
        char *keyword;
        keyword = build_keyword("SO->BF transmat", key2);

        num_so = rd_nso(key2);
        usotbf = matrix<double>(num_so,num_so);

        ptr = PSIO_ZERO;
        for(i=0;i<num_so;i++)
                psio->read(PSIF_CHKPT, keyword, (char *) usotbf[i], (int) num_so*sizeof(double), ptr, &ptr);

        free(keyword);
        return usotbf;
  }
  else
    return rd_usotao(key2);
}

void Chkpt::wt_usotbf(double **usotbf, const char *key2)
{
        int num_so, i;
        psio_address ptr;
        char *keyword;
        keyword = build_keyword("SO->BF transmat", key2);

        num_so = rd_nso(key2);

        ptr = PSIO_ZERO;
        for(i=0;i<num_so;i++)
                psio->write(PSIF_CHKPT, keyword, (char *) usotbf[i], (int) num_so*sizeof(double), ptr, &ptr);

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_usotbf(): Reads in the SO to basis functions transformation matrix
**
** takes no arguments.
**
** returns: usotbf = Read in a num_so by num_so matrix of doubles
**
** \ingroup CHKPT
*/
        double **chkpt_rd_usotbf(void)
        {
                return _default_chkpt_lib_->rd_usotbf();
        }

/*!
** chkpt_wt_usotbf(): Writes out the SO to basis functions transformation
**                    matrix
**
** \param usotbf = A num_so by num_so matrix of doubles
**
** returns: none
**
** \ingroup CHKPT
*/
        void chkpt_wt_usotbf(double **usotbf, const char *key2)
        {
                _default_chkpt_lib_->wt_usotbf(usotbf, key2);
        }
}
