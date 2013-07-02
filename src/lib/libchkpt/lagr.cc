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

double **Chkpt::rd_lagr(void)
{
        int nmo;
        double **lagr;
        char *keyword;
        keyword = build_keyword("MO Lagrangian");

        nmo = rd_nmo();

        lagr = matrix<double>(nmo,nmo);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) lagr[0],
                nmo*nmo*sizeof(double));

        free(keyword);
        return lagr;
}

void Chkpt::wt_lagr(double **lagr)
{
        int nmo;
        char *keyword;
        keyword = build_keyword("MO Lagrangian");

        nmo = rd_nmo();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) lagr[0],
                nmo*nmo*sizeof(double));

        free(keyword);
}

double **Chkpt::rd_alpha_lagr(void)
{
        int nmo;
        double **lagr;
        char *keyword;
        keyword = build_keyword("Alpha MO Lagrangian");

        nmo = rd_nmo();

        lagr = matrix<double>(nmo,nmo);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) lagr[0],
                nmo*nmo*sizeof(double));

        free(keyword);
        return lagr;
}

void Chkpt::wt_alpha_lagr(double **lagr)
{
        int nmo;
        char *keyword;
        keyword = build_keyword("Alpha MO Lagrangian");

        nmo = rd_nmo();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) lagr[0],
                nmo*nmo*sizeof(double));

        free(keyword);
}

double **Chkpt::rd_beta_lagr(void)
{
        int nmo;
        double **lagr;
        char *keyword;
        keyword = build_keyword("Beta MO Lagrangian");

        nmo = rd_nmo();

        lagr = matrix<double>(nmo,nmo);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) lagr[0],
                nmo*nmo*sizeof(double));

        free(keyword);
        return lagr;
}

void Chkpt::wt_beta_lagr(double **lagr)
{
        int nmo;
        char *keyword;
        keyword = build_keyword("Beta MO Lagrangian");

        nmo = rd_nmo();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) lagr[0],
                nmo*nmo*sizeof(double));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_lagr():  Reads in the MO lagrangian matrix for RHF/ROHF.
**
** Parameters: none
**
** Returns:
**   double **lagr, a matrix nmo by nmo.
**
** \ingroup CHKPT
*/
        double **chkpt_rd_lagr(void)
        {
                double **lagr;
                lagr = _default_chkpt_lib_->rd_lagr();
                return lagr;
        }


/*!
** chkpt_wt_lagr():  Writes the MO lagrangian matrix for RHF/ROHF.
**
** \param lagr = Lagrangian matrix with dimensions nmo by nmo.
**
** Returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_lagr(double **lagr)
        {
                _default_chkpt_lib_->wt_lagr(lagr);
        }

/*!
** chkpt_rd_alpha_lagr():  Reads in the alpha MO lagrangian matrix for UHF.
**
** Parameters: none
**
** Returns:
**   double **lagr, a matrix nmo by nmo.
** \ingroup CHKPT
*/
        double **chkpt_rd_alpha_lagr(void)
        {
                double **lagr;
                lagr = _default_chkpt_lib_->rd_alpha_lagr();
                return lagr;
        }

/*!
** chkpt_wt_alpha_lagr():  Writes the alpha MO lagrangian matrix for UHF.
**
** \param lagr = Lagrangian matrix of size nmo by nmo.
**
** returns: none
** \ingroup CHKPT
*/
void chkpt_wt_alpha_lagr(double **lagr)
{
        _default_chkpt_lib_->wt_alpha_lagr(lagr);
}


/*!
** chkpt_rd_beta_lagr():  Reads in the beta MO lagrangian matrix for UHF.
**
** takes no arguments.
**
** returns:
**	double **lagr	a matrix nmo by nmo.
** \ingroup CHKPT
*/
        double **chkpt_rd_beta_lagr(void)
        {
                double **lagr;
                lagr = _default_chkpt_lib_->rd_beta_lagr();
                return lagr;
        }

/*!
** chkpt_wt_beta_lagr():  Writes the beta MO lagrangian matrix for UHF.
**
** \param lagr = Lagrangian matrix of size nmo by nmo.
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_beta_lagr(double **lagr)
        {
                _default_chkpt_lib_->wt_beta_lagr(lagr);
        }
}
