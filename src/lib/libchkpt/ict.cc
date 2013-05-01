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

int **Chkpt::rd_ict(void)
{
        int i, natom, nirreps;
        int **ict;
        psio_address ptr;
        char *keyword;
        keyword = build_keyword("ICT Table");

        nirreps = rd_nirreps();
        natom = rd_natom();

        ptr = PSIO_ZERO;
        ict = (int **) malloc(sizeof(char *) * nirreps);
        for(i=0; i < nirreps; i++) {
                ict[i] = (int *) malloc(sizeof(int) * natom);
                psio->read(PSIF_CHKPT, keyword, (char *) ict[i],
                  natom*sizeof(int), ptr, &ptr);
        }

        free(keyword);
        return ict;
}

void Chkpt::wt_ict(int **ict)
{
        int i, natom, nirreps;
        psio_address ptr;
        char *keyword;
        keyword = build_keyword("ICT Table");

        nirreps = rd_nirreps();
        natom = rd_natom();

        ptr = PSIO_ZERO;
        for(i=0; i < nirreps; i++)
                psio->write(PSIF_CHKPT, keyword, (char *) ict[i],
                  natom*sizeof(int), ptr, &ptr);

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_ict():  Reads the transformation properties of the nuclei
**     under the operations allowed for the particular symmetry point group
**     in which the molecule is considered.
**
**   takes no arguments.
**
**   returns: ict = a matrix of integers. Each row corresponds
**     to a particular symmetry operation, while each column corresponds to
**     a particular atom.  The value of ict[2][1], then, should be interpreted
**     in the following manner: under the third symmetry operation of the
**     relavant point group, the second atom is placed in the location
**     originally occupied by the atom with the index ict[2][1].
** \ingroup CHKPT
*/
        int **chkpt_rd_ict(void)
        {
                return _default_chkpt_lib_->rd_ict();
        }

/*!
** chkpt_wt_ict():  Reads the transformation properties of the nuclei
**     under the operations allowed for the particular symmetry point group
**     in which the molecule is considered.
**
**   arguments:
**   \param ict = a matrix of integers. Each row corresponds
**     to a particular symmetry operation, while each column corresponds to
**     a particular atom.  The value of ict[2][1], then, should be interpreted
**     in the following manner: under the third symmetry operation of the
**     relavant point group, the second atom is placed in the location
**     originally occupied by the atom with the index ict[2][1].
**
**   returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_ict(int **ict)
        {
                _default_chkpt_lib_->wt_ict(ict);
        }
}
