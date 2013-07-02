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

#include <cstdlib>
#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

int *Chkpt::rd_symoper(void)
{
        int *symoper;
        int nirreps;
        char *keyword;
        keyword = build_keyword("Cotton -> local map");

        nirreps = rd_nirreps();
        symoper = array<int>(nirreps);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) symoper, nirreps*sizeof(int));

        free(keyword);
        return symoper;
}

void Chkpt::wt_symoper(int *symoper)
{
        int nirreps;
        char *keyword;
        keyword = build_keyword("Cotton -> local map");

        nirreps = rd_nirreps();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) symoper, nirreps*sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int *chkpt_rd_symoper()
** Reads in the mapping array between "canonical" ordering of symmetry
** operations in the point group and the one defined in symmetry.h
**
** returns: symoper = Array nirrep long
**
** \ingroup CHKPT
*/
        int *chkpt_rd_symoper(void)
        {
                return _default_chkpt_lib_->rd_symoper();
        }

/*!
** void chkpt_wt_symoper(int *)
** Writes out the mapping array between "canonical" ordering of symmetry
** operations in the point group and the one defined in symmetry.h
**
** \param symoper = Array nirrep long
**
** returns: none
**
** \ingroup CHKPT
*/
        void chkpt_wt_symoper(int *symoper)
        {
                _default_chkpt_lib_->wt_symoper(symoper);
        }
}
