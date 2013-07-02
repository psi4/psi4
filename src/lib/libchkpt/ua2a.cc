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

int *Chkpt::rd_ua2a(void)
{
        int *ua2a;
        int num_unique_atoms;
        char *keyword;
        keyword = build_keyword("Unique atom -> full atom map");

        num_unique_atoms = rd_num_unique_atom();
        ua2a = array<int>(num_unique_atoms);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) ua2a, num_unique_atoms*sizeof(int));

        free(keyword);
        return ua2a;
}

void Chkpt::wt_ua2a(int *ua2a)
{
        int num_unique_atoms;
        char *keyword;
        keyword = build_keyword("Unique atom -> full atom map");

        num_unique_atoms = rd_num_unique_atom();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) ua2a, num_unique_atoms*sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int *chkpt_rd_ua2a()
** Reads in a mapping array from the symmetry-unique atom
** list to the full atom list
**
** returns: ua2a = Read in an array num_unique_atom long
**
** \ingroup CHKPT
*/
        int *chkpt_rd_ua2a(void)
        {
                return _default_chkpt_lib_->rd_ua2a();
        }

/*!
** void chkpt_wt_ua2a(int *)
** Writes out a mapping array from the symmetry-unique atom
** list to the full atom list
**
** \param ua2a = An array num_unique_atom long
**
** returns: none
*/
        void chkpt_wt_ua2a(int *ua2a)
        {
                _default_chkpt_lib_->wt_ua2a(ua2a);
        }
}
