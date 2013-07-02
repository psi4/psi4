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

int Chkpt::rd_nprim(const char *key2)
{
    int nprim;
    char *keyword;
    keyword = build_keyword("Num. primitives", key2);

    psio->read_entry(PSIF_CHKPT, keyword, (char *) &nprim, sizeof(int));

    free(keyword);
    return nprim;
}

void Chkpt::wt_nprim(int nprim, const char *key2)
{
    char *keyword;
    keyword = build_keyword("Num. primitives", key2);

    psio->write_entry(PSIF_CHKPT, keyword, (char *) &nprim, sizeof(int));

    free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_nprim()
** Reads in the total number of primitive Gaussian functions
** (only primitives of symmetry independent atoms are taken into account!).
**
** returns: nprim = total number of primitive Gaussian functions.
** \ingroup CHKPT
*/
int chkpt_rd_nprim(void)
{
    return _default_chkpt_lib_->rd_nprim();
}

/*!
** void chkpt_wt_nprim(int)
** Writes out the total number of primitive Gaussian functions
** (only primitives of symmetry independent atoms are taken into account!).
**
** \param nprim = total number of primitive Gaussian functions.
**
** returns: none
** \ingroup CHKPT
*/
void chkpt_wt_nprim(int nprim, const char *key2)
{
    _default_chkpt_lib_->wt_nprim(nprim, key2);
}
}
