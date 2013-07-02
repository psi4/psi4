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

int Chkpt::rd_rottype(void)
{
        int rottype;
        char *keyword;
        keyword = build_keyword("Rotor type");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &rottype, sizeof(int));

        free(keyword);
        return rottype;
}

void Chkpt::wt_rottype(int rottype)
{
        char *keyword;
        keyword = build_keyword("Rotor type");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &rottype, sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_rottype()
** Reads in type of the rigid rotor molecule represents.
**
** returns: rottype = type of rigid rotor. Allowed values are:
**            0 - asymmetric top
**            1 - symmetric top
**            2 - spherical top
**            3 - linear molecule
**            6 - atom
** \ingroup CHKPT
*/
        int chkpt_rd_rottype(void)
        {
                return _default_chkpt_lib_->rd_rottype();
        }

/*!
** void chkpt_wt_rottype(int)
** Reads in type of the rigid rotor molecule represents.
**
** \param rottype = type of rigid rotor. Allowed values are:
**            0 - asymmetric top
**            1 - symmetric top
**            2 - spherical top
**            3 - linear molecule
**            6 - atom
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_rottype(int rottype)
        {
                _default_chkpt_lib_->wt_rottype(rottype);
        }
}
