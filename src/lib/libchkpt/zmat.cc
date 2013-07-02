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

struct z_entry *Chkpt::rd_zmat(void)
{
        int nallatom;
        struct z_entry *z_geom;
        char *keyword;
        keyword = build_keyword("Z-matrix");

        nallatom = rd_nallatom();
        z_geom = (struct z_entry *) malloc(nallatom*(sizeof(struct z_entry)));

        psio->read_entry(PSIF_CHKPT, keyword, (char *) z_geom,
                sizeof(struct z_entry)*nallatom);

        free(keyword);
        return z_geom;
}

void Chkpt::wt_zmat(struct z_entry *z_geom)
{
        int nallatom;
        char *keyword;
        keyword = build_keyword("Z-matrix");

        nallatom = rd_nallatom();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) z_geom,
                sizeof(struct z_entry)*nallatom);

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_zmat():  Reads in the z_matrix.
**
**   takes no arguments.
**
**   returns: z_geom = An array natom long which contains
**     a z_entry struct for each atom
**
** \ingroup CHKPT
*/
        struct z_entry *chkpt_rd_zmat(void)
        {
                return _default_chkpt_lib_->rd_zmat();
        }

/*!
** chkpt_wt_zmat():  Writes out the z_matrix.
**
**  \param z_geom = An array natom long which contains
**     a z_entry struct for each atom
**
** returns: none
**
** \ingroup CHKPT
*/
        void chkpt_wt_zmat(struct z_entry *z_geom)
        {
                _default_chkpt_lib_->wt_zmat(z_geom);
        }
}
