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

int *Chkpt::rd_atom_position(void)
{
        int *atom_position, natom;
        char *keyword;
        keyword = build_keyword("Atomic symm positions");

        natom = rd_natom();
        atom_position = array<int>(natom);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) atom_position,
                natom*sizeof(int));

        free(keyword);
        return atom_position;
}

void Chkpt::wt_atom_position(int *atom_position)
{
        int natom;
        char *keyword;
        keyword = build_keyword("Atomic symm positions");

        natom = rd_natom();

        psio->write_entry(PSIF_CHKPT, keyword,
                (char *) atom_position, natom*sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int *chkpt_rd_atom_position()
**
** Reads in symmetry positions of atoms:
**    Possible values are as follows:
**	1   - atom in general position
**      2   - atom on c2z axis
**	4   - atom on c2y axis
**	8   - atom on c2x axis
**	16  - atom in the inversion center
**	32  - atom in the sigma_xy plane
**	64  - atom in the sigma_xz plane
**	128 - atom in the sigma_yz plane
**	This data is sufficient to define stabilizers of the nuclei.
**
** Returns: int *atom_position, an array of symmetry positions of atoms
** \ingroup CHKPT
*/
        int *chkpt_rd_atom_position(void)
        {
                return _default_chkpt_lib_->rd_atom_position();
        }


/*!
** chkpt_wt_atom_position()
**
** Writes out symmetry positions of atoms:
**    Possible values are as follows:
**	1   - atom in general position
**      2   - atom on c2z axis
**	4   - atom on c2y axis
**	8   - atom on c2x axis
**	16  - atom in the inversion center
**	32  - atom in the sigma_xy plane
**	64  - atom in the sigma_xz plane
**	128 - atom in the sigma_yz plane
**	This data is sufficient to define stabilizers of the nuclei.
**
** \param atom_position = an array of symmetry positions of atoms
**
** Returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_atom_position(int *atom_position)
        {
                _default_chkpt_lib_->wt_atom_position(atom_position);
        }
}

