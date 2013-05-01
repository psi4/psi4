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
  \file read and write the rotational symmetry number
  \ingroup CHKPT
*/

#include <cstdlib>
#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

int Chkpt::rd_rot_symm_num(void)
{
  int rot_symm_num;
  char *keyword;
  keyword = build_keyword("Rotational Symmetry Number");

  psio->read_entry(PSIF_CHKPT, keyword, (char *) &rot_symm_num, sizeof(int));

  free(keyword);
  return rot_symm_num;
}

void Chkpt::wt_rot_symm_num(int rot_symm_num)
{
  char *keyword;
  keyword = build_keyword("Rotational Symmetry Number");

  psio->write_entry(PSIF_CHKPT, keyword, (char *) &rot_symm_num, sizeof(int));

  free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_rot_symm_num()
** Reads the rotational symmetry number.
**
** returns: rot_symm_num = rotational symmetry number
** \ingroup CHKPT
*/
  int chkpt_rd_rot_symm_num(void)
  {
  return _default_chkpt_lib_->rd_rot_symm_num();
  }

/*!
** void chkpt_wt_rot_symm_num(int)
** Writes the rotational symmetry number.
**
** \param rot_symm_num = rotational symmetry number
**
** \ingroup CHKPT
*/
  void chkpt_wt_rot_symm_num(int rot_symm_num)
  {
  _default_chkpt_lib_->wt_rot_symm_num(rot_symm_num);
  }
}
