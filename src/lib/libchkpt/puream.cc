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

bool Chkpt::rd_puream(const char *key2)
{
        int puream;
        char *keyword;
        keyword = build_keyword("Pure Harmonics?", key2);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &puream, sizeof(int));

        free(keyword);
        return (puream == 1);
}

void Chkpt::wt_puream(bool puream, const char *key2)
{
        char *keyword;
        keyword = build_keyword("Pure Harmonics?", key2);

        int p = puream ? 1 : 0;
        psio->write_entry(PSIF_CHKPT, keyword, (char *) &p, sizeof(int));

        free(keyword);
}

extern "C" {
/*!
**  int chkpt_rd_puream()
**  Reads whether cartesian or spherical harmonics are used (Psi is currently
**  limited to only using one type of functions at a time)
**
**  returns: 1 (harmonics) or 0 (cartesian)
**  \ingroup CHKPT
*/
        int chkpt_rd_puream(void)
        {
                return _default_chkpt_lib_->rd_puream() ? 1 : 0;
        }

/*!
**  void chkpt_wt_puream(int)
**  Writes whether cartesian or spherical harmonics are used (Psi is currently
**  limited to only using one type of functions at a time)
**
**  \param 1 (harmonics) or 0 (cartesian)
**
**  returns: none
**  \ingroup CHKPT
*/
        void chkpt_wt_puream(int puream, const char *key2)
        {
                _default_chkpt_lib_->wt_puream(puream == 1, key2);
        }
}
