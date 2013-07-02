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

Chkpt::~Chkpt()
{
    // The chkpt might be closed...check
    if (psio->open_check(PSIF_CHKPT))
        psio->close(PSIF_CHKPT, 1);
    psio = NULL;
}

/*!
**  chkpt_close()  closes up the checkpoint file.
**
**  Parameters: none, but chkpt_init must already have been called for
**    this to work.
**
**  Returns: none
**  \ingroup CHKPT
*/
extern "C" {
    int chkpt_close(void)
    {
        if (_default_chkpt_lib_) {
            _default_chkpt_lib_.reset();
        }
        return 0;
    }
}
