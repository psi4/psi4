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

int Chkpt::exist(const char *keyword)
{
        int exists=0;

        if (psio->tocscan(PSIF_CHKPT, keyword) != NULL)
                exists=1;

        return exists;
}

extern "C" {
/*!
** chkpt_exist(): Checks to see if entry already exists in chkpt file. Note
** this function should be called only by functions in the chkpt library, as
** the calling function prepends the prefix.
**
**   takes no arguments.
**
**   returns: 1 if entry exists, 0 otherwise
**
** \ingroup CHKPT
*/
        int chkpt_exist(const char *keyword)
        {
                return(_default_chkpt_lib_->exist(keyword));
        }
}
