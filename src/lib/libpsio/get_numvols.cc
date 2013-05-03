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
 \ingroup PSIO
 */

#include <boost/shared_ptr.hpp>
#include <cstdio>
#include <cstdlib>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

namespace psi {

unsigned int PSIO::get_numvols(unsigned int unit) {
  std::string charnum;
  charnum = filecfg_kwd("PSI", "NVOLUME", unit);
  if (!charnum.empty())
    return ((unsigned int)atoi(charnum.c_str()));
  charnum = filecfg_kwd("PSI", "NVOLUME", -1);
  if (!charnum.empty())
    return ((unsigned int)atoi(charnum.c_str()));
  charnum = filecfg_kwd("DEFAULT", "NVOLUME", unit);
  if (!charnum.empty())
    return ((unsigned int)atoi(charnum.c_str()));
  charnum = filecfg_kwd("DEFAULT", "NVOLUME", -1);
  if (!charnum.empty())
    return ((unsigned int)atoi(charnum.c_str()));

  // assume that the default has been provided already
  abort();
}

  unsigned int psio_get_numvols_default(void) {
    std::string charnum;

    charnum = _default_psio_lib_->filecfg_kwd("PSI", "NVOLUME", -1);
    if (!charnum.empty())
      return ((unsigned int)atoi(charnum.c_str()));
    charnum = _default_psio_lib_->filecfg_kwd("DEFAULT", "NVOLUME", -1);
    if (!charnum.empty())
      return ((unsigned int)atoi(charnum.c_str()));

    // assume that the default has been provided already
    abort();
  }

}

