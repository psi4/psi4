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
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

namespace psi {

void PSIO::write_entry(unsigned int unit, const char *key, char *buffer, ULI size) {
  psio_address end = PSIO_ZERO;
  write(unit, key, buffer, size, PSIO_ZERO, &end);
}

  /*!
   ** PSIO_WRITE_ENTRY()
   **
   ** \ingroup PSIO
   */

  int psio_write_entry(unsigned int unit, const char *key, char *buffer, ULI size) {
    _default_psio_lib_->write_entry(unit, key, buffer, size);
    return 1;
  }

}

