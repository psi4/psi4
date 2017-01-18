/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

/*!
 \file
 \ingroup PSIO
 */

#include <cstdlib>
#include <unistd.h>
 #include "psi4/pragma.h"
 PRAGMA_WARNING_PUSH
 PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
 #include <memory>
 PRAGMA_WARNING_POP
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"

namespace psi {

void PSIO::tocwrite(unsigned int unit) {
  unsigned int i;
  psio_ud *this_unit;
  psio_tocentry *this_entry;
  ULI entry_size;
  psio_address address;

  this_unit = &(psio_unit[unit]);
  entry_size = sizeof(psio_tocentry) - 2*sizeof(psio_tocentry *);

  if (!open_check(unit))
    return;

  wt_toclen(unit, this_unit->toclen);

  this_entry = this_unit->toc;
  address = psio_get_address(PSIO_ZERO, sizeof(ULI));
  for (i=0; i < this_unit->toclen; i++) {
    rw(unit, (char *) this_entry, address, entry_size, 1);
    this_entry = this_entry->next;
    if (this_entry != NULL)
      address = this_entry->sadd;
  }
}

  /*!
   ** PSIO_TOCWRITE(): Write the table of contents for file number 'unit'.
   **
   ** \param unit  = The PSI unit to which we will write the TOC.
   **
   ** NB: This function should NOT call psio_error because the latter calls it!
   **
   ** \ingroup PSIO
   */

  int psio_tocwrite(unsigned int unit) {
    _default_psio_lib_->tocwrite(unit);
    return 1;
  }

}
