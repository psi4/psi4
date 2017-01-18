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

#include <cstdio>
#include <unistd.h>
#include <cstdlib>
 #include "psi4/pragma.h"
 PRAGMA_WARNING_PUSH
 PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
 #include <memory>
 PRAGMA_WARNING_POP
#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"

namespace psi {

bool PSIO::tocdel(unsigned int unit, const char *key) {
  psio_tocentry *this_entry = tocscan(unit, key);

  if (this_entry == NULL) return false;

  psio_tocentry *last_entry = this_entry->last;
  psio_tocentry *next_entry = this_entry->next;

  if(next_entry == NULL) last_entry->next = NULL;
  else {
    last_entry->next = next_entry;
    next_entry->last = last_entry;
  }

  free(this_entry);
  psio_ud *this_unit = &(psio_unit[unit]);
  this_unit->toclen--;

  return true;
}

bool psio_tocdel(unsigned int unit, const char *key) {
  return _default_psio_lib_->tocdel(unit,key);
}

}
