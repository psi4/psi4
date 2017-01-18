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
 #include "psi4/pragma.h"
 PRAGMA_WARNING_PUSH
 PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
 #include <memory>
 PRAGMA_WARNING_POP
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/psi4-dec.h"
namespace psi {



void PSIO::tocprint(unsigned int unit) {
  psio_tocentry *this_entry;

  bool already_open = open_check(unit);
  if(!already_open) open(unit, PSIO_OPEN_OLD);

  this_entry = psio_unit[unit].toc;

  outfile->Printf( "\nTable of Contents for Unit %5u\n", unit);
  outfile->Printf(
          "----------------------------------------------------------------------------\n");
  outfile->Printf(
          "Key                                   Spage    Soffset      Epage    Eoffset\n");
  outfile->Printf(
          "----------------------------------------------------------------------------\n");

  while (this_entry != NULL) {
    outfile->Printf("%-32s %10lu %10lu %10lu %10lu\n", this_entry->key,
            this_entry->sadd.page, this_entry->sadd.offset,
            this_entry->eadd.page, this_entry->eadd.offset);
    this_entry = this_entry->next;
  }
  outfile->Printf("\n");


  if(!already_open) close(unit, 1); // keep
}

  /*!
   ** PSIO_TOCPRINT(): Print the table of contents for the given unit
   **
   ** \ingroup PSIO
   */

  void psio_tocprint(unsigned int unit) {
    return _default_psio_lib_->tocprint(unit);
  }


}
