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
#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/psi4-dec.h"
namespace psi {

unsigned int PSIO::toclen(unsigned int unit) {
  unsigned int len=0;
  psio_tocentry *this_entry;

  this_entry = psio_unit[unit].toc;

  while (this_entry != NULL) {
    ++len;
    this_entry = this_entry->next;
  }

  return (len);
}

ULI PSIO::rd_toclen(unsigned int unit) {
  int errcod, stream;
  psio_ud *this_unit;
  ULI len;

  this_unit = &(psio_unit[unit]);

  /* Seek vol[0] to its beginning */
  stream = this_unit->vol[0].stream;

    errcod = ::lseek(stream, 0L, SEEK_SET);

  if (errcod == -1)
    psio_error(unit, PSIO_ERROR_LSEEK);

  /* Read the value */

    errcod = ::read(stream, (char *) &len, sizeof(ULI));


  if(errcod != sizeof(ULI)) return(0); /* assume that all is well (see comments above) */

  return(len);
}

void PSIO::wt_toclen(unsigned int unit, ULI len) {
  int errcod, stream;
  psio_ud *this_unit;

  this_unit = &(psio_unit[unit]);

  /* Seek vol[0] to its beginning */
  stream = this_unit->vol[0].stream;

    errcod = ::lseek(stream, 0L, SEEK_SET);

  if (errcod == -1) {
    ::fprintf(stderr, "Error in PSIO_WT_TOCLEN()!\n");
    exit(_error_exit_code_);
  }

  /* Write the value */

    errcod = ::write(stream, (char *) &len, sizeof(ULI));

  if(errcod != sizeof(ULI)) {
    ::fprintf(stderr, "PSIO_ERROR: Failed to write toclen to unit %d.\n", unit);
    fflush(stderr);
    throw PSIEXCEPTION("PSIO Error");
  }
}

unsigned long int psio_rd_toclen(unsigned int unit) {
  return _default_psio_lib_->rd_toclen(unit);
}

}
