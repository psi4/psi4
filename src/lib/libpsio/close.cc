/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

/*! \defgroup PSIO libpsio: The PSI I/O Library */

/*!
 ** \file
 ** \ingroup PSIO
 */

#include <unistd.h>
#include <cstring>
#include <cstdlib>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include "psi4-dec.h"
#include "../libparallel2/Communicator.h"
#include "../libparallel2/ParallelEnvironment.h"
namespace psi {

void PSIO::close(unsigned int unit, int keep) {
  unsigned int i;
  psio_ud *this_unit;
  psio_tocentry *this_entry, *next_entry;

  this_unit = &(psio_unit[unit]);

  /* First check to see if this unit is already closed */
  if (this_unit->vol[0].stream == -1)
    psio_error(unit, PSIO_ERROR_RECLOSE);

  /* Dump the current TOC back out to disk */
  tocwrite(unit);

  /* Free the TOC */
  this_entry = this_unit->toc;
  for (i=0; i < this_unit->toclen; i++) {
    next_entry = this_entry->next;
    free(this_entry);
    this_entry = next_entry;
  }

  /* Close each volume (remove if necessary) and free the path */
  for (i=0; i < this_unit->numvols; i++) {
    int errcod;
    boost::shared_ptr<const LibParallel::Communicator> Comm=
          WorldComm->GetComm();
    if (Comm->Me() == 0) {
      errcod = ::close(this_unit->vol[i].stream);
    }
    Comm->Bcast(&errcod, 1, 0);
    if (errcod == -1)
      psio_error(unit,PSIO_ERROR_CLOSE);
    /* Delete the file completely if requested */
    if(!keep) unlink(this_unit->vol[i].path);
    PSIOManager::shared_object()->close_file(std::string(this_unit->vol[i].path), unit, (keep ? true : false));

    free(this_unit->vol[i].path);
    this_unit->vol[i].path = NULL;
    this_unit->vol[i].stream = -1;
  }

  /* Reset the global page stats to zero */
  this_unit->numvols = 0;
  this_unit->toclen = 0;
  this_unit->toc = NULL;
}

int psio_close(unsigned int unit, int keep) {
  _default_psio_lib_->close(unit, keep);
  return 0;
}

}