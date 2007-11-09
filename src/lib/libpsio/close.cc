/*! \defgroup PSIO libpsio: The PSI I/O Library */

/*!
 ** \file close.cc
 ** \ingroup (PSIO)
 */

#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

using namespace psi;

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
    
    if (::close(this_unit->vol[i].stream) == -1)
    psio_error(unit,PSIO_ERROR_CLOSE);

    /* Delete the file completely if requested */
    if(!keep) unlink(this_unit->vol[i].path);

    free(this_unit->vol[i].path);
    this_unit->vol[i].path = NULL;
    this_unit->vol[i].stream = -1;
  }

  /* Reset the global page stats to zero */
  this_unit->numvols = 0;
  this_unit->toclen = 0;
}

extern "C" {
  int psio_close(unsigned int unit, int keep) {
    _default_psio_lib_->close(unit, keep);
    return 0;
  }
}

