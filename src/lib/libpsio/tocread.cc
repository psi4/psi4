/*!
 \file
 \ingroup PSIO
 */

#include <unistd.h>
#include <cstdlib>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

namespace psi {

void PSIO::tocread(unsigned int unit) {
  unsigned int i;
  int entry_size;
  psio_ud *this_unit;
  psio_tocentry *last_entry, *this_entry;
  psio_address address;
  
  this_unit = &(psio_unit[unit]);
  entry_size = sizeof(psio_tocentry) - 2*sizeof(psio_tocentry *);
  
  /* This wasn't doing anything. */
  //if (!open_check(unit))
  //  ;
  
  /* grab the number of records */
  this_unit->toclen = rd_toclen(unit);
  
  /* Malloc room for the TOC */
  if (this_unit->toclen) {
    this_unit->toc = (psio_tocentry *) malloc(sizeof(psio_tocentry));
    this_entry = this_unit->toc;
    this_entry->last = NULL;
    for (i=1; i < this_unit->toclen; i++) {
      last_entry = this_entry;
      this_entry = (psio_tocentry *) malloc(sizeof(psio_tocentry));
      last_entry->next = this_entry;
      this_entry->last = last_entry;
    }
    this_entry->next = NULL;
  }
  
  /* Read the TOC entry-by-entry */
  this_entry = this_unit->toc;
  address = psio_get_address(PSIO_ZERO, sizeof(ULI)); /* start one ULI after the top of the file */
  for (i=0; i < this_unit->toclen; i++) {
    rw(unit, (char *) this_entry, address, entry_size, 0);
    address = this_entry->eadd;
    this_entry = this_entry->next;
  }
}

}

