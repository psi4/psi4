/*!
 \file
 \ingroup PSIO
 */

#include <cstdio>
#include <unistd.h>
#include <cstdlib>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

namespace psi {

unsigned int PSIO::toclen(unsigned int unit) {
  unsigned int toclen=0;
  psio_tocentry *this_entry;
  
  this_entry = psio_unit[unit].toc;
  
  while (this_entry != NULL) {
    ++toclen;
    this_entry = this_entry->next;
  }
  
  return (toclen);
}

ULI PSIO::rd_toclen(unsigned int unit) {
  int errcod, stream;
  psio_ud *this_unit;
  ULI toclen;
  
  this_unit = &(psio_unit[unit]);
  
  /* Seek vol[0] to its beginning */
  stream = this_unit->vol[0].stream;
  errcod = lseek(stream, 0L, SEEK_SET);
  if (errcod == -1)
    psio_error(unit, PSIO_ERROR_LSEEK);
  
  /* Read the value */
  errcod =:: read(stream, (char *) &toclen, sizeof(ULI));
  if(errcod != sizeof(ULI)) return(0); /* assume that all is well (see comments above) */

  return(toclen);
}

void PSIO::wt_toclen(unsigned int unit, ULI toclen) {
  int errcod, stream;
  psio_ud *this_unit;
  
  this_unit = &(psio_unit[unit]);
  
  /* Seek vol[0] to its beginning */
  stream = this_unit->vol[0].stream;
  errcod = lseek(stream, 0L, SEEK_SET);
  if (errcod == -1) {
    fprintf(stderr, "Error in PSIO_WT_TOCLEN()!\n");
    exit(_error_exit_code_);
  }
  
  /* Write the value */
  errcod =:: write(stream, (char *) &toclen, sizeof(ULI));
  if(errcod != sizeof(ULI)) {
    fprintf(stderr, "PSIO_ERROR: Failed to write toclen to unit %d.\n", unit);
    psio_error(unit,PSIO_ERROR_WRITE);
  }
}

}

