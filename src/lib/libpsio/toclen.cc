/*!
 \file
 \ingroup PSIO
 */

#include <cstdio>
#include <unistd.h>
#include <cstdlib>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <libparallel/parallel.h>

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
  if (Communicator::world->me() == 0) {
    errcod = ::lseek(stream, 0L, SEEK_SET);
  }
  Communicator::world->bcast(&(errcod), 1, 0);
  //Communicator::world->raw_bcast(&(errcod), sizeof(int), 0);
  if (errcod == -1)
    psio_error(unit, PSIO_ERROR_LSEEK);
  
  /* Read the value */
  if (Communicator::world->me() == 0) {
    errcod = ::read(stream, (char *) &len, sizeof(ULI));
  }

  Communicator::world->bcast(&(errcod), 1, 0);
  Communicator::world->bcast(&(len), 1, 0);
  //Communicator::world->raw_bcast(&(errcod), sizeof(int), 0);
  //Communicator::world->raw_bcast(&(len), sizeof(ULI), 0);

  if(errcod != sizeof(ULI)) return(0); /* assume that all is well (see comments above) */

  return(len);
}

void PSIO::wt_toclen(unsigned int unit, ULI len) {
  int errcod, stream;
  psio_ud *this_unit;
  
  this_unit = &(psio_unit[unit]);
  
  /* Seek vol[0] to its beginning */
  stream = this_unit->vol[0].stream;
  if (Communicator::world->me() == 0) {
    errcod = ::lseek(stream, 0L, SEEK_SET);
  }
  Communicator::world->bcast(&(errcod), 1, 0);
  //Communicator::world->raw_bcast(&(errcod), sizeof(int), 0);
  if (errcod == -1) {
    fprintf(stderr, "Error in PSIO_WT_TOCLEN()!\n");
    exit(_error_exit_code_);
  }
  
  /* Write the value */
  if (Communicator::world->me() == 0) {
    errcod = ::write(stream, (char *) &len, sizeof(ULI));
  }
  Communicator::world->bcast(&(errcod), 1, 0);
  //Communicator::world->raw_bcast(&(errcod), sizeof(int), 0);
  if(errcod != sizeof(ULI)) {
    fprintf(stderr, "PSIO_ERROR: Failed to write toclen to unit %d.\n", unit);
    psio_error(unit,PSIO_ERROR_WRITE);
  }
}

unsigned long int psio_rd_toclen(unsigned int unit) {
  return _default_psio_lib_->rd_toclen(unit);
}

}

