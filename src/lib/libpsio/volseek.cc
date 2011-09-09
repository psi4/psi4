/*!
 \file
 \ingroup PSIO
 */

#include <unistd.h>
#include <libpsio/psio.h>
#include <libparallel/parallel.h>

/* This is strictly used to avoid overflow errors on lseek() calls */
#define PSIO_BIGNUM 10000

namespace psi {
  /*!
   ** PSIO_VOLSEEK()
   **
   ** \ingroup PSIO
   */
  int psio_volseek(psio_vol *vol, ULI page, ULI offset, ULI numvols) {
    int stream, errcod;
    ULI bignum, total_offset;
    
    bignum = PSIO_BIGNUM*numvols;
    
    stream = vol->stream;
    
    /* Set file pointer to beginning of file */
    if (Communicator::world->me() == 0)
        errcod = lseek(stream, (ULI) 0, SEEK_SET);
    Communicator::world->bcast(&(errcod), 1, 0);
    //Communicator::world->raw_bcast(&errcod, sizeof(int), 0);
    if (errcod == -1)
      return (errcod);
    
    /* lseek() through large chunks of the file to avoid offset overflows */
    for (; page > bignum; page -= bignum) {
      total_offset = PSIO_BIGNUM * PSIO_PAGELEN;
      if (Communicator::world->me() == 0)
          errcod = lseek(stream, total_offset, SEEK_CUR);
      Communicator::world->bcast(&(errcod), 1, 0);
      //Communicator::world->raw_bcast(&errcod, sizeof(int), 0);
      if (errcod == -1)
        return (errcod);
    }
    
    /* Now compute the final offset including the page-relative term */
    total_offset = (ULI) page/numvols; /* This should truncate */
    total_offset *= PSIO_PAGELEN;
    total_offset += offset; /* Add the page-relative term */
    if (Communicator::world->me() == 0)
        errcod = lseek(stream, total_offset, SEEK_CUR);
    Communicator::world->bcast(&(errcod), 1, 0);
    //Communicator::world->raw_bcast(&errcod, sizeof(int), 0);
    if (errcod == -1)
      return (errcod);
    
    return (0);
  }

}

