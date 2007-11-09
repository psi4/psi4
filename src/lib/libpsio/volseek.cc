/*!
 \file volseek.cc
 \ingroup (PSIO)
 */

#include <unistd.h>
#include <libpsio/psio.h>

/* This is strictly used to avoid overflow errors on lseek() calls */
#define PSIO_BIGNUM 10000

extern "C" {
  /*!
   ** PSIO_VOLSEEK()
   **
   ** \ingroup (PSIO)
   */
  int psio_volseek(psio_vol *vol, ULI page, ULI offset, ULI numvols) {
    int stream, errcod;
    ULI bignum, total_offset;
    
    bignum = PSIO_BIGNUM*numvols;
    
    stream = vol->stream;
    
    /* Set file pointer to beginning of file */
    errcod = lseek(stream, (ULI) 0, SEEK_SET);
    if (errcod == -1)
      return (errcod);
    
    /* lseek() through large chunks of the file to avoid offset overflows */
    for (; page > bignum; page -= bignum) {
      total_offset = PSIO_BIGNUM * PSIO_PAGELEN;
      errcod = lseek(stream, total_offset, SEEK_CUR);
      if (errcod == -1)
        return (errcod);
    }
    
    /* Now compute the final offset including the page-relative term */
    total_offset = (ULI) page/numvols; /* This should truncate */
    total_offset *= PSIO_PAGELEN;
    total_offset += offset; /* Add the page-relative term */
    errcod = lseek(stream, total_offset, SEEK_CUR);
    if (errcod == -1)
      return (errcod);
    
    return (0);
  }
}
