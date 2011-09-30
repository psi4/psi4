/*!
 ** \file
 ** \ingroup PSIO
 */

#include <cstdio>
#include <cstdlib>
#include <exception.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>

namespace psi {

  /*!
   ** \ingroup PSIO
   **
   ** PSIO_ERROR(): Print out an error message for libpsio.
   **
   ** \param unit   = file number
   ** \param errval = error code (defined symbolically, PSIO_ERROR_XXX)
   **
   */
  void psio_error(unsigned int unit, unsigned int errval) {
    int i;

    fprintf(stderr, "PSIO_ERROR: unit = %d\n", unit);
    /* Try to save the TOCs for all open units */
    /* psio_tocwrite() does not call psio_error() so this is OK */
    for (i=0; i < PSIO_MAXUNIT; i++)
      psio_tocwrite(i);

    switch (errval) {
      case PSIO_ERROR_INIT:
        fprintf(stderr, "PSIO_ERROR: %d (I/O inititalization failed)\n", PSIO_ERROR_INIT);
        break;
      case PSIO_ERROR_DONE:
        fprintf(stderr, "PSIO_ERROR: %d (I/O cleanup failed)\n", PSIO_ERROR_DONE);
        break;
      case PSIO_ERROR_MAXVOL:
        fprintf(stderr, "PSIO_ERROR: %d (maximum number of volumes exceeded)\n", PSIO_ERROR_MAXVOL);
        break;
      case PSIO_ERROR_NOVOLPATH:
        fprintf(stderr, "PSIO_ERROR: %d (no volume path given)\n", PSIO_ERROR_NOVOLPATH);
        break;
      case PSIO_ERROR_IDENTVOLPATH:
        fprintf(stderr, "PSIO_ERROR: %d (two identical volume paths)\n", PSIO_ERROR_IDENTVOLPATH);
        break;
      case PSIO_ERROR_OPEN:
        fprintf(stderr, "PSIO_ERROR: %d (file not open or open call failed)\n", PSIO_ERROR_OPEN);
        break;
      case PSIO_ERROR_REOPEN:
        fprintf(stderr, "PSIO_ERROR: %d (file is already open)\n", PSIO_ERROR_REOPEN);
        break;
      case PSIO_ERROR_CLOSE:
        fprintf(stderr, "PSIO_ERROR: %d (file close failed)\n", PSIO_ERROR_CLOSE);
        break;
      case PSIO_ERROR_RECLOSE:
        fprintf(stderr, "PSIO_ERROR: %d (file is already closed)\n", PSIO_ERROR_RECLOSE);
        break;
      case PSIO_ERROR_OSTAT:
        fprintf(stderr, "PSIO_ERROR: %d (invalid status flag for file open)\n", PSIO_ERROR_OSTAT);
        break;
      case PSIO_ERROR_LSEEK:
        fprintf(stderr, "PSIO_ERROR: %d (lseek failed)\n", PSIO_ERROR_LSEEK);
        break;
      case PSIO_ERROR_NOTOCENT:
        fprintf(stderr, "PSIO_ERROR: %d (no such TOC entry)\n", PSIO_ERROR_NOTOCENT);
        break;
      case PSIO_ERROR_TOCENTSZ:
        fprintf(stderr, "PSIO_ERROR: %d (TOC entry size mismatch)\n", PSIO_ERROR_TOCENTSZ);
        break;
      case PSIO_ERROR_KEYLEN:
        fprintf(stderr, "PSIO_ERROR: %d (TOC key too long)\n", PSIO_ERROR_KEYLEN);
        break;
      case PSIO_ERROR_BLKSIZ:
        fprintf(stderr, "PSIO_ERROR: %d (Requested blocksize invalid)\n", PSIO_ERROR_BLKSIZ);
        break;
      case PSIO_ERROR_BLKSTART:
        fprintf(stderr, "PSIO_ERROR: %d (Incorrect block start address)\n", PSIO_ERROR_BLKSTART);
        break;
      case PSIO_ERROR_BLKEND:
        fprintf(stderr, "PSIO_ERROR: %d (Incorrect block end address)\n", PSIO_ERROR_BLKEND);
        break;
      case PSIO_ERROR_MAXUNIT:
        fprintf(stderr, "PSIO_ERROR: %d (Maximum unit number exceeded)\n", PSIO_ERROR_MAXUNIT);
        fprintf(stderr, "Open failed because unit %d exceeds ", unit);
        fprintf(stderr, "PSIO_MAXUNIT = %d.\n", PSIO_MAXUNIT);
        break;
    }
    fflush(stderr);
    throw PSIEXCEPTION("PSIO Error");
    //exit(PSIO::_error_exit_code_);
  }

}

