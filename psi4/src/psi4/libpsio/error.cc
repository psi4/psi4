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
 ** \file
 ** \ingroup PSIO
 */

#include <cstdio>
#include <cstdlib>
#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"

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

    fprintf(stderr, "PSIO_ERROR: unit = %d, errval = %d\n", unit, errval);
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
        fprintf(stderr, "PSIO_ERROR:\n");
        fprintf(stderr, "PSIO_ERROR: Check the location of your scratch directory which can be\n");
        fprintf(stderr, "PSIO_ERROR: specified via the $PSI_SCRATCH environment variable or in\n");
        fprintf(stderr, "PSIO_ERROR: the $HOME/.psi4rc file.\n");
        fprintf(stderr, "PSIO_ERROR:\n");
        fprintf(stderr, "PSIO_ERROR: Please note that the scratch directory must exist and be\n");
        fprintf(stderr, "PSIO_ERROR: writable by Psi4\n");
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
      case PSIO_ERROR_WRITE:
        fprintf(stderr, "PSIO_ERROR: %d (error writing to file)\n", PSIO_ERROR_WRITE);
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
