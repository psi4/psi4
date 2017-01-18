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

#include <unistd.h>
#include "psi4/libpsio/psio.h"
#include "psi4/psi4-dec.h"
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
