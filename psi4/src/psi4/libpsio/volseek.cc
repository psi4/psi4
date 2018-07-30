/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*!
 \file
 \ingroup PSIO
 */

#ifdef _MSC_VER
#include <io.h>
#define SYSTEM_LSEEK ::_lseek
#else
#include <unistd.h>
#define SYSTEM_LSEEK ::lseek
#endif
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
  int psio_volseek(psio_vol *vol, size_t page, size_t offset, size_t numvols) {
    int stream, errcod;
    size_t bignum, total_offset;

    bignum = PSIO_BIGNUM*numvols;

    stream = vol->stream;

    /* Set file pointer to beginning of file */
        errcod = SYSTEM_LSEEK(stream, (size_t) 0, SEEK_SET);
    if (errcod == -1)
      return (errcod);

    /* lseek() through large chunks of the file to avoid offset overflows */
    for (; page > bignum; page -= bignum) {
      total_offset = PSIO_BIGNUM * PSIO_PAGELEN;
          errcod = SYSTEM_LSEEK(stream, total_offset, SEEK_CUR);
      if (errcod == -1)
        return (errcod);
    }

    /* Now compute the final offset including the page-relative term */
    total_offset = (size_t) page/numvols; /* This should truncate */
    total_offset *= PSIO_PAGELEN;
    total_offset += offset; /* Add the page-relative term */
        errcod = SYSTEM_LSEEK(stream, total_offset, SEEK_CUR);
    if (errcod == -1)
      return (errcod);

    return (0);
  }

}
