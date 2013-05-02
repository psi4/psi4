/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

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
    if (WorldComm->me() == 0)
        errcod = lseek(stream, (ULI) 0, SEEK_SET);
    WorldComm->bcast(&(errcod), 1, 0);
    //WorldComm->raw_bcast(&errcod, sizeof(int), 0);
    if (errcod == -1)
      return (errcod);
    
    /* lseek() through large chunks of the file to avoid offset overflows */
    for (; page > bignum; page -= bignum) {
      total_offset = PSIO_BIGNUM * PSIO_PAGELEN;
      if (WorldComm->me() == 0)
          errcod = lseek(stream, total_offset, SEEK_CUR);
      WorldComm->bcast(&(errcod), 1, 0);
      //WorldComm->raw_bcast(&errcod, sizeof(int), 0);
      if (errcod == -1)
        return (errcod);
    }
    
    /* Now compute the final offset including the page-relative term */
    total_offset = (ULI) page/numvols; /* This should truncate */
    total_offset *= PSIO_PAGELEN;
    total_offset += offset; /* Add the page-relative term */
    if (WorldComm->me() == 0)
        errcod = lseek(stream, total_offset, SEEK_CUR);
    WorldComm->bcast(&(errcod), 1, 0);
    //WorldComm->raw_bcast(&errcod, sizeof(int), 0);
    if (errcod == -1)
      return (errcod);
    
    return (0);
  }

}

