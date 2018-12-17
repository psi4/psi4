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

#include "psi4/libpsio/psio.h"
#include "psi4/psi4-dec.h"

#ifdef _MSC_VER
#include <io.h>
#define SYSTEM_LSEEK ::_lseeki64
#else
#include <unistd.h>
#define SYSTEM_LSEEK ::lseek
#endif

/* This is strictly used to avoid overflow errors on lseek() calls */
#define PSIO_BIGNUM 10000

namespace psi {
/*!
 ** PSIO_VOLSEEK()
 **
 ** \ingroup PSIO
 */
int psio_volseek(psio_vol *vol, size_t page, size_t offset, size_t numvols) {

    int stream = vol->stream;

    /* Set file pointer to beginning of file */
    if (SYSTEM_LSEEK(stream, 0, SEEK_SET) == -1)
        return -1;

    /* lseek() through large chunks of the file to avoid offset overflows */
    size_t bignum = PSIO_BIGNUM * numvols;
    for (; page > bignum; page -= bignum)
        if (SYSTEM_LSEEK(stream, PSIO_BIGNUM * PSIO_PAGELEN, SEEK_CUR) == -1)
            return -1;

    /* Now compute the final offset including the page-relative term */
    size_t final_offset = (page / numvols) * PSIO_PAGELEN + offset;
    if (SYSTEM_LSEEK(stream, final_offset, SEEK_CUR) == -1)
        return -1;

    return 0;
}

}  // namespace psi
