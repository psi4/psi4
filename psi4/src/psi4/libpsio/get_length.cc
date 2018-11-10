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

namespace psi {
/*!
 ** PSIO_GET_LENGTH(): Given a start page and offset for two data sets,
 ** compute the number of bytes between them.  Note that eadd denotes the
 ** beginning of the next entry and not the end of the current entry.
 **
 ** \ingroup PSIO
 */

size_t psio_get_length(psio_address sadd, psio_address eadd) {
    size_t full_page_bytes;

    /* Number of bytes on fullpages */
    full_page_bytes = (eadd.page - sadd.page - 1) * PSIO_PAGELEN;

    /* Uh, size_t will NEVER be less than 0 */
    // if (full_page_bytes < 0) { /* We're on a single page */
    //  return (eadd.offset - sadd.offset);
    //} else if (full_page_bytes == 0) { /* We're on the next page */
    if (full_page_bytes == 0) { /* We're on the next page */
        return ((PSIO_PAGELEN - sadd.offset) + eadd.offset);
    } else {
        return ((PSIO_PAGELEN - sadd.offset) + full_page_bytes + eadd.offset);
    }
}

}  // namespace psi
