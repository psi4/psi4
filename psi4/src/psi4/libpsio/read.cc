/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
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

#include <cstdlib>
#include <cstring>
#include "psi4/pragma.h"
PRAGMA_WARNING_PUSH
PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
#include <memory>
PRAGMA_WARNING_POP
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"

namespace psi {

void PSIO::read(size_t unit, const char *key, char *buffer, size_t size, psio_address start, psio_address *end) {
    psio_ud *this_unit;
    psio_tocentry *this_entry;
    psio_address start_toc, start_data, end_data; /* global addresses */
    size_t tocentry_size;

    if (!open_check(unit)) {
        fprintf(stderr, "PSIO_ERROR: Must open file %zu before reading it\n", unit);
        psio_error(unit, PSIO_ERROR_UNOPENED);
    }

    this_unit = &(psio_unit[unit]);

    /* Find the entry in the TOC */
    this_entry = tocscan(unit, key);

    tocentry_size = sizeof(psio_tocentry) - 2 * sizeof(psio_tocentry *);

    if (this_entry == nullptr) {
        fprintf(stderr, "PSIO_ERROR: Can't find TOC Entry %s\n", key);
        psio_error(unit, PSIO_ERROR_NOTOCENT);
    } else {
        /* Compute the global starting page and offset for the data */
        start_toc = this_entry->sadd;
        start_data = psio_get_address(start_toc, tocentry_size);
        start_data = psio_get_global_address(start_data, start);

        /* Make sure the block starts and ends within the entry */
        if (start_data.page > this_entry->eadd.page) {
            fprintf(stderr, "PSIO_ERROR: Start page %zu > this entry end page %zu\n", start_data.page,
                    this_entry->eadd.page);
            psio_error(unit, PSIO_ERROR_BLKSTART);
        } else if ((start_data.page == this_entry->eadd.page) && (start_data.offset > this_entry->eadd.offset)) {
            fprintf(stderr, "PSIO_ERROR: Start data offset %zu > this entry end address offset %zu\n",
                    start_data.offset, this_entry->eadd.offset);
            psio_error(unit, PSIO_ERROR_BLKSTART);
        }

        end_data = psio_get_address(start_data, size);
        if ((end_data.page > this_entry->eadd.page))
            psio_error(unit, PSIO_ERROR_BLKEND);
        else if ((end_data.page == this_entry->eadd.page) && (end_data.offset > this_entry->eadd.offset))
            psio_error(unit, PSIO_ERROR_BLKEND);

        /* Update end (an entry-relative address) for the caller */
        *end = psio_get_address(start, size);
    }

    /* Now read the actual data from the unit */
    rw(unit, buffer, start_data, size, 0);

#ifdef PSIO_STATS
    psio_readlen[unit] += size;
#endif
}

/*!
 ** PSIO_READ(): Reads data from within a TOC entry from a PSI file.
 **
 **  \param unit   = The PSI unit number used to identify the file to all
 **                  read and write functions.
 **  \param key    = The TOC keyword identifying the desired entry.
 **  \param buffer = The buffer to store the data as it is read.
 **  \param size   = The number of bytes to read.
 **  \param start  = The entry-relative starting page/offset of the desired data.
 **  \param end    = A pointer to the entry-relative page/offset for the next
 **                  byte after the end of the read request.
 **
 ** \ingroup PSIO
 */

int psio_read(size_t unit, const char *key, char *buffer, size_t size, psio_address start, psio_address *end) {
    _default_psio_lib_->read(unit, key, buffer, size, start, end);
    return 1;
}

}  // namespace psi
