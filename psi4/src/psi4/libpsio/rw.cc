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

#include <cstdio>
#ifdef _MSC_VER
#include <io.h>
#define SYSTEM_READ ::_read
#define SYSTEM_WRITE ::_write
#else
#include <unistd.h>
#define SYSTEM_READ ::read
#define SYSTEM_WRITE ::write
#endif
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/psi4-dec.h"

namespace psi {

void PSIO::rw(size_t unit, char *buffer, psio_address address, size_t size, int wrt) {
    int errcod;
    size_t i;
    size_t errcod_uli;
    size_t page, offset;
    size_t buf_offset;
    size_t this_page, this_page_max, this_page_total;
    size_t first_vol, this_vol, numvols;
    size_t bytes_left, num_full_pages;
    psio_ud *this_unit;

    this_unit = &(psio_unit[unit]);
    numvols = this_unit->numvols;
    page = address.page;
    offset = address.offset;

    /* Seek all volumes to correct starting positions */
    first_vol = page % numvols;
    errcod = psio_volseek(&(this_unit->vol[first_vol]), page, offset, numvols);
    if (errcod == -1) psio_error(unit, PSIO_ERROR_LSEEK);
    for (i = 1, this_page = page + 1; i < numvols; i++, this_page++) {
        this_vol = this_page % numvols;
        errcod = psio_volseek(&(this_unit->vol[this_vol]), this_page, (size_t)0, numvols);
        if (errcod == -1) psio_error(unit, PSIO_ERROR_LSEEK);
    }

    /* Number of bytes left on the first page */
    this_page_max = PSIO_PAGELEN - offset;

    /* If we have enough room on this page, use it */
    if (size <= this_page_max)
        this_page_total = size;
    else
        this_page_total = this_page_max;
    buf_offset = 0;
    if (wrt) {
        errcod_uli = SYSTEM_WRITE(this_unit->vol[first_vol].stream, &(buffer[buf_offset]), this_page_total);
        if (errcod_uli != this_page_total) psio_error(unit, PSIO_ERROR_WRITE);
    } else {
        errcod_uli = SYSTEM_READ(this_unit->vol[first_vol].stream, &(buffer[buf_offset]), this_page_total);
        if (errcod_uli != this_page_total) psio_error(unit, PSIO_ERROR_READ);
    }

    /* Total number of bytes remaining to be read/written */
    bytes_left = size - this_page_total;

    /* Read/Write all the full pages */
    num_full_pages = bytes_left / PSIO_PAGELEN;
    buf_offset += this_page_total;
    for (i = 0, this_page = page + 1; i < num_full_pages; i++, this_page++) {
        this_vol = this_page % numvols;
        this_page_total = PSIO_PAGELEN;
        if (wrt) {
            errcod_uli = SYSTEM_WRITE(this_unit->vol[this_vol].stream, &(buffer[buf_offset]), this_page_total);
            if (errcod_uli != this_page_total) psio_error(unit, PSIO_ERROR_WRITE);
        } else {
            errcod_uli = SYSTEM_READ(this_unit->vol[this_vol].stream, &(buffer[buf_offset]), this_page_total);
            if (errcod_uli != this_page_total) psio_error(unit, PSIO_ERROR_READ);
        }
        buf_offset += this_page_total;
    }

    /* Read/Write the final partial page */
    bytes_left -= num_full_pages * PSIO_PAGELEN;
    this_vol = this_page % numvols;
    if (bytes_left) {
        if (wrt) {
            errcod_uli = SYSTEM_WRITE(this_unit->vol[this_vol].stream, &(buffer[buf_offset]), bytes_left);
            if (errcod_uli != bytes_left) psio_error(unit, PSIO_ERROR_WRITE);
        } else {
            errcod_uli = SYSTEM_READ(this_unit->vol[this_vol].stream, &(buffer[buf_offset]), bytes_left);
            if (errcod_uli != bytes_left) psio_error(unit, PSIO_ERROR_READ);
        }
    }
}

/*!
 ** PSIO_RW(): Central function for all reads and writes on a PSIO unit.
 **
 ** \params unit    = The PSI unit number.
 ** \params buffer  = The buffer containing the bytes for the read/write event.
 ** \params address = the PSIO global address for the start of the read/write.
 ** \params size    = The number of bytes to read/write.
 ** \params         = Indicates if the call is to read (0) or write (0) the input data.
 **
 ** \ingroup PSIO
 */

int psio_rw(size_t unit, char *buffer, psio_address address, size_t size, int wrt) {
    _default_psio_lib_->rw(unit, buffer, address, size, wrt);
    return 1;
}

}  // namespace psi
