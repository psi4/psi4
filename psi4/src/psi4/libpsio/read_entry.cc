/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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
#include "psi4/libpsio/psio.hpp"

namespace psi {

void PSIO::read_entry(size_t unit, const char *key, char *buffer, size_t size) {
    psio_address end;
    read(unit, key, buffer, size, PSIO_ZERO, &end);
}

/*!
 ** PSIO_READ_ENTRY(): Reads an entire TOC entry from a PSI file.
 **
 **  \param unit   = The PSI unit number used to identify the file to all read
 **                  and write functions.
 **  \param key    = The TOC keyword identifying the desired entry.
 **  \param buffer = The buffer to store the data as it is read.
 **  \param size   = The number of bytes to read.
 **
 ** Note that the value of size is not directly compared to the actual
 ** size of the entry, but care is taken to ensure that the end of the
 ** entry is not surpassed.
 **
 ** \ingroup PSIO
 */

int psio_read_entry(size_t unit, const char *key, char *buffer, size_t size) {
    psio_address end;
    return psio_read(unit, key, buffer, size, PSIO_ZERO, &end);
}

}  // namespace psi
