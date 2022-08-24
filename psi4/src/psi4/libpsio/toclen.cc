/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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
#define SYSTEM_LSEEK ::_lseek
#define SYSTEM_READ ::_read
#define SYSTEM_WRITE ::_write
#else
#include <unistd.h>
#define SYSTEM_LSEEK ::lseek
#define SYSTEM_READ ::read
#define SYSTEM_WRITE ::write
#endif
#include <cstdlib>
#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/psi4-dec.h"
namespace psi {

size_t PSIO::toclen(size_t unit) {
    size_t len = 0;
    psio_tocentry *this_entry;

    this_entry = psio_unit[unit].toc;

    while (this_entry != nullptr) {
        ++len;
        this_entry = this_entry->next;
    }

    return (len);
}

size_t PSIO::rd_toclen(size_t unit) {
    int errcod, stream;
    psio_ud *this_unit;
    size_t len;

    this_unit = &(psio_unit[unit]);

    /* Seek vol[0] to its beginning */
    stream = this_unit->vol[0].stream;

    errcod = SYSTEM_LSEEK(stream, 0L, SEEK_SET);

    if (errcod == -1){
        perror("LSEEK failed. Error description from the OS: ");
        fflush(stderr);
        psio_error(unit, PSIO_ERROR_LSEEK);
    }

    /* Read the value */

    errcod = SYSTEM_READ(stream, (char *)&len, sizeof(size_t));

    if (errcod != sizeof(size_t)){
        perror("READ failed in rd_toclen(). Error description from the OS: ");
        fflush(stderr);
        return (0); /* assume that all is well (see comments above) */
    }

    return (len);
}

void PSIO::wt_toclen(size_t unit, size_t len) {
    int errcod, stream;
    psio_ud *this_unit;

    this_unit = &(psio_unit[unit]);

    /* Seek vol[0] to its beginning */
    stream = this_unit->vol[0].stream;

    errcod = SYSTEM_LSEEK(stream, 0L, SEEK_SET);

    if (errcod == -1) {
        perror("LSEEK failed. Error description from the OS: ");
        ::fprintf(stderr, "Error in PSIO_WT_TOCLEN()! Cannot seek vol[0] to its beginning, unit %zu.\n", unit);
        fflush(stderr);
        throw PSIEXCEPTION("PSIO Error");
    }

    /* Write the value */

    errcod = SYSTEM_WRITE(stream, (char *)&len, sizeof(size_t));

    if (errcod != sizeof(size_t)) {
        perror("WRITE failed. Error description from the OS: ");
        ::fprintf(stderr, "PSIO_ERROR: Failed to write toclen to unit %zu.\n", unit);
        fflush(stderr);
        throw PSIEXCEPTION("PSIO Error");
    }
}

size_t psio_rd_toclen(size_t unit) { return _default_psio_lib_->rd_toclen(unit); }

}  // namespace psi
