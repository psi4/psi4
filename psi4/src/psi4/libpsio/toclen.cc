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
#include "psi4/libpsi4util/PsiOutStream.h"
namespace psi {

size_t PSIO::toclen(const size_t unit) {
    size_t len = 0;
    psio_tocentry *this_entry = psio_unit[unit].toc;

    while (this_entry != nullptr) {
        ++len;
        this_entry = this_entry->next;
    }

    return (len);
}

/// @brief Seek the stream of the vol[0] of a unit to its beginning
/// @param unit : file unit number to rewind
void PSIO::rewind_toclen(const size_t unit) {
    const auto stream = psio_unit[unit].vol[0].stream;
    const auto errcod = SYSTEM_LSEEK(stream, 0L, SEEK_SET);
    const auto sys_errno = errno;
    if (errcod == -1) {
        std::string errmsg = "LSEEK failed. Error description from the OS: " + decode_errno(sys_errno);
        errmsg += "\nCannot seek vol[0] to its beginning, unit ";
        errmsg += std::to_string(unit) + ".\n";
        psio_error(unit, PSIO_ERROR_LSEEK, errmsg);
    }
}

size_t PSIO::rd_toclen(const size_t unit) {
    // Seek to the beginning
    rewind_toclen(unit);

    // Read the value
    size_t len;
    const auto stream = psio_unit[unit].vol[0].stream;
    const auto errcod = SYSTEM_READ(stream, (char *)&len, sizeof(size_t));
    const auto sys_errno = errno;
    if (errcod != sizeof(size_t)) {
#ifdef DEBUG
        if (errcod == -1) {
            std::string errmsg = "READ failed. Error description from the OS: " + decode_errno(sys_errno);
            errmsg += "\nError in PSIO::rd_toclen()! Cannot read TOC length, unit ";
            errmsg += std::to_string(unit) + ".\n";
            errmsg += "Assuming a length of zero and continuing...\n";
            outfile->Printf(errmsg);
        }
#endif
        return (0);  // assume that all is well (see comments in psio.hpp)
    }
    return (len);
}

void PSIO::wt_toclen(const size_t unit, const size_t len) {
    // Seek to the beginning
    rewind_toclen(unit);

    // Write the value
    const auto stream = psio_unit[unit].vol[0].stream;
    const auto errcod = SYSTEM_WRITE(stream, (char *)&len, sizeof(size_t));
    const auto sys_errno = errno;
    if (errcod != sizeof(size_t)) {
        std::string errmsg = "WRITE failed. Error description from the OS: " + decode_errno(sys_errno);
        errmsg += "\nError in PSIO::wt_toclen()! Cannot write TOC length, unit ";
        errmsg += std::to_string(unit) + ".\n";
        psio_error(unit, PSIO_ERROR_WRITE, errmsg);
    }
}

size_t psio_rd_toclen(size_t unit) { return _default_psio_lib_->rd_toclen(unit); }

}  // namespace psi
