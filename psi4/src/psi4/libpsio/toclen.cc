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
#include <cstdlib>
#include <cerrno>
#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/psi4-dec.h"
namespace psi {

/// @brief Compute the length of the TOC for a given unit using the in-core TOC list
/// @param unit : file unit number
/// @return length of the TOC for a given unit
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
    if (!open_check(unit)) psio_error(unit, PSIO_ERROR_UNOPENED);
    const auto stream = psio_unit[unit].vol[0].stream;
    const auto errcod = SYSTEM_LSEEK(stream, 0L, SEEK_SET);
    const int saved_errno = errno;
    if (errcod == -1) {
        const std::string errmsg =
            psio_compose_err_msg("LSEEK failed.", "Cannot seek vol[0] to its beginning", unit, saved_errno);
        psio_error(unit, PSIO_ERROR_LSEEK, errmsg);
    }
}

/// @brief Read the length of the TOC for a given unit directly from the file. Note that we do not exit if the read
/// request of the toclen from the file fails to read the expected number of bytes. This is because the request may be
/// to an new file for which the toclen has not yet been written. (We allow the user open files with status
/// PSIO_OPEN_OLD even if they don't exist, because sometimes you can't know this in advance.)
/// @param unit : file unit number to read TOC length from
/// @return length of the TOC for a given unit
size_t PSIO::rd_toclen(const size_t unit) {
    if (!open_check(unit)) psio_error(unit, PSIO_ERROR_UNOPENED);
    // Seek to the beginning
    rewind_toclen(unit);
    // Read the value
    size_t len;
    const auto stream = psio_unit[unit].vol[0].stream;
    const auto errcod = SYSTEM_READ(stream, (char *)&len, sizeof(size_t));
    const int saved_errno = errno;
    if (errcod != sizeof(size_t)) {
        if (errcod == -1) {
            const std::string errmsg = psio_compose_err_msg(
                "READ failed.", "Error in PSIO::rd_toclen()! Cannot read TOC length", unit, saved_errno);
            psio_error(unit, PSIO_ERROR_READ, errmsg);
        }
        return (0);  // assume that all is well (see comments above)
    }
    return (len);
}

/// @brief Write the length of the TOC for a given unit directly to the file.
/// @param unit : file unit number to write TOC length to
/// @param len  : length value to write
void PSIO::wt_toclen(const size_t unit, const size_t len) {
    if (!open_check(unit)) psio_error(unit, PSIO_ERROR_UNOPENED);
    // Seek to the beginning
    rewind_toclen(unit);
    // Write the value
    const auto stream = psio_unit[unit].vol[0].stream;
    const auto errcod = SYSTEM_WRITE(stream, (char *)&len, sizeof(size_t));
    const int saved_errno = errno;
    if (errcod != sizeof(size_t)) {
        const std::string errmsg = psio_compose_err_msg(
            "WRITE failed.", "Error in PSIO::wt_toclen()! Cannot write TOC length", unit, saved_errno);
        psio_error(unit, PSIO_ERROR_WRITE, errmsg);
    }
}

size_t psio_rd_toclen(size_t unit) { return _default_psio_lib_->rd_toclen(unit); }

}  // namespace psi
