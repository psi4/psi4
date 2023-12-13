/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

#include <cstring>
#include "psi4/pragma.h"
#include <memory>
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"

namespace psi {

psio_tocentry *PSIO::tocscan(size_t unit, const char *key) {
    psio_tocentry *this_entry;

    if (key == nullptr) return (nullptr);

    if ((strlen(key) + 1) > PSIO_KEYLEN) psio_error(unit, PSIO_ERROR_KEYLEN);

    bool already_open = open_check(unit);
    if (!already_open) open(unit, PSIO_OPEN_OLD);

    this_entry = psio_unit[unit].toc;

    while (this_entry != nullptr) {
        if (!strcmp(this_entry->key, key)) {
            if (!already_open) close(unit, 1);  // keep
            return (this_entry);
        }
        this_entry = this_entry->next;
    }

    if (!already_open) close(unit, 1);  // keep
    return (nullptr);
}

/*!
 ** PSIO_TOCSCAN(): Scans the TOC for a particular keyword and returns either
 ** a pointer to the entry or nullptr to the caller.
 **
 ** \ingroup PSIO
 */

psio_tocentry *psio_tocscan(size_t unit, const char *key) { return _default_psio_lib_->tocscan(unit, key); }

bool PSIO::tocentry_exists(size_t unit, const char *key) {
    psio_tocentry *this_entry;

    if (key == nullptr) return (true);

    if ((strlen(key) + 1) > PSIO_KEYLEN) psio_error(unit, PSIO_ERROR_KEYLEN);

    bool already_open = open_check(unit);
    if (!already_open) open(unit, PSIO_OPEN_OLD);

    this_entry = psio_unit[unit].toc;

    while (this_entry != nullptr) {
        if (!strcmp(this_entry->key, key)) {
            if (!already_open) close(unit, 1);  // keep
            return (true);
        }
        this_entry = this_entry->next;
    }

    if (!already_open) close(unit, 1);  // keep
    return (false);
}

/*!
 ** PSIO_TOCSCAN(): Scans the TOC for a particular keyword and returns either
 ** a pointer to the entry or nullptr to the caller.
 **
 ** \ingroup PSIO
 */

bool psio_tocentry_exists(size_t unit, const char *key) { return _default_psio_lib_->tocentry_exists(unit, key); }

}  // namespace psi
