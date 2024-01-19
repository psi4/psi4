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

#include <cstring>
#include <cstdlib>
#include "psi4/pragma.h"
#include <memory>
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
namespace psi {

/*!
 ** PSIO_TOCCLEAN(): Delete all TOC entries after the given key.
 ** If a blank key is given, the entire TOC will be wiped.
 **
 ** \ingroup PSIO
 */
void PSIO::tocclean(size_t unit, const char *key) {
    psio_tocentry *this_entry, *last_entry, *prev_entry;
    psio_ud *this_unit;

    this_unit = &(psio_unit[unit]);

    this_entry = tocscan(unit, key);
    if (this_entry == nullptr) {
        if (!strcmp(key, ""))
            this_entry = this_unit->toc;
        else {
            fprintf(stderr, "PSIO_ERROR: Can't find TOC Entry %s in unit %zu\n", key, unit);
            psio_error(unit, PSIO_ERROR_NOTOCENT);
        }
    } else
        this_entry = this_entry->next;

    /* Get the end of the TOC and work backwards */
    last_entry = toclast(unit);

    while ((last_entry != this_entry) && (last_entry != nullptr)) {
        /* Now free all the remaining members */
        prev_entry = last_entry->last;
        free(last_entry);
        last_entry = prev_entry;
        this_unit->toclen--;
    }

    /* Update on disk */
    wt_toclen(unit, this_unit->toclen);
    tocwrite(unit);
}
}  // namespace psi
