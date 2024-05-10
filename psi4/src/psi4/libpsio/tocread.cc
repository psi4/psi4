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

#include <cstdlib>
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"

namespace psi {

void PSIO::tocread(size_t unit) {
    size_t i;
    int entry_size;
    psio_ud *this_unit;
    psio_tocentry *last_entry, *this_entry;
    psio_address address;

    this_unit = &(psio_unit[unit]);
    entry_size = sizeof(psio_tocentry) - 2 * sizeof(psio_tocentry *);

    /* This wasn't doing anything. */
    // if (!open_check(unit))
    //  ;

    /* grab the number of records */
    this_unit->toclen = rd_toclen(unit);

    /* Malloc room for the TOC */
    if (this_unit->toclen) {
        this_unit->toc = (psio_tocentry *)malloc(sizeof(psio_tocentry));
        this_entry = this_unit->toc;
        this_entry->last = nullptr;
        for (i = 1; i < this_unit->toclen; i++) {
            last_entry = this_entry;
            this_entry = (psio_tocentry *)malloc(sizeof(psio_tocentry));
            last_entry->next = this_entry;
            this_entry->last = last_entry;
        }
        this_entry->next = nullptr;
    } else {
        this_unit->toclen = 0;
        this_unit->toc = nullptr;
    }

    /* Read the TOC entry-by-entry */
    this_entry = this_unit->toc;
    address = psio_get_address(PSIO_ZERO, sizeof(size_t)); /* start one size_t after the top of the file */
    for (i = 0; i < this_unit->toclen; i++) {
        rw(unit, (char *)this_entry, address, entry_size, 0);
        address = this_entry->eadd;
        this_entry = this_entry->next;
    }
}

}  // namespace psi
