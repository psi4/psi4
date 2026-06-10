/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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

// Standalone round-trip test for libpsio. Writes several TOC entries to a new
// unit (including a multi-page entry that exercises the paging path in
// PSIO::rw), closes and reopens the file, reads everything back, and checks the
// table-of-contents lookups. This is the behavior gate for any internal
// refactor of the file/TOC/paging machinery.

#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"

namespace {

using namespace psi;

bool g_failed = false;

void check(bool cond, const std::string &what) {
    if (!cond) {
        std::cerr << "  FAIL: " << what << "\n";
        g_failed = true;
    }
}

// A TOC entry's worth of test data of a given byte size, deterministically
// filled so we can verify it after a round trip.
std::vector<char> make_payload(std::size_t nbytes, unsigned seed) {
    std::vector<char> v(nbytes);
    for (std::size_t i = 0; i < nbytes; ++i) v[i] = static_cast<char>((i * 31 + seed) & 0xff);
    return v;
}

}  // namespace

int main() {
    psio_init();
    auto psio = PSIO::shared_object();

    const size_t unit = 84;

    // Sizes chosen to span the single-partial-page, multi-full-page, and
    // trailing-partial-page branches of PSIO::rw (PSIO_PAGELEN = 65536).
    struct Entry {
        std::string key;
        std::size_t size;
        unsigned seed;
    };
    std::vector<Entry> entries = {
        {"Alpha small", 64, 1},
        {"Beta one page", PSIO_PAGELEN, 2},
        {"Gamma multi-page", 3 * PSIO_PAGELEN + 137, 3},
        {"Delta empty-ish", 8, 4},
    };

    std::vector<std::vector<char>> payloads;
    for (auto &e : entries) payloads.push_back(make_payload(e.size, e.seed));

    // --- Write pass ---
    psio->open(unit, PSIO_OPEN_NEW);
    for (std::size_t i = 0; i < entries.size(); ++i) {
        psio->write_entry(unit, entries[i].key.c_str(), payloads[i].data(), entries[i].size);
    }
    psio->close(unit, /*keep*/ 1);
    check(!psio->open_check(unit), "unit closed after write pass");

    // --- Read pass (reopen, verify each entry + TOC lookups) ---
    psio->open(unit, PSIO_OPEN_OLD);
    check(psio->open_check(unit) != 0, "unit open after reopen");

    for (std::size_t i = 0; i < entries.size(); ++i) {
        check(psio->tocentry_exists(unit, entries[i].key.c_str()), "tocentry_exists: " + entries[i].key);
        std::vector<char> got(entries[i].size, 0);
        psio->read_entry(unit, entries[i].key.c_str(), got.data(), entries[i].size);
        check(got == payloads[i], "payload round-trips: " + entries[i].key);
    }

    check(!psio->tocentry_exists(unit, "No Such Key"), "absent key reported absent");

    psio->close(unit, /*keep*/ 0);  // erase the scratch file

    if (g_failed) {
        std::cerr << "psio_roundtrip: FAILED\n";
        return 1;
    }
    std::cout << "psio_roundtrip: all checks passed (" << entries.size() << " entries)\n";
    return 0;
}
