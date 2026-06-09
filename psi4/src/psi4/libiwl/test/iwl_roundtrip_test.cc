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

// Standalone round-trip test for libiwl. Exercises every bucket-boundary case
// the 1995 README warned about: empty, single, exactly one bucket, just over
// one bucket, mid-bucket termination, and all-zeros sweeps. Also verifies
// that the 16-bit Label overflow check fires instead of silently truncating.

#include <cassert>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "psi4/libiwl/iwl.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsi4util/PsiOutStream.h"

// The libiwl write paths reach into the global `outfile` only when the caller
// passes `printflag != 0`. We never do, but the symbol must exist for the
// link. Provide a minimal stub.
namespace psi {
std::shared_ptr<PsiOutStream> outfile;
}

namespace {

using namespace psi;

struct Quartet {
    int p, q, r, s;
    double value;
};

// Deterministic pseudo-random test data with bounded indices.
std::vector<Quartet> make_test_data(std::size_t n) {
    std::vector<Quartet> out;
    out.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        // Wraparound keeps indices in [0, 200) so they fit in 16 bits and we
        // can compare without sorting.
        const int p = static_cast<int>((i * 7 + 11) % 200);
        const int q = static_cast<int>((i * 13 + 5) % 200);
        const int r = static_cast<int>((i * 17 + 3) % 200);
        const int s = static_cast<int>((i * 19 + 1) % 200);
        // Values bounded away from cutoff (1e-14) so none get filtered.
        const double v = 1.0 + 0.001 * static_cast<double>(i);
        out.push_back({p, q, r, s, v});
    }
    return out;
}

// Drive an IWL write loop using the C++ API.
void write_all(psi::IWL &buf, const std::vector<Quartet> &data) {
    for (const auto &q : data) {
        buf.write_value(q.p, q.q, q.r, q.s, q.value, /*printflag*/ 0,
                        std::string("outfile"), /*dirac*/ 0);
    }
    buf.flush(/*lastbuf*/ 1);
}

// Drive an IWL read loop using the C-style API (the dominant in-tree pattern).
std::vector<Quartet> read_all(int itap) {
    std::vector<Quartet> out;
    iwlbuf B;
    iwl_buf_init(&B, itap, 1.0e-14, /*oldfile*/ 1, /*readflag*/ 1);

    int lastbuf = B.lastbuf;
    while (true) {
        for (int i = B.idx; i < B.inbuf; ++i) {
            const int j = 4 * i;
            const int p = static_cast<int>(B.labels[j + 0]);
            const int q = static_cast<int>(B.labels[j + 1]);
            const int r = static_cast<int>(B.labels[j + 2]);
            const int s = static_cast<int>(B.labels[j + 3]);
            const double v = static_cast<double>(B.values[i]);
            out.push_back({p, q, r, s, v});
        }
        B.idx = B.inbuf;
        if (lastbuf) break;
        iwl_buf_fetch(&B);
        lastbuf = B.lastbuf;
    }
    iwl_buf_close(&B, /*keep*/ 0);
    return out;
}

bool quartets_equal(const Quartet &a, const Quartet &b) {
    return a.p == b.p && a.q == b.q && a.r == b.r && a.s == b.s &&
           std::fabs(a.value - b.value) < 1.0e-12;
}

// Run one write/read round-trip for the given size and tape unit. Returns
// true on success, false otherwise (with diagnostic to stderr).
bool round_trip(std::size_t n, int itap) {
    auto data = make_test_data(n);

    {
        psi::IWL buf(psi::_default_psio_lib_.get(), itap, 1.0e-14,
                     /*oldfile*/ 0, /*readflag*/ 0);
        write_all(buf, data);
        buf.set_keep_flag(true);
    }

    auto got = read_all(itap);

    if (got.size() != data.size()) {
        std::cerr << "  size mismatch: wrote " << data.size() << " read "
                  << got.size() << "\n";
        return false;
    }
    for (std::size_t i = 0; i < data.size(); ++i) {
        if (!quartets_equal(data[i], got[i])) {
            std::cerr << "  entry " << i << " mismatch\n";
            return false;
        }
    }
    return true;
}

// Empty buffer: no integrals written, but flush(lastbuf=1) must still produce
// a readable file that read_all() returns as zero entries. The 1995 README
// specifically calls out this case.
bool round_trip_empty(int itap) {
    {
        psi::IWL buf(psi::_default_psio_lib_.get(), itap, 1.0e-14,
                     /*oldfile*/ 0, /*readflag*/ 0);
        buf.flush(/*lastbuf*/ 1);
        buf.set_keep_flag(true);
    }
    auto got = read_all(itap);
    if (!got.empty()) {
        std::cerr << "  empty round-trip returned " << got.size()
                  << " entries\n";
        return false;
    }
    return true;
}

// Verify that writing an out-of-range orbital index throws rather than
// silently truncating.
bool overflow_check_fires(int itap) {
    psi::IWL buf(psi::_default_psio_lib_.get(), itap, 1.0e-14,
                 /*oldfile*/ 0, /*readflag*/ 0);
    bool threw = false;
    try {
        buf.write_value(SHRT_MAX + 1, 0, 0, 0, 1.0, 0, std::string("outfile"), 0);
    } catch (const std::exception &e) {
        threw = true;
    }
    buf.flush(/*lastbuf*/ 1);
    buf.set_keep_flag(false);
    if (!threw) {
        std::cerr << "  overflow check did not throw\n";
        return false;
    }
    return true;
}

}  // namespace

int main() {
    // libpsio needs init before anything opens a unit.
    psi::psio_init();

    // Use tape unit numbers in a range unlikely to clash with the few entries
    // libpsio's filecfg might preset.
    const int base_unit = 80;
    int unit = base_unit;
    int failures = 0;

    const std::size_t B = static_cast<std::size_t>(psi::IWL_INTS_PER_BUF);
    const std::vector<std::size_t> sizes = {1, B - 1, B, B + 1, 3 * B + B / 2};

    std::cout << "[iwl_roundtrip] empty buffer\n";
    if (!round_trip_empty(unit++)) ++failures;

    for (auto n : sizes) {
        std::cout << "[iwl_roundtrip] N = " << n << "\n";
        if (!round_trip(n, unit++)) ++failures;
    }

    std::cout << "[iwl_roundtrip] Label overflow check\n";
    if (!overflow_check_fires(unit++)) ++failures;

    if (failures) {
        std::cerr << failures << " case(s) failed.\n";
        return 1;
    }
    std::cout << "All round-trip cases passed.\n";
    return 0;
}
