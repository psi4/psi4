/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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
// one bucket, and mid-bucket termination.
//
// Every case is run through all four write/read API combinations -- {C, C++}
// write x {C, C++} read -- since libiwl exposes both a C-style (iwl_buf_*) and
// a C++ (psi::IWL) interface and both are in use in-tree; a byte written by one
// must read back identically through either.
//
// Also covers cutoff filtering. write_value() keeps an integral only when
// |value| > cutoff, so sub-cutoff entries are *dropped*, not stored as zeros --
// nothing is read back for them at all. That makes filtering a repacking
// operation: survivors slide forward into earlier bucket slots, so which
// bucket the dropped run falls in changes the packing downstream of it. The
// positional cases below (first/second/second-to-last/last bucket fully
// filtered, plus all-filtered and interleaved) pin that behavior down.

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "psi4/libiwl/iwl.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"

// Hosting globals (`outfile`, `restart_id`, `psi_file_prefix`, `global_dpd_`)
// live in test/stubs.cc; the link line pulls that TU in.

namespace {

using namespace psi;

// The IWL cutoff used throughout: write_value() keeps an integral only when
// |value| is strictly greater than this.
constexpr double CUTOFF = 1.0e-14;

struct Quartet {
    int p, q, r, s;
    double value;
};

// Deterministic pseudo-random test data with bounded indices.
std::vector<Quartet> make_test_data(const std::size_t n) {
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

// Sub-cutoff stand-ins. Exactly +0.0 is the common case in a real integrals
// sweep, but a nonzero magnitude below the cutoff must be dropped just the
// same, and |value| == cutoff must *also* drop because write_value() tests
// strictly greater-than. Cycling all three keeps that boundary honest.
double subcutoff_value(const std::size_t i) {
    switch (i % 3) {
        case 0:
            return 0.0;
        case 1:
            return (i % 6 == 1) ? 1.0e-20 : -1.0e-20;
        default:
            return 1.0e-14;  // == cutoff, so filtered by the strict >
    }
}

// A writer emits `data` to tape `itap` and leaves the file on disk for the
// reader; a reader pulls it back. The two write and two read functions below
// are the C++ (psi::IWL) and C-style (iwl_buf_*) sides of the same operation.
using Writer = void (*)(const int itap, const std::vector<Quartet>& data);
using Reader = std::vector<Quartet> (*)(const int itap);

// Write via the C++ API (psi::IWL).
void write_all_cpp(const int itap, const std::vector<Quartet>& data) {
    IWL buf(_default_psio_lib_.get(), itap, CUTOFF, /*oldfile*/ 0, /*readflag*/ 0);
    for (const auto& q : data) {
        buf.write_value(q.p, q.q, q.r, q.s, q.value, /*printflag*/ 0, std::string("outfile"), /*dirac*/ 0);
    }
    buf.flush(/*lastbuf*/ 1);
    buf.set_keep_flag(true);  // keep the file on disk for the reader
}

// Write via the C-style API (iwl_buf_*).
void write_all_c(const int itap, const std::vector<Quartet>& data) {
    iwlbuf B;
    iwl_buf_init(&B, itap, CUTOFF, /*oldfile*/ 0, /*readflag*/ 0);
    for (const auto& q : data) {
        iwl_buf_wrt_val(&B, q.p, q.q, q.r, q.s, q.value, /*printflag*/ 0, std::string("outfile"), /*dirac*/ 0);
    }
    iwl_buf_flush(&B, /*lastbuf*/ 1);
    iwl_buf_close(&B, /*keep*/ 1);  // keep the file on disk for the reader
}

// Read via the C-style API (the dominant in-tree pattern).
std::vector<Quartet> read_all_c(const int itap) {
    std::vector<Quartet> out;
    iwlbuf B;
    iwl_buf_init(&B, itap, CUTOFF, /*oldfile*/ 1, /*readflag*/ 1);

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
    iwl_buf_close(&B, /*keep*/ 0);  // done reading; delete the scratch file
    return out;
}

// Read via the C++ API (psi::IWL). fetch() refills the buffer each pass; the
// label/value pointers are stable across fetches.
std::vector<Quartet> read_all_cpp(const int itap) {
    std::vector<Quartet> out;
    IWL buf(_default_psio_lib_.get(), itap, CUTOFF, /*oldfile*/ 1, /*readflag*/ 0);
    const Label* labels = buf.labels();
    const Value* values = buf.values();

    int lastbuf;
    do {
        buf.fetch();
        lastbuf = buf.last_buffer();
        for (int i = 0; i < buf.buffer_count(); ++i) {
            const int j = 4 * i;
            const int p = static_cast<int>(labels[j + 0]);
            const int q = static_cast<int>(labels[j + 1]);
            const int r = static_cast<int>(labels[j + 2]);
            const int s = static_cast<int>(labels[j + 3]);
            const double v = static_cast<double>(values[i]);
            out.push_back({p, q, r, s, v});
        }
    } while (!lastbuf);
    buf.set_keep_flag(false);  // done reading; delete the scratch file
    return out;
}

// Exact comparison, deliberately. IWL neither scales nor rounds: write_value()
// casts to Value (= double, so a no-op) and put()/fetch() move the raw bytes
// through libpsio. A surviving integral must therefore come back bit-identical,
// and any tolerance here would only hide a real corruption.
bool quartets_equal(const Quartet& a, const Quartet& b) {
    return a.p == b.p && a.q == b.q && a.r == b.r && a.s == b.s && a.value == b.value;
}

// Run one write/read round-trip for the given size and tape unit through the
// supplied write/read APIs. Returns true on success, false otherwise (with
// diagnostic to stderr).
bool round_trip(const std::size_t n, const int itap, const Writer write, const Reader read, const std::string& api) {
    const auto data = make_test_data(n);
    write(itap, data);
    const auto got = read(itap);

    if (got.size() != data.size()) {
        std::cerr << "  " << api << ": size mismatch: wrote " << data.size() << " read " << got.size() << "\n";
        return false;
    }
    for (std::size_t i = 0; i < data.size(); ++i) {
        if (!quartets_equal(data[i], got[i])) {
            std::cerr << "  " << api << ": entry " << i << " mismatch\n";
            return false;
        }
    }
    return true;
}

// Empty buffer: no integrals written, but flush(lastbuf=1) must still produce a
// readable file that the reader returns as zero entries. The 1995 README
// specifically calls out this case.
bool round_trip_empty(const int itap, const Writer write, const Reader read, const std::string& api) {
    write(itap, {});
    const auto got = read(itap);
    if (!got.empty()) {
        std::cerr << "  " << api << ": empty round-trip returned " << got.size() << " entries\n";
        return false;
    }
    return true;
}

// Round-trip N integrals of which a contiguous run [lo, hi) is pushed below the
// cutoff. Only the survivors should come back, in order and bit-exact; the
// dropped run must leave no trace (no zero-valued placeholder, no gap).
bool round_trip_filtered(const std::size_t n, const std::size_t lo, const std::size_t hi, const int itap,
                         const std::string& what, const Writer write, const Reader read, const std::string& api) {
    auto data = make_test_data(n);
    std::vector<Quartet> expected;
    for (std::size_t i = 0; i < n; ++i) {
        if (i >= lo && i < hi) {
            data[i].value = subcutoff_value(i);
        } else {
            expected.push_back(data[i]);
        }
    }

    write(itap, data);
    const auto got = read(itap);

    if (got.size() != expected.size()) {
        std::cerr << "  " << api << " " << what << ": expected " << expected.size() << " survivors, read " << got.size()
                  << "\n";
        return false;
    }
    for (std::size_t i = 0; i < expected.size(); ++i) {
        if (!quartets_equal(expected[i], got[i])) {
            std::cerr << "  " << api << " " << what << ": survivor " << i << " mismatch\n";
            return false;
        }
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
    int failures = 0;

    constexpr std::size_t B = IWL_INTS_PER_BUF;
    const std::vector<std::size_t> sizes = {1, B - 1, B, B + 1, 3 * B + B / 2};

    const std::size_t N4 = 4 * B;
    struct FilterCase {
        const char* what;
        std::size_t lo, hi;
    };
    // Cutoff filtering, positioned bucket by bucket over a four-bucket write.
    // Each case drops a full bucket's worth of input, so the survivors have to
    // repack across the remaining buckets -- the first and last cases also
    // exercise "the very first thing written is dropped" and "the file ends on
    // a dropped run", where an off-by-one in flush/lastbuf would show up.
    const std::vector<FilterCase> filtered = {
        {"first bucket below cutoff", 0, B},
        {"second bucket below cutoff", B, 2 * B},
        {"second-to-last bucket below cutoff", 2 * B, 3 * B},
        {"last bucket below cutoff", 3 * B, N4},
        {"all below cutoff", 0, N4},
    };

    // Every case runs through all four write/read API combinations. Each reader
    // deletes its file when done, so the unit range is free to reuse per pass.
    struct Api {
        const char* name;
        Writer write;
        Reader read;
    };
    const std::vector<Api> apis = {
        {"C++w/Cr", write_all_cpp, read_all_c},
        {"Cw/C++r", write_all_c, read_all_cpp},
        {"C++w/C++r", write_all_cpp, read_all_cpp},
        {"Cw/Cr", write_all_c, read_all_c},
    };

    for (const auto& api : apis) {
        int unit = base_unit;
        std::cout << "[iwl_roundtrip] ===== API: " << api.name << " =====\n";

        std::cout << "[iwl_roundtrip] empty buffer\n";
        if (!round_trip_empty(unit++, api.write, api.read, api.name)) ++failures;

        for (const auto n : sizes) {
            std::cout << "[iwl_roundtrip] N = " << n << "\n";
            if (!round_trip(n, unit++, api.write, api.read, api.name)) ++failures;
        }

        for (const auto& fc : filtered) {
            std::cout << "[iwl_roundtrip] " << fc.what << "\n";
            if (!round_trip_filtered(N4, fc.lo, fc.hi, unit++, fc.what, api.write, api.read, api.name)) ++failures;
        }

        // A sub-cutoff run that straddles a bucket boundary, so survivors on
        // both sides of it repack into the same bucket.
        std::cout << "[iwl_roundtrip] run straddling a bucket boundary\n";
        if (!round_trip_filtered(3 * B, B - 7, B + 11, unit++, "straddling run", api.write, api.read, api.name))
            ++failures;
    }

    if (failures) {
        std::cerr << failures << " case(s) failed.\n";
        return 1;
    }
    std::cout << "All round-trip cases passed.\n";
    return 0;
}
