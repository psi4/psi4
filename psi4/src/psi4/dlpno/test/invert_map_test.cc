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


// Standalone unit test for psi::dlpno::invert_map. Checks hand-worked
// inversions (including empty rows, y values with no preimage, and the
// ascending-x ordering of the output lists), a round trip, and that every
// out-of-range y -- above, at, and below the valid range [0, ny) -- is
// rejected with a PsiException instead of indexing out of bounds.
//
// invert_map is header-only, so this executable links no psi4 library; the
// only extra source it compiles is libpsi4util/exception.cc, for PsiException.

#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "psi4/dlpno/invert_map.h"
#include "psi4/libpsi4util/exception.h"

namespace {

using psi::dlpno::invert_map;
using psi::dlpno::SparseMap;

int g_checks = 0;
bool g_failed = false;

std::string str(const SparseMap& m) {
    std::ostringstream os;
    os << "{";
    for (std::size_t i = 0; i < m.size(); ++i) {
        os << (i ? ", " : "") << "{";
        for (std::size_t j = 0; j < m[i].size(); ++j) os << (j ? ", " : "") << m[i][j];
        os << "}";
    }
    os << "}";
    return os.str();
}

void check(bool cond, const std::string& what) {
    ++g_checks;
    if (!cond) {
        std::cerr << "  FAIL: " << what << "\n";
        g_failed = true;
    }
}

// Checks that inverting x_to_y gives expected. A throw is reported as a failure,
// so one bad case does not hide the ones after it.
void check_inverts_to(const SparseMap& x_to_y, int ny, const SparseMap& expected, const std::string& what) {
    ++g_checks;
    try {
        const SparseMap got = invert_map(x_to_y, ny);
        if (got != expected) {
            std::cerr << "  FAIL: " << what << ": invert_map(" << str(x_to_y) << ", " << ny << ") gave " << str(got)
                      << ", expected " << str(expected) << "\n";
            g_failed = true;
        }
    } catch (const std::exception& e) {
        std::cerr << "  FAIL: " << what << ": invert_map(" << str(x_to_y) << ", " << ny << ") threw: " << e.what()
                  << "\n";
        g_failed = true;
    }
}

// Checks that inverting x_to_y is rejected with a PsiException.
void check_rejects(const SparseMap& x_to_y, int ny, const std::string& what) {
    ++g_checks;
    try {
        const SparseMap got = invert_map(x_to_y, ny);
        std::cerr << "  FAIL: " << what << ": invert_map(" << str(x_to_y) << ", " << ny << ") returned " << str(got)
                  << " instead of throwing\n";
        g_failed = true;
    } catch (const psi::PsiException&) {
        // expected
    } catch (const std::exception& e) {
        std::cerr << "  FAIL: " << what << ": threw something other than PsiException: " << e.what() << "\n";
        g_failed = true;
    }
}

void test_hand_worked() {
    // x=2 has no y; y=3 has no x; y=0 and y=2 have two preimages each, which
    // must come out in ascending x.
    check_inverts_to({{0, 2}, {1}, {}, {0, 1, 2}}, 4, {{0, 3}, {1, 3}, {0, 3}, {}}, "general map");

    // The order of y within an x's list does not affect the result.
    check_inverts_to({{2, 0}, {0, 2}}, 3, {{0, 1}, {}, {0, 1}}, "unsorted input lists");

    // A one-to-one map inverts to its inverse permutation.
    check_inverts_to({{2}, {0}, {1}}, 3, {{1}, {2}, {0}}, "permutation");

    // A DLPNO-shaped case, like lmo_to_riatoms_ext_ -> riatom_to_lmos_ext_:
    // 5 LMOs over 4 atoms, where atom 2 is in no LMO's domain.
    check_inverts_to({{0, 1}, {1}, {1, 3}, {3}, {0, 3}}, 4, {{0, 4}, {0, 1, 2}, {}, {2, 3, 4}},
                     "LMO -> atom domains");
}

void test_edge_cases() {
    // Both ends of the valid range are accepted.
    check_inverts_to({{0, 4}}, 5, {{0}, {}, {}, {}, {0}}, "y = 0 and y = ny - 1 accepted");

    check_inverts_to({}, 3, {{}, {}, {}}, "no x values: ny empty lists");
    check_inverts_to({{}, {}}, 2, {{}, {}}, "only empty rows: ny empty lists");
    check_inverts_to({}, 0, {}, "no x values and ny = 0");
    check_inverts_to({{}, {}}, 0, {}, "only empty rows and ny = 0");
}

void test_round_trip() {
    // Inverting twice gives the original map back when every list is sorted and
    // free of duplicates. Build such a map deterministically.
    const int nx = 13, ny = 11;
    SparseMap x_to_y(nx);
    std::size_t nentries = 0;
    for (int x = 0; x < nx; ++x) {
        std::set<int> ys;
        for (int k = 0; k < x % 4; ++k) ys.insert((3 * x + 5 * k) % ny);
        x_to_y[x].assign(ys.begin(), ys.end());
        nentries += ys.size();
    }

    try {
        const SparseMap y_to_x = invert_map(x_to_y, ny);
        check(static_cast<int>(y_to_x.size()) == ny, "round trip: inverse has ny lists");

        std::size_t ninverse = 0;
        for (const auto& xs : y_to_x) ninverse += xs.size();
        check(ninverse == nentries, "round trip: inverse has as many entries as the original");

        check_inverts_to(y_to_x, nx, x_to_y, "round trip: inverting twice gives the original back");
    } catch (const std::exception& e) {
        check(false, std::string("round trip: invert_map threw: ") + e.what());
    }
}

void test_rejects_out_of_range() {
    check_rejects({{0, 3}}, 3, "y == ny");
    check_rejects({{7}}, 3, "y > ny");
    check_rejects({{-1}}, 3, "y == -1 (the sentinel used by the *_dense_ maps)");
    check_rejects({{0}}, 0, "any y with ny == 0");
    // Out-of-range entries are found anywhere in the map, not just in the first row.
    check_rejects({{0, 1}, {}, {2}, {1, 5}}, 3, "out-of-range y in the last row");
}

}  // namespace

int main() {
    test_hand_worked();
    test_edge_cases();
    test_round_trip();
    test_rejects_out_of_range();

    if (g_failed) {
        std::cerr << "dlpno_invert_map: FAILED\n";
        return 1;
    }
    std::cout << "dlpno_invert_map: all checks passed (" << g_checks << " checks)\n";
    return 0;
}
