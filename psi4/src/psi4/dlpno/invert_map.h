/*
 * @BEGIN LICENSE
 * Psi4: an open-source quantum chemistry software package
 * Copyright (c) 2007-2026 The Psi4 Developers.
 * The copyrights for code used from other parties are included in the corresponding files.
 * This file is part of Psi4.
 * Psi4 is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with Psi4; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 * @END LICENSE
 */
#ifndef PSI4_SRC_DLPNO_INVERT_MAP_H_
#define PSI4_SRC_DLPNO_INVERT_MAP_H_

// invert_map lives in its own header, rather than in sparse.h, for two reasons: translation units that include sparse.h
// do not have to parse it (its only caller is dlpno.cc), and it can be unit-tested (test/invert_map_test.cc) without
// linking any psi4 library. That is also why this header does not include sparse.h, which would pull in
// libmints/matrix.h.
#include "psi4/libpsi4util/exception.h"
#include <iterator>
#include <vector>

namespace psi {
namespace dlpno {

// Same alias as in sparse.h. Declaring it again with the same type is legal, and
// a mismatch between the two is a compile error in any TU that includes both.
using SparseMap = std::vector<std::vector<int>>;

/* Args: SparseMap from x to y, number of possible y indices
 * Returns: SparseMap from y to x
 */
[[nodiscard]] inline SparseMap invert_map(const SparseMap& x_to_y, int ny) {
    std::vector<int> count(ny, 0);
    for (const auto& ys : x_to_y)
        for (int y : ys) {
            if (y < 0 || y >= ny)
                throw PSIEXCEPTION("Cannot invert an invalid map! ny = " + std::to_string(ny) +
                                   ", y = " + std::to_string(y));
            ++count[y];
        }

    SparseMap y_to_x(ny);
    for (int y = 0; y < ny; ++y) y_to_x[y].reserve(count[y]);

    for (int x = 0; x < std::ssize(x_to_y); ++x)  // ascending x keeps each list sorted
        for (int y : x_to_y[x]) y_to_x[y].push_back(x);

    return y_to_x;
}
}  // namespace dlpno
}  // namespace psi
#endif  // PSI4_SRC_DLPNO_INVERT_MAP_H_
