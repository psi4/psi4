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

#include <algorithm>

#include "libpsi4util.h"

namespace psi {

/**
 * Generate combinations of 0,1,...,(n-1) taken k at a time
 * @param n
 * @param k
 * @param combinations a vector<vector<int> > that will store all the combinations
 */
void generate_combinations(int n, int k, std::vector<std::vector<int>>& combinations) {
    if ((n > 0) && (k > 0)) {
        std::vector<int> combination;
        auto* a = new bool[n];
        for (int i = 0; i < n - k; ++i) a[i] = false;
        for (int i = n - k; i < n; ++i) a[i] = true;
        do {
            combination.clear();
            for (int i = 0; i < n; ++i) {
                if (a[i]) combination.push_back(i);
            }
            combinations.push_back(combination);
        } while (std::next_permutation(a, a + n));
        delete[] a;
    }
}
}  // namespace psi
