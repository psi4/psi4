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

#include <vector>
#include "psi4/libmints/dimension.h"

namespace psi {
namespace cctransort {

std::vector<int> pitzer2qt(std::vector<Dimension> &spaces) {
    int nirreps = spaces[0].n();

    Dimension total(nirreps);
    for (int h = 0; h < nirreps; h++)
        for (int i = 0; i < spaces.size(); i++) total[h] += spaces[i][h];
    int nmo = total.sum();

    std::vector<int> order(nmo);
    order.assign(nmo, 0);

    Dimension offset(nirreps);
    offset[0] = 0;
    for (int h = 1; h < nirreps; h++) offset[h] = offset[h - 1] + total[h - 1];

    int count = 0;

    for (int j = 0; j < spaces.size(); j++)
        for (int h = 0; h < nirreps; h++) {
            int this_offset = offset[h];
            for (int k = 0; k < j; k++) this_offset += spaces[k][h];
            for (int i = 0; i < spaces[j][h]; i++) order[this_offset + i] = count++;
        }

    return order;
}

}  // namespace cctransort
}  // namespace psi
