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

#include <cmath>
#include "scf.h"

namespace psi {
namespace mcscf {

void SCF::construct_S_inverse_sqrt() {
    SBlockVector lambda("lambda", nirreps, sopi);
    SBlockMatrix L("L", nirreps, sopi, sopi);
    SBlockMatrix Lambda("Lambda", nirreps, sopi, sopi);

    S.diagonalize(L, lambda);

    //   lambda->print();
    //   L->print();

    for (int h = 0; h < nirreps; ++h) {
        for (int i = 0; i < sopi[h]; ++i) {
            Lambda->set(h, i, i, 1.0 / sqrt(lambda->get(h, i)));
        }
    }

    T.multiply(false, true, Lambda, L);
    S_sqrt_inv.multiply(false, false, L, T);

    for (int h = 0; h < nirreps; ++h) {
        for (int i = 0; i < sopi[h]; ++i) {
            Lambda->set(h, i, i, sqrt(lambda->get(h, i)));
        }
    }

    T.multiply(false, true, Lambda, L);
    S_sqrt.multiply(false, false, L, T);
}

}  // namespace mcscf
}  // namespace psi
