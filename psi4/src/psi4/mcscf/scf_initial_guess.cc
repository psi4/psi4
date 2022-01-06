/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include "psi4/liboptions/liboptions.h"

#include "scf.h"

namespace psi {
namespace mcscf {

void SCF::initial_guess() {
    using namespace psi;

    SBlockMatrix H_t("H_t", nirreps, sopi, sopi);
    SBlockVector eigenvectors("H_t", nirreps, sopi);

    transform(H, H_t, S_sqrt_inv);

    H_t.diagonalize(C_t, eigenvectors);

    C.multiply(false, false, S_sqrt_inv, C_t);

    epsilon = eigenvectors;

    guess_occupation();
}

}  // namespace mcscf
}  // namespace psi
