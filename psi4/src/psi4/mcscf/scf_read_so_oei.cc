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

#define MAX(i, j) ((i > j) ? i : j)
#define MIN(i, j) ((i > j) ? j : i)
#define INDEX(i, j) ((i > j) ? (ioff[(i)] + (j)) : (ioff[(j)] + (i)))

#include <iostream>

#include "psi4/psifiles.h"
#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libpsi4util/libpsi4util.h"

#include "scf.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/mintshelper.h"

extern FILE* outfile;

namespace psi {
namespace mcscf {

void SCF::read_so_oei() {

    S_ = SharedMatrix(mintshelper()->so_overlap()->clone());
    H_ = SharedMatrix(mintshelper()->so_kinetic()->clone());
    H_->add(mintshelper()->so_potential());

    // Grab the overlap integrals in Pitzer order
    for (int h = 0; h < nirreps; h++) {
        for (int i = 0; i < S->get_rows(h); i++) {
            for (int j = 0; j < S->get_rows(h); j++) {
                S->set(h, i, j, S_->get(h,i,j));
            }
        }
    }
    // Grab the one-electron integrals in Pitzer order
    for (int h = 0; h < nirreps; h++) {
        for (int i = 0; i < H->get_rows(h); i++) {
            for (int j = 0; j < H->get_cols(h); j++) {
                int ij = INDEX(H->get_abs_row(h, i), H->get_abs_col(h, j));
                H->set(h, i, j, H_->get(h, i, j));
            }
        }
    }
    if (options_.get_int("DEBUG") > 4) {
        S->print();
        H->print();
    }
}

}  // namespace mcscf
}  // namespace psi
