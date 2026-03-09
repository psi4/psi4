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

#include "psi4/libqt/qt.h"
#include "defines.h"
#include "dfocc.h"

using namespace psi;

namespace psi {
namespace dfoccwave {

void DFOCC::ccd_F_intr_low() {
    // defs
    SharedTensor2d K, T, U, Tau;

    // Read ints
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
    K->read(psio_, PSIF_DFOCC_INTS);
    Tau = std::make_shared<Tensor2d>("T2 (Q|IA)", nQ, naoccA, navirA);
    Tau->read(psio_, PSIF_DFOCC_AMPS);

    // F_mi +=  \sum_{Q,e} T_ie^Q b_me^Q
    FijA->zero();
    FijA->contract332(false, true, navirA, K, Tau, 1.0, 1.0);

    // F_ae -=  \sum_{Q,m} T_ma^Q b_me^Q
    FabA->contract(true, false, navirA, navirA, nQ * naoccA, Tau, K, -1.0, 0.0);
    K.reset();
    Tau.reset();

    // outfile->Printf("\tF int done.\n");
}  // end ccd_F_intr_low

}  // namespace dfoccwave
}  // namespace psi
