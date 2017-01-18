/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/libqt/qt.h"
#include "defines.h"
#include "dfocc.h"

using namespace psi;
using namespace std;


namespace psi{ namespace dfoccwave{

void DFOCC::ccd_3index_intr_low()
{

    // defs
    SharedTensor2d K, T, U, Tau;

    // T(Q,ia) = \sum_{jb} b_jb^Q u_ij^ab = \sum_{jb} b(Q,jb) U(jb,ia)
    U = SharedTensor2d(new Tensor2d("U2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    U->read_symm(psio_, PSIF_DFOCC_AMPS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA));
    K->read(psio_, PSIF_DFOCC_INTS);
    T = SharedTensor2d(new Tensor2d("T2 (Q|IA)", nQ, naoccA, navirA));
    T->gemm(false, false, K, U, 1.0, 0.0);
    U.reset();
    K.reset();
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

    //outfile->Printf("\t3indices done.\n");

}// end ccd_3index_intr_low
}} // End Namespaces
