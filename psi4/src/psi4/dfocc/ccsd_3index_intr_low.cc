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

void DFOCC::ccsd_3index_intr_low()
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
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

    // t(Q) = 2\sum_{mf} t_m^f b_mf^Q
    T1c->gemv(false, K, t1A, 2.0, 0.0);

    // t(Q,ij) = \sum_{e} t_i^e b_je^Q
    T = SharedTensor2d(new Tensor2d("T1 (Q|IJ)", nQ, naoccA, naoccA));
    T->contract233(false, true, naoccA, naoccA, t1A, K, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);

    // Ttilde(Q,ia) = \sum_{m} t_m^a t_im^Q
    U = SharedTensor2d(new Tensor2d("T1tilde (Q|IA)", nQ, naoccA, navirA));
    U->contract(false, false, nQ * naoccA, navirA, naoccA, T, t1A, 1.0, 0.0);
    T.reset();
    U->write(psio_, PSIF_DFOCC_AMPS);
    U.reset();

    // Tau(Q,ia) = \sum_{mf} (2*Tau_im^af - Tau_mi^af) b_mf^Q
    T = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    T->read_symm(psio_, PSIF_DFOCC_AMPS);
    Tau = SharedTensor2d(new Tensor2d("Tau2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tau->dirprd224(t1A, t1A, 0.5, 0.0);
    Tau->add(T);
    T.reset();
    U = SharedTensor2d(new Tensor2d("2*Tau2(IA,JB) - Tau2(IB,JA)", naoccA, navirA, naoccA, navirA));
    U->sort(1432, Tau, 1.0, 0.0);
    U->scale(-1.0);
    U->axpy(Tau, 2.0);
    Tau.reset();
    T = SharedTensor2d(new Tensor2d("Tau2 (Q|IA)", nQ, naoccA, navirA));
    T->gemm(false, true, K, U, 1.0, 0.0);
    U.reset();
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

    // t(Q,ab) = \sum_{m} t_m^a b_mb^Q
    T = SharedTensor2d(new Tensor2d("T1 (Q|AB)", nQ, navirA, navirA));
    T->contract233(true, false, navirA, navirA, t1A, K, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    K.reset();
    T.reset();

    // t(Q,ia) = \sum_{f} t_i^f b_fa^Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    T = SharedTensor2d(new Tensor2d("T1 (Q|IA)", nQ, naoccA, navirA));
    T->contract233(false, false, naoccA, navirA, t1A, K, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();
    K.reset();

    // t(Q,ai) = \sum_{m} t_m^a b_mi^Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    T = SharedTensor2d(new Tensor2d("T1 (Q|AI)", nQ, navirA, naoccA));
    T->contract233(true, false, navirA, naoccA, t1A, K, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();
    K.reset();

    // Read intermediates
    T = SharedTensor2d(new Tensor2d("T1 (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    U = SharedTensor2d(new Tensor2d("T1 (Q|AI)", nQ, navirA, naoccA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    K = SharedTensor2d(new Tensor2d("T1tilde (Q|IA)", nQ, naoccA, navirA));
    K->read(psio_, PSIF_DFOCC_AMPS);
    // t'(Q,ia) = t(Q,ia) - t(Q,ai) - tt(Q,ia)
    Tau = SharedTensor2d(new Tensor2d("T1p (Q|IA)", nQ, naoccA, navirA));
    Tau->swap_3index_col(U);
    Tau->add(K);
    K.reset();
    Tau->scale(-1.0);
    Tau->add(T);
    Tau->write(psio_, PSIF_DFOCC_AMPS);
    Tau.reset();
    // Tau'(Q,ia) = t(Q,ia) + Tau(Q,ai)
    K = SharedTensor2d(new Tensor2d("Tau2 (Q|IA)", nQ, naoccA, navirA));
    K->read(psio_, PSIF_DFOCC_AMPS);
    Tau = SharedTensor2d(new Tensor2d("Tau2p (Q|IA)", nQ, naoccA, navirA));
    Tau->copy(T);
    Tau->add(K);
    T.reset();
    Tau->write(psio_, PSIF_DFOCC_AMPS);
    Tau.reset();
    // Tau''(Q,ia) = -t(Q,ai) + Tau(Q,ai)
    Tau = SharedTensor2d(new Tensor2d("Tau2pp (Q|IA)", nQ, naoccA, navirA));
    Tau->swap_3index_col(U);
    U.reset();
    Tau->scale(-1.0);
    Tau->add(K);
    K.reset();
    Tau->write(psio_, PSIF_DFOCC_AMPS);
    Tau.reset();
    //outfile->Printf("\t3indices done.\n");

}// end ccsd_3index_intr_low
}} // End Namespaces
