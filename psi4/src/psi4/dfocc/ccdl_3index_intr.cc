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

void DFOCC::ccdl_3index_intr() {
    // defs
    SharedTensor2d K, L, T, U, Tau, V, V2, Vij, Vai, X, Y, Z;

    // L(Q,ia) = \sum_{jb} b_jb^Q Ut_ij^ab = \sum_{jb} b(Q,jb) Ut(jb,ia)
    U = std::make_shared<Tensor2d>("Ut2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    ccsd_u2_amps(U, l2);
    T = std::make_shared<Tensor2d>("L2 (Q|IA)", nQ, naoccA, navirA);
    T->gemm(false, false, bQiaA, U, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

    // Build G_mi and G_ae
    // G_mi = \sum_{n,e,f} U_mn^ef L_in^ef = U(mn,ef) L(in,ef)
    U = std::make_shared<Tensor2d>("U2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    U->read_symm(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T->sort(1324, U, 1.0, 0.0);
    U.reset();
    L = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    L->sort(1324, l2, 1.0, 0.0);
    // GijA->contract442(1, 1, U, l2, 1.0, 0.0);
    GijA->contract(false, true, naoccA, naoccA, naoccA * navirA * navirA, T, L, 1.0, 0.0);

    // G_ae = -\sum_{m,n,f} U_mn^ef L_mn^af = L(mn,fa) U(mn,fe)
    // GabA->contract442(2, 2, l2, U, -1.0, 0.0);
    GabA->contract(true, false, navirA, navirA, naoccA * naoccA * navirA, L, T, -1.0, 0.0);
    T.reset();
    L.reset();

    // G_Q = 2\sum_{ef} G_ef b_ef^Q
    gQ->gemv(false, bQabA, GabA, 2.0, 0.0);

    // G(Q,ia) = \sum_{m} G_mi b_ma^Q
    T = std::make_shared<Tensor2d>("G (Q|IA)", nQ, naoccA, navirA);
    T->contract233(true, false, naoccA, navirA, GijA, bQiaA, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

    // G(Q,ai) = \sum_{e} G_ae b_ie^Q
    T = std::make_shared<Tensor2d>("G (Q|AI)", nQ, navirA, naoccA);
    T->contract233(false, true, navirA, naoccA, GabA, bQiaA, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

    // Build V_ijkl
    // V_ijkl = \sum_{ef} T_ij^ef L_kl^ef
    t2 = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    t2->read_symm(psio_, PSIF_DFOCC_AMPS);
    // t2.reset();
    U = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    U->sort(1324, t2, 1.0, 0.0);
    L = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    L->sort(1324, l2, 1.0, 0.0);
    V = std::make_shared<Tensor2d>("V <IJ|KL>", naoccA, naoccA, naoccA, naoccA);
    V->gemm(false, true, U, L, 1.0, 0.0);
    L.reset();
    U.reset();
    V->write(psio_, PSIF_DFOCC_AMPS);

    // V_ij^Q = \sum_{mn} (2*V_imjn - V_imnj) b_mn^Q = B(Q,mn) Y(mn,ij)
    X = std::make_shared<Tensor2d>("X <IJ|KL>", naoccA, naoccA, naoccA, naoccA);
    // X_imjn = 2*V_imjn - V_imnj
    X->tei_cs1_anti_symm(V, V);
    V.reset();
    // Y_mnij = X_imjn
    Y = std::make_shared<Tensor2d>("Y <IJ|KL>", naoccA, naoccA, naoccA, naoccA);
    Y->sort(2413, X, 1.0, 0.0);
    X.reset();
    T = std::make_shared<Tensor2d>("V (Q|IJ)", nQ, naoccA, naoccA);
    T->gemm(false, false, bQijA, Y, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

    // Build V_iajb
    // V_iajb = 1/2 \sum_{me} T'(ib,me) L'(me,ja)
    T = std::make_shared<Tensor2d>("T2p (IB|JA)", naoccA, navirA, naoccA, navirA);
    ccsd_t2_prime_amps(T, t2);
    t2.reset();
    L = std::make_shared<Tensor2d>("L2p (IB|JA)", naoccA, navirA, naoccA, navirA);
    ccsd_t2_prime_amps(L, l2);
    X = std::make_shared<Tensor2d>("X (IB|JA)", naoccA, navirA, naoccA, navirA);
    X->gemm(false, false, T, L, 0.5, 0.0);
    T.reset();
    L.reset();
    V = std::make_shared<Tensor2d>("V (IA|JB)", naoccA, navirA, naoccA, navirA);
    V->sort(1432, X, 1.0, 0.0);
    X.reset();

    // V_ij^Q' = 2\sum_{ef} V_iejf b_ef^Q
    Vij = std::make_shared<Tensor2d>("Vp (Q|IJ)", nQ, naoccA, naoccA);
    X = std::make_shared<Tensor2d>("X (AB|IJ)", navirA, navirA, naoccA, naoccA);
    X->sort(2413, V, 1.0, 0.0);
    Vij->gemm(false, false, bQabA, X, 2.0, 0.0);
    X.reset();

    // V_ai^Q = \sum_{me} V_maie b_me^Q
    Vai = std::make_shared<Tensor2d>("V (Q|AI)", nQ, navirA, naoccA);
    X = std::make_shared<Tensor2d>("X (IB|AJ)", naoccA, navirA, navirA, naoccA);
    X->sort(1423, V, 1.0, 0.0);
    V.reset();
    Vai->gemm(false, false, bQiaA, X, 1.0, 0.0);
    X.reset();

    // Build V_iabj
    // V_iabj = -1/2 \sum_{me} U(ib,me) L(me,ja)
    U = std::make_shared<Tensor2d>("U2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    U->read_symm(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (IB|JA)", naoccA, navirA, naoccA, navirA);
    X->gemm(false, false, U, l2, -0.5, 0.0);
    U.reset();
    // V_iabj += 1/2 \sum_{me} T(ib,me) L'(me,ja)
    t2 = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    t2->read_symm(psio_, PSIF_DFOCC_AMPS);
    L = std::make_shared<Tensor2d>("L2p (IB|JA)", naoccA, navirA, naoccA, navirA);
    ccsd_t2_prime_amps(L, l2);
    X->gemm(false, false, t2, L, 0.5, 1.0);
    t2.reset();
    L.reset();
    V = std::make_shared<Tensor2d>("V (IA|BJ)", naoccA, navirA, navirA, naoccA);
    V->sort(1423, X, 1.0, 0.0);
    X.reset();

    // V_ij^Q' -= \sum_{ef} V_iefj b_ef^Q
    X = std::make_shared<Tensor2d>("X (AB|IJ)", navirA, navirA, naoccA, naoccA);
    X->sort(2314, V, 1.0, 0.0);
    Vij->gemm(false, false, bQabA, X, -1.0, 1.0);
    X.reset();
    Vij->write(psio_, PSIF_DFOCC_AMPS);
    Vij.reset();

    // V_ai^Q -= 2\sum_{me} V_maei b_me^Q
    X = std::make_shared<Tensor2d>("X (IB|AJ)", naoccA, navirA, navirA, naoccA);
    X->sort(1324, V, 1.0, 0.0);
    V.reset();
    Vai->gemm(false, false, bQiaA, X, -2.0, 1.0);
    X.reset();
    Vai->write(psio_, PSIF_DFOCC_AMPS);
    Vai.reset();

    // outfile->Printf("\t3indices done.\n");

}  // end ccdl_3index_intr
}  // namespace dfoccwave
}  // namespace psi
