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

void DFOCC::ccsdl_3index_intr()
{

    // defs
    SharedTensor2d K, L, T, U, Tau, V, V2, Vij, Vai, X, Y, Z;

    // L(Q,ia) = \sum_{jb} b_jb^Q Ut_ij^ab = \sum_{jb} b(Q,jb) Ut(jb,ia)
    U = SharedTensor2d(new Tensor2d("Ut2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_u2_amps(U,l2);
    T = SharedTensor2d(new Tensor2d("L2 (Q|IA)", nQ, naoccA, navirA));
    T->gemm(false, false, bQiaA, U, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

    // Lt(Q,ia) = \sum_{jb} t_jb^Q Ut_ij^ab = \sum_{jb} t(Q,jb) Ut(jb,ia)
    K = SharedTensor2d(new Tensor2d("T1 (Q|IA)", nQ, naoccA, navirA));
    K->read(psio_, PSIF_DFOCC_AMPS);
    T = SharedTensor2d(new Tensor2d("L2t (Q|IA)", nQ, naoccA, navirA));
    T->gemm(false, false, K, U, 1.0, 0.0);
    U.reset();
    K.reset();
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

    // l(Q,ia) = \sum_{e} l_i^e b_ea^Q
    T = SharedTensor2d(new Tensor2d("L1 (Q|IA)", nQ, naoccA, navirA));
    T->contract233(false, false, naoccA, navirA, l1A, bQabA, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

    // lt(Q,ia) = \sum_{e} l_i^e t_ea^Q
    T = SharedTensor2d(new Tensor2d("L1t (Q|IA)", nQ, naoccA, navirA));
    U = SharedTensor2d(new Tensor2d("T1 (Q|AB)", nQ, navirA, navirA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->contract233(false, false, naoccA, navirA, l1A, U, 1.0, 0.0);
    U.reset();
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

    // Build G_mi and G_ae
    // G_mi = \sum_{n,e,f} U_mn^ef L_in^ef = U(mn,ef) L(in,ef)
    U = SharedTensor2d(new Tensor2d("U2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    U->read_symm(psio_, PSIF_DFOCC_AMPS);
    T = SharedTensor2d(new Tensor2d("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    T->sort(1324, U, 1.0, 0.0);
    U.reset();
    L = SharedTensor2d(new Tensor2d("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    L->sort(1324, l2, 1.0, 0.0);
    GijA->contract(false, true, naoccA, naoccA, naoccA*navirA*navirA, T, L, 1.0, 0.0);

    // G_ae = -\sum_{m,n,f} U_mn^ef L_mn^af = L(mn,fa) U(mn,fe)
    GabA->contract(true, false, navirA, navirA, naoccA*naoccA*navirA, L, T, -1.0, 0.0);
    T.reset();
    L.reset();

    // G_Q = 2\sum_{ef} G_ef b_ef^Q
    gQ->gemv(false, bQabA, GabA, 2.0, 0.0);

    // G(Q,ia) = \sum_{m} G_mi b_ma^Q
    T = SharedTensor2d(new Tensor2d("G (Q|IA)", nQ, naoccA, navirA));
    T->contract233(true, false, naoccA, navirA, GijA, bQiaA, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

    // G(Q,ai) = \sum_{e} G_ae b_ie^Q
    T = SharedTensor2d(new Tensor2d("G (Q|AI)", nQ, navirA, naoccA));
    T->contract233(false, true, navirA, naoccA, GabA, bQiaA, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);

    // G_Q' = 2\sum_{me} t_m^e G_em^Q
    U = SharedTensor2d(new Tensor2d("T1 <A|I>", navirA, naoccA));
    U = t1A->transpose();
    gQp->gemv(false, T, U, 2.0, 0.0);
    T.reset();
    U.reset();

    // Build V_ijkl
    // V_ijkl = \sum_{ef} Tau_ij^ef L_kl^ef
    t2 = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    t2->read_symm(psio_, PSIF_DFOCC_AMPS);
    Tau = SharedTensor2d(new Tensor2d("Tau (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_tau_amps(Tau,t2);
    //t2.reset();
    U = SharedTensor2d(new Tensor2d("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA));
    U->sort(1324, Tau, 1.0, 0.0);
    Tau.reset();
    L = SharedTensor2d(new Tensor2d("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    L->sort(1324, l2, 1.0, 0.0);
    V = SharedTensor2d(new Tensor2d("V <IJ|KL>", naoccA, naoccA, naoccA, naoccA));
    V->gemm(false, true, U, L, 1.0, 0.0);
    L.reset();
    U.reset();
    V->write(psio_, PSIF_DFOCC_AMPS);

    // V_ij^Q = \sum_{mn} (2*V_imjn - V_imnj) b_mn^Q = B(Q,mn) Y(mn,ij)
    X = SharedTensor2d(new Tensor2d("X <IJ|KL>", naoccA, naoccA, naoccA, naoccA));
    // X_imjn = 2*V_imjn - V_imnj
    X->tei_cs1_anti_symm(V,V);
    V.reset();
    // Y_mnij = X_imjn
    Y = SharedTensor2d(new Tensor2d("Y <IJ|KL>", naoccA, naoccA, naoccA, naoccA));
    Y->sort(2413, X, 1.0, 0.0);
    X.reset();
    T = SharedTensor2d(new Tensor2d("V (Q|IJ)", nQ, naoccA, naoccA));
    T->gemm(false, false, bQijA, Y, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();
    // Vt_ij^Q = \sum_{mn} (2*V_imjn - V_imnj) t_nm^Q
    T = SharedTensor2d(new Tensor2d("Vt (Q|IJ)", nQ, naoccA, naoccA));
    U = SharedTensor2d(new Tensor2d("T1 (Q|IJ)", nQ, naoccA, naoccA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    L = SharedTensor2d(new Tensor2d("T1 (Q|JI)", nQ, naoccA, naoccA));
    L->swap_3index_col(U);
    U.reset();
    T->gemm(false, false, L, Y, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();
    L.reset();
    Y.reset();

    // Build V_iajb
    // V_iajb = 1/2 \sum_{me} T'(ib,me) L'(me,ja)
    T = SharedTensor2d(new Tensor2d("T2p (IB|JA)", naoccA, navirA, naoccA, navirA));
    ccsd_t2_prime_amps(T,t2);
    t2.reset();
    L = SharedTensor2d(new Tensor2d("L2p (IB|JA)", naoccA, navirA, naoccA, navirA));
    ccsd_t2_prime_amps(L,l2);
    X = SharedTensor2d(new Tensor2d("X (IB|JA)", naoccA, navirA, naoccA, navirA));
    X->gemm(false, false, T, L, 0.5, 0.0);
    T.reset();
    L.reset();
    V = SharedTensor2d(new Tensor2d("V (IA|JB)", naoccA, navirA, naoccA, navirA));
    V->sort(1432, X, 1.0, 0.0);
    X.reset();
    //V->print();

    // V_ij^Q' = 2\sum_{ef} V_iejf b_ef^Q
    Vij = SharedTensor2d(new Tensor2d("Vp (Q|IJ)", nQ, naoccA, naoccA));
    X = SharedTensor2d(new Tensor2d("X (AB|IJ)", navirA, navirA, naoccA, naoccA));
    X->sort(2413, V, 1.0, 0.0);
    Vij->gemm(false, false, bQabA, X, 2.0, 0.0);
    X.reset();

    // V_ai^Q = \sum_{me} V_maie b_me^Q
    Vai = SharedTensor2d(new Tensor2d("V (Q|AI)", nQ, navirA, naoccA));
    X = SharedTensor2d(new Tensor2d("X (IB|AJ)", naoccA, navirA, navirA, naoccA));
    X->sort(1423, V, 1.0, 0.0);
    V.reset();
    Vai->gemm(false, false, bQiaA, X, 1.0, 0.0);
    X.reset();

    // Build V_iabj
    // V_iabj = -1/2 \sum_{me} U(ib,me) L(me,ja)
    U = SharedTensor2d(new Tensor2d("U2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    U->read_symm(psio_, PSIF_DFOCC_AMPS);
    X = SharedTensor2d(new Tensor2d("X (IB|JA)", naoccA, navirA, naoccA, navirA));
    X->gemm(false, false, U, l2, -0.5, 0.0);
    U.reset();
    // V_iabj += 1/2 \sum_{me} T(ib,me) L'(me,ja)
    t2 = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    t2->read_symm(psio_, PSIF_DFOCC_AMPS);
    L = SharedTensor2d(new Tensor2d("L2p (IB|JA)", naoccA, navirA, naoccA, navirA));
    ccsd_t2_prime_amps(L,l2);
    X->gemm(false, false, t2, L, 0.5, 1.0);
    t2.reset();
    L.reset();
    V = SharedTensor2d(new Tensor2d("V (IA|BJ)", naoccA, navirA, navirA, naoccA));
    V->sort(1423, X, 1.0, 0.0);
    X.reset();

    // V_ij^Q' -= \sum_{ef} V_iefj b_ef^Q
    X = SharedTensor2d(new Tensor2d("X (AB|IJ)", navirA, navirA, naoccA, naoccA));
    X->sort(2314, V, 1.0, 0.0);
    Vij->gemm(false, false, bQabA, X, -1.0, 1.0);
    X.reset();
    Vij->write(psio_, PSIF_DFOCC_AMPS);
    Vij.reset();

    // V_ai^Q -= 2\sum_{me} V_maei b_me^Q
    X = SharedTensor2d(new Tensor2d("X (IB|AJ)", naoccA, navirA, navirA, naoccA));
    X->sort(1324, V, 1.0, 0.0);
    V.reset();
    Vai->gemm(false, false, bQiaA, X, -2.0, 1.0);
    X.reset();
    Vai->write(psio_, PSIF_DFOCC_AMPS);
    //Vai->print();
    Vai.reset();

    // Build L_ijka
    // L_ijka = \sum{e} L(ja,ie) T(k,e)
    X = SharedTensor2d(new Tensor2d("X <JA|IK>", naoccA, navirA, naoccA, naoccA));
    X->contract(false, true, naoccA*navirA*naoccA, naoccA, navirA, l2, t1A, 1.0, 0.0);
    L = SharedTensor2d(new Tensor2d("L <IJ|KA>", naoccA, naoccA, naoccA, navirA));
    L->sort(3142, X, 1.0, 0.0);
    X.reset();
    L->write(psio_, PSIF_DFOCC_AMPS);

    // Z_ij^Q = \sum_{me} (2*L_imje - L_mije) t_me^Q
    X = SharedTensor2d(new Tensor2d("X <IM|JE>", naoccA, naoccA, naoccA, navirA));
    X->tei_cs2_anti_symm(L,L);
    L.reset();
    // Y_meij = X_imje
    Y = SharedTensor2d(new Tensor2d("Y <ME|IJ>", naoccA, navirA, naoccA, naoccA));
    Y->sort(2413, X, 1.0, 0.0);
    T = SharedTensor2d(new Tensor2d("T1 (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    Z = SharedTensor2d(new Tensor2d("Zeta (Q|IJ)", nQ, naoccA, naoccA));
    Z->gemm(false, false, T, Y, 1.0, 0.0);
    Y.reset();
    T.reset();
    Z->write(psio_, PSIF_DFOCC_AMPS);
    Z.reset();

    // Z_ia^Q = -\sum_{mn} (2*L_mina - L_imna) t_mn^Q
    // Y_mnia = X_mina
    Y = SharedTensor2d(new Tensor2d("Y <MN|IA>", naoccA, naoccA, naoccA, navirA));
    Y->sort(1324, X, 1.0, 0.0);
    X.reset();
    T = SharedTensor2d(new Tensor2d("T1 (Q|IJ)", nQ, naoccA, naoccA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    Z = SharedTensor2d(new Tensor2d("Zeta (Q|IA)", nQ, naoccA, navirA));
    Z->gemm(false, false, T, Y, -1.0, 0.0);
    T.reset();
    Y.reset();
    Z->write(psio_, PSIF_DFOCC_AMPS);
    Z.reset();

    //outfile->Printf("\t3indices done.\n");

}// end ccsdl_3index_intr
}} // End Namespaces
