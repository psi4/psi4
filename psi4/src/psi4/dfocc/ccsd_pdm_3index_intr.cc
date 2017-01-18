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

void DFOCC::ccsd_pdm_3index_intr()
{

    // defs
    SharedTensor2d K, L, T, U, Tau, V, V2, Vij, Vai, Vab, X, Y, Z;
    SharedTensor2d Vijka, Vijak;
    SharedTensor2d Vs, Va, Ts, Ta, S, A;

    // Tau(Q,ia) = \sum_{jb} b_jb^Q Ubar_ij^ab = \sum_{jb} b(Q,jb) Ubar(jb,ia)
    U = SharedTensor2d(new Tensor2d("2*Tau(ia,jb) - Tau(ib,ja)", naoccA, navirA, naoccA, navirA));
    U->read_symm(psio_, PSIF_DFOCC_AMPS);
    T = SharedTensor2d(new Tensor2d("Tau (Q|IA)", nQ, naoccA, navirA));
    T->gemm(false, false, bQiaA, U, 1.0, 0.0);
    U.reset();
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

    // l(Q) = 2\sum_{me} l_m^e b_me^Q
    L1c->gemv(false, bQiaA, l1A, 2.0, 0.0);

    // L(Q,ia) = \sum_{jb} b_jb^Q Ut_ij^ab = \sum_{jb} b(Q,jb) Ut(jb,ia)
    U = SharedTensor2d(new Tensor2d("Ut2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_u2_amps(U,l2);
    T = SharedTensor2d(new Tensor2d("L2 (Q|IA)", nQ, naoccA, navirA));
    T->gemm(false, false, bQiaA, U, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

    // l(Q,ia) = \sum_{e} l_i^e b_ea^Q
    T = SharedTensor2d(new Tensor2d("L1 (Q|IA)", nQ, naoccA, navirA));
    T->contract233(false, false, naoccA, navirA, l1A, bQabA, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

    // l(Q,ai) = \sum_{m} l_m^a b_mi^Q
    T = SharedTensor2d(new Tensor2d("L1 (Q|AI)", nQ, navirA, naoccA));
    T->contract233(true, false, navirA, naoccA, l1A, bQijA, 1.0, 0.0);
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

    // Build Gt_mi and Gt_ae
    // Gt_mi = G_mi + \sum_{e} t_m^e l_i^e
    GtijA->copy(GijA);
    GtijA->gemm(false, true, t1A, l1A, 1.0, 1.0);

    // Gt_ae = G_ae - \sum_{m} t_m^e l_m^a
    GtabA->copy(GabA);
    GtabA->gemm(true, false, l1A, t1A, -1.0, 1.0);

    // G_Q = 2\sum_{ef} G_ef b_ef^Q
    gQ->gemv(false, bQabA, GabA, 2.0, 0.0);

    // Gt_Q = 2\sum_{mn} G_mn b_mn^Q
    gQt->gemv(false, bQijA, GijA, 2.0, 0.0);

    // G(Q,ij) = \sum_{m} G_im b_mj^Q
    T = SharedTensor2d(new Tensor2d("G (Q|IJ)", nQ, naoccA, naoccA));
    T->contract233(false, false, naoccA, naoccA, GijA, bQijA, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

    // G(Q,ai) = \sum_{e} G_ae b_ie^Q
    T = SharedTensor2d(new Tensor2d("G (Q|AI)", nQ, navirA, naoccA));
    T->contract233(false, true, navirA, naoccA, GabA, bQiaA, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

    // Gt(Q,ia) = \sum_{m} Gt_mi b_ma^Q
    T = SharedTensor2d(new Tensor2d("Gt (Q|IA)", nQ, naoccA, navirA));
    T->contract233(true, false, naoccA, navirA, GtijA, bQiaA, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

    // Gt(Q,ai) = \sum_{e} Gt_ae b_ie^Q
    T = SharedTensor2d(new Tensor2d("Gt (Q|AI)", nQ, navirA, naoccA));
    T->contract233(false, true, navirA, naoccA, GtabA, bQiaA, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

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

    // Y_ijab = \sum(mn) Tau_mn^ab V_ijmn
    Y = SharedTensor2d(new Tensor2d("Y <IJ|AB>", naoccA, naoccA, navirA, navirA));
    Y->gemm(false, false, V, U, 1.0, 0.0);
    U.reset();
    Y->write(psio_, PSIF_DFOCC_AMPS);
    Y.reset();

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
    L.reset();
    Y.reset();
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

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
    V->write(psio_, PSIF_DFOCC_AMPS);

    // Build V_ijka
    // V_ijka = \sum_{e} t_j^e V_ieka
    // Xeika = Vieka
    X = SharedTensor2d(new Tensor2d("X (EI|KA)", navirA, naoccA, naoccA, navirA));
    X->sort(2134, V, 1.0, 0.0);
    Y = SharedTensor2d(new Tensor2d("V (JI|KA)", naoccA, naoccA, naoccA, navirA));
    Y->contract(false, false, naoccA, naoccA*naoccA*navirA, navirA, t1A, X, 1.0, 0.0);
    X.reset();
    Vijka = SharedTensor2d(new Tensor2d("V (IJ|KA)", naoccA, naoccA, naoccA, navirA));
    Vijka->sort(2134, Y, 1.0, 0.0);
    Y.reset();

    // V_ij^Q' <= 2\sum_{ef} V_iejf b_ef^Q
    Vij = SharedTensor2d(new Tensor2d("Vp (Q|IJ)", nQ, naoccA, naoccA));
    X = SharedTensor2d(new Tensor2d("X (AB|IJ)", navirA, navirA, naoccA, naoccA));
    X->sort(2413, V, 1.0, 0.0);
    Vij->gemm(false, false, bQabA, X, 2.0, 0.0);
    X.reset();

    // V_ai^Q <= \sum_{me} V_maie b_me^Q
    Vai = SharedTensor2d(new Tensor2d("V (Q|AI)", nQ, navirA, naoccA));
    X = SharedTensor2d(new Tensor2d("X (IB|AJ)", naoccA, navirA, navirA, naoccA));
    X->sort(1423, V, 1.0, 0.0);
    Vai->gemm(false, false, bQiaA, X, 1.0, 0.0);
    X.reset();

    // V_ab^Q = 2\sum_{mn} V_manb b_mn^Q
    X = SharedTensor2d(new Tensor2d("X (MN|AB)", naoccA, naoccA, navirA, navirA));
    X->sort(1324, V, 1.0, 0.0);
    V.reset();
    Vab = SharedTensor2d(new Tensor2d("V (Q|AB)", nQ, navirA, navirA));
    Vab->gemm(false, false, bQijA, X, 2.0, 0.0);
    Vab->write(psio_, PSIF_DFOCC_AMPS);
    Vab.reset();

    // Vt_ab^Q = 2\sum_{mn} V_manb t_nm^Q
    Vab = SharedTensor2d(new Tensor2d("Vt (Q|AB)", nQ, navirA, navirA));
    U = SharedTensor2d(new Tensor2d("T1 (Q|IJ)", nQ, naoccA, naoccA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    L = SharedTensor2d(new Tensor2d("T1 (Q|JI)", nQ, naoccA, naoccA));
    L->swap_3index_col(U);
    U.reset();
    Vab->gemm(false, false, L, X, 2.0, 0.0);
    X.reset();
    L.reset();
    Vab->write(psio_, PSIF_DFOCC_AMPS);
    Vab.reset();


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
    V->write(psio_, PSIF_DFOCC_AMPS);
    X.reset();

    // Build V_ijak
    // V_ijak = \sum_{e} t_j^e V_ieak
    // Xeiak = Vieak
    X = SharedTensor2d(new Tensor2d("X (EI|AK)", navirA, naoccA, navirA, naoccA));
    X->sort(2134, V, 1.0, 0.0);
    Y = SharedTensor2d(new Tensor2d("V (JI|AK)", naoccA, naoccA, navirA, naoccA));
    Y->contract(false, false, naoccA, naoccA*naoccA*navirA, navirA, t1A, X, 1.0, 0.0);
    X.reset();
    Vijak = SharedTensor2d(new Tensor2d("V (IJ|AK)", naoccA, naoccA, navirA, naoccA));
    Vijak->sort(2134, Y, 1.0, 0.0);
    Y.reset();

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
    Vai->gemm(false, false, bQiaA, X, -2.0, 1.0);
    X.reset();
    Vai->write(psio_, PSIF_DFOCC_AMPS);
    Vai.reset();

    // V_ab^Q -= \sum_{mn} V_mabn b_mn^Q
    X = SharedTensor2d(new Tensor2d("X (MN|AB)", naoccA, naoccA, navirA, navirA));
    X->sort(1423, V, 1.0, 0.0);
    V.reset();
    Vab = SharedTensor2d(new Tensor2d("V (Q|AB)", nQ, navirA, navirA));
    Vab->read(psio_, PSIF_DFOCC_AMPS);
    Vab->gemm(false, false, bQijA, X, -1.0, 1.0);
    Vab->write(psio_, PSIF_DFOCC_AMPS);
    Vab.reset();

    // Vt_ab^Q -= \sum_{mn} V_mabn t_nm^Q
    Vab = SharedTensor2d(new Tensor2d("Vt (Q|AB)", nQ, navirA, navirA));
    Vab->read(psio_, PSIF_DFOCC_AMPS);
    U = SharedTensor2d(new Tensor2d("T1 (Q|IJ)", nQ, naoccA, naoccA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    L = SharedTensor2d(new Tensor2d("T1 (Q|JI)", nQ, naoccA, naoccA));
    L->swap_3index_col(U);
    U.reset();
    Vab->gemm(false, false, L, X, -1.0, 1.0);
    X.reset();
    L.reset();
    Vab->write(psio_, PSIF_DFOCC_AMPS);
    Vab.reset();

    // cal{V}_ij^Q = \sum_{me} (2V_imje - V_imej) b_me^Q
    X = SharedTensor2d(new Tensor2d("X (IM|JE)", naoccA, naoccA, naoccA, navirA));
    X->tei_cs1_anti_symm(Vijka, Vijak);
    Vijka.reset();
    Vijak.reset();
    Y = SharedTensor2d(new Tensor2d("Y (ME|IJ)", naoccA, navirA, naoccA, naoccA));
    Y->sort(2413, X, 1.0, 0.0);
    X.reset();
    Vij = SharedTensor2d(new Tensor2d("calV (Q|IJ)", nQ, naoccA, naoccA));
    Vij->gemm(false, false, bQiaA, Y, 1.0, 0.0);
    Y.reset();
    Vij->write(psio_, PSIF_DFOCC_AMPS);
    Vij.reset();

    // Build L_ijka
    // L_ijka = \sum{e} L(ja,ie) T(k,e)
    X = SharedTensor2d(new Tensor2d("X <JA|IK>", naoccA, navirA, naoccA, naoccA));
    X->contract(false, true, naoccA*navirA*naoccA, naoccA, navirA, l2, t1A, 1.0, 0.0);
    L = SharedTensor2d(new Tensor2d("L <IJ|KA>", naoccA, naoccA, naoccA, navirA));
    L->sort(3142, X, 1.0, 0.0);
    X.reset();
    //L->write(psio_, PSIF_DFOCC_AMPS);


    /*
    // Build tEta_ia^Q = \sum(ef) L'_ieaf b_ef^Q
    // Build L'_ieaf = 2*Lt_ieaf - Lt_iefa
    // Build L'_ieaf = \sum(m,n) (2*Tau_mn^af - Tau_mn^fa)  L_mnie
    Tau = SharedTensor2d(new Tensor2d("2*Tau(ia,jb) - Tau(ib,ja)", naoccA, navirA, naoccA, navirA));
    Tau->read_symm(psio_, PSIF_DFOCC_AMPS);
    U = SharedTensor2d(new Tensor2d("2*Tau(ij,ab) - Tau(ji,ab)", naoccA, naoccA, navirA, navirA));
    U->sort(1324, Tau, 1.0, 0.0);
    Tau.reset();
    X = SharedTensor2d(new Tensor2d("Lambda <IE|AF>", naoccA, navirA, navirA, navirA));
    X->gemm(true, false, L, U, 1.0, 0.0);
    U.reset();
    Y = SharedTensor2d(new Tensor2d("Temp <EF|IA>", navirA, navirA, naoccA, navirA));
    Y->sort(2413, X, 1.0, 0.0);
    Z = SharedTensor2d(new Tensor2d("Eta2 (Q|IA)", nQ, naoccA, navirA));
    Z->gemm(false, false, bQabA, Y, 1.0, 0.0);
    Y.reset();
    Z->write(psio_, PSIF_DFOCC_AMPS);
    Z.reset();

    // Build tEta_ab^Q = \sum(me) L'_maeb b_me^Q
    Y = SharedTensor2d(new Tensor2d("Temp <ME|AB>", naoccA, navirA, navirA, navirA));
    Y->sort(1324, X, 1.0, 0.0);
    X.reset();
    Z = SharedTensor2d(new Tensor2d("Eta2 (Q|AB)", nQ, navirA, navirA));
    Z->gemm(false, false, bQiaA, Y, 1.0, 0.0);
    Y.reset();
    Z->write(psio_, PSIF_DFOCC_AMPS);
    Z.reset();
    */

    // Build tEta_ia^Q = \sum(ef) L'_ieaf b_ef^Q
    // Build L'_ieaf = 2*Lt_ieaf - Lt_iefa
    // Build L'_ieaf = \sum(m,n) (2*Tau_mn^af - Tau_mn^fa)  L_mnie
    // (+)U(ij, ab) = 1/2 (U'_ij^ab + U'_ji^ab) * (2 - \delta_{ij})
    // (-)U(ij, ab) = 1/2 (U'_ij^ab - U'_ji^ab) * (2 - \delta_{ij})
    Tau = SharedTensor2d(new Tensor2d("2*Tau(ia,jb) - Tau(ib,ja)", naoccA, navirA, naoccA, navirA));
    Tau->read_symm(psio_, PSIF_DFOCC_AMPS);
    U = SharedTensor2d(new Tensor2d("2*Tau(ij,ab) - Tau(ji,ab)", naoccA, naoccA, navirA, navirA));
    U->sort(1324, Tau, 1.0, 0.0);
    Tau.reset();
    Ts = SharedTensor2d(new Tensor2d("(+)U' [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    Ta = SharedTensor2d(new Tensor2d("(-)U' [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    Ts->symm_row_packed4(U);
    Ta->antisymm_row_packed4(U);
    U.reset();

    // Symmetric & Anti-symmetric contributions
    Vs = SharedTensor2d(new Tensor2d("(+)L[I] (E, M>=N)", navirA, ntri_ijAA));
    Va = SharedTensor2d(new Tensor2d("(-)L[I] (E, M>=N)", navirA, ntri_ijAA));
    S = SharedTensor2d(new Tensor2d("S[I] (E, A>=F)", navirA, ntri_abAA));
    A = SharedTensor2d(new Tensor2d("A[I] (E, A>=F)", navirA, ntri_abAA));
    X = SharedTensor2d(new Tensor2d("X[I] (A,EF)", navirA, navirA, navirA));
    Y = SharedTensor2d(new Tensor2d("Eta2[I] (Q|A)", nQ, navirA));
    Z = SharedTensor2d(new Tensor2d("Eta2 (Q|IA)", nQ, naoccA, navirA));
    // Main loop
    for(int i = 0 ; i < naoccA; ++i){

            // Form (+)L[i](e, m>=n)
            #pragma omp parallel for
            for(int e = 0 ; e < navirA; ++e){
                int ie = ia_idxAA->get(i,e);
                for(int m = 0 ; m < naoccA; ++m){
                    for(int n = 0 ; n <= m; ++n){
                        int mn = ij_idxAA->get(m,n);
                        int nm = ij_idxAA->get(n,m);
                        int mn2 = index2(m,n);
                        double value1 = 0.5 * ( L->get(mn, ie) + L->get(nm, ie) );
                        double value2 = 0.5 * ( L->get(mn, ie) - L->get(nm, ie) );
                        Vs->set(e, mn2, value1);
                        Va->set(e, mn2, value2);
                    }
                }
            }

            // Form S[i](e, a>=f) = \sum_{m>=n} U'(m>=n,a>=f) L[i](e, m>=n)
            S->gemm(false, false, Vs, Ts, 1.0, 0.0);
            A->gemm(false, false, Va, Ta, 1.0, 0.0);

            // Form S(ie,af) & A(ie,af)-->X[i](a,ef)
            #pragma omp parallel for
            for(int e = 0 ; e < navirA; ++e){
                for(int a = 0 ; a < navirA; ++a){
                    for(int f = 0 ; f < navirA; ++f){
                        int ef = ab_idxAA->get(e,f);
                        int af = index2(a,f);
                        int perm = ( a > f ) ? 1 : -1;
                        double value = S->get(e,af) + (perm * A->get(e,af));
                        X->set(a, ef, value);
                    }
                }
            }

	    // Eta2[i](Q,a) = \sum(ef) b(Q,ef) * X[i](a,ef)
            Y->gemm(false, true, bQabA, X, 1.0, 0.0);

	    // Form Eta2
            #pragma omp parallel for
            for(int Q = 0 ; Q < nQ; ++Q){
                for(int a = 0 ; a < navirA; ++a){
                    int ia = ia_idxAA->get(i,a);
		    Z->set(Q, ia, Y->get(Q,a));
                }
            }
    }
    Vs.reset();
    Va.reset();
    S.reset();
    A.reset();
    X.reset();
    Y.reset();
    Z->write(psio_, PSIF_DFOCC_AMPS);
    Z.reset();

    // Build tEta_ab^Q = \sum(me) L'_maeb b_me^Q
    // Symmetric & Anti-symmetric contributions
    Vs = SharedTensor2d(new Tensor2d("(+)L[A] (M, K>=N)", naoccA, ntri_ijAA));
    Va = SharedTensor2d(new Tensor2d("(-)L[A] (M, K>=N)", naoccA, ntri_ijAA));
    S = SharedTensor2d(new Tensor2d("S[A] (M, E>=B)", naoccA, ntri_abAA));
    A = SharedTensor2d(new Tensor2d("A[A] (M, E>=B)", naoccA, ntri_abAA));
    X = SharedTensor2d(new Tensor2d("X[A] (B,ME)", navirA, naoccA, navirA));
    Y = SharedTensor2d(new Tensor2d("Eta2[A] (Q|B)", nQ, navirA));
    Z = SharedTensor2d(new Tensor2d("Eta2 (Q|AB)", nQ, navirA, navirA));
    // Main loop
    for(int a = 0 ; a < navirA; ++a){

            // Form (+)L[a](m, k>=n)
            #pragma omp parallel for
            for(int m = 0 ; m < naoccA; ++m){
                int ma = ia_idxAA->get(m,a);
                for(int k = 0 ; k < naoccA; ++k){
                    for(int n = 0 ; n <= k; ++n){
                        int kn = ij_idxAA->get(k,n);
                        int nk = ij_idxAA->get(n,k);
                        int kn2 = index2(k,n);
                        double value1 = 0.5 * ( L->get(kn, ma) + L->get(nk, ma) );
                        double value2 = 0.5 * ( L->get(kn, ma) - L->get(nk, ma) );
                        Vs->set(m, kn2, value1);
                        Va->set(m, kn2, value2);
                    }
                }
            }

            // Form S[a](m, e>=b) = \sum_{k>=n} U'(k>=n,e>=b) L[a](m, k>=n)
            S->gemm(false, false, Vs, Ts, 1.0, 0.0);
            A->gemm(false, false, Va, Ta, 1.0, 0.0);

            // Form S(ma,eb) & A(ma,eb)-->X[a](b,me)
            #pragma omp parallel for
            for(int m = 0 ; m < naoccA; ++m){
                for(int e = 0 ; e < navirA; ++e){
                    int me = ia_idxAA->get(m,e);
                    for(int b = 0 ; b < navirA; ++b){
                        int eb = index2(e,b);
                        int perm = ( e > b ) ? 1 : -1;
                        double value = S->get(m,eb) + (perm * A->get(m,eb));
                        X->set(b, me, value);
                    }
                }
            }

	    // Eta2[a](Q,b) = \sum(me) b(Q,me) * X[a](b,me)
            Y->gemm(false, true, bQiaA, X, 1.0, 0.0);

	    // Form Eta2
            #pragma omp parallel for
            for(int Q = 0 ; Q < nQ; ++Q){
                for(int b = 0 ; b < navirA; ++b){
                    int ab = ab_idxAA->get(a,b);
		    Z->set(Q, ab, Y->get(Q,b));
                }
            }
    }
    Vs.reset();
    Va.reset();
    Ts.reset();
    Ta.reset();
    S.reset();
    A.reset();
    X.reset();
    Y.reset();
    Z->write(psio_, PSIF_DFOCC_AMPS);
    Z.reset();

    // Z_ij^Q = \sum_{me} (2*L_imje - L_mije) t_me^Q
    X = SharedTensor2d(new Tensor2d("X <IM|JE>", naoccA, naoccA, naoccA, navirA));
    X->tei_cs2_anti_symm(L,L);
    L.reset();
    // Y_meij = X_imje
    Y = SharedTensor2d(new Tensor2d("Y <ME|IJ>", naoccA, navirA, naoccA, naoccA));
    Y->sort(2413, X, 1.0, 0.0);
    // Build Z_ij^Q
    T = SharedTensor2d(new Tensor2d("T1 (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    Z = SharedTensor2d(new Tensor2d("Zeta (Q|IJ)", nQ, naoccA, naoccA));
    Z->gemm(false, false, T, Y, 1.0, 0.0);
    T.reset();
    Z->write(psio_, PSIF_DFOCC_AMPS);
    Z.reset();

    // n_ij^Q = \sum_{me} (2*L_imje - L_mije) b_me^Q
    Z = SharedTensor2d(new Tensor2d("Eta (Q|IJ)", nQ, naoccA, naoccA));
    Z->gemm(false, false, bQiaA, Y, 1.0, 0.0);
    Y.reset();
    Z->write(psio_, PSIF_DFOCC_AMPS);
    Z.reset();

    // Z_ia^Q = -\sum_{mn} (2*L_mina - L_imna) t_mn^Q
    // Y_mnia = X_mina
    Y = SharedTensor2d(new Tensor2d("Y <MN|IA>", naoccA, naoccA, naoccA, navirA));
    Y->sort(1324, X, 1.0, 0.0);
    X.reset();
    // Build Z_ia^Q
    T = SharedTensor2d(new Tensor2d("T1 (Q|IJ)", nQ, naoccA, naoccA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    Z = SharedTensor2d(new Tensor2d("Zeta (Q|IA)", nQ, naoccA, navirA));
    Z->gemm(false, false, T, Y, -1.0, 0.0);
    T.reset();
    Z->write(psio_, PSIF_DFOCC_AMPS);
    Z.reset();

    // n_ia^Q = -\sum_{mn} (2*L_mina - L_imna) b_mn^Q
    Z = SharedTensor2d(new Tensor2d("Eta (Q|IA)", nQ, naoccA, navirA));
    Z->gemm(false, false, bQijA, Y, -1.0, 0.0);
    Y.reset();
    Z->write(psio_, PSIF_DFOCC_AMPS);
    Z.reset();

    //outfile->Printf("\t3indices done.\n");

}// end ccsd_pdm_3index_intr

//======================================================================
//    Build y_ia^Q
//======================================================================
void DFOCC::ccsd_pdm_yQia()
{
    // defs
    SharedTensor2d K, L, T, U, Tau, V, V2, X, Y, Y2, Z;
    SharedTensor2d Yt, Yp;

    // Read Viajb
    V = SharedTensor2d(new Tensor2d("V (IA|JB)", naoccA, navirA, naoccA, navirA));
    V->read(psio_, PSIF_DFOCC_AMPS);

    // Read T2
    t2 = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    t2->read_symm(psio_, PSIF_DFOCC_AMPS);

    // Y_iajb = -\sum(me) T(ia,me) X(me,jb)
    // X(me,jb) = V(je,mb)
    X = SharedTensor2d(new Tensor2d("X (ME|JB)", naoccA, navirA, naoccA, navirA));
    X->sort(3214, V, 1.0, 0.0);
    Y = SharedTensor2d(new Tensor2d("Y (IA|JB)", naoccA, navirA, naoccA, navirA));
    Y->gemm(false, false, t2, X, -1.0, 0.0);
    X.reset();

    // build T_mi^ae + 2 * T_m^a * T_i^e
    T = SharedTensor2d(new Tensor2d("T2p (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_t2_prime_amps(T,t2);
    t2.reset();
    #pragma omp parallel for
    for(int i = 0 ; i < naoccA; ++i){
        for(int j = 0 ; j < naoccA; ++j){
            for(int a = 0 ; a < navirA; ++a){
                int ia = ia_idxAA->get(i,a);
                for(int b = 0 ; b < navirA; ++b){
                    int jb = ia_idxAA->get(j,b);
                    T->add(ia, jb, 2.0 * t1A->get(i,b) * t1A->get(j,a) );
                }
            }
        }
    }

    // Y_iabj = -\sum(me) T'(ia,me) X(me,bj)
    // X(me,bj) = V(je,mb)
    X = SharedTensor2d(new Tensor2d("X (ME|BJ)", naoccA, navirA, navirA, naoccA));
    X->sort(3241, V, 1.0, 0.0);
    V.reset();
    Y2 = SharedTensor2d(new Tensor2d("Y (IA|BJ)", naoccA, navirA, navirA, naoccA));
    Y2->gemm(false, false, T, X, -1.0, 0.0);
    T.reset();
    X.reset();
    Y2->write(psio_, PSIF_DFOCC_AMPS);
    Y2.reset();

    // Read Viabj
    V = SharedTensor2d(new Tensor2d("V (IA|BJ)", naoccA, navirA, navirA, naoccA));
    V->read(psio_, PSIF_DFOCC_AMPS);

    // Read U2
    U = SharedTensor2d(new Tensor2d("U2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    U->read_symm(psio_, PSIF_DFOCC_AMPS);

    // Build Vt_iabj = V_iabj - T_i^b * L_j^a
    // build U_im^ae - 2 * T_m^a * T_i^e
    #pragma omp parallel for
    for(int i = 0 ; i < naoccA; ++i){
        for(int j = 0 ; j < naoccA; ++j){
            for(int a = 0 ; a < navirA; ++a){
                int ia = ia_idxAA->get(i,a);
                for(int b = 0 ; b < navirA; ++b){
                    int jb = ia_idxAA->get(j,b);
                    int bj = ai_idxAA->get(b,j);
                    U->subtract(ia, jb, 2.0 * t1A->get(i,b) * t1A->get(j,a) );
                    V->subtract(ia, bj, t1A->get(i,b) * l1A->get(j,a) );
                }
            }
        }
    }

    // Y_iajb += \sum(me) [U(ia,me) + 2 * T_m^a * T_i^e] X(me,jb)
    // X(me,jb) = Vt_(je,bm)
    X = SharedTensor2d(new Tensor2d("X (ME|JB)", naoccA, navirA, naoccA, navirA));
    X->sort(4213, V, 1.0, 0.0);
    V.reset();
    Y->gemm(false, false, U, X, 1.0, 1.0);
    U.reset();
    X.reset();

    //===================================
    // Build Yt_ijab
    //===================================

    // Yt_ijab = -1/2 (Y_iajb + Y_jbia)
    Yt = SharedTensor2d(new Tensor2d("Y2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    Yt->sort(1324, Y, -0.5, 0.0);
    Yt->sort(3142, Y, -0.5, 1.0);
    Y.reset();

    // Yt_ijab += 1/2 Y_ijab
    U = SharedTensor2d(new Tensor2d("Y <IJ|AB>", naoccA, naoccA, navirA, navirA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    Yt->axpy(U, 0.5);
    U.reset();

    // Yt_ijab -= 1/2 (Y_jabi + Y_ibaj)
    Y = SharedTensor2d(new Tensor2d("Y (IA|BJ)", naoccA, navirA, navirA, naoccA));
    Y->read(psio_, PSIF_DFOCC_AMPS);
    Yt->sort(4123, Y, -0.5, 1.0);
    Yt->sort(1432, Y, -0.5, 1.0);
    Y.reset();

    // Build S_i^a = \sum(m) T_m^a \sum(e) T_i^e * L_m^e
    SharedTensor2d Sij = SharedTensor2d(new Tensor2d("S <I|J>", naoccA, naoccA));
    SharedTensor2d s1 = SharedTensor2d(new Tensor2d("S <I|A>", naoccA, navirA));
    Sij->gemm(false, true, t1A, l1A, 1.0, 0.0);
    s1->gemm(false, false, Sij, t1A, 1.0, 0.0);
    Sij.reset();

    // Yt_ijab += 3/2 (t_i^a * s_j^b + t_j^b * t_i^a)
    #pragma omp parallel for
    for(int i = 0 ; i < naoccA; ++i){
        for(int j = 0 ; j < naoccA; ++j){
            int ij = ij_idxAA->get(i,j);
            for(int a = 0 ; a < navirA; ++a){
                for(int b = 0 ; b < navirA; ++b){
                    int ab = ab_idxAA->get(a,b);
		    double value = ( t1A->get(i,a) * s1->get(j,b) ) + ( t1A->get(j,b) * s1->get(i,a) );
                    Yt->add(ij, ab, 1.5 * value);
                }
            }
        }
    }
    s1.reset();

    // Yt_ijab -= V_jabi + V_ibaj
    V = SharedTensor2d(new Tensor2d("V (IA|BJ)", naoccA, navirA, navirA, naoccA));
    V->read(psio_, PSIF_DFOCC_AMPS);
    Yt->sort(4123, V, -1.0, 1.0);
    Yt->sort(1432, V, -1.0, 1.0);

    //===================================
    // Build Y'_ijab
    //===================================

    // Y'_ijab = Yt_ijab + V_jabi + V_ibaj
    Yp = SharedTensor2d(new Tensor2d("Y2p <IJ|AB>", naoccA, naoccA, navirA, navirA));
    Yp->copy(Yt);
    Yp->sort(4123, V, 1.0, 1.0);
    Yp->sort(1432, V, 1.0, 1.0);
    V.reset();

    // Y'_ijab -= V_iajb + V_jbia
    V = SharedTensor2d(new Tensor2d("V (IA|JB)", naoccA, navirA, naoccA, navirA));
    V->read(psio_, PSIF_DFOCC_AMPS);
    Yp->sort(1324, V, -1.0, 1.0);
    Yp->sort(3142, V, -1.0, 1.0);
    V.reset();

    //===================================
    // Build y_ia^Q
    //===================================

    // X_imae = 2*Yt_imae - Yp_imea
    X = SharedTensor2d(new Tensor2d("X <IM|AE>", naoccA, naoccA, navirA, navirA));
    X->tei_cs1_anti_symm(Yt,Yp);
    Yt.reset();
    Yp.reset();
    // Y(me,ia) = X(im,ae)
    Y = SharedTensor2d(new Tensor2d("Y (ME|IA)", naoccA, navirA, naoccA, navirA));
    Y->sort(2413, X, 1.0, 0.0);
    X.reset();
    // y_ia^Q = \sum(me) b(Q,me) * Y(me,ia)
    Z = SharedTensor2d(new Tensor2d("Y (Q|IA)", nQ, naoccA, navirA));
    Z->gemm(false, false, bQiaA, Y, 1.0, 0.0);
    Y.reset();
    Z->write(psio_, PSIF_DFOCC_AMPS);
    Z.reset();

}// end ccsd_pdm_yQia

}} // End Namespaces
