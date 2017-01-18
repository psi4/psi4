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

void DFOCC::ccd_pdm_3index_intr()
{

    // defs
    SharedTensor2d K, L, T, U, Tau, V, V2, Vij, Vai, Vab, X, Y, Z;
    SharedTensor2d Vijka, Vijak;
    SharedTensor2d Vs, Va, Ts, Ta, S, A;

    // L(Q,ia) = \sum_{jb} b_jb^Q Ut_ij^ab = \sum_{jb} b(Q,jb) Ut(jb,ia)
    U = SharedTensor2d(new Tensor2d("Ut2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_u2_amps(U,l2);
    T = SharedTensor2d(new Tensor2d("L2 (Q|IA)", nQ, naoccA, navirA));
    T->gemm(false, false, bQiaA, U, 1.0, 0.0);
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

    // Gt_Q = 2\sum_{mn} G_mn b_mn^Q
    gQt->gemv(false, bQijA, GijA, 2.0, 0.0);

    // G(Q,ij) = \sum_{m} G_im b_mj^Q
    T = SharedTensor2d(new Tensor2d("G (Q|IJ)", nQ, naoccA, naoccA));
    T->contract233(false, false, naoccA, naoccA, GijA, bQijA, 1.0, 0.0);
    //T->cont233("IJ", "IM", "MJ", GijA, bQijA, 1.0, 0.0); // it works
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

    // G(Q,ai) = \sum_{e} G_ae b_ie^Q
    T = SharedTensor2d(new Tensor2d("G (Q|AI)", nQ, navirA, naoccA));
    T->contract233(false, true, navirA, naoccA, GabA, bQiaA, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

    // G(Q,ia) = \sum_{m} G_mi b_ma^Q
    T = SharedTensor2d(new Tensor2d("G (Q|IA)", nQ, naoccA, navirA));
    T->contract233(true, false, naoccA, navirA, GijA, bQiaA, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

    // Build V_ijkl
    // V_ijkl = \sum_{ef} T_ij^ef L_kl^ef
    t2 = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    t2->read_symm(psio_, PSIF_DFOCC_AMPS);
    U = SharedTensor2d(new Tensor2d("T <IJ|AB>", naoccA, naoccA, navirA, navirA));
    U->sort(1324, t2, 1.0, 0.0);
    //t2.reset();
    L = SharedTensor2d(new Tensor2d("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    L->sort(1324, l2, 1.0, 0.0);
    V = SharedTensor2d(new Tensor2d("V <IJ|KL>", naoccA, naoccA, naoccA, naoccA));
    V->gemm(false, true, U, L, 1.0, 0.0);
    L.reset();

    // Y_ijab = \sum(mn) T_mn^ab V_ijmn
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
    X.reset();
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
    X.reset();
    Vab->write(psio_, PSIF_DFOCC_AMPS);
    Vab.reset();

    //outfile->Printf("\t3indices done.\n");

}// end ccd_pdm_3index_intr

//======================================================================
//    Build y_ia^Q
//======================================================================
void DFOCC::ccd_pdm_yQia()
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

    // Y_iabj = -\sum(me) T'(ia,me) X(me,bj)
    T = SharedTensor2d(new Tensor2d("T2p (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_t2_prime_amps(T,t2);
    t2.reset();
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

    // Y_iajb += \sum(me) U(ia,me) X(me,jb)
    // X(me,jb) = V_(je,bm)
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

}// end ccd_pdm_yQia

}} // End Namespaces
