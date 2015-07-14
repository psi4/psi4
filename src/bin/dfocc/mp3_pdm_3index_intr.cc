/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#include <libqt/qt.h>
#include "defines.h"
#include "dfocc.h"

using namespace psi;
using namespace std;


namespace psi{ namespace dfoccwave{
  
void DFOCC::mp3_pdm_3index_intr()
{

    // defs
    SharedTensor2d K, L, T, T1, T2, U, Tau, V, V2, Vij, Vai, Vab, X, Y, Z;
    SharedTensor2d Vijka, Vijak;
    SharedTensor2d Vs, Va, Ts, Ta, S, A;
    SharedTensor2d Yt, Yp;

    // Read DF integrals
    bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA));
    bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA));
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQijA->read(psio_, PSIF_DFOCC_INTS);
    bQiaA->read(psio_, PSIF_DFOCC_INTS);
    bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);


    // Read amps
    U = SharedTensor2d(new Tensor2d("U2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    U->read_symm(psio_, PSIF_DFOCC_AMPS);
    T1 = SharedTensor2d(new Tensor2d("T2_1 (IA|JB)", naoccA, navirA, naoccA, navirA));
    T1->read_symm(psio_, PSIF_DFOCC_AMPS);

    // T(Q,ia) = \sum_{jb} b_jb^Q U_ij^ab = \sum_{jb} b(Q,jb) U(jb,ia)
    T = SharedTensor2d(new Tensor2d("T2 (Q|IA)", nQ, naoccA, navirA));
    T->gemm(false, false, bQiaA, U, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

    // Build G_mi and G_ae
    // G_mi = \sum_{n,e,f} U_mn^ef L_in^ef(1) = U(mn,ef) L(in,ef)(1)
    T = SharedTensor2d(new Tensor2d("U2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    T->sort(1324, U, 1.0, 0.0);
    U.reset();
    L = SharedTensor2d(new Tensor2d("T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    L->sort(1324, T1, 1.0, 0.0);
    GijA->contract(false, true, naoccA, naoccA, naoccA*navirA*navirA, T, L, 1.0, 0.0);

    // G_ae = -\sum_{m,n,f} U_mn^ef L_mn^af(1) = L(mn,fa)(1) U(mn,fe) 
    GabA->contract(true, false, navirA, navirA, naoccA*naoccA*navirA, L, T, -1.0, 0.0);
    T.reset();
    L.reset();

    // Read amps
    T2 = SharedTensor2d(new Tensor2d("T2_2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    T2->read_symm(psio_, PSIF_DFOCC_AMPS);

    // G_mi += \sum_{n,e,f} U_mn^ef(1) L_in^ef(2) = U(mn,ef)(1) L(in,ef)(2)
    U = SharedTensor2d(new Tensor2d("U2_1 (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_u2_amps(U,T1);
    T = SharedTensor2d(new Tensor2d("U2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    T->sort(1324, U, 1.0, 0.0);
    U.reset();
    L = SharedTensor2d(new Tensor2d("T2_2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    L->sort(1324, T2, 1.0, 0.0);
    T2.reset();
    GijA->contract(false, true, naoccA, naoccA, naoccA*navirA*navirA, T, L, 1.0, 1.0);

    // G_ae += -\sum_{m,n,f} U_mn^ef(1) L_mn^af(2) = L(mn,fa)(2) U(mn,fe)(1) 
    GabA->contract(true, false, navirA, navirA, naoccA*naoccA*navirA, L, T, -1.0, 1.0);
    T.reset();
    L.reset();


    // Build V_ijkl
    // V_ijkl = \sum_{ef} T_ij^ef(1) L_kl^ef(1)
    U = SharedTensor2d(new Tensor2d("T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    U->sort(1324, T1, 1.0, 0.0);
    V = SharedTensor2d(new Tensor2d("V <IJ|KL>", naoccA, naoccA, naoccA, naoccA));
    V->gemm(false, true, U, U, 1.0, 0.0);
    U.reset();
    //V->write(psio_, PSIF_DFOCC_AMPS);

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
    // V_iajb = 1/2 \sum_{me} T'(ib,me)(1) L'(me,ja)(1)
    T = SharedTensor2d(new Tensor2d("T2p (IB|JA)", naoccA, navirA, naoccA, navirA));
    ccsd_t2_prime_amps(T,T1);
    X = SharedTensor2d(new Tensor2d("X (IB|JA)", naoccA, navirA, naoccA, navirA));
    X->gemm(false, false, T, T, 0.5, 0.0);
    T.reset();
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
    // V_iabj = -1/2 \sum_{me} U(ib,me)(1) L(me,ja)(1)
    U = SharedTensor2d(new Tensor2d("U2_1 (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_u2_amps(U,T1);
    X = SharedTensor2d(new Tensor2d("X (IB|JA)", naoccA, navirA, naoccA, navirA));
    X->gemm(false, false, U, T1, -0.5, 0.0);
    U.reset();
    // V_iabj += 1/2 \sum_{me} T(ib,me)(1) L'(me,ja)(1)
    L = SharedTensor2d(new Tensor2d("L2p (IB|JA)", naoccA, navirA, naoccA, navirA));
    ccsd_t2_prime_amps(L,T1);
    X->gemm(false, false, T1, L, 0.5, 1.0);
    T1.reset();
    L.reset();
    V = SharedTensor2d(new Tensor2d("V (IA|BJ)", naoccA, navirA, navirA, naoccA));
    V->sort(1423, X, 1.0, 0.0);
    X.reset();
    V->write(psio_, PSIF_DFOCC_AMPS);

    // V_ij^Q' -= \sum_{ef} V_iefj b_ef^Q
    X = SharedTensor2d(new Tensor2d("X (AB|IJ)", navirA, navirA, naoccA, naoccA));
    X->sort(2314, V, 1.0, 0.0);
    Vij->gemm(false, false, bQabA, X, -1.0, 1.0); 
    X.reset();
    Vij->write(psio_, PSIF_DFOCC_AMPS);
    Vij.reset();

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

    //===================================
    // Build Yt_ijab
    //===================================

    // Yt_ijab -= V_jabi + V_ibaj
    Yt = SharedTensor2d(new Tensor2d("Y2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    V = SharedTensor2d(new Tensor2d("V (IA|BJ)", naoccA, navirA, navirA, naoccA));
    V->read(psio_, PSIF_DFOCC_AMPS);
    Yt->sort(4123, V, -1.0, 0.0);
    Yt->sort(1432, V, -1.0, 1.0);
    V.reset();

    //===================================
    // Build Y'_ijab
    //===================================

    // Y'_ijab -= V_iajb + V_jbia
    Yp = SharedTensor2d(new Tensor2d("Y2p <IJ|AB>", naoccA, naoccA, navirA, navirA));
    V = SharedTensor2d(new Tensor2d("V (IA|JB)", naoccA, navirA, naoccA, navirA));
    V->read(psio_, PSIF_DFOCC_AMPS);
    Yp->sort(1324, V, -1.0, 0.0);
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

    // Free ints
    bQijA.reset();
    bQiaA.reset();
    bQabA.reset();
    //outfile->Printf("\t3indices done.\n");

}// end mp3_pdm_3index_intr


}} // End Namespaces

