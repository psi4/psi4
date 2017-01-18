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

void DFOCC::ccsdl_l1_amps()
{

    // defs
    SharedTensor2d K, L, T1, T, U, Tau, X, Y, Z, W, W2, V;

    // l_i^a <= Ftia
    //FiaA->print();
    l1newA->copy(FiaA);

    // l_i^a <= \sum_{e} l_i^e Ft_ea
    l1newA->gemm(false, false, l1A, FtabA, 1.0, 1.0);

    // l_i^a <= -\sum_{m} l_m^a Ft_im
    l1newA->gemm(false, false, FtijA, l1A, -1.0, 1.0);

    // l_i^a <= -\sum_{m} G_mi F_ma
    l1newA->gemm(true, false, GijA, FiaA, -1.0, 1.0);

    // l_i^a <= \sum_{me} l_m^e (2*W_ieam - W_iema)
    // l_i^a <= \sum_{me} l_m^e [2*W(ia,me) - W'(ia,me)]
    X = SharedTensor2d(new Tensor2d("X (ME|JB)", naoccA, navirA, naoccA, navirA));
    W = SharedTensor2d(new Tensor2d("WL (ME|JB)", naoccA, navirA, naoccA, navirA));
    W->read(psio_, PSIF_DFOCC_AMPS);
    X->axpy(W, 2.0);
    W.reset();
    W = SharedTensor2d(new Tensor2d("WLp (ME|JB)", naoccA, navirA, naoccA, navirA));
    W->read(psio_, PSIF_DFOCC_AMPS);
    X->axpy(W, -1.0);
    W.reset();
    l1newA->gemv(false, X, l1A, 1.0, 1.0);
    X.reset();

    // l_i^a <= -\sum_{mn} G_mn (2*W_mina - W_imna)
    W = SharedTensor2d(new Tensor2d("WL (MN|IE)", naoccA, naoccA, naoccA, navirA));
    //W->read(psio_, PSIF_DFOCC_AMPS);
    ccsdl_Wmnie_direct(W);
    X = SharedTensor2d(new Tensor2d("X (MN|IE)", naoccA, naoccA, naoccA, navirA));
    X->tei_cs2_anti_symm(W,W);
    W.reset();
    Y = SharedTensor2d(new Tensor2d("Y (IA|MN)", naoccA, navirA, naoccA, naoccA));
    // Y_iamn = X_mina
    Y->sort(2413, X, 1.0, 0.0);
    X.reset();
    l1newA->gemv(false, Y, GijA, -1.0, 1.0);
    Y.reset();

    // l_i^a <= -\sum_{mne} Ut_mn^ae W_iemn
    U = SharedTensor2d(new Tensor2d("Ut2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_u2_amps(U,l2);
    X = SharedTensor2d(new Tensor2d("X (BI|JA)", navirA, naoccA, naoccA, navirA));
    // X_emna = Ut(ma,ne)
    X->sort(4132, U, 1.0, 0.0);
    U.reset();
    W = SharedTensor2d(new Tensor2d("WL (MB|IJ)", naoccA, navirA, naoccA, naoccA));
    W->read(psio_, PSIF_DFOCC_AMPS);
    l1newA->contract(false, false, naoccA, navirA, navirA*naoccA*naoccA, W, X, -1.0, 1.0);
    W.reset();
    X.reset();

    // l_i^a <= -\sum_{mne} (2*L_imne - L_mine) Z_neam
    // l_i^a <= -\sum_{mne} (2*L_imne - L_mine) Z(na,me)
    L = SharedTensor2d(new Tensor2d("L <IJ|KA>", naoccA, naoccA, naoccA, navirA));
    L->read(psio_, PSIF_DFOCC_AMPS);
    K = SharedTensor2d(new Tensor2d("K <IJ|KA>", naoccA, naoccA, naoccA, navirA));
    K->tei_cs2_anti_symm(L,L);
    L.reset();
    Z = SharedTensor2d(new Tensor2d("Z (ME|JB)", naoccA, navirA, naoccA, navirA));
    Z->read(psio_, PSIF_DFOCC_AMPS);
    X = SharedTensor2d(new Tensor2d("X (IJ|AB)", naoccA, naoccA, navirA, navirA));
    // X_mnea = Z(na,me)
    X->sort(3142, Z, 1.0, 0.0);
    Z.reset();
    l1newA->contract(false, false, naoccA, navirA, navirA*naoccA*naoccA, K, X, -1.0, 1.0);
    X.reset();

    // l_i^a <= -\sum_{mne} (2*L_mine - L_imne) Z_nema
    // l_i^a <= -\sum_{mne} (2*L_mine - L_imne) Z'(na,me)
    Y = SharedTensor2d(new Tensor2d("Y <IJ|KA>", naoccA, naoccA, naoccA, navirA));
    Y->sort(2134, K, 1.0, 0.0);
    K.reset();
    Z = SharedTensor2d(new Tensor2d("Zp (ME|JB)", naoccA, navirA, naoccA, navirA));
    Z->read(psio_, PSIF_DFOCC_AMPS);
    X = SharedTensor2d(new Tensor2d("X (IJ|AB)", naoccA, naoccA, navirA, navirA));
    // X_mnea = Z(na,me)
    X->sort(3142, Z, 1.0, 0.0);
    Z.reset();
    l1newA->contract(false, false, naoccA, navirA, navirA*naoccA*naoccA, Y, X, -1.0, 1.0);
    X.reset();
    Y.reset();

    // l_i^a <= \sum_{Q,e} (L_ie^Q + Lt_ie^Q + 2V_ei^Q + Z_ie^Q) b_ea^Q
    T = SharedTensor2d(new Tensor2d("Temp (Q|IA)", nQ, naoccA, navirA));
    V = SharedTensor2d(new Tensor2d("V (Q|AI)", nQ, navirA, naoccA));
    V->read(psio_, PSIF_DFOCC_AMPS);
    T->swap_3index_col(V);
    V.reset();
    T->scale(2.0);
    U = SharedTensor2d(new Tensor2d("L2 (Q|IA)", nQ, naoccA, navirA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->add(U);
    U.reset();
    U = SharedTensor2d(new Tensor2d("L2t (Q|IA)", nQ, naoccA, navirA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->add(U);
    U.reset();
    U = SharedTensor2d(new Tensor2d("Zeta (Q|IA)", nQ, naoccA, navirA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->add(U);
    U.reset();
    Tau = SharedTensor2d(new Tensor2d("Temp (Q|AI)", nQ, navirA, naoccA));
    Tau->swap_3index_col(T);
    T.reset();
    l1newA->contract(true, false, naoccA, navirA, nQ * navirA, Tau, bQabA, 1.0, 1.0);
    Tau.reset();

    // l_i^a <= \sum_{Q,m} (V_mi^Q + Vt_mi^Q - 2V'_mi^Q - Z_im^Q) b_ma^Q
    T = SharedTensor2d(new Tensor2d("Temp (Q|IJ)", nQ, naoccA, naoccA));
    U = SharedTensor2d(new Tensor2d("V (Q|IJ)", nQ, naoccA, naoccA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->copy(U);
    U.reset();
    U = SharedTensor2d(new Tensor2d("Vt (Q|IJ)", nQ, naoccA, naoccA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->add(U);
    U.reset();
    U = SharedTensor2d(new Tensor2d("Vp (Q|IJ)", nQ, naoccA, naoccA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, -2.0);
    U.reset();
    U = SharedTensor2d(new Tensor2d("Zeta (Q|IJ)", nQ, naoccA, naoccA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    Tau = SharedTensor2d(new Tensor2d("Tau (Q|IJ)", nQ, naoccA, naoccA));
    Tau->swap_3index_col(U);
    U.reset();
    T->axpy(Tau, -1.0);
    Tau.reset();
    l1newA->contract(true, false, naoccA, navirA, nQ * naoccA, T, bQiaA, 1.0, 1.0);
    T.reset();

    // l_i^a <= \sum_{Q} (Gp_Q - G_Q) b_ia^Q
    SharedTensor1d gQp2 = SharedTensor1d(new Tensor1d("CCSDL G_Qp - G_Q", nQ));
    gQp2->copy(gQp);
    gQp2->axpy(gQ, -1.0);
    l1newA->gemv(true, bQiaA, gQp2, 1.0, 1.0);
    gQp2.reset();

    // l_i^a <= \sum_{Q,e} G_ei^Q (b_ea^Q - t_ea^Q)
    T = SharedTensor2d(new Tensor2d("T1 (Q|AB)", nQ, navirA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    T->scale(-1.0);
    T->add(bQabA);
    U = SharedTensor2d(new Tensor2d("G (Q|AI)", nQ, navirA, naoccA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    l1newA->contract(true, false, naoccA, navirA, nQ * navirA, U, T, 1.0, 1.0);
    T.reset();
    U.reset();

    // Denom
    for(int i = 0 ; i < naoccA; ++i){
        for(int a = 0 ; a < navirA; ++a){
            double value = FockA->get(i + nfrzc, i + nfrzc) - FockA->get(a + noccA, a + noccA);
            l1newA->set(i, a, l1newA->get(i, a) / value);
        }
    }
    //l1newA->print();

}// end ccsdl_l1_amps
}} // End Namespaces
