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

void DFOCC::ccsdl_l1_amps() {

    // RHF
    if (reference_ == "RESTRICTED") {
        // defs
        SharedTensor2d K, L, T1, T, U, Tau, X, Y, Z, W, W2, V, tL1;

        // l_i^a <= Ftia
        // FiaA->print();
        l1newA->copy(FiaA);

        // l_i^a <= \sum_{e} l_i^e Ft_ea
        l1newA->gemm(false, false, l1A, FtabA, 1.0, 1.0);

        // l_i^a <= -\sum_{m} l_m^a Ft_im
        l1newA->gemm(false, false, FtijA, l1A, -1.0, 1.0);

        // l_i^a <= -\sum_{m} G_mi F_ma
        l1newA->gemm(true, false, GijA, FiaA, -1.0, 1.0);

        // l_i^a <= \sum_{me} l_m^e (2*W_ieam - W_iema)
        // l_i^a <= \sum_{me} l_m^e [2*W(ia,me) - W'(ia,me)]
        X = std::make_shared<Tensor2d>("X (ME|JB)", naoccA, navirA, naoccA, navirA);
        W = std::make_shared<Tensor2d>("WL (ME|JB)", naoccA, navirA, naoccA, navirA);
        W->read(psio_, PSIF_DFOCC_AMPS);
        X->axpy(W, 2.0);
        W.reset();
        W = std::make_shared<Tensor2d>("WLp (ME|JB)", naoccA, navirA, naoccA, navirA);
        W->read(psio_, PSIF_DFOCC_AMPS);
        X->axpy(W, -1.0);
        W.reset();
        l1newA->gemv(false, X, l1A, 1.0, 1.0);
        X.reset();

        // l_i^a <= -\sum_{mn} G_mn (2*W_mina - W_imna)
        W = std::make_shared<Tensor2d>("WL (MN|IE)", naoccA, naoccA, naoccA, navirA);
        // W->read(psio_, PSIF_DFOCC_AMPS);
        ccsdl_Wmnie_direct(W);
        X = std::make_shared<Tensor2d>("X (MN|IE)", naoccA, naoccA, naoccA, navirA);
        X->tei_cs2_anti_symm(W, W);
        W.reset();
        Y = std::make_shared<Tensor2d>("Y (IA|MN)", naoccA, navirA, naoccA, naoccA);
        // Y_iamn = X_mina
        Y->sort(2413, X, 1.0, 0.0);
        X.reset();
        l1newA->gemv(false, Y, GijA, -1.0, 1.0);
        Y.reset();

        // l_i^a <= -\sum_{mne} Ut_mn^ae W_iemn
        U = std::make_shared<Tensor2d>("Ut2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        ccsd_u2_amps(U, l2);
        X = std::make_shared<Tensor2d>("X (BI|JA)", navirA, naoccA, naoccA, navirA);
        // X_emna = Ut(ma,ne)
        X->sort(4132, U, 1.0, 0.0);
        U.reset();
        W = std::make_shared<Tensor2d>("WL (MB|IJ)", naoccA, navirA, naoccA, naoccA);
        W->read(psio_, PSIF_DFOCC_AMPS);
        l1newA->contract(false, false, naoccA, navirA, navirA * naoccA * naoccA, W, X, -1.0, 1.0);
        W.reset();
        X.reset();

        // l_i^a <= -\sum_{mne} (2*L_imne - L_mine) Z_neam
        // l_i^a <= -\sum_{mne} (2*L_imne - L_mine) Z(na,me)
        L = std::make_shared<Tensor2d>("L <IJ|KA>", naoccA, naoccA, naoccA, navirA);
        L->read(psio_, PSIF_DFOCC_AMPS);
        K = std::make_shared<Tensor2d>("K <IJ|KA>", naoccA, naoccA, naoccA, navirA);
        K->tei_cs2_anti_symm(L, L);
        L.reset();
        Z = std::make_shared<Tensor2d>("Z (ME|JB)", naoccA, navirA, naoccA, navirA);
        Z->read(psio_, PSIF_DFOCC_AMPS);
        X = std::make_shared<Tensor2d>("X (IJ|AB)", naoccA, naoccA, navirA, navirA);
        // X_mnea = Z(na,me)
        X->sort(3142, Z, 1.0, 0.0);
        Z.reset();
        l1newA->contract(false, false, naoccA, navirA, navirA * naoccA * naoccA, K, X, -1.0, 1.0);
        X.reset();

        // l_i^a <= -\sum_{mne} (2*L_mine - L_imne) Z_nema
        // l_i^a <= -\sum_{mne} (2*L_mine - L_imne) Z'(na,me)
        Y = std::make_shared<Tensor2d>("Y <IJ|KA>", naoccA, naoccA, naoccA, navirA);
        Y->sort(2134, K, 1.0, 0.0);
        K.reset();
        Z = std::make_shared<Tensor2d>("Zp (ME|JB)", naoccA, navirA, naoccA, navirA);
        Z->read(psio_, PSIF_DFOCC_AMPS);
        X = std::make_shared<Tensor2d>("X (IJ|AB)", naoccA, naoccA, navirA, navirA);
        // X_mnea = Z(na,me)
        X->sort(3142, Z, 1.0, 0.0);
        Z.reset();
        l1newA->contract(false, false, naoccA, navirA, navirA * naoccA * naoccA, Y, X, -1.0, 1.0);
        X.reset();
        Y.reset();

        // l_i^a <= \sum_{Q,e} (L_ie^Q + Lt_ie^Q + 2V_ei^Q + Z_ie^Q) b_ea^Q
        T = std::make_shared<Tensor2d>("Temp (Q|IA)", nQ, naoccA, navirA);
        V = std::make_shared<Tensor2d>("V (Q|AI)", nQ, navirA, naoccA);
        V->read(psio_, PSIF_DFOCC_AMPS);
        T->swap_3index_col(V);
        V.reset();
        T->scale(2.0);
        U = std::make_shared<Tensor2d>("L2 (Q|IA)", nQ, naoccA, navirA);
        U->read(psio_, PSIF_DFOCC_AMPS);
        T->add(U);
        U.reset();
        U = std::make_shared<Tensor2d>("L2t (Q|IA)", nQ, naoccA, navirA);
        U->read(psio_, PSIF_DFOCC_AMPS);
        T->add(U);
        U.reset();
        U = std::make_shared<Tensor2d>("Zeta (Q|IA)", nQ, naoccA, navirA);
        U->read(psio_, PSIF_DFOCC_AMPS);
        T->add(U);
        U.reset();
        Tau = std::make_shared<Tensor2d>("Temp (Q|AI)", nQ, navirA, naoccA);
        Tau->swap_3index_col(T);
        T.reset();
        l1newA->contract(true, false, naoccA, navirA, nQ * navirA, Tau, bQabA, 1.0, 1.0);
        Tau.reset();

        // l_i^a <= \sum_{Q,m} (V_mi^Q + Vt_mi^Q - 2V'_mi^Q - Z_im^Q) b_ma^Q
        T = std::make_shared<Tensor2d>("Temp (Q|IJ)", nQ, naoccA, naoccA);
        U = std::make_shared<Tensor2d>("V (Q|IJ)", nQ, naoccA, naoccA);
        U->read(psio_, PSIF_DFOCC_AMPS);
        T->copy(U);
        U.reset();
        U = std::make_shared<Tensor2d>("Vt (Q|IJ)", nQ, naoccA, naoccA);
        U->read(psio_, PSIF_DFOCC_AMPS);
        T->add(U);
        U.reset();
        U = std::make_shared<Tensor2d>("Vp (Q|IJ)", nQ, naoccA, naoccA);
        U->read(psio_, PSIF_DFOCC_AMPS);
        T->axpy(U, -2.0);
        U.reset();
        U = std::make_shared<Tensor2d>("Zeta (Q|IJ)", nQ, naoccA, naoccA);
        U->read(psio_, PSIF_DFOCC_AMPS);
        Tau = std::make_shared<Tensor2d>("Tau (Q|IJ)", nQ, naoccA, naoccA);
        Tau->swap_3index_col(U);
        U.reset();
        T->axpy(Tau, -1.0);
        Tau.reset();
        l1newA->contract(true, false, naoccA, navirA, nQ * naoccA, T, bQiaA, 1.0, 1.0);
        T.reset();

        // l_i^a <= \sum_{Q} (Gp_Q - G_Q) b_ia^Q
        SharedTensor1d gQp2 = std::make_shared<Tensor1d>("CCSDL G_Qp - G_Q", nQ);
        gQp2->copy(gQp);
        gQp2->axpy(gQ, -1.0);
        l1newA->gemv(true, bQiaA, gQp2, 1.0, 1.0);
        gQp2.reset();

        // l_i^a <= \sum_{Q,e} G_ei^Q (b_ea^Q - t_ea^Q)
        T = std::make_shared<Tensor2d>("T1 (Q|AB)", nQ, navirA, navirA);
        T->read(psio_, PSIF_DFOCC_AMPS);
        T->scale(-1.0);
        T->add(bQabA);
        U = std::make_shared<Tensor2d>("G (Q|AI)", nQ, navirA, naoccA);
        U->read(psio_, PSIF_DFOCC_AMPS);
        l1newA->contract(true, false, naoccA, navirA, nQ * navirA, U, T, 1.0, 1.0);
        T.reset();
        U.reset();

        if (wfn_type_ == "DF-CCSD(T)") {
            tL1 = std::make_shared<Tensor2d>("(T)L <I|A>", naoccA, navirA);
            tL1->read(psio_, PSIF_DFOCC_AMPS);
            l1newA->axpy(tL1, 1.0);
            tL1.reset();
        }

        // Denom
        for (int i = 0; i < naoccA; ++i) {
            for (int a = 0; a < navirA; ++a) {
                double value = FockA->get(i + nfrzc, i + nfrzc) - FockA->get(a + noccA, a + noccA);
                l1newA->set(i, a, l1newA->get(i, a) / value);
            }
        }
        // l1newA->print();
    }  // if (reference_ == "RESTRICTED")

    // UHF
    else if (reference_ == "UNRESTRICTED") {
        // Alpha Part
        SharedTensor2d J, W, I, K, X, Y, T, Z, L, U, V, Tau;
        // l_I^A  = Ft_IA + \sum_{E} l_I^E Ft_EA - \sum_{M} l_M^A Ft_IM - \sum_{M} G_MI Ft_MA + \sum_{ME} l_M^E W_IEAM + \sum_{me} l_m^e W_IeAm
        //        - \sum_{MN} G_MN W_MINA - \sum_{mn} G_mn W_mInA - 0.5 * \sum_{MN} \sum_{E} L_MN^AE W_IEMN - \sum_{Mn} \sum_{e} L_Mn^Ae W_IeMn
        //        - \sum_{MN} \sum_{E} L_IMNE Z_NEAM - \sum_{mN} \sum_{e} L_ImNe Z_NeAm + \sum_{mn} \sum_{E} L_mInE Z_nEAm
        //        + \sum_{Q,E} (L_IE^Q + Lt_IE^Q + 2V_EI^Q + Z_IE^Q) b_EA^Q + \sum_{Q,M} (V_MI^Q + Vt_MI^Q - 2Vp_MI^Q - Z_IM^Q) b_MA^Q
        //        + \sum_{Q} (Gp_Q - G_Q) b_IA^Q + \sum_{Q,E} G_EI^Q (b_EA^Q - t_EA^Q)                                                              (115)

        // l_I^A += Ft_IA  (1)
        l1newA->copy(FiaA);
        // l_I^A += \sum_{E} l_I^E Ft_EA  (2)
        l1newA->gemm(false, false, l1A, FtabA, 1.0, 1.0);
        // l_I^A -= \sum_{M} l_M^A Ft_IM  (3)
        l1newA->gemm(false, false, FtijA, l1A, -1.0, 1.0);
        // l_I^A -= \sum_{M} G_MI Ft_MA   (4)
        l1newA->gemm(true, false, GijA, FiaA, -1.0, 1.0);

        // l_I^A += \sum_{ME} l_M^E W_IEAM  (5)
        W = std::make_shared<Tensor2d>("WL (ME|JB)", naoccA, navirA, naoccA, navirA);
        W->read(psio_, PSIF_DFOCC_AMPS);
        l1newA->gemv(false, W, l1A, 1.0, 1.0);
        W.reset();

        // l_I^A += \sum_{me} l_m^e W_IeAm  (6)
        W = std::make_shared<Tensor2d>("WL (ME|jb)", naoccA, navirA, naoccB, navirB);
        W->read(psio_, PSIF_DFOCC_AMPS);
        l1newA->gemv(false, W, l1B, 1.0, 1.0);
        W.reset();

        // l_I^A -= \sum_{MN} G_MN W_MINA   (7)
        X = std::make_shared<Tensor2d>("WL (MN|IE)", naoccA, naoccA, naoccA, navirA);
        X->read(psio_, PSIF_DFOCC_AMPS);
        W = std::make_shared<Tensor2d>("X (IA|MN)", naoccA, navirA, naoccA, naoccA);
        W->sort(2413, X, 1.0, 0.0);
        X.reset();
        l1newA->gemv(false, W, GijA, -1.0, 1.0);
        W.reset();

        // l_I^A -= \sum_{mn} G_mn W_mInA   (8)
        X = std::make_shared<Tensor2d>("WL (mN|iE)", naoccB, naoccA, naoccB, navirA);
        X->read(psio_, PSIF_DFOCC_AMPS);
        W = std::make_shared<Tensor2d>("X (IA|mn)", naoccA, navirA, naoccB, naoccB);
        W->sort(2413, X, 1.0, 0.0);
        X.reset();
        l1newA->gemv(false, W, GijB, -1.0, 1.0);
        W.reset();

        // l_I^A -= 0.5 * \sum_{MN} \sum_{E} L_MN^AE W_IEMN   (9)
        Y = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        Y->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        L = std::make_shared<Tensor2d>("L <EM|NA>", navirA, naoccA, naoccA, navirA);
        L->sort(4123, Y, 1.0, 0.0);
        Y.reset();
        W = std::make_shared<Tensor2d>("WL (MB|IJ)", naoccA, navirA, naoccA, naoccA);
        W->read(psio_, PSIF_DFOCC_AMPS);
        l1newA->contract(false, false, naoccA, navirA, navirA * naoccA * naoccA, W, L, -0.5, 1.0);
        W.reset();
        L.reset();

        // l_I^A -= \sum_{Mn} \sum_{e} L_Mn^Ae W_IeMn         (10)
        Y = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        Y->read(psio_, PSIF_DFOCC_AMPS);
        L = std::make_shared<Tensor2d>("L <eM|nA>", navirB, naoccA, naoccB, navirA);
        L->sort(4123, Y, 1.0, 0.0);
        Y.reset();
        W = std::make_shared<Tensor2d>("WL (Mb|Ij)", naoccA, navirB, naoccA, naoccB);
        W->read(psio_, PSIF_DFOCC_AMPS);
        l1newA->contract(false, false, naoccA, navirA, navirB * naoccA * naoccB, W, L, -1.0, 1.0);
        W.reset();
        L.reset();

        // l_I^A -= \sum_{MN} \sum_{E} L_IMNE Z_NEAM          (11)
        X = std::make_shared<Tensor2d>("Z (ME|JB)", naoccA, navirA, naoccA, navirA);
        X->read(psio_, PSIF_DFOCC_AMPS);
        Z = std::make_shared<Tensor2d>("Z (MN|EA)", naoccA, naoccA, navirA, navirA);
        Z->sort(3142, X, 1.0, 0.0);
        X.reset();
        L = std::make_shared<Tensor2d>("L <IJ|KA>", naoccA, naoccA, naoccA, navirA);
        L->read(psio_, PSIF_DFOCC_AMPS);
        l1newA->contract(false, false, naoccA, navirA, naoccA * naoccA * navirA, L, Z, -1.0, 1.0);
        L.reset();
        Z.reset();

        // l_I^A -= \sum_{mN} \sum_{e} L_ImNe Z_NeAm          (12)
        Y = std::make_shared<Tensor2d>("Z (ME|jb)", naoccA, navirA, naoccB, navirB);
        Y->read(psio_, PSIF_DFOCC_AMPS);
        Z = std::make_shared<Tensor2d>("Z (mN|eA)", naoccB, naoccA, navirB, navirA);
        Z->sort(3142, Y, 1.0, 0.0);
        Y.reset();
        L = std::make_shared<Tensor2d>("L <Ij|Ka>", naoccA, naoccB, naoccA, navirB);
        L->read(psio_, PSIF_DFOCC_AMPS);
        l1newA->contract(false, false, naoccA, navirA, naoccB * naoccA * navirB, L, Z, -1.0, 1.0);
        L.reset();
        Z.reset();

        // l_I^A += \sum_{mn} \sum_{E} L_mInE Z_nEAm          (13)
        Y = std::make_shared<Tensor2d>("Z (mE|jB)", naoccB, navirA, naoccB, navirA);
        Y->read(psio_, PSIF_DFOCC_AMPS);
        Z = std::make_shared<Tensor2d>("Z (mn|EA)", naoccB, naoccB, navirA, navirA);
        Z->sort(3142, Y, 1.0, 0.0);
        Y.reset();
        X = std::make_shared<Tensor2d>("L <iJ|kA>", naoccB, naoccA, naoccB, navirA);
        X->read(psio_, PSIF_DFOCC_AMPS);
        L = std::make_shared<Tensor2d>("L <Im|nE>", naoccA, naoccB, naoccB, navirA);
        L->sort(2134, X, 1.0, 0.0);
        X.reset();
        l1newA->contract(false, false, naoccA, navirA, naoccB * naoccB * navirA, L, Z, 1.0, 1.0);
        L.reset();
        Z.reset();

        // l_I^A += \sum_{Q,E} (L_IE^Q + Lt_IE^Q + 2V_EI^Q + Z_IE^Q) b_EA^Q    (14)
        T = std::make_shared<Tensor2d>("Temp (Q|IA)", nQ, naoccA, navirA);
        V = std::make_shared<Tensor2d>("V (Q|AI)", nQ, navirA, naoccA);
        V->read(psio_, PSIF_DFOCC_AMPS);
        T->swap_3index_col(V);
        V.reset();
        T->scale(2.0);
        U = std::make_shared<Tensor2d>("L2 (Q|IA)", nQ, naoccA, navirA);
        U->read(psio_, PSIF_DFOCC_AMPS);
        T->add(U);
        U.reset();
        U = std::make_shared<Tensor2d>("L2t (Q|IA)", nQ, naoccA, navirA);
        U->read(psio_, PSIF_DFOCC_AMPS);
        T->add(U);
        U.reset();
        U = std::make_shared<Tensor2d>("Zeta (Q|IA)", nQ, naoccA, navirA);
        U->read(psio_, PSIF_DFOCC_AMPS);
        T->add(U);
        U.reset();
        X = std::make_shared<Tensor2d>("Temp (Q|AI)", nQ, navirA, naoccA);
        X->swap_3index_col(T);
        T.reset();

        l1newA->contract(true, false, naoccA, navirA, nQ * navirA, X, bQabA, 1.0, 1.0);
        X.reset();

        // l_I^A += \sum_{Q,M} (V_MI^Q + Vt_MI^Q - 2Vp_MI^Q - Z_IM^Q) b_MA^Q   (15)
        T = std::make_shared<Tensor2d>("Temp (Q|IJ)", nQ, naoccA, naoccA);
        U = std::make_shared<Tensor2d>("V (Q|IJ)", nQ, naoccA, naoccA);
        U->read(psio_, PSIF_DFOCC_AMPS);
        T->copy(U);
        U.reset();
        U = std::make_shared<Tensor2d>("Vt (Q|IJ)", nQ, naoccA, naoccA);
        U->read(psio_, PSIF_DFOCC_AMPS);
        T->add(U);
        U.reset();
        U = std::make_shared<Tensor2d>("Vp (Q|IJ)", nQ, naoccA, naoccA);
        U->read(psio_, PSIF_DFOCC_AMPS);
        T->axpy(U, -2.0);
        U.reset();
        U = std::make_shared<Tensor2d>("Zeta (Q|IJ)", nQ, naoccA, naoccA);
        U->read(psio_, PSIF_DFOCC_AMPS);
        Tau = std::make_shared<Tensor2d>("Tau (Q|IJ)", nQ, naoccA, naoccA);
        Tau->swap_3index_col(U);
        U.reset();
        T->axpy(Tau, -1.0);
        Tau.reset();

        l1newA->contract(true, false, naoccA, navirA, nQ * naoccA, T, bQiaA, 1.0, 1.0);
        T.reset();

        // l_I^A += \sum_{Q} (Gp_Q - G_Q) b_IA^Q              (16)
        SharedTensor1d gQp2 = std::make_shared<Tensor1d>("CCSDL G_Qp - G_Q", nQ);
        gQp2->copy(gQp);
        gQp2->axpy(gQ, -1.0);
        l1newA->gemv(true, bQiaA, gQp2, 1.0, 1.0);
        //gQp2.reset();

        // l_I^A +=  \sum_{Q,E} G_EI^Q (b_EA^Q - t_EA^Q)      (17)
        T = std::make_shared<Tensor2d>("T1 (Q|AB)", nQ, navirA, navirA);
        T->read(psio_, PSIF_DFOCC_AMPS);
        T->scale(-1.0);

        T->add(bQabA);
        U = std::make_shared<Tensor2d>("G (Q|AI)", nQ, navirA, naoccA);
        U->read(psio_, PSIF_DFOCC_AMPS);
        l1newA->contract(true, false, naoccA, navirA, nQ * navirA, U, T, 1.0, 1.0);
        T.reset();
        U.reset();

        if (wfn_type_ == "DF-CCSD(T)") {
            SharedTensor2d tL1 = std::make_shared<Tensor2d>("(T)L <I|A>", naoccA, navirA);
            tL1->read(psio_, PSIF_DFOCC_AMPS);
            l1newA->axpy(tL1, 1.0);
            tL1.reset();
        }

        // Denom
        for (int i = 0; i < naoccA; ++i) {
            for (int a = 0; a < navirA; ++a) {
                double value = FockA->get(i + nfrzc, i + nfrzc) - FockA->get(a + noccA, a + noccA);
                l1newA->set(i, a, l1newA->get(i, a) / value);
            }
        }

        // Beta Part
        // l_i^a = Ft_IA + \sum_{e} l_i^e Ft_ea - \sum_{e} l_m^a Ft_im - \sum_{m} G_mi Ft_ma  + \sum_{me} l_m^e W_ieam + \sum_{ME} l_M^E W_iEaM
        //       - \sum_{mn} G_mn W_mina - \sum_{MN} G_MN W_MiNa - 0.5 * \sum_{mn} \sum_{e} L_mn^ae W_iemn - \sum_{Mn} \sum_{E} L_Mn^Ea W_iEnM
        //       - \sum_{mn} \sum_{e} L_imne Z_neam  + \sum_{MN} \sum_{e} L_MiNe Z_NeaM - \sum_{Mn} \sum_{E} L_iMnE Z_nEaM
        //       + \sum_{Q,e} (L_ie^Q + Lt_ie^Q + 2V_ei^Q + Z_ie^Q) b_ea^Q  + \sum_{Q,m} (V_mi^Q + Vt_mi^Q - 2Vp_mi^Q - Z_im^Q) b_ma^Q
        //       + \sum_{Q} (Gp_Q - G_Q) b_ia^Q + \sum_{Q,e} G_ei^Q (b_ea^Q - t_ea^Q)                                                              (116)

        // l_i^a += Ft_ia  (1)
        l1newB->copy(FiaB);

        // l_i^a += \sum_{e} l_i^e Ft_ea (2)
        l1newB->gemm(false, false, l1B, FtabB, 1.0, 1.0);

        // l_i^a -= \sum_{e} l_m^a Ft_im (3)
        l1newB->gemm(false, false, FtijB, l1B, -1.0, 1.0);

        // l_i^a -= \sum_{m} G_mi Ft_ma  (4)
        l1newB->gemm(true, false, GijB, FiaB, -1.0, 1.0);

        // l_i^a += \sum_{me} l_m^e W_ieam (5)
        W = std::make_shared<Tensor2d>("WL (me|jb)", naoccB, navirB, naoccB, navirB);
        W->read(psio_, PSIF_DFOCC_AMPS);
        l1newB->gemv(false, W, l1B, 1.0, 1.0);
        W.reset();

        // l_i^a += \sum_{ME} l_M^E W_iEaM (6)
        W = std::make_shared<Tensor2d>("WL (me|JB)", naoccB, navirB, naoccA, navirA);
        W->read(psio_, PSIF_DFOCC_AMPS);
        l1newB->gemv(false, W, l1A, 1.0, 1.0);
        W.reset();

        // l_i^a -= \sum_{mn} G_mn W_mina  (7)
        X = std::make_shared<Tensor2d>("WL (mn|ie)", naoccB, naoccB, naoccB, navirB);
        X->read(psio_, PSIF_DFOCC_AMPS);
        W = std::make_shared<Tensor2d>("X (ia|mn)", naoccB, navirB, naoccB, naoccB);
        W->sort(2413, X, 1.0, 0.0);
        X.reset();
        l1newB->gemv(false, W, GijB, -1.0, 1.0);
        W.reset();

        // l_i^a -= \sum_{MN} G_MN W_MiNa  (8)
        X = std::make_shared<Tensor2d>("WL (Mn|Ie)", naoccA, naoccB, naoccA, navirB);
        X->read(psio_, PSIF_DFOCC_AMPS);
        W = std::make_shared<Tensor2d>("X (ia|MN)", naoccB, navirB, naoccA, naoccA);
        W->sort(2413, X, 1.0, 0.0);
        X.reset();
        l1newB->gemv(false, W, GijA, -1.0, 1.0);
        W.reset();

        // l_i^a -= 0.5 * \sum_{mn} \sum_{e} L_mn^ae W_iemn  (9)
        Y = std::make_shared<Tensor2d>("L2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        Y->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        L = std::make_shared<Tensor2d>("L <em|na>", navirB, naoccB, naoccB, navirB);
        L->sort(4123, Y, 1.0, 0.0);
        Y.reset();
        W = std::make_shared<Tensor2d>("WL (mb|ij)", naoccB, navirB, naoccB, naoccB);
        W->read(psio_, PSIF_DFOCC_AMPS);
        l1newB->contract(false, false, naoccB, navirB, navirB * naoccB * naoccB, W, L, -0.5, 1.0);
        W.reset();
        L.reset();

        // l_i^a -= \sum_{Mn} \sum_{E} L_Mn^Ea W_iEnM        (10)
        Y = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        Y->read(psio_, PSIF_DFOCC_AMPS);
        L = std::make_shared<Tensor2d>("L (En|Ma)", navirA, naoccB, naoccA, navirB);
        L->sort(3214, Y, 1.0, 0.0);
        Y.reset();
        W = std::make_shared<Tensor2d>("WL (mB|iJ)", naoccB, navirA, naoccB, naoccA);
        W->read(psio_, PSIF_DFOCC_AMPS);
        l1newB->contract(false, false, naoccB, navirB, navirA * naoccB * naoccA, W, L, -1.0, 1.0);
        W.reset();
        L.reset();

        // l_i^a -= \sum_{mn} \sum_{e} L_imne Z_neam         (11)
        X = std::make_shared<Tensor2d>("Z (me|jb)", naoccB, navirB, naoccB, navirB);
        X->read(psio_, PSIF_DFOCC_AMPS);
        Z = std::make_shared<Tensor2d>("Z (mn|ea)", naoccB, naoccB, navirB, navirB);
        Z->sort(3142, X, 1.0, 0.0);
        X.reset();
        L = std::make_shared<Tensor2d>("L <ij|ka>", naoccB, naoccB, naoccB, navirB);
        L->read(psio_, PSIF_DFOCC_AMPS);
        l1newB->contract(false, false, naoccB, navirB, naoccB * naoccB * navirB, L, Z, -1.0, 1.0);
        L.reset();
        Z.reset();

        // l_i^a += \sum_{MN} \sum_{e} L_MiNe Z_NeaM         (12)
        X = std::make_shared<Tensor2d>("Z (Me|Jb)", naoccA, navirB, naoccA, navirB);
        X->read(psio_, PSIF_DFOCC_AMPS);
        Z = std::make_shared<Tensor2d>("Z (MN|ea)", naoccA, naoccA, navirB, navirB);
        Z->sort(3142, X, 1.0, 0.0);
        X.reset();
        Y = std::make_shared<Tensor2d>("L <Ij|Ka>", naoccA, naoccB, naoccA, navirB);
        Y->read(psio_, PSIF_DFOCC_AMPS);
        L = std::make_shared<Tensor2d>("L <iM|Ne>", naoccB, naoccA, naoccA, navirB);
        L->sort(2134, Y, 1.0, 0.0);
        Y.reset();
        l1newB->contract(false, false, naoccB, navirB, naoccA * naoccA * navirB, L, Z, 1.0, 1.0);
        L.reset();
        Z.reset();

        // l_i^a -= \sum_{Mn} \sum_{E} L_iMnE Z_nEaM         (13)
        X = std::make_shared<Tensor2d>("Z (me|JB)", naoccB, navirB, naoccA, navirA);
        X->read(psio_, PSIF_DFOCC_AMPS);
        Z = std::make_shared<Tensor2d>("Z (Mn|Ea)", naoccA, naoccB, navirA, navirB);
        Z->sort(3142, X, 1.0, 0.0);
        X.reset();
        L = std::make_shared<Tensor2d>("L <iJ|kA>", naoccB, naoccA, naoccB, navirA);
        L->read(psio_, PSIF_DFOCC_AMPS);
        l1newB->contract(false, false, naoccB, navirB, naoccA * naoccB * navirA, L, Z, -1.0, 1.0);
        L.reset();
        Z.reset();

        // l_i^a += \sum_{Q,e} (L_ie^Q + Lt_ie^Q + 2V_ei^Q + Z_ie^Q) b_ea^Q    (14)
        T = std::make_shared<Tensor2d>("Temp (Q|ia)", nQ, naoccB, navirB);
        V = std::make_shared<Tensor2d>("V (Q|ai)", nQ, navirB, naoccB);
        V->read(psio_, PSIF_DFOCC_AMPS);
        T->swap_3index_col(V);
        V.reset();
        T->scale(2.0);
        U = std::make_shared<Tensor2d>("L2 (Q|ia)", nQ, naoccB, navirB);
        U->read(psio_, PSIF_DFOCC_AMPS);
        T->add(U);
        U.reset();
        U = std::make_shared<Tensor2d>("L2t (Q|ia)", nQ, naoccB, navirB);
        U->read(psio_, PSIF_DFOCC_AMPS);
        T->add(U);
        U.reset();
        U = std::make_shared<Tensor2d>("Zeta (Q|ia)", nQ, naoccB, navirB);
        U->read(psio_, PSIF_DFOCC_AMPS);
        T->add(U);
        U.reset();
        X = std::make_shared<Tensor2d>("Temp (Q|ai)", nQ, navirB, naoccB);
        X->swap_3index_col(T);
        T.reset();

        l1newB->contract(true, false, naoccB, navirB, nQ * navirB, X, bQabB, 1.0, 1.0);
        X.reset();

        // l_i^a += \sum_{Q,m} (V_mi^Q + Vt_mi^Q - 2Vp_mi^Q - Z_im^Q) b_ma^Q   (15)
        T = std::make_shared<Tensor2d>("Temp (Q|ij)", nQ, naoccB, naoccB);
        U = std::make_shared<Tensor2d>("V (Q|ij)", nQ, naoccB, naoccB);
        U->read(psio_, PSIF_DFOCC_AMPS);
        T->copy(U);
        U.reset();
        U = std::make_shared<Tensor2d>("Vt (Q|ij)", nQ, naoccB, naoccB);
        U->read(psio_, PSIF_DFOCC_AMPS);
        T->add(U);
        U.reset();
        U = std::make_shared<Tensor2d>("Vp (Q|ij)", nQ, naoccB, naoccB);
        U->read(psio_, PSIF_DFOCC_AMPS);
        T->axpy(U, -2.0);
        U.reset();
        U = std::make_shared<Tensor2d>("Zeta (Q|ij)", nQ, naoccB, naoccB);
        U->read(psio_, PSIF_DFOCC_AMPS);
        Tau = std::make_shared<Tensor2d>("Tau (Q|ij)", nQ, naoccB, naoccB);
        Tau->swap_3index_col(U);
        U.reset();
        T->axpy(Tau, -1.0);
        Tau.reset();

        l1newB->contract(true, false, naoccB, navirB, nQ * naoccB, T, bQiaB, 1.0, 1.0);
        T.reset();

        // l_i^a += \sum_{Q} (Gp_Q - G_Q) b_ia^Q                               (16)
        //SharedTensor1d gQp2 = std::make_shared<Tensor1d>("CCSDL G_Qp - G_Q", nQ);
        //gQp2->copy(gQp);
        //gQp2->axpy(gQ, -1.0);
        l1newB->gemv(true, bQiaB, gQp2, 1.0, 1.0);
        gQp2.reset();

        // l_i^a += \sum_{Q,e} G_ei^Q (b_ea^Q - t_ea^Q)                        (17)
        T = std::make_shared<Tensor2d>("T1 (Q|ab)", nQ, navirB, navirB);
        T->read(psio_, PSIF_DFOCC_AMPS);
        T->scale(-1.0);
        T->add(bQabB);
        U = std::make_shared<Tensor2d>("G (Q|ai)", nQ, navirB, naoccB);
        U->read(psio_, PSIF_DFOCC_AMPS);
        l1newB->contract(true, false, naoccB, navirB, nQ * navirB, U, T, 1.0, 1.0);
        T.reset();
        U.reset();

        if (wfn_type_ == "DF-CCSD(T)") {
            SharedTensor2d tL1 = std::make_shared<Tensor2d>("(T)L <i|a>", naoccB, navirB);
            tL1->read(psio_, PSIF_DFOCC_AMPS);
            l1newB->axpy(tL1, 1.0);
            tL1.reset();
        }

        // Denom
        for (int i = 0; i < naoccB; ++i) {
            for (int a = 0; a < navirB; ++a) {
                double value = FockB->get(i + nfrzc, i + nfrzc) - FockB->get(a + noccB, a + noccB);
                l1newB->set(i, a, l1newB->get(i, a) / value);
            }
        }
        //l1newB->print();

    }  // else if (reference_ == "UNRESTRICTED")

}  // end ccsdl_l1_amps
}  // namespace dfoccwave
}  // namespace psi
