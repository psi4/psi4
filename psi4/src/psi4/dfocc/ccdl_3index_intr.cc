/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

    // RHF
    if (reference_ == "RESTRICTED") {
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
    }// if (reference_ == "RESTRICTED")

    // UHF
    else if (reference_ == "UNRESTRICTED") {
        SharedTensor2d K, L, L2, T, T2, U, Tau, V, V2, Vij, Vai, X, Y, Z;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // V_ijkl Intermediates ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // V_IJKL = 0.5 * \sum_{EF} T_IJ^EF L_KL^EF (41)
        Tau = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        Tau->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        // Contraction
        L = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        L->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        V = std::make_shared<Tensor2d>("V <IJ|KL>", naoccA, naoccA, naoccA, naoccA);
        V->gemm(false, true, Tau, L, 0.5, 0.0);
        L.reset();
        Tau.reset();
        V->write(psio_, PSIF_DFOCC_AMPS);
        V.reset();

        // V_ijkl = \sum_{ef} T_ij^ef L_kl^ef (42)
        Tau = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        Tau->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        // Contraction
        L = std::make_shared<Tensor2d>("L2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        L->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        V = std::make_shared<Tensor2d>("V <ij|kl>", naoccB, naoccB, naoccB, naoccB);
        V->gemm(false, true, Tau, L, 0.5, 0.0);
        V->write(psio_, PSIF_DFOCC_AMPS);
        V.reset();
        L.reset();
        Tau.reset();

        // V_IjKl = \sum_{Ef} T_Ij^Ef L_Kl^Ef (43)
        Tau = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        Tau->read(psio_, PSIF_DFOCC_AMPS);
        L = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        L->read(psio_, PSIF_DFOCC_AMPS);
        V = std::make_shared<Tensor2d>("V <Ij|Kl>", naoccA, naoccB, naoccA, naoccB);
        V->gemm(false, true, Tau, L, 1.0, 0.0);
        V->write(psio_, PSIF_DFOCC_AMPS);

        /*
        SharedTensor2d C = std::make_shared<Tensor2d>("V <IJ|KL>", naoccA, naoccA, naoccA, naoccA);
        C->read(psio_, PSIF_DFOCC_AMPS);
        C->axpy(V, 1.0);
        C->print();
        C.reset();
        */

        V.reset();
        L.reset();
        Tau.reset();

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // V(Q,ij) Intermediates /////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Alpha Block
        // V_IJ^Q <= \sum_{MN} V_IMJN b_MN^Q = B(Q,MN) Y(MN,IJ) (63)
        // Y_MNIJ = V_IMJN
        V = std::make_shared<Tensor2d>("V <IJ|KL>", naoccA, naoccA, naoccA, naoccA);
        V->read(psio_, PSIF_DFOCC_AMPS);
        Y = std::make_shared<Tensor2d>("Y <MN|IJ>", naoccA, naoccA, naoccA, naoccA);
        Y->sort(2413, V, 1.0, 0.0);
        V.reset();

        T = std::make_shared<Tensor2d>("V (Q|IJ)", nQ, naoccA, naoccA);
        T->gemm(false, false, bQijA, Y, 1.0, 0.0);
        Y.reset();
        // V_IJ^Q += \sum_{mn} V_ImJn b_mn^Q = B(Q,mn) Y(mn,IJ)
        // Y_mnIJ = V_ImJn
        X = std::make_shared<Tensor2d>("V <Ij|Kl>", naoccA, naoccB, naoccA, naoccB);
        X->read(psio_, PSIF_DFOCC_AMPS);
        Y = std::make_shared<Tensor2d>("Y <mn|IJ>", naoccB, naoccB, naoccA, naoccA);
        Y->sort(2413, X, 1.0, 0.0);
        X.reset();

        T->gemm(false, false, bQijB, Y, 1.0, 1.0);
        Y.reset();
        T->write(psio_, PSIF_DFOCC_AMPS);
        T.reset();

        // Beta Block
        // V_ij^Q <= \sum_{mn} V_imjn b_mn^Q = B(Q,mn) Y(mn,ij) (64)
        // Y_mnij = V_imjn
        V = std::make_shared<Tensor2d>("V <ij|kl>", naoccB, naoccB, naoccB, naoccB);
        V->read(psio_, PSIF_DFOCC_AMPS);
        Y = std::make_shared<Tensor2d>("Y <mn|ij>", naoccB, naoccB, naoccB, naoccB);
        Y->sort(2413, V, 1.0, 0.0);
        V.reset();

        T = std::make_shared<Tensor2d>("V (Q|ij)", nQ, naoccB, naoccB);
        T->gemm(false, false, bQijB, Y, 1.0, 0.0);
        Y.reset();
        // V_ij^Q += \sum_{MN} V_MiNj b_MN^Q = B(Q,MN) Y(MN,ij)
        // Y_MNij = V_MiNj
        X = std::make_shared<Tensor2d>("V <Ij|Kl>", naoccA, naoccB, naoccA, naoccB);
        X->read(psio_, PSIF_DFOCC_AMPS);
        Y = std::make_shared<Tensor2d>("Y <MN|ij>", naoccA, naoccA, naoccB, naoccB);
        Y->sort(1324, X, 1.0, 0.0);
        X.reset();

        T->gemm(false, false, bQijA, Y, 1.0, 1.0);
        T->write(psio_, PSIF_DFOCC_AMPS);
        T.reset();
        Y.reset();

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // V_iajb Intermediates ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // V_IAJB = 1/2 \sum_{ME} T(IM,BE) L(JM,AE) + 1/2 \sum_{me} T(Im,Be) L(Jm,Ae) (47)
        // V_IAJB = 1/2 \sum_{ME} T(IM,BE) L(JM,AE)
        T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        T = std::make_shared<Tensor2d>("T2p (IB|ME)", naoccA, navirA, naoccA, navirA);
        T->sort(1324, T2, 1.0, 0.0);
        T2.reset();
        L2 = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        L2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        L = std::make_shared<Tensor2d>("L2p (ME|JA)", naoccA, navirA, naoccA, navirA);
        L->sort(2413, L2, 1.0, 0.0);
        L2.reset();
        X = std::make_shared<Tensor2d>("X (IB|JA)", naoccA, navirA, naoccA, navirA);
        X->gemm(false, false, T, L, 0.5, 0.0);
        T.reset();
        L.reset();
        // V_IAJB +=  1/2 \sum_{me} T(Im,Be) L(Jm,Ae)
        T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T2->read(psio_, PSIF_DFOCC_AMPS);
        T = std::make_shared<Tensor2d>("T2p (IB|me)", naoccA, navirA, naoccB, navirB);
        T->sort(1324, T2, 1.0, 0.0);
        T2.reset();
        L2 = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        L2->read(psio_, PSIF_DFOCC_AMPS);
        L = std::make_shared<Tensor2d>("L2p (me|JA)", naoccB, navirB, naoccA, navirA);
        L->sort(2413, L2, 1.0, 0.0);
        L2.reset();
        X->gemm(false, false, T, L, 0.5, 1.0);
        T.reset();
        L.reset();
        V = std::make_shared<Tensor2d>("V (IA|JB)", naoccA, navirA, naoccA, navirA);
        V->sort(1432, X, 1.0, 0.0);
        X.reset();
        V->write(psio_, PSIF_DFOCC_AMPS);
        V.reset();

        // V_iajb = 1/2 \sum_{me} T(im,be) L(jm,ae) + 1/2 \sum_{ME} T(Mi,Eb) L(Mj,Ea) (48)
        // V_iajb = 1/2 \sum_{me} T(im,be) L(jm,ae)
        T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        T = std::make_shared<Tensor2d>("T2 (ib|me)", naoccB, navirB, naoccB, navirB);
        T->sort(1324, T2, 1.0, 0.0);
        T2.reset();
        L2 = std::make_shared<Tensor2d>("L2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        L2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        L = std::make_shared<Tensor2d>("L2p (me|ja)", naoccB, navirB, naoccB, navirB);
        L->sort(2413, L2, 1.0, 0.0);
        L2.reset();
        X = std::make_shared<Tensor2d>("X (ib|ja)", naoccB, navirB, naoccB, navirB);
        X->gemm(false, false, T, L, 0.5, 0.0);
        T.reset();
        L.reset();
        // V_iajb += 1/2 \sum_{ME} T(Mi,Eb) L(Mj,Ea)
        T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T2->read(psio_, PSIF_DFOCC_AMPS);
        T = std::make_shared<Tensor2d>("T2p (ib|ME)", naoccB, navirB, naoccA, navirA);
        T->sort(2413, T2, 1.0, 0.0);
        T2.reset();
        L2 = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        L2->read(psio_, PSIF_DFOCC_AMPS);
        L = std::make_shared<Tensor2d>("L2p (ME|ja)", naoccA, navirA, naoccB, navirB);
        L->sort(1324, L2, 1.0, 0.0);
        L2.reset();
        X->gemm(false, false, T, L, 0.5, 1.0);
        T.reset();
        L.reset();
        V = std::make_shared<Tensor2d>("V (ia|jb)", naoccB, navirB, naoccB, navirB);
        V->sort(1432, X, 1.0, 0.0);
        X.reset();
        V->write(psio_, PSIF_DFOCC_AMPS);
        V.reset();

        // V_IaJb = 1/2 \sum_{mE} T(Im,Eb) L(Jm,Ea) (49)
        T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T2->read(psio_, PSIF_DFOCC_AMPS);
        T = std::make_shared<Tensor2d>("T2p (Ib|mE)", naoccA, navirB, naoccB, navirA);
        T->sort(1423, T2, 1.0, 0.0);
        T2.reset();
        L2 = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        L2->read(psio_, PSIF_DFOCC_AMPS);
        L = std::make_shared<Tensor2d>("L2p (mE|Ja)", naoccB, navirA, naoccA, navirB);
        L->sort(2314, L2, 1.0, 0.0);
        L2.reset();
        X = std::make_shared<Tensor2d>("X (Ib|Ja)", naoccA, navirB, naoccA, navirB);
        X->gemm(false, false, T, L, 0.5, 0.0);
        T.reset();
        L.reset();
        V = std::make_shared<Tensor2d>("V (Ia|Jb)", naoccA, navirB, naoccA, navirB);
        V->sort(1432, X, 1.0, 0.0);
        X.reset();
        V->write(psio_, PSIF_DFOCC_AMPS);
        V.reset();

        // V_iAjB = 1/2 \sum_{Me} T(Mi,Be) L(Mj,Ae) (50)
        T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T2->read(psio_, PSIF_DFOCC_AMPS);
        T = std::make_shared<Tensor2d>("T2p (iB|Me)", naoccB, navirA, naoccA, navirB);
        T->sort(2314, T2, 1.0, 0.0);
        T2.reset();
        L2 = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        L2->read(psio_, PSIF_DFOCC_AMPS);
        L = std::make_shared<Tensor2d>("L2p (Me|jA)", naoccA, navirB, naoccB, navirA);
        L->sort(1423, L2, 1.0, 0.0);
        L2.reset();
        X = std::make_shared<Tensor2d>("X (iB|jA)", naoccB, navirA, naoccB, navirA);
        X->gemm(false, false, T, L, 0.5, 0.0);
        T.reset();
        L.reset();
        V = std::make_shared<Tensor2d>("V (iA|jB)", naoccB, navirA, naoccB, navirA);
        V->sort(1432, X, 1.0, 0.0);
        X.reset();
        V->write(psio_, PSIF_DFOCC_AMPS);
        V.reset();

        // V_IajB = 1/2 \sum_{ME} T(IM,BE) L(Mj,Ea) + 1/2 \sum_{me} T(Im,Be) L(jm,ae) (51)
        // V_IajB = 1/2 \sum_{ME} T(IM,BE) L(Mj,Ea)
        T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        T = std::make_shared<Tensor2d>("T2p (IA|JB)", naoccA, navirA, naoccA, navirA);
        T->sort(1324, T2, 1.0, 0.0);
        T2.reset();
        L2 = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        L2->read(psio_, PSIF_DFOCC_AMPS);
        L = std::make_shared<Tensor2d>("L2p (IA|jb)", naoccA, navirA, naoccB, navirB);
        L->sort(1324, L2, 1.0, 0.0);
        L2.reset();
        X = std::make_shared<Tensor2d>("X (IB|ja)", naoccA, navirA, naoccB, navirB);
        X->gemm(false, false, T, L, 0.5, 0.0);
        T.reset();
        L.reset();
        // V_IajB += 1/2 \sum_{me} T(IB,me) L(me,ja)
        T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T2->read(psio_, PSIF_DFOCC_AMPS);
        T = std::make_shared<Tensor2d>("T2p (IA|jb)", naoccA, navirA, naoccB, navirB);
        T->sort(1324, T2, 1.0, 0.0);
        T2.reset();
        L2 = std::make_shared<Tensor2d>("L2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        L2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        L = std::make_shared<Tensor2d>("L2p (ia|jb)", naoccB, navirB, naoccB, navirB);
        L->sort(2413, L2, 1.0, 0.0);
        L2.reset();
        X->gemm(false, false, T, L, 0.5, 1.0);
        T.reset();
        L.reset();
        V = std::make_shared<Tensor2d>("V (Ia|jB)", naoccA, navirB, naoccB, navirA);
        V->sort(1432, X, 1.0, 0.0);
        X.reset();
        V->write(psio_, PSIF_DFOCC_AMPS);
        V.reset();

        // V_iAJb = 1/2 \sum_{me} T(im,be) L(Jm,Ae) + 1/2 \sum_{ME} T(Mi,Eb) L(JM,AE) (52)
        // V_iAJb = 1/2 \sum_{me} T(im,be) L(Jm,Ae)
        T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        T = std::make_shared<Tensor2d>("T2 (ia|jb)", naoccB, navirB, naoccB, navirB);
        T->sort(1324, T2, 1.0, 0.0);
        T2.reset();
        U = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        U->read(psio_, PSIF_DFOCC_AMPS);
        L = std::make_shared<Tensor2d>("L2 (me|JA)", naoccB, navirB, naoccA, navirA);
        L->sort(2413, U, 1.0, 0.0);
        U.reset();
        X = std::make_shared<Tensor2d>("X (ib|JA)", naoccB, navirB, naoccA, navirA);
        X->gemm(false, false, T, L, 0.5, 0.0);
        T.reset();
        L.reset();
        // V_iAJb += 1/2 \sum_{ME} T(Mi,Eb) L(JM,AE)
        T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T2->read(psio_, PSIF_DFOCC_AMPS);
        T = std::make_shared<Tensor2d>("T2p (ib|ME)", naoccB, navirB, naoccA, navirA);
        T->sort(2413, T2, 1.0, 0.0);
        T2.reset();
        L2 = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        L2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        L = std::make_shared<Tensor2d>("L2p (IA|JB)", naoccA, navirA, naoccA, navirA);
        L->sort(2413, L2, 1.0, 0.0);
        L2.reset();
        X->gemm(false, false, T, L, 0.5, 1.0);
        T.reset();
        L.reset();
        V = std::make_shared<Tensor2d>("V (iA|Jb)", naoccB, navirA, naoccA, navirB);
        V->sort(1432, X, 1.0, 0.0);
        X.reset();
        V->write(psio_, PSIF_DFOCC_AMPS);
        V.reset();

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // V(Qp,ij) Intermediates ////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Alpha Block
        // V_IJ^Q' <= \sum_{EF} V_IEJF b_EF^Q  (67)
        X = std::make_shared<Tensor2d>("V (IA|JB)", naoccA, navirA, naoccA, navirA);
        X->read(psio_, PSIF_DFOCC_AMPS);
        Y = std::make_shared<Tensor2d>("Y (AB|IJ)", navirA, navirA, naoccA, naoccA);
        Y->sort(2413, X, 1.0, 0.0);
        X.reset();

        V = std::make_shared<Tensor2d>("Vp (Q|IJ)", nQ, naoccA, naoccA);
        V->gemm(false, false, bQabA, Y, 1.0, 0.0);
        Y.reset();
        // V_IJ^Q' += \sum_{ef} V_IeJf b_ef^Q
        X = std::make_shared<Tensor2d>("V (Ia|Jb)", naoccA, navirB, naoccA, navirB);
        X->read(psio_, PSIF_DFOCC_AMPS);
        Y = std::make_shared<Tensor2d>("Y (ab|IJ)", navirB, navirB, naoccA, naoccA);
        Y->sort(2413, X, 1.0, 0.0);
        X.reset();

        V->gemm(false, false, bQabB, Y, 1.0, 1.0);
        V->write(psio_, PSIF_DFOCC_AMPS);
        V.reset();
        Y.reset();

        // Beta Block
        // V_ij^Q' <= \sum_{ef} V_iejf b_ef^Q  (68)
        X = std::make_shared<Tensor2d>("V (ia|jb)", naoccB, navirB, naoccB, navirB);
        X->read(psio_, PSIF_DFOCC_AMPS);
        Y = std::make_shared<Tensor2d>("Y (ab|ij)", navirB, navirB, naoccB, naoccB);
        Y->sort(2413, X, 1.0, 0.0);
        X.reset();

        V = std::make_shared<Tensor2d>("Vp (Q|ij)", nQ, naoccB, naoccB);
        V->gemm(false, false, bQabB, Y, 1.0, 0.0);
        Y.reset();
        // V_ij^Q' += \sum_{EF} V_iEjF b_EF^Q
        X = std::make_shared<Tensor2d>("V (iA|jB)", naoccB, navirA, naoccB, navirA);
        X->read(psio_, PSIF_DFOCC_AMPS);
        Y = std::make_shared<Tensor2d>("Y (AB|ij)", navirA, navirA, naoccB, naoccB);
        Y->sort(2413, X, 1.0, 0.0);
        X.reset();

        V->gemm(false, false, bQabA, Y, 1.0, 1.0);
        V->write(psio_, PSIF_DFOCC_AMPS);
        V.reset();
        Y.reset();


        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // V(Q,ai) Intermediates /////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Alpha Block
        // V_AI^Q = \sum_{ME} V_MAIE b_ME^Q   (69)
        X = std::make_shared<Tensor2d>("V (IA|JB)", naoccA, navirA, naoccA, navirA);
        X->read(psio_, PSIF_DFOCC_AMPS);
        Y = std::make_shared<Tensor2d>("X (IB|AJ)", naoccA, navirA, navirA, naoccA);
        Y->sort(1423, X, 1.0, 0.0);
        X.reset();

        V = std::make_shared<Tensor2d>("V (Q|AI)", nQ, navirA, naoccA);
        V->gemm(false, false, bQiaA, Y, 1.0, 0.0);
        Y.reset();
        // V_AI^Q += \sum_{me} V_mAIe b_me^Q
        X = std::make_shared<Tensor2d>("V (iA|Jb)", naoccB, navirA, naoccA, navirB);
        X->read(psio_, PSIF_DFOCC_AMPS);
        Y = std::make_shared<Tensor2d>("X (ib|AJ)", naoccB, navirB, navirA, naoccA);
        Y->sort(1423, X, 1.0, 0.0);
        X.reset();

        V->gemm(false, false, bQiaB, Y, 1.0, 1.0);
        V->write(psio_, PSIF_DFOCC_AMPS);
        V.reset();
        Y.reset();

        // Beta Block
        // V_ai^Q = \sum_{me} V_maie b_me^Q  (70)
        X = std::make_shared<Tensor2d>("V (ia|jb)", naoccB, navirB, naoccB, navirB);
        X->read(psio_, PSIF_DFOCC_AMPS);
        Y = std::make_shared<Tensor2d>("V (ib|aj)", naoccB, navirB, navirB, naoccB);
        Y->sort(1423, X, 1.0, 0.0);
        X.reset();
        V = std::make_shared<Tensor2d>("V (Q|ai)", nQ, navirB, naoccB);

        V->gemm(false, false, bQiaB, Y, 1.0, 0.0);
        Y.reset();
        // V_ai^Q += \sum_{ME} V_MaiE b_ME^Q
        X = std::make_shared<Tensor2d>("V (Ia|jB)", naoccA, navirB, naoccB, navirA);
        X->read(psio_, PSIF_DFOCC_AMPS);
        Y = std::make_shared<Tensor2d>("V (IB|aj)", naoccA, navirA, navirB, naoccB);
        Y->sort(1423, X, 1.0, 0.0);
        X.reset();

        V->gemm(false, false, bQiaA, Y, 1.0, 1.0);
        V->write(psio_, PSIF_DFOCC_AMPS);
        V.reset();
        Y.reset();

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // G Intermediates /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // ALPHA BLOCK
        // G_MI <= 1/2 \sum_{N,E,F} T(MN,EF) L(IN,EF) (59)
        T = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        T->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        L = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        L->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        GijA->contract(false, true, naoccA, naoccA, naoccA * navirA * navirA, T, L, 0.5, 0.0);

        // G_AE <= -1/2 \sum_{MNF} L(MN,FA) T(MN,FE) (61)
        GabA->contract(true, false, navirA, navirA, naoccA * naoccA * navirA, L, T, -0.5, 0.0);
        T.reset();
        L.reset();

        // G_MI += \sum_{n,E,f} T(Mn,Ef) L(In,Ef) (59)
        T = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T->read(psio_, PSIF_DFOCC_AMPS);
        L = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        L->read(psio_, PSIF_DFOCC_AMPS);
        GijA->contract(false, true, naoccA, naoccA, naoccB * navirA * navirB, T, L, 1.0, 1.0);

        // G_AE += - \sum_{Mnf} L(Mn,Af) T(Mn,Ef) (61)
        X = std::make_shared<Tensor2d>("X <Ij|bA>", naoccA, naoccB, navirB, navirA);
        X->sort(1243, L, 1.0, 0.0);
        L.reset();
        Y = std::make_shared<Tensor2d>("Y <Ij|bA>", naoccA, naoccB, navirB, navirA);
        Y->sort(1243, T, 1.0, 0.0);
        T.reset();
        GabA->contract(true, false, navirA, navirA, naoccA * naoccB * navirB, X, Y, -1.0, 1.0);
        X.reset();
        Y.reset();

        // BETA BLOCK
        // G_mi <= 1/2 \sum_{nef} T(mn,ef) L(in,ef) (60)
        T = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        T->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        L = std::make_shared<Tensor2d>("L2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        L->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        GijB->contract(false, true, naoccB, naoccB, naoccB * navirB * navirB, T, L, 0.5, 0.0);

        // G_ae <= -1/2 \sum_mnf} L(mn,fa) T(mn,fe) (62)
        GabB->contract(true, false, navirB, navirB, naoccB * naoccB * navirB, L, T, -0.5, 0.0);
        T.reset();
        L.reset();

        // G_ae += - \sum_{mNF} L(Nm,Fa) T(Nm,Fe) (62)
        T = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T->read(psio_, PSIF_DFOCC_AMPS);
        L = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        L->read(psio_, PSIF_DFOCC_AMPS);
        GabB->contract(true, false, navirB, navirB, naoccB * naoccA * navirA, L, T, -1.0, 1.0);

        // G_mi += \sum_{N,e,F} T(Nm,Fe) L(Ni,Fe) (60)
        X = std::make_shared<Tensor2d>("X <jI|Ab>", naoccB, naoccA, navirA, navirB);
        X->sort(2134, T, 1.0, 0.0);
        T.reset();
        Y = std::make_shared<Tensor2d>("Y <jI|Ab>", naoccB, naoccA, navirA, navirB);
        Y->sort(2134, L, 1.0, 0.0);
        L.reset();
        GijB->contract(false, true, naoccB, naoccB, naoccA * navirB * navirA, X, Y, 1.0, 1.0);
        X.reset();
        Y.reset();

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // G(Q) Intermediates/////////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // G_Q = \sum_{EF} G_EF b_EF^Q + \sum_{ef} G_ef b_ef^Q  (75)
        gQ->gemv(false, bQabA, GabA, 1.0, 0.0);
        gQ->gemv(false, bQabB, GabB, 1.0, 1.0);

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // G(Q,ia) Intermediates//////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // G(Q,IA) = \sum_{M} G_MI b_MA^Q  (78)
        T = std::make_shared<Tensor2d>("G (Q|IA)", nQ, naoccA, navirA);
        T->contract233(true, false, naoccA, navirA, GijA, bQiaA, 1.0, 0.0);
        T->write(psio_, PSIF_DFOCC_AMPS);
        T.reset();

        // G(Q,ia) = \sum_{m} G_mi b_ma^Q  (79)
        T = std::make_shared<Tensor2d>("G (Q|ia)", nQ, naoccB, navirB);
        T->contract233(true, false, naoccB, navirB, GijB, bQiaB, 1.0, 0.0);
        T->write(psio_, PSIF_DFOCC_AMPS);
        T.reset();

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // G(Q,ai) Intermediates//////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // G(Q,AI) = \sum_{E} G_AE b_IE^Q  (76)
        X = std::make_shared<Tensor2d>("G (Q|AI)", nQ, navirA, naoccA);
        X->contract233(false, true, navirA, naoccA, GabA, bQiaA, 1.0, 0.0);
        X->write(psio_, PSIF_DFOCC_AMPS);

        // G(Q,ai) = \sum_{e} G_ae b_ie^Q  (77)
        Y = std::make_shared<Tensor2d>("G (Q|ai)", nQ, navirB, naoccB);
        Y->contract233(false, true, navirB, naoccB, GabB, bQiaB, 1.0, 0.0);
        Y->write(psio_, PSIF_DFOCC_AMPS);

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // L(Q,ia) Intermediates//////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Alpha Block
        // L(Q,IA) = \sum_{JB} b_JB^Q L_IJ^AB + \sum_{jb} b_jb^Q L_Ij^Ab   (85)
        // L(Q,IA) = \sum_{JB} b_JB^Q L_IJ^AB
        L2 = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        L2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        U = std::make_shared<Tensor2d>("L2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        U->sort(1324, L2, 1.0, 0.0);
        L2.reset();

        T = std::make_shared<Tensor2d>("L2 (Q|IA)", nQ, naoccA, navirA);
        T->gemm(false, true, bQiaA, U, 1.0, 0.0);
        U.reset();
        // L(Q,IA) += \sum_{jb} b_jb^Q L_Ij^Ab
        L2 = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        L2->read(psio_, PSIF_DFOCC_AMPS);
        U = std::make_shared<Tensor2d>("L2 (IA|jb)", naoccA, navirA, naoccB, navirB);
        U->sort(1324, L2, 1.0, 0.0);
        L2.reset();

        T->gemm(false, true, bQiaB, U, 1.0, 1.0);
        U.reset();
        T->write(psio_, PSIF_DFOCC_AMPS);
        T.reset();

        // Beta Block
        // L(Q,ia) = \sum_{jb} b_jb^Q L_ij^ab + \sum_{JB} b_JB^Q L_JiBa    (86)
        // L(Q,ia) = \sum_{jb} b_jb^Q L_ij^ab
        L2 = std::make_shared<Tensor2d>("L2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        L2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        U = std::make_shared<Tensor2d>("L2 (ia|jb)", naoccB, navirB, naoccB, navirB);
        U->sort(1324, L2, 1.0, 0.0);
        L2.reset();
        T = std::make_shared<Tensor2d>("L2 (Q|ia)", nQ, naoccB, navirB);
        T->gemm(false, true, bQiaB, U, 1.0, 0.0);
        U.reset();
        // L(Q,ia) += \sum_{JB} b_JB^Q L_JiBa
        L2 = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        L2->read(psio_, PSIF_DFOCC_AMPS);
        U = std::make_shared<Tensor2d>("L2 (IA|jb)", naoccA, navirA, naoccB, navirB);
        U->sort(1324, L2, 1.0, 0.0);
        L2.reset();

        T->gemm(false, false, bQiaA, U, 1.0, 1.0);
        U.reset();
        T->write(psio_, PSIF_DFOCC_AMPS);
    }// else if (reference_ == "UNRESTRICTED")

    // outfile->Printf("\t3indices done.\n");

}  // end ccdl_3index_intr
}  // namespace dfoccwave
}  // namespace psi
