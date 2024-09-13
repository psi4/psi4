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

void DFOCC::ccd_pdm_3index_intr() {

    // RHF
    if (reference_ == "RESTRICTED") {
        // defs
        SharedTensor2d K, L, T, U, Tau, V, V2, Vij, Vai, Vab, X, Y, Z;
        SharedTensor2d Vijka, Vijak;
        SharedTensor2d Vs, Va, Ts, Ta, S, A;

        // L(Q,ia) = \sum_{jb} b_jb^Q Ut_ij^ab = \sum_{jb} b(Q,jb) Ut(jb,ia)
        U = std::make_shared<Tensor2d>("Ut2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        ccsd_u2_amps(U, l2);
        U->write_symm(psio_, PSIF_DFOCC_AMPS);
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
        GijA->contract(false, true, naoccA, naoccA, naoccA * navirA * navirA, T, L, 1.0, 0.0);

        // G_ae = -\sum_{m,n,f} U_mn^ef L_mn^af = L(mn,fa) U(mn,fe)
        GabA->contract(true, false, navirA, navirA, naoccA * naoccA * navirA, L, T, -1.0, 0.0);
        T.reset();
        L.reset();

        // G_Q = 2\sum_{ef} G_ef b_ef^Q
        gQ->gemv(false, bQabA, GabA, 2.0, 0.0);

        // Gt_Q = 2\sum_{mn} G_mn b_mn^Q
        gQt->gemv(false, bQijA, GijA, 2.0, 0.0);

        // G(Q,ij) = \sum_{m} G_im b_mj^Q
        T = std::make_shared<Tensor2d>("G (Q|IJ)", nQ, naoccA, naoccA);
        T->contract233(false, false, naoccA, naoccA, GijA, bQijA, 1.0, 0.0);
        // T->cont233("IJ", "IM", "MJ", GijA, bQijA, 1.0, 0.0); // it works
        T->write(psio_, PSIF_DFOCC_AMPS);
        T.reset();

        // G(Q,ai) = \sum_{e} G_ae b_ie^Q
        T = std::make_shared<Tensor2d>("G (Q|AI)", nQ, navirA, naoccA);
        T->contract233(false, true, navirA, naoccA, GabA, bQiaA, 1.0, 0.0);
        T->write(psio_, PSIF_DFOCC_AMPS);
        T.reset();

        // G(Q,ia) = \sum_{m} G_mi b_ma^Q
        T = std::make_shared<Tensor2d>("G (Q|IA)", nQ, naoccA, navirA);
        T->contract233(true, false, naoccA, navirA, GijA, bQiaA, 1.0, 0.0);
        T->write(psio_, PSIF_DFOCC_AMPS);
        T.reset();

        // Build V_ijkl
        // V_ijkl = \sum_{ef} T_ij^ef L_kl^ef
        t2 = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        t2->read_symm(psio_, PSIF_DFOCC_AMPS);
        U = std::make_shared<Tensor2d>("T <IJ|AB>", naoccA, naoccA, navirA, navirA);
        U->sort(1324, t2, 1.0, 0.0);
        // t2.reset();
        L = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        L->sort(1324, l2, 1.0, 0.0);
        V = std::make_shared<Tensor2d>("V <IJ|KL>", naoccA, naoccA, naoccA, naoccA);
        V->gemm(false, true, U, L, 1.0, 0.0);
        L.reset();

        // Y_ijab = \sum(mn) T_mn^ab V_ijmn
        Y = std::make_shared<Tensor2d>("Y <IJ|AB>", naoccA, naoccA, navirA, navirA);
        Y->gemm(false, false, V, U, 1.0, 0.0);
        U.reset();
        Y->write(psio_, PSIF_DFOCC_AMPS);
        Y.reset();

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
        Y.reset();
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
        V->write(psio_, PSIF_DFOCC_AMPS);

        // V_ij^Q' <= 2\sum_{ef} V_iejf b_ef^Q
        Vij = std::make_shared<Tensor2d>("Vp (Q|IJ)", nQ, naoccA, naoccA);
        X = std::make_shared<Tensor2d>("X (AB|IJ)", navirA, navirA, naoccA, naoccA);
        X->sort(2413, V, 1.0, 0.0);
        Vij->gemm(false, false, bQabA, X, 2.0, 0.0);
        X.reset();

        // V_ai^Q <= \sum_{me} V_maie b_me^Q
        Vai = std::make_shared<Tensor2d>("V (Q|AI)", nQ, navirA, naoccA);
        X = std::make_shared<Tensor2d>("X (IB|AJ)", naoccA, navirA, navirA, naoccA);
        X->sort(1423, V, 1.0, 0.0);
        Vai->gemm(false, false, bQiaA, X, 1.0, 0.0);
        X.reset();

        // V_ab^Q = 2\sum_{mn} V_manb b_mn^Q
        X = std::make_shared<Tensor2d>("X (MN|AB)", naoccA, naoccA, navirA, navirA);
        X->sort(1324, V, 1.0, 0.0);
        V.reset();
        Vab = std::make_shared<Tensor2d>("V (Q|AB)", nQ, navirA, navirA);
        Vab->gemm(false, false, bQijA, X, 2.0, 0.0);
        X.reset();
        Vab->write(psio_, PSIF_DFOCC_AMPS);
        Vab.reset();

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
        V->write(psio_, PSIF_DFOCC_AMPS);
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
        Vai->gemm(false, false, bQiaA, X, -2.0, 1.0);
        X.reset();
        Vai->write(psio_, PSIF_DFOCC_AMPS);
        Vai.reset();

        // V_ab^Q -= \sum_{mn} V_mabn b_mn^Q
        X = std::make_shared<Tensor2d>("X (MN|AB)", naoccA, naoccA, navirA, navirA);
        X->sort(1423, V, 1.0, 0.0);
        V.reset();
        Vab = std::make_shared<Tensor2d>("V (Q|AB)", nQ, navirA, navirA);
        Vab->read(psio_, PSIF_DFOCC_AMPS);
        Vab->gemm(false, false, bQijA, X, -1.0, 1.0);
        X.reset();
        Vab->write(psio_, PSIF_DFOCC_AMPS);
    }// end if (reference_ == "RESTRICTED")

    // UHF
    else if (reference_ == "UNRESTRICTED") {
        SharedTensor2d T, T2, L2, Tau, X, Y, Z, V, Vt, U, L, Lt, Vab, Vtab, Vij;

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Gtilde(Q) Intermediate ////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Gtilde(Q) = \sum{M,N} G(M,N) b(Q,MN) + \sum{m,n} G(m,n) b(Q,mn)       (81)
        gQt = std::make_shared<Tensor1d>("CCSD PDM G_Qt", nQ);
        gQt->gemv(false, bQijA, GijA, 1.0, 0.0);
        gQt->gemv(false, bQijB, GijB, 1.0, 1.0);

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // G(Q,IJ) Intermediates /////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // G(Q,IJ) = \sum{M} G(I,M) b(Q,JM)        (83)
        T = std::make_shared<Tensor2d>("G (Q|IJ)", nQ, naoccA, naoccA);
        T->contract233(false, false, naoccA, naoccA, GijA, bQijA, 1.0, 0.0);
        T->write(psio_, PSIF_DFOCC_AMPS);
        T.reset();
        // G(Q,ij) = \sum{m} G(i,m) b(Q,jm)        (84)
        T = std::make_shared<Tensor2d>("G (Q|ij)", nQ, naoccB, naoccB);
        T->contract233(false, false, naoccB, naoccB, GijB, bQijB, 1.0, 0.0);
        T->write(psio_, PSIF_DFOCC_AMPS);
        T.reset();
        //std::cout << "I am here \n";

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // G(Q,AI) Intermediates ////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // G(Q,AI) = \sum{E} G(A,E) b(Q,IE)        (89)
        T = std::make_shared<Tensor2d>("Gt (Q|AI)", nQ, navirA, naoccA);
        T->contract233(false, true, navirA, naoccA, GabA, bQiaA, 1.0, 0.0);
        T->write(psio_, PSIF_DFOCC_AMPS);
        T.reset();
        // G(Q,ai) = \sum{e} G(a,e) b(Q,ie)        (90)
        T = std::make_shared<Tensor2d>("Gt (Q|ai)", nQ, navirB, naoccB);
        T->contract233(false, true, navirB, naoccB, GabB, bQiaB, 1.0, 0.0);
        T->write(psio_, PSIF_DFOCC_AMPS);
        T.reset();

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // V(Q,AB) and Vtilde(Q,AB) Intermediates ////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Alpha Block

        // V(Q,AB) = \sum{M,N} V(MA,NB) b(Q,MN)
        V = std::make_shared<Tensor2d>("V (IA|JB)", naoccA, navirA, naoccA, navirA);
        V->read(psio_, PSIF_DFOCC_AMPS);
        X = std::make_shared<Tensor2d>("X (MN|AB)", naoccA, naoccA, navirA, navirA);
        X->sort(1324, V, 1.0, 0.0);
        V.reset();
        Vab = std::make_shared<Tensor2d>("V (Q|AB)", nQ, navirA, navirA);
        Vab->gemm(false, false, bQijA, X, 1.0, 0.0);
        X.reset();

        // V(Q,AB) += \sum{m,n} V(mA,nB) b(Q,mn)
        V = std::make_shared<Tensor2d>("V (iA|jB)", naoccB, navirA, naoccB, navirA);
        V->read(psio_, PSIF_DFOCC_AMPS);
        X = std::make_shared<Tensor2d>("X (mn|AB)", naoccB, naoccB, navirA, navirA);
        X->sort(1324, V, 1.0, 0.0);
        V.reset();
        Vab->gemm(false, false, bQijB, X, 1.0, 1.0);
        X.reset();
        Vab->write(psio_, PSIF_DFOCC_AMPS);
        Vab.reset();

        // Beta Block
        // V(Q,ab) = \sum{m,n} V(ma,nb) b(Q,mn)
        V = std::make_shared<Tensor2d>("V (ia|jb)", naoccB, navirB, naoccB, navirB);
        V->read(psio_, PSIF_DFOCC_AMPS);
        X = std::make_shared<Tensor2d>("X (mn|ab)", naoccB, naoccB, navirB, navirB);
        X->sort(1324, V, 1.0, 0.0);
        V.reset();
        Vab = std::make_shared<Tensor2d>("V (Q|ab)", nQ, navirB, navirB);
        Vab->gemm(false, false, bQijB, X, 1.0, 0.0);
        X.reset();

        // V(Q,ab) += \sum{M,N} V(Ma,Nb) b(Q,MN)
        V = std::make_shared<Tensor2d>("V (Ia|Jb)", naoccA, navirB, naoccA, navirB);
        V->read(psio_, PSIF_DFOCC_AMPS);
        X = std::make_shared<Tensor2d>("X (MN|ab)", naoccA, naoccA, navirB, navirB);
        X->sort(1324, V, 1.0, 0.0);
        V.reset();
        Vab->gemm(false, false, bQijA, X, 1.0, 1.0);
        X.reset();
        Vab->write(psio_, PSIF_DFOCC_AMPS);
        Vab.reset();

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Nu(Q,IJ) Intermediates ////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //// Nu(Q,IJ) = \sum{M,E} V(IM,JE) b(Q,ME) + \sum{m,e} V(Im,Je) b(Q,me)       (98)
        //V = std::make_shared<Tensor2d>("V <IJ|KA>", naoccA, naoccA, naoccA, navirA);
        //V->read(psio_, PSIF_DFOCC_AMPS);
        //X = std::make_shared<Tensor2d>("X (ME|IJ)", naoccA, navirA, naoccA, naoccA);
        //X->sort(2413, V, 1.0, 0.0);
        //V.reset();
        //Vij = std::make_shared<Tensor2d>("calV (Q|IJ)", nQ, naoccA, naoccA);
        //Vij->gemm(false, false, bQiaA, X, 1.0, 0.0);
        //X.reset();
        //V = std::make_shared<Tensor2d>("V <Ij|Ka>", naoccA, naoccB, naoccA, navirB);
        //V->read(psio_, PSIF_DFOCC_AMPS);
        //X = std::make_shared<Tensor2d>("X (me|IJ)", naoccB, navirB, naoccA, naoccA);
        //X->sort(2413, V, 1.0, 0.0);
        //V.reset();
        //Vij->gemm(false, false, bQiaB, X, 1.0, 1.0);
        //X.reset();
        //Vij->write(psio_, PSIF_DFOCC_AMPS);
        //Vij.reset();

        //// Nu(Q,ij) = \sum{m,e} V(im,je) b(Q,me) + \sum{M,E} V(iM,jE) b(Q,ME)       (99)
        //V = std::make_shared<Tensor2d>("V <ij|ka>", naoccB, naoccB, naoccB, navirB);
        //V->read(psio_, PSIF_DFOCC_AMPS);
        //X = std::make_shared<Tensor2d>("X (me|ij)", naoccB, navirB, naoccB, naoccB);
        //X->sort(2413, V, 1.0, 0.0);
        //V.reset();
        //Vij = std::make_shared<Tensor2d>("calV (Q|ij)", nQ, naoccB, naoccB);
        //Vij->gemm(false, false, bQiaB, X, 1.0, 0.0);
        //X.reset();
        //V = std::make_shared<Tensor2d>("V <iJ|kA>", naoccB, naoccA, naoccB, navirA);
        //V->read(psio_, PSIF_DFOCC_AMPS);
        //X = std::make_shared<Tensor2d>("X (ME|ij)", naoccA, navirA, naoccB, naoccB);
        //X->sort(2413, V, 1.0, 0.0);
        //V.reset();
        //Vij->gemm(false, false, bQiaA, X, 1.0, 1.0);
        //X.reset();
        //Vij->write(psio_, PSIF_DFOCC_AMPS);
        //Vij.reset();
     }// end else if (reference_ == "UNRESTRICTED")

    // outfile->Printf("\t3indices done.\n");

}  // end ccd_pdm_3index_intr

//======================================================================
//    Build y_ia^Q
//======================================================================
void DFOCC::ccd_pdm_yQia() {

    // RHF
    if (reference_ == "RESTRICTED") {
        // defs
        SharedTensor2d K, L, T, U, Tau, V, V2, X, Y, Y2, Z;
        SharedTensor2d Yt, Yp;

        // Read Viajb
        V = std::make_shared<Tensor2d>("V (IA|JB)", naoccA, navirA, naoccA, navirA);
        V->read(psio_, PSIF_DFOCC_AMPS);

        // Read T2
        t2 = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        t2->read_symm(psio_, PSIF_DFOCC_AMPS);

        // Y_iajb = -\sum(me) T(ia,me) X(me,jb)
        // X(me,jb) = V(je,mb)
        X = std::make_shared<Tensor2d>("X (ME|JB)", naoccA, navirA, naoccA, navirA);
        X->sort(3214, V, 1.0, 0.0);
        Y = std::make_shared<Tensor2d>("Y (IA|JB)", naoccA, navirA, naoccA, navirA);
        Y->gemm(false, false, t2, X, -1.0, 0.0);
        X.reset();

        // Y_iabj = -\sum(me) T'(ia,me) X(me,bj)
        T = std::make_shared<Tensor2d>("T2p (IA|JB)", naoccA, navirA, naoccA, navirA);
        ccsd_t2_prime_amps(T, t2);
        t2.reset();
        // X(me,bj) = V(je,mb)
        X = std::make_shared<Tensor2d>("X (ME|BJ)", naoccA, navirA, navirA, naoccA);
        X->sort(3241, V, 1.0, 0.0);
        V.reset();
        Y2 = std::make_shared<Tensor2d>("Y (IA|BJ)", naoccA, navirA, navirA, naoccA);
        Y2->gemm(false, false, T, X, -1.0, 0.0);
        T.reset();
        X.reset();
        Y2->write(psio_, PSIF_DFOCC_AMPS);
        Y2.reset();

        // Read Viabj
        V = std::make_shared<Tensor2d>("V (IA|BJ)", naoccA, navirA, navirA, naoccA);
        V->read(psio_, PSIF_DFOCC_AMPS);

        // Read U2
        U = std::make_shared<Tensor2d>("U2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        U->read_symm(psio_, PSIF_DFOCC_AMPS);

        // Y_iajb += \sum(me) U(ia,me) X(me,jb)
        // X(me,jb) = V_(je,bm)
        X = std::make_shared<Tensor2d>("X (ME|JB)", naoccA, navirA, naoccA, navirA);
        X->sort(4213, V, 1.0, 0.0);
        V.reset();
        Y->gemm(false, false, U, X, 1.0, 1.0);
        U.reset();
        X.reset();

        //===================================
        // Build Yt_ijab
        //===================================

        // Yt_ijab = -1/2 (Y_iajb + Y_jbia)
        Yt = std::make_shared<Tensor2d>("Y2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        Yt->sort(1324, Y, -0.5, 0.0);
        Yt->sort(3142, Y, -0.5, 1.0);
        Y.reset();

        // Yt_ijab += 1/2 Y_ijab
        U = std::make_shared<Tensor2d>("Y <IJ|AB>", naoccA, naoccA, navirA, navirA);
        U->read(psio_, PSIF_DFOCC_AMPS);
        Yt->axpy(U, 0.5);
        U.reset();

        // Yt_ijab -= 1/2 (Y_jabi + Y_ibaj)
        Y = std::make_shared<Tensor2d>("Y (IA|BJ)", naoccA, navirA, navirA, naoccA);
        Y->read(psio_, PSIF_DFOCC_AMPS);
        Yt->sort(4123, Y, -0.5, 1.0);
        Yt->sort(1432, Y, -0.5, 1.0);
        Y.reset();

        // Yt_ijab -= V_jabi + V_ibaj
        V = std::make_shared<Tensor2d>("V (IA|BJ)", naoccA, navirA, navirA, naoccA);
        V->read(psio_, PSIF_DFOCC_AMPS);
        Yt->sort(4123, V, -1.0, 1.0);
        Yt->sort(1432, V, -1.0, 1.0);

        //===================================
        // Build Y'_ijab
        //===================================

        // Y'_ijab = Yt_ijab + V_jabi + V_ibaj
        Yp = std::make_shared<Tensor2d>("Y2p <IJ|AB>", naoccA, naoccA, navirA, navirA);
        Yp->copy(Yt);
        Yp->sort(4123, V, 1.0, 1.0);
        Yp->sort(1432, V, 1.0, 1.0);
        V.reset();

        // Y'_ijab -= V_iajb + V_jbia
        V = std::make_shared<Tensor2d>("V (IA|JB)", naoccA, navirA, naoccA, navirA);
        V->read(psio_, PSIF_DFOCC_AMPS);
        Yp->sort(1324, V, -1.0, 1.0);
        Yp->sort(3142, V, -1.0, 1.0);
        V.reset();

        //===================================
        // Build y_ia^Q
        //===================================

        // X_imae = 2*Yt_imae - Yp_imea
        X = std::make_shared<Tensor2d>("X <IM|AE>", naoccA, naoccA, navirA, navirA);
        X->tei_cs1_anti_symm(Yt, Yp);
        Yt.reset();
        Yp.reset();
        // Y(me,ia) = X(im,ae)
        Y = std::make_shared<Tensor2d>("Y (ME|IA)", naoccA, navirA, naoccA, navirA);
        Y->sort(2413, X, 1.0, 0.0);
        X.reset();
        // y_ia^Q = \sum(me) b(Q,me) * Y(me,ia)
        Z = std::make_shared<Tensor2d>("Y (Q|IA)", nQ, naoccA, navirA);
        Z->gemm(false, false, bQiaA, Y, 1.0, 0.0);
        Y.reset();
        Z->write(psio_, PSIF_DFOCC_AMPS);
    }// end if (reference_ == "RESTRICTED")

    // UHF
    else if (reference_ == "UNRESTRICTED") {
        SharedTensor2d T, T2, L2, Tau, X, Y, Z, V, Vt, U, L, Lt, Yt;
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Y (IA,JB)  Intermediates //////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // AAAA Block
        // Y(IA,JB) = \sum{M,E} [t(MI,AE)] * V(JE,MB) -  \sum{m,e} t(Im,Ae) * V(Je,mB)       (55)
        //X = std::make_shared<Tensor2d>("X (MA|IE)", naoccA, navirA, naoccA, navirA);
        //X->sort(1324, T2, 1.0, 0.0);
        //U->sort(3214, X, 1.0, 0.0);
        //X.reset();
        U = std::make_shared<Tensor2d>("X (IA|ME)", naoccA, navirA, naoccA, navirA);
        T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        U->sort(2314, T2, 1.0, 0.0);
        T2.reset();
        Vt = std::make_shared<Tensor2d>("V (IA|JB)", naoccA, navirA, naoccA, navirA);
        Vt->read(psio_, PSIF_DFOCC_AMPS);
        V = std::make_shared<Tensor2d>("Vtx (ME|JB)", naoccA, navirA, naoccA, navirA);
        V->sort(3214, Vt, 1.0, 0.0);
        Vt.reset();
        Y = std::make_shared<Tensor2d>("Y (IA|JB)", naoccA, navirA, naoccA, navirA);
        Y->gemm(false, false, U, V, 1.0, 0.0);
        U.reset();
        V.reset();

        T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T2->read(psio_, PSIF_DFOCC_AMPS);
        X = std::make_shared<Tensor2d>("X (IA|me)", naoccA, navirA, naoccB, navirB);
        X->sort(1324, T2, 1.0, 0.0);
        T2.reset();
        Vt = std::make_shared<Tensor2d>("V (Ia|jB)", naoccA, navirB, naoccB, navirA);
        Vt->read(psio_, PSIF_DFOCC_AMPS);
        V = std::make_shared<Tensor2d>("Vtx (me|JB)", naoccB, navirB, naoccA, navirA);
        V->sort(3214, Vt, 1.0, 0.0);
        Vt.reset();
        Y->gemm(false, false, X, V, -1.0, 1.0);
        V.reset();
        X.reset();
        Y->write(psio_, PSIF_DFOCC_AMPS);
        Y.reset();

        // BBBB Block
        // Y(ia,jb) = \sum{m,e} [t(mi,ae) ] * V(je,mb) -  \sum{M,E} t(Mi,Ea) * V(jE,Mb)       (56)
        //X = std::make_shared<Tensor2d>("X (ma|ie)", naoccB, navirB, naoccB, navirB);
        //X->sort(1324, T2, 1.0, 1.0);
        //X.reset();
        //U->sort(3214, X, 1.0, 0.0);
        T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        U = std::make_shared<Tensor2d>("X (ia|me)", naoccB, navirB, naoccB, navirB);
        U->sort(2314, T2, 1.0, 0.0);
        T2.reset();
        Vt = std::make_shared<Tensor2d>("V (ia|jb)", naoccB, navirB, naoccB, navirB);
        Vt->read(psio_, PSIF_DFOCC_AMPS);
        V = std::make_shared<Tensor2d>("Vtx (me|jb)", naoccB, navirB, naoccB, navirB);
        V->sort(3214, Vt, 1.0, 0.0);
        Vt.reset();
        Y = std::make_shared<Tensor2d>("Y (ia|jb)", naoccB, navirB, naoccB, navirB);
        Y->gemm(false, false, U, V, 1.0, 0.0);
        U.reset();
        V.reset();

        T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T2->read(psio_, PSIF_DFOCC_AMPS);
        X = std::make_shared<Tensor2d>("X (ia|ME)", naoccB, navirB, naoccA, navirA);
        X->sort(2413, T2, 1.0, 0.0);
        T2.reset();
        Vt = std::make_shared<Tensor2d>("V (iA|Jb)", naoccB, navirA, naoccA, navirB);
        Vt->read(psio_, PSIF_DFOCC_AMPS);
        V = std::make_shared<Tensor2d>("Vtx (ME|jb)", naoccA, navirA, naoccB, navirB);
        V->sort(3214, Vt, 1.0, 0.0);
        Vt.reset();
        Y->gemm(false, false, X, V, -1.0, 1.0);
        V.reset();
        X.reset();
        Y->write(psio_, PSIF_DFOCC_AMPS);
        Y.reset();

        // AABB Block
        // Y_{IAjb} = \sum_{M,E}(t_{MI,AE}) V_{jEMb} - \sum_{m,e} t_{Im}^{Ae} V_{jemb}
        //X = std::make_shared<Tensor2d>("X (MA|IE)", naoccA, navirA, naoccA, navirA);
        //X->sort(1324, T2, 1.0, 1.0);
        //X.reset();
        //U->sort(3214, X, 1.0, 0.0);
        T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        U = std::make_shared<Tensor2d>("X (IA|ME)", naoccA, navirA, naoccA, navirA);
        U->sort(2314, T2, 1.0, 0.0);
        T2.reset();
        Vt = std::make_shared<Tensor2d>("V (iA|Jb)", naoccB, navirA, naoccA, navirB);
        Vt->read(psio_, PSIF_DFOCC_AMPS);
        V = std::make_shared<Tensor2d>("Vtx (ME|jb)", naoccA, navirA, naoccB, navirB);
        V->sort(3214, Vt, 1.0, 0.0);
        Vt.reset();
        Y = std::make_shared<Tensor2d>("Y (IA|jb)", naoccA, navirA, naoccB, navirB);
        Y->gemm(false, false, U, V, 1.0, 0.0);
        U.reset();
        V.reset();

        T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T2->read(psio_, PSIF_DFOCC_AMPS);

        X = std::make_shared<Tensor2d>("X (IA|me)", naoccA, navirA, naoccB, navirB);
        X->sort(1324, T2, 1.0, 0.0);
        T2.reset();
        Vt = std::make_shared<Tensor2d>("V (ia|jb)", naoccB, navirB, naoccB, navirB);
        Vt->read(psio_, PSIF_DFOCC_AMPS);
        V = std::make_shared<Tensor2d>("Vtx (me|jb)", naoccB, navirB, naoccB, navirB);
        V->sort(3214, Vt, 1.0, 0.0);
        Vt.reset();
        Y->gemm(false, false, X, V, -1.0, 1.0);
        V.reset();
        X.reset();
        Y->write(psio_, PSIF_DFOCC_AMPS);
        Y.reset();

        // BBAA Block
        // Y_{iaJB} = \sum_{m,e} t_{mi,ae} {V}_{JemB} - \sum_{M,E}(t_{Mi}^{Ea} {V}_{JEMB}
        //X = std::make_shared<Tensor2d>("X (ma|ie)", naoccB, navirB, naoccB, navirB);
        //X->sort(1324, T2, 1.0, 1.0);
        //U->sort(3214, X, 1.0, 0.0);
        //X.reset();
        T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        U = std::make_shared<Tensor2d>("X (ia|me)", naoccB, navirB, naoccB, navirB);
        U->sort(2314, T2, 1.0, 0.0);
        T2.reset();
        Vt = std::make_shared<Tensor2d>("V (Ia|jB)", naoccA, navirB, naoccB, navirA);
        Vt->read(psio_, PSIF_DFOCC_AMPS);
        V = std::make_shared<Tensor2d>("Vtx (me|JB)", naoccB, navirB, naoccA, navirA);
        V->sort(3214, Vt, 1.0, 0.0);
        Vt.reset();
        Y = std::make_shared<Tensor2d>("Y (ia|JB)", naoccB, navirB, naoccA, navirA);
        Y->gemm(false, false, U, V, 1.0, 0.0);
        U.reset();
        V.reset();

        T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T2->read(psio_, PSIF_DFOCC_AMPS);
        X = std::make_shared<Tensor2d>("X (ia|ME)", naoccB, navirB, naoccA, navirA);
        X->sort(2413, T2, 1.0, 0.0);
        T2.reset();
        Vt = std::make_shared<Tensor2d>("V (IA|JB)", naoccA, navirA, naoccA, navirA);
        Vt->read(psio_, PSIF_DFOCC_AMPS);
        V = std::make_shared<Tensor2d>("Vtx (ME|JB)", naoccA, navirA, naoccA, navirA);
        V->sort(3214, Vt, 1.0, 0.0);
        Vt.reset();
        Y->gemm(false, false, X, V, -1.0, 1.0);
        V.reset();
        X.reset();
        Y->write(psio_, PSIF_DFOCC_AMPS);
        Y.reset();

        // ABBA Block
        // Y_{IajB} = \sum_{m,E} (t_{Im,Ea}) {V}_{jEmB}
        //X = std::make_shared<Tensor2d>("X (ma|IE)", naoccB, navirB, naoccA, navirA);
        //X->sort(2413, T2, 1.0, 1.0);
        //X.reset();
        //U->sort(3214, X, 1.0, 0.0);
        T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T2->read(psio_, PSIF_DFOCC_AMPS);
        U = std::make_shared<Tensor2d>("X (Ia|mE)", naoccA, navirB, naoccB, navirA);
        U->sort(1423, T2, 1.0, 0.0);
        T2.reset();
        Vt = std::make_shared<Tensor2d>("V (iA|jB)", naoccB, navirA, naoccB, navirA);
        Vt->read(psio_, PSIF_DFOCC_AMPS);
        V = std::make_shared<Tensor2d>("Vtx (mE|jB)", naoccB, navirA, naoccB, navirA);
        V->sort(3214, Vt, 1.0, 0.0);
        Vt.reset();
        Y = std::make_shared<Tensor2d>("Y (Ia|jB)", naoccA, navirB, naoccB, navirA);
        Y->gemm(false, false, U, V, 1.0, 0.0);
        U.reset();
        V.reset();
        Y->write(psio_, PSIF_DFOCC_AMPS);
        Y.reset();

        // BAAB Block
        // Y_{iAJb} = \sum_{M,e} (t_{Mi,Ae}) {V}_{JeMb}
        //X = std::make_shared<Tensor2d>("X (MA|ie)", naoccA, navirA, naoccB, navirB);
        //X->sort(1324, T2, 1.0, 1.0);
        //U->sort(3214, X, 1.0, 0.0);
        //X.reset();
        T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T2->read(psio_, PSIF_DFOCC_AMPS);
        U = std::make_shared<Tensor2d>("X (iA|Me)", naoccB, navirA, naoccA, navirB);
        U->sort(2314, T2, 1.0, 0.0);
        T2.reset();
        Vt = std::make_shared<Tensor2d>("V (Ia|Jb)", naoccA, navirB, naoccA, navirB);
        Vt->read(psio_, PSIF_DFOCC_AMPS);
        V = std::make_shared<Tensor2d>("Vtx (Me|Jb)", naoccA, navirB, naoccA, navirB);
        V->sort(3214, Vt, 1.0, 0.0);
        Vt.reset();
        Y = std::make_shared<Tensor2d>("Y (iA|Jb)", naoccB, navirA, naoccA, navirB);
        Y->gemm(false, false, U, V, 1.0, 0.0);
        U.reset();
        V.reset();
        Y->write(psio_, PSIF_DFOCC_AMPS);
        Y.reset();

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Y (IJ,AB)  Intermediates //////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Y(IJ,AB) = 0.5 * \sum{M,N} T(MN,AB) * V(IJ,MN)       (60)
        Tau = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        Tau->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        V = std::make_shared<Tensor2d>("V <IJ|KL>", naoccA, naoccA, naoccA, naoccA);
        V->read(psio_, PSIF_DFOCC_AMPS);
        Y = std::make_shared<Tensor2d>("Y <IJ|AB>", naoccA, naoccA, navirA, navirA);
        Y->gemm(false, false, V, Tau, 0.5, 0.0);
        Tau.reset();
        V.reset();
        Y->write(psio_, PSIF_DFOCC_AMPS);
        Y.reset();

        // Y(ij,ab) = 0.5 * \sum{m,n} T(mn,ab) * V(ij,mn)       (61)
        Tau = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        Tau->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        V = std::make_shared<Tensor2d>("V <ij|kl>", naoccB, naoccB, naoccB, naoccB);
        V->read(psio_, PSIF_DFOCC_AMPS);
        Y = std::make_shared<Tensor2d>("Y <ij|ab>", naoccB, naoccB, navirB, navirB);
        Y->gemm(false, false, V, Tau, 0.5, 0.0);
        Tau.reset();
        V.reset();
        Y->write(psio_, PSIF_DFOCC_AMPS);
        Y.reset();

        // Y(Ij,Ab) = \sum{M,n} T(Mn,Ab) * V(Ij,Mn)        (62)
        Tau = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        Tau->read(psio_, PSIF_DFOCC_AMPS);
        V = std::make_shared<Tensor2d>("V <Ij|Kl>", naoccA, naoccB, naoccA, naoccB);
        V->read(psio_, PSIF_DFOCC_AMPS);
        Y = std::make_shared<Tensor2d>("Y <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        Y->gemm(false, false, V, Tau, 1.0, 0.0);
        Tau.reset();
        V.reset();
        Y->write(psio_, PSIF_DFOCC_AMPS);
        Y.reset();

        // Y(iJ,aB) = \sum{m,N} T(Nm,Ba) * V(iJ,mN) - 0.5 * \sum{M,n} T(Mn,Ba) * V(iJ,Mn)       (63)   // hatali sonucta kontrol et.
        Tau = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        Tau->read(psio_, PSIF_DFOCC_AMPS);
        X = std::make_shared<Tensor2d>("X <mN|aB>", naoccB, naoccA, navirB, navirA);
        X->sort(2143, Tau, 1.0, 0.0);
        Tau.reset();
        V = std::make_shared<Tensor2d>("V <Ij|Kl>", naoccA, naoccB, naoccA, naoccB);
        V->read(psio_, PSIF_DFOCC_AMPS);
        U = std::make_shared<Tensor2d>("X <iJ|mN>", naoccB, naoccA, naoccB, naoccA);
        U->sort(2143, V, 1.0, 0.0);
        V.reset();
        Y = std::make_shared<Tensor2d>("Y <iJ|aB>", naoccB, naoccA, navirB, navirA);
        Y->gemm(false, false, U, X, 1.0, 0.0);
        X.reset();
        U.reset();
        Y->write(psio_, PSIF_DFOCC_AMPS);
        Y.reset();

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Ytilde (IJ,AB)  Intermediates /////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Ytilde(IJ,AB) = 0.5 * [Y(IJ,AB) - Y(IA,JB) + Y(JA,IB) + Y(IB,JA) - Y(JB,IA)]
                        // + V(IB,JA) + V(JA,IB)    (65)

        Y = std::make_shared<Tensor2d>("Y <IJ|AB>", naoccA, naoccA, navirA, navirA);
        Y->read(psio_, PSIF_DFOCC_AMPS);
        Yt = std::make_shared<Tensor2d>("Y2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        Yt->axpy(Y, 0.5);
        Y.reset();
        U = std::make_shared<Tensor2d>("Y (IA|JB)", naoccA, navirA, naoccA, navirA);
        U->read(psio_, PSIF_DFOCC_AMPS);
        Yt->sort(1324, U, -0.5, 1.0);
        Yt->sort(3124, U, 0.5, 1.0);
        Yt->sort(1342, U, 0.5, 1.0);
        Yt->sort(3142, U, -0.5, 1.0);
        U.reset();

        V = std::make_shared<Tensor2d>("V (IA|JB)", naoccA, navirA, naoccA, navirA);
        V->read(psio_, PSIF_DFOCC_AMPS);
        Yt->sort(1342, V, 1.0, 1.0);
        Yt->sort(3124, V, 1.0, 1.0);
        V.reset();
        Yt->write(psio_, PSIF_DFOCC_AMPS);
        Yt.reset();

        // Ytilde(ij,ab) = 0.5 * [Y(ij,ab) - Y(ia,jb) + Y(ja,ib) + Y(ib,ja) - Y(jb,ia)]
                        // + V(ib,ja) + V(ja,ib)    (66)
        Y = std::make_shared<Tensor2d>("Y <ij|ab>", naoccB, naoccB, navirB, navirB);
        Y->read(psio_, PSIF_DFOCC_AMPS);
        Yt = std::make_shared<Tensor2d>("Y2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        Yt->axpy(Y, 0.5);
        Y.reset();
        U = std::make_shared<Tensor2d>("Y (ia|jb)", naoccB, navirB, naoccB, navirB);
        U->read(psio_, PSIF_DFOCC_AMPS);
        Yt->sort(1324, U, -0.5, 1.0);
        Yt->sort(3124, U, 0.5, 1.0);
        Yt->sort(1342, U, 0.5, 1.0);
        Yt->sort(3142, U, -0.5, 1.0);
        U.reset();

        V  = std::make_shared<Tensor2d>("V (ia|jb)", naoccB, navirB, naoccB, navirB);
        V->read(psio_, PSIF_DFOCC_AMPS);
        Yt->sort(1342, V, 1.0, 1.0);
        Yt->sort(3124, V, 1.0, 1.0);
        V.reset();
        Yt->write(psio_, PSIF_DFOCC_AMPS);
        Yt.reset();

        // Ytilde(Ij,Ab) = 0.5 * [Y(Ij,Ab) - Y(IA,jb) + Y(jA,Ib) + Y(Ib,jA) - Y(jb,IA)]
                        // + V(Ib,jA) + V(jA,Ib)     (67)
        Y = std::make_shared<Tensor2d>("Y <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        Y->read(psio_, PSIF_DFOCC_AMPS);
        Yt = std::make_shared<Tensor2d>("Y2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        Yt->axpy(Y, 0.5);
        Y.reset();
        Y = std::make_shared<Tensor2d>("Y (IA|jb)", naoccA, navirA, naoccB, navirB);
        Y->read(psio_, PSIF_DFOCC_AMPS);
        Yt->sort(1324, Y, -0.5, 1.0);
        Y.reset();
        Y = std::make_shared<Tensor2d>("Y (iA|Jb)", naoccB, navirA, naoccA, navirB);
        Y->read(psio_, PSIF_DFOCC_AMPS);
        Yt->sort(3124, Y, 0.5, 1.0);
        Y.reset();
        Y = std::make_shared<Tensor2d>("Y (Ia|jB)", naoccA, navirB, naoccB, navirA);
        Y->read(psio_, PSIF_DFOCC_AMPS);
        Yt->sort(1342, Y, 0.5, 1.0);
        Y.reset();
        Y = std::make_shared<Tensor2d>("Y (ia|JB)", naoccB, navirB, naoccA, navirA);
        Y->read(psio_, PSIF_DFOCC_AMPS);
        Yt->sort(3142, Y, -0.5, 1.0);
        Y.reset();

        V = std::make_shared<Tensor2d>("V (Ia|jB)", naoccA, navirB, naoccB, navirA);
        V->read(psio_, PSIF_DFOCC_AMPS);
        Yt->sort(1342, V, 1.0, 1.0);
        V.reset();
        V = std::make_shared<Tensor2d>("V (iA|Jb)", naoccB, navirA, naoccA, navirB);
        V->read(psio_, PSIF_DFOCC_AMPS);
        Yt->sort(3124, V, 1.0, 1.0);
        V.reset();
        Yt->write(psio_, PSIF_DFOCC_AMPS);
        Yt.reset();

        // Ytilde(iJ,aB) = 0.5 * [Y(iJ,aB) - Y(ia,JB) + Y(Ja,iB) + Y(iB,Ja) - Y(JB,ia)]
                        // + V(iB,Ja) + V(Ja,iB)    (65)
        Y = std::make_shared<Tensor2d>("Y <iJ|aB>", naoccB, naoccA, navirB, navirA);
        Y->read(psio_, PSIF_DFOCC_AMPS);
        Yt = std::make_shared<Tensor2d>("Y2 <iJ|aB>", naoccB, naoccA, navirB, navirA);
        Yt->axpy(Y, 0.5);
        Y.reset();

        Y = std::make_shared<Tensor2d>("Y (ia|JB)", naoccB, navirB, naoccA, navirA);
        Y->read(psio_, PSIF_DFOCC_AMPS);
        Yt->sort(1324, Y, -0.5, 1.0);
        Y.reset();

        Y = std::make_shared<Tensor2d>("Y (Ia|jB)", naoccA, navirB, naoccB, navirA);
        Y->read(psio_, PSIF_DFOCC_AMPS);
        Yt->sort(3124, Y, 0.5, 1.0);
        Y.reset();

        Y = std::make_shared<Tensor2d>("Y (iA|Jb)", naoccB, navirA, naoccA, navirB);
        Y->read(psio_, PSIF_DFOCC_AMPS);
        Yt->sort(1342, Y, 0.5, 1.0);
        Y.reset();

        Y = std::make_shared<Tensor2d>("Y (IA|jb)", naoccA, navirA, naoccB, navirB);
        Y->read(psio_, PSIF_DFOCC_AMPS);
        Yt->sort(3142, Y, -0.5, 1.0);
        Y.reset();

        V = std::make_shared<Tensor2d>("V (iA|Jb)", naoccB, navirA, naoccA, navirB);
        V->read(psio_, PSIF_DFOCC_AMPS);
        Yt->sort(1342, V, 1.0, 1.0);
        V.reset();
        V = std::make_shared<Tensor2d>("V (Ia|jB)", naoccA, navirB, naoccB, navirA);
        V->read(psio_, PSIF_DFOCC_AMPS);
        Yt->sort(3124, V, 1.0, 1.0);
        V.reset();
        Yt->write(psio_, PSIF_DFOCC_AMPS);
        Yt.reset();

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // y(Q,IA) Intermediates /////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // y(Q,IA) = \sum{M,E} Ytilde(IM,AE) b(Q,ME) + \sum{m,e} Ytilde(Im,Ae) b(Q,me)       (113)
        Yt = std::make_shared<Tensor2d>("Y2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        Yt->read(psio_, PSIF_DFOCC_AMPS);
        X = std::make_shared<Tensor2d>("X (ME|IA)", naoccA, navirA, naoccA, navirA);
        X->sort(2413, Yt, 1.0, 0.0);
        Yt.reset();
        Z = std::make_shared<Tensor2d>("Y (Q|IA)", nQ, naoccA, navirA);
        Z->gemm(false, false, bQiaA, X, 1.0, 0.0);
        X.reset();
        Yt = std::make_shared<Tensor2d>("Y2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        Yt->read(psio_, PSIF_DFOCC_AMPS);
        X = std::make_shared<Tensor2d>("X (me|IA)", naoccB, navirB, naoccA, navirA);
        X->sort(2413, Yt, 1.0, 0.0);
        Yt.reset();
        Z->gemm(false, false, bQiaB, X, 1.0, 1.0);
        X.reset();
        Z->write(psio_, PSIF_DFOCC_AMPS);
        Z.reset();

        // y(Q,ia) = \sum{m,e} Ytilde(im,ae) b(Q,me) + \sum{M,E} Ytilde(iM,aE) b(Q,ME)       (114)
        Yt = std::make_shared<Tensor2d>("Y2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        Yt->read(psio_, PSIF_DFOCC_AMPS);
        X = std::make_shared<Tensor2d>("X (me|ia)", naoccB, navirB, naoccB, navirB);
        X->sort(2413, Yt, 1.0, 0.0);
        Yt.reset();
        Z = std::make_shared<Tensor2d>("Y (Q|ia)", nQ, naoccB, navirB);
        Z->gemm(false, false, bQiaB, X, 1.0, 0.0);
        X.reset();
        Yt = std::make_shared<Tensor2d>("Y2 <iJ|aB>", naoccB, naoccA, navirB, navirA);
        Yt->read(psio_, PSIF_DFOCC_AMPS);
        X = std::make_shared<Tensor2d>("X (ME|ia)", naoccA, navirA, naoccB, navirB);
        X->sort(2413, Yt, 1.0, 0.0);
        Yt.reset();
        Z->gemm(false, false, bQiaA, X, 1.0, 1.0);
        X.reset();
        Z->write(psio_, PSIF_DFOCC_AMPS);
    }// end else if (reference_ == "UNRESTRICTED")
}  // end ccd_pdm_yQia

}  // namespace dfoccwave
}  // namespace psi
