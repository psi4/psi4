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

void DFOCC::ccsd_3index_intr() {

    // RHF
    if (reference_ == "RESTRICTED") {
        // defs
        SharedTensor2d K, T, U, Tau;

        // T(Q,ia) = \sum_{jb} b_jb^Q u_ij^ab = \sum_{jb} b(Q,jb) U(jb,ia)
        U = std::make_shared<Tensor2d>("U2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        ccsd_u2_amps(U, t2);
        T = std::make_shared<Tensor2d>("T2 (Q|IA)", nQ, naoccA, navirA);
        T->gemm(false, false, bQiaA, U, 1.0, 0.0);
        U.reset();
        T->write(psio_, PSIF_DFOCC_AMPS);
        T.reset();

        // t(Q) = 2\sum_{mf} t_m^f b_mf^Q
        T1c->gemv(false, bQiaA, t1A, 2.0, 0.0);

        // t(Q,ij) = \sum_{e} t_i^e b_je^Q
        T = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
        T->contract233(false, true, naoccA, naoccA, t1A, bQiaA, 1.0, 0.0);
        T->write(psio_, PSIF_DFOCC_AMPS);

        // Ttilde(Q,ia) = \sum_{m} t_m^a t_im^Q
        U = std::make_shared<Tensor2d>("T1tilde (Q|IA)", nQ, naoccA, navirA);
        U->contract(false, false, nQ * naoccA, navirA, naoccA, T, t1A, 1.0, 0.0);
        T.reset();
        U->write(psio_, PSIF_DFOCC_AMPS);
        U.reset();

        // Tau(Q,ia) = \sum_{mf} (2*Tau_im^af - Tau_mi^af) b_mf^Q
        Tau = std::make_shared<Tensor2d>("Tau2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        ccsd_tau_tilde_amps(Tau, t2);
        U = std::make_shared<Tensor2d>("2*Tau2(IA,JB) - Tau2(IB,JA)", naoccA, navirA, naoccA, navirA);
        U->sort(1432, Tau, 1.0, 0.0);
        U->scale(-1.0);
        U->axpy(Tau, 2.0);
        Tau.reset();
        T = std::make_shared<Tensor2d>("Tau2 (Q|IA)", nQ, naoccA, navirA);
        T->gemm(false, true, bQiaA, U, 1.0, 0.0);
        U.reset();
        T->write(psio_, PSIF_DFOCC_AMPS);
        T.reset();

        // t(Q,ab) = \sum_{m} t_m^a b_mb^Q
        T = std::make_shared<Tensor2d>("T1 (Q|AB)", nQ, navirA, navirA);
        T->contract233(true, false, navirA, navirA, t1A, bQiaA, 1.0, 0.0);
        T->write(psio_, PSIF_DFOCC_AMPS);
        T.reset();

        // t(Q,ia) = \sum_{f} t_i^f b_fa^Q
        T = std::make_shared<Tensor2d>("T1 (Q|IA)", nQ, naoccA, navirA);
        T->contract233(false, false, naoccA, navirA, t1A, bQabA, 1.0, 0.0);
        T->write(psio_, PSIF_DFOCC_AMPS);
        T.reset();

        // t(Q,ai) = \sum_{m} t_m^a b_mi^Q
        T = std::make_shared<Tensor2d>("T1 (Q|AI)", nQ, navirA, naoccA);
        T->contract233(true, false, navirA, naoccA, t1A, bQijA, 1.0, 0.0);
        T->write(psio_, PSIF_DFOCC_AMPS);
        T.reset();

        // Read intermediates
        T = std::make_shared<Tensor2d>("T1 (Q|IA)", nQ, naoccA, navirA);
        T->read(psio_, PSIF_DFOCC_AMPS);
        U = std::make_shared<Tensor2d>("T1 (Q|AI)", nQ, navirA, naoccA);
        U->read(psio_, PSIF_DFOCC_AMPS);
        K = std::make_shared<Tensor2d>("T1tilde (Q|IA)", nQ, naoccA, navirA);
        K->read(psio_, PSIF_DFOCC_AMPS);
        // t'(Q,ia) = t(Q,ia) - t(Q,ai) - tt(Q,ia)
        Tau = std::make_shared<Tensor2d>("T1p (Q|IA)", nQ, naoccA, navirA);
        Tau->swap_3index_col(U);
        Tau->add(K);
        K.reset();
        Tau->scale(-1.0);
        Tau->add(T);
        Tau->write(psio_, PSIF_DFOCC_AMPS);
        Tau.reset();
        // Tau'(Q,ia) = t(Q,ia) + Tau(Q,ia)
        K = std::make_shared<Tensor2d>("Tau2 (Q|IA)", nQ, naoccA, navirA);
        K->read(psio_, PSIF_DFOCC_AMPS);
        Tau = std::make_shared<Tensor2d>("Tau2p (Q|IA)", nQ, naoccA, navirA);
        Tau->copy(T);
        Tau->add(K);
        T.reset();
        Tau->write(psio_, PSIF_DFOCC_AMPS);
        Tau.reset();
        // Tau''(Q,ia) = -t(Q,ai) + Tau(Q,ia)
        Tau = std::make_shared<Tensor2d>("Tau2pp (Q|IA)", nQ, naoccA, navirA);
        Tau->swap_3index_col(U);
        U.reset();
        Tau->scale(-1.0);
        Tau->add(K);
        K.reset();
        Tau->write(psio_, PSIF_DFOCC_AMPS);
    }  // if (reference_ == "RESTRICTED")

    // UHF
    else if (reference_ == "UNRESTRICTED") {
        SharedTensor2d TQiaA, TQiaB, TautAA, TautBB, TautAB, TautQA, TautQB;
        SharedTensor2d tQovA, tQovB, tQooA, tQooB, tQvoA, tQvoB, tQvvA, tQvvB, ttQA, ttQB;
        SharedTensor2d L, M, N, X, Y, T, K, U, Tau, T2;

        // T(Q,IA) = \sum(J,B) (T(IJ,AB) * b(Q,JB)) + \sum(j,b) (T2AB(Ij,Ab) * b(Q,jb)) (45)
        // T(Q,IA) += \sum(J,B) (T(IJ,AB) * b(Q,JB))
        L = std::make_shared<Tensor2d>("L (IA|JB)", naoccA, navirA, naoccA, navirA);
        T = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        T->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        L->sort(1324, T, 1.0, 0.0);
        T.reset();
        TQiaA = std::make_shared<Tensor2d>("T2 (Q|IA)", nQ, naoccA, navirA);
        TQiaA->gemm(false, true, bQiaA, L, 1.0, 0.0);
        L.reset();
        // T(Q,IA) += \sum(j,b) (T2AB(Ij,Ab) * b(Q,jb))
        M = std::make_shared<Tensor2d>("M (IA|jb)", naoccA, navirA, naoccB, navirB);
        T = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T->read(psio_, PSIF_DFOCC_AMPS);
        M->sort(1324, T, 1.0, 0.0);
        T.reset();
        TQiaA->gemm(false, true, bQiaB, M, 1.0, 1.0);
        M.reset();
        TQiaA->write(psio_, PSIF_DFOCC_AMPS);
        TQiaA.reset();

        // T(Q,ia) = \sum(j,b) (T2BB(ij,ab) * b(Q,jb)) + \sum(J,B) (T2AB(Ji,Ba) * b(Q,JB)) (46)
        // T(Q,ia) +=  \sum(j,b) (T2BB(ij,ab) * b(Q,jb))
        N = std::make_shared<Tensor2d>("N (ia|jb)", naoccB, navirB, naoccB, navirB);
        T = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        T->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        N->sort(1324, T, 1.0, 0.0);
        T.reset();
        TQiaB = std::make_shared<Tensor2d>("T2 (Q|ia)", nQ, naoccB, navirB);
        TQiaB->gemm(false, true, bQiaB, N, 1.0, 0.0);
        N.reset();
        // T(Q,ia) +=  \sum(J,B) (T2AB(Ji,Ba) * b(Q,JB))
        M = std::make_shared<Tensor2d>("M (IA|jb)", naoccA, navirA, naoccB, navirB);
        T = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T->read(psio_, PSIF_DFOCC_AMPS);
        M->sort(1324, T, 1.0, 0.0);
        T.reset();
        TQiaB->gemm(false, false, bQiaA, M, 1.0, 1.0);
        M.reset();
        TQiaB->write(psio_, PSIF_DFOCC_AMPS);
        TQiaB.reset();

        // tQ = \sum(M,F) (t1A(M,F) * b(Q,MF)) + \sum(m,f) (t1A(m,f) * b(Q,mf)) (48)
        // tQ += \sum(M,F) (b(Q,MF) * t1A(M,F))
        T1c->gemv(false, bQiaA, t1A, 1.0, 0.0);
        // tQ += \sum(m,f) (b(Q,mf) * t1A(m,f))
        T1c->gemv(false, bQiaB, t1B, 1.0, 1.0);
        //t(Q,IA) = \sum(F) t1A(I,F) * b(Q,AF)  (50)
        tQovA = std::make_shared<Tensor2d>("T1 (Q|IA)", nQ, naoccA, navirA);
        tQovA->contract233(false, false, naoccA, navirA, t1A, bQabA, 1.0, 0.0);
        tQovA->write(psio_, PSIF_DFOCC_AMPS);
        tQovA.reset();

        //t(Q,ia) = \sum(f) t1B(i,f) * b(Q,af)  (51)
        tQovB = std::make_shared<Tensor2d>("T1 (Q|ia)", nQ, naoccB, navirB);
        tQovB->contract233(false, false, naoccB, navirB, t1B, bQabB, 1.0, 0.0);
        tQovB->write(psio_, PSIF_DFOCC_AMPS);
        tQovB.reset();

        //t(Q,IJ) = \sum(E) t1A(I,E) * b(Q,JE)   (53)
        tQooA = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
        tQooA->contract233(false, true, naoccA, naoccA, t1A, bQiaA, 1.0, 0.0);
        tQooA->write(psio_, PSIF_DFOCC_AMPS);
        tQooA.reset();

        //t(Q,ij) = \sum(e) t1B(i,e) * b(Q,je)   (54)
        tQooB = std::make_shared<Tensor2d>("T1 (Q|ij)", nQ, naoccB, naoccB);
        tQooB->contract233(false, true, naoccB, naoccB, t1B, bQiaB, 1.0, 0.0);
        tQooB->write(psio_, PSIF_DFOCC_AMPS);
        tQooB.reset();

        //t(Q,AI) = \sum(M) t1A(M,A) * b(Q,MI)  (56)
        tQvoA = std::make_shared<Tensor2d>("T1 (Q|AI)", nQ, navirA, naoccA);
        tQvoA->contract233(true, false, navirA, naoccA, t1A, bQijA, 1.0, 0.0);
        tQvoA->write(psio_, PSIF_DFOCC_AMPS);
        tQvoA.reset();

        //t(Q,ai) = \sum(m) t1A(m,a) * b(Q,mi)  (57)
        tQvoB = std::make_shared<Tensor2d>("T1 (Q|ai)", nQ, navirB, naoccB);
        tQvoB->contract233(true, false, navirB, naoccB, t1B, bQijB, 1.0, 0.0);
        tQvoB->write(psio_, PSIF_DFOCC_AMPS);
        tQvoB.reset();

        //t(Q,AB) = \sum(M) t1A(M,A) * b(Q,MB)  (59)
        tQvvA = std::make_shared<Tensor2d>("T1 (Q|AB)", nQ, navirA, navirA);
        tQvvA->contract233(true, false, navirA, navirA, t1A, bQiaA, 1.0, 0.0);
        tQvvA->write(psio_, PSIF_DFOCC_AMPS);
        tQvvA.reset();

        //t(Q,ab) = \sum(m) t1A(m,a) * b(Q,mb)  (60)
        tQvvB = std::make_shared<Tensor2d>("T1 (Q|ab)", nQ, navirB, navirB);
        tQvvB->contract233(true, false, navirB, navirB, t1B, bQiaB, 1.0, 0.0);
        tQvvB->write(psio_, PSIF_DFOCC_AMPS);
        tQvvB.reset();

        // Tautilde(Q,IA) = \sum(M,F) Taut(IM,AF) * b(Q,MF) + \sum(m,f) Taut(Im,Af) * b(Q,mf) (62)
        // Tautilde(Q,IA) += \sum(M,F) Taut(IM,AF) * b(Q,MF)
        T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        Tau = std::make_shared<Tensor2d>("Tautilde <IJ|AB>", naoccA, naoccA, navirA, navirA);
        uccsd_tau_tilde_amps(naoccA, naoccA, navirA, navirA, Tau, T2, t1A, t1A);
        T2.reset();
        L = std::make_shared<Tensor2d>("L (IA|JB)", naoccA, navirA, naoccA, navirA);
        L->sort(1324, Tau, 1.0, 0.0);
        Tau.reset();
        TautQA = std::make_shared<Tensor2d>("Tau2 (Q|IA)", nQ, naoccA, navirA);
        TautQA->gemm(false, true, bQiaA, L, 1.0, 0.0);
        L.reset();

        // Tautilde(Q,IA) += \sum(m,f) Taut(Im,Af) * b(Q,mf)
        Tau = std::make_shared<Tensor2d>("Tautilde <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T2->read(psio_, PSIF_DFOCC_AMPS);
        uccsd_tau_tilde_amps_OS(naoccA, naoccB, navirA, navirB, Tau, T2, t1A, t1B);
        T2.reset();
        M = std::make_shared<Tensor2d>("M (IA|jb)", naoccA, navirA, naoccB, navirB);
        M->sort(1324, Tau, 1.0, 0.0);
        Tau.reset();
        TautQA->gemm(false, true, bQiaB, M, 1.0, 1.0);
        M.reset();
        TautQA->write(psio_, PSIF_DFOCC_AMPS);
        TautQA.reset();

        // Tautilde(Q,ia) = \sum(m,f) Taut(im,af) * b(Q,mf) + \sum(M,F) Taut(Mi,Fa) * b(Q,MF)  (63)
        // Tautilde(Q,ia) +=  \sum(m,f) Taut(im,af) * b(Q,mf)
        T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        Tau = std::make_shared<Tensor2d>("Tautilde <ij|ab>", naoccB, naoccB, navirB, navirB);
        uccsd_tau_tilde_amps(naoccB, naoccB, navirB, navirB, Tau, T2, t1B, t1B);
        T2.reset();
        N = std::make_shared<Tensor2d>("N (ia|jb)", naoccB, navirB, naoccB, navirB);
        N->sort(1324, Tau, 1.0, 0.0);
        Tau.reset();
        TautQB = std::make_shared<Tensor2d>("Tau2 (Q|ia)", nQ, naoccB, navirB);
        TautQB->gemm(false, true, bQiaB, N, 1.0, 0.0);
        N.reset();
        // Tautilde(Q,ia) += \sum(M,F) Taut(Mi,Fa) * b(Q,MF)
        Tau = std::make_shared<Tensor2d>("Tautilde <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T2->read(psio_, PSIF_DFOCC_AMPS);
        uccsd_tau_tilde_amps_OS(naoccA, naoccB, navirA, navirB, Tau, T2, t1A, t1B);
        T2.reset();
        M = std::make_shared<Tensor2d>("M (IA|jb)", naoccA, navirA, naoccB, navirB);
        M->sort(1324, Tau, 1.0, 0.0);
        Tau.reset();
        TautQB->gemm(false, false, bQiaA, M, 1.0, 1.0);
        M.reset();
        TautQB->write(psio_, PSIF_DFOCC_AMPS);
        TautQB.reset();

        // ttilde(Q,IA) = \sum(M) t(Q,IM) * t(M,A) (65)
        tQooA = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
        tQooA->read(psio_, PSIF_DFOCC_AMPS);
        ttQA = std::make_shared<Tensor2d>("T1tilde (Q|IA)", nQ, naoccA, navirA);
        ttQA->contract(false, false, nQ*naoccA, navirA, naoccA, tQooA, t1A, 1.0, 0.0);
        tQooA.reset();
        ttQA->write(psio_, PSIF_DFOCC_AMPS);
        ttQA.reset();

        // ttilde(Q,ia) = \sum(m) t(Q,im) * t(m,a) (66)
        tQooB = std::make_shared<Tensor2d>("T1 (Q|ij)", nQ, naoccB, naoccB);
        tQooB->read(psio_, PSIF_DFOCC_AMPS);
        ttQB = std::make_shared<Tensor2d>("T1tilde (Q|ia)", nQ, naoccB, navirB);
        ttQB->contract(false, false, nQ*naoccB, navirB, naoccB, tQooB, t1B, 1.0, 0.0);
        tQooB.reset();
        ttQB->write(psio_, PSIF_DFOCC_AMPS);
        ttQB.reset();

        // ttilde(Q,AI) = \sum(E) t(Q,AE) * t(I,E) (89)
        T = std::make_shared<Tensor2d>("T1 (Q|AB)", nQ, navirA, navirA);
        T->read(psio_, PSIF_DFOCC_AMPS);
        X = std::make_shared<Tensor2d>("Temp (Q|IA)", nQ, naoccA, navirA);
        X->contract233(false, true, naoccA, navirA, t1A, T, 1.0, 0.0);
        T.reset();
        U = std::make_shared<Tensor2d>("T1tilde (Q|AI)", nQ, navirA, naoccA);
        U->swap_3index_col(X);
        X.reset();
        U->write(psio_, PSIF_DFOCC_AMPS);
        U.reset();

        // ttilde(Q,ia) = \sum(m) t(Q,im) * t(m,a) (90)
        T = std::make_shared<Tensor2d>("T1 (Q|ab)", nQ, navirB, navirB);
        T->read(psio_, PSIF_DFOCC_AMPS);
        X = std::make_shared<Tensor2d>("Temp (Q|ia)", nQ, naoccB, navirB);
        X->contract233(false, true, naoccB, navirB, t1B, T, 1.0, 0.0);
        T.reset();
        U = std::make_shared<Tensor2d>("T1tilde (Q|ai)", nQ, navirB, naoccB);
        U->swap_3index_col(X);
        X.reset();
        U->write(psio_, PSIF_DFOCC_AMPS);
        U.reset();

        // Alpha Block
        // Read intermediates for Alpha Blocks
        T = std::make_shared<Tensor2d>("T1 (Q|IA)", nQ, naoccA, navirA);
        T->read(psio_, PSIF_DFOCC_AMPS);
        U = std::make_shared<Tensor2d>("T1 (Q|AI)", nQ, navirA, naoccA);
        U->read(psio_, PSIF_DFOCC_AMPS);
        K = std::make_shared<Tensor2d>("T1tilde (Q|IA)", nQ, naoccA, navirA);
        K->read(psio_, PSIF_DFOCC_AMPS);
        // t'(Q,IA) = t(Q,IA) - t(Q,AI) - tt(Q,IA)   (68)
        Tau = std::make_shared<Tensor2d>("T1p (Q|IA)", nQ, naoccA, navirA);
        Tau->swap_3index_col(U);
        Tau->add(K);
        K.reset();
        Tau->scale(-1.0);
        Tau->add(T);
        Tau->write(psio_, PSIF_DFOCC_AMPS);
        Tau.reset();
        // Tau'(Q,IA) = t(Q,IA) + Tau(Q,IA)   (71)
        K = std::make_shared<Tensor2d>("Tau2 (Q|IA)", nQ, naoccA, navirA);
        K->read(psio_, PSIF_DFOCC_AMPS);
        Tau = std::make_shared<Tensor2d>("Tau2p (Q|IA)", nQ, naoccA, navirA);
        Tau->copy(T);
        Tau->axpy(K, 1.0);
        T.reset();
        Tau->write(psio_, PSIF_DFOCC_AMPS);
        Tau.reset();
        // Tau''(Q,IA) = -t(Q,AI) + Tau(Q,IA)  (74)
        Tau = std::make_shared<Tensor2d>("Tau2pp (Q|IA)", nQ, naoccA, navirA);
        Tau->swap_3index_col(U);
        U.reset();
        Tau->scale(-1.0);
        Tau->add(K);
        K.reset();
        Tau->write(psio_, PSIF_DFOCC_AMPS);
        Tau.reset();

        // Beta Block
        // Read intermediates for Beta Blocks
        T = std::make_shared<Tensor2d>("T1 (Q|ia)", nQ, naoccB, navirB);
        T->read(psio_, PSIF_DFOCC_AMPS);
        U = std::make_shared<Tensor2d>("T1 (Q|ai)", nQ, navirB, naoccB);
        U->read(psio_, PSIF_DFOCC_AMPS);
        K = std::make_shared<Tensor2d>("T1tilde (Q|ia)", nQ, naoccB, navirB);
        K->read(psio_, PSIF_DFOCC_AMPS);

        // t'(Q,ia) = t(Q,ia) - t(Q,ai) - tt(Q,ia)  (69)
        Tau = std::make_shared<Tensor2d>("T1p (Q|ia)", nQ, naoccB, navirB);
        Tau->swap_3index_col(U);
        Tau->add(K);
        K.reset();
        Tau->scale(-1.0);
        Tau->add(T);
        Tau->write(psio_, PSIF_DFOCC_AMPS);
        Tau.reset();

        // Tau'(Q,ia) = t(Q,ia) + Tau(Q,ia)    (72)
        K = std::make_shared<Tensor2d>("Tau2 (Q|ia)", nQ, naoccB, navirB);
        K->read(psio_, PSIF_DFOCC_AMPS);
        Tau = std::make_shared<Tensor2d>("Tau2p (Q|ia)", nQ, naoccB, navirB);
        Tau->copy(T);
        Tau->axpy(K, 1.0);
        T.reset();
        Tau->write(psio_, PSIF_DFOCC_AMPS);
        Tau.reset();

        // Tau''(Q,ia) = -t(Q,ai) + Tau(Q,ia)  (75)
        Tau = std::make_shared<Tensor2d>("Tau2pp (Q|ia)", nQ, naoccB, navirB);
        Tau->swap_3index_col(U);
        U.reset();
        Tau->scale(-1.0);
        Tau->add(K);
        K.reset();
        Tau->write(psio_, PSIF_DFOCC_AMPS);
        Tau.reset();

    }  // else if (reference_ == "UNRESTRICTED")

}  // end ccsd_3index_intr

//============================================================================
//    DF-UCCSD: Tau_IJ^AB = T_IJ^AB + t_I^A t_J^B - t_I^B t_J^A
//============================================================================
void DFOCC::uccsd_tau_amps(int occ1, int occ2, int vir1, int vir2, SharedTensor2d &Tau, SharedTensor2d &T2, SharedTensor2d &T1a, SharedTensor2d &T1b) {
     for(int i=0; i<occ1; ++i) {
         for(int j=0; j<occ2; ++j) {
             int ij = (i * occ2) + j;
             for(int a=0; a<vir1; ++a) {
                 for(int b=0; b<vir2; ++b) {
                     int ab = (a * vir2) + b;
                     double value = T2->get(ij,ab) + (T1a->get(i,a) * T1b->get(j,b)) - (T1a->get(i,b) * T1b->get(j,a));
                     Tau->set(ij,ab,value);
                 }
             }
         }
      }
}  // end ccsd_tau_amps

//============================================================================
//    DF-UCCSD: Tau_Ij^Ab = T_Ij^Ab + t_I^A t_j^b
//============================================================================
void DFOCC::uccsd_tau_amps_OS(int occ1, int occ2, int vir1, int vir2, SharedTensor2d &Tau, SharedTensor2d &T2, SharedTensor2d &T1a, SharedTensor2d &T1b) {
     for(int i=0; i<occ1; ++i) {
         for(int j=0; j<occ2; ++j) {
             int ij = (i * occ2) + j;
             for(int a=0; a<vir1; ++a) {
                 for(int b=0; b<vir2; ++b) {
                     int ab = (a * vir2) + b;
                     double value = T2->get(ij,ab) + (T1a->get(i,a) * T1b->get(j,b));
                     Tau->set(ij,ab,value);
                 }
             }
         }
      }
}  // end ccsd_tau_amps_OS

//=============================================================================
//    DF-UCCSD: Tautilde_IJ^AB = T_IJ^AB + 1/2 (t_I^A t_J^B - t_I^B t_J^A)
//=============================================================================
void DFOCC::uccsd_tau_tilde_amps(int occ1, int occ2, int vir1, int vir2, SharedTensor2d &Tautilde, SharedTensor2d &T2, SharedTensor2d &T1a, SharedTensor2d &T1b) {
     for(int i=0; i<occ1; ++i) {
         for(int j=0; j<occ2; ++j) {
             int ij = (i * occ2) + j;
             for(int a=0; a<vir1; ++a) {
                 for(int b=0; b<vir2; ++b) {
                     int ab = (a * vir2) + b;
                     double value = T2->get(ij,ab) + 0.5 * ((T1a->get(i,a) * T1b->get(j,b)) - (T1a->get(i,b) * T1b->get(j,a)));
                     Tautilde->set(ij,ab,value);
                 }
             }
         }
      }
}  // end ccsd_tau_tilde_amps


//=============================================================================
//    DF-UCCSD: Tautilde_Ij^Ab = T_Ij^Ab + 1/2 (t_I^A t_j^b)
//=============================================================================
void DFOCC::uccsd_tau_tilde_amps_OS(int occ1, int occ2, int vir1, int vir2, SharedTensor2d &Tautilde, SharedTensor2d &T2, SharedTensor2d &T1a, SharedTensor2d &T1b) {
     for(int i=0; i<occ1; ++i) {
         for(int j=0; j<occ2; ++j) {
             int ij = (i * occ2) + j;
             for(int a=0; a<vir1; ++a) {
                 for(int b=0; b<vir2; ++b) {
                     int ab = (a * vir2) + b;
                     double value = T2->get(ij,ab) + 0.5 * ((T1a->get(i,a) * T1b->get(j,b)));
                     Tautilde->set(ij,ab,value);
                 }
             }
         }
      }
}  // end ccsd_tau_tilde_amps


}  // namespace dfoccwave
}  // namespace psi
