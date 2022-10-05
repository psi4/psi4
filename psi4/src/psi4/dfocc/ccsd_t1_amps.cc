/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

void DFOCC::ccsd_t1_amps() {

    // RHF
    if (reference_ == "RESTRICTED") {
        // defs
        SharedTensor2d K, T1, T, U, Tau;

        // t_i^a <= \sum_{e} t_i^e Fae
        t1newA->gemm(false, true, t1A, FabA, 1.0, 0.0);

        // t_i^a <= -\sum_{m} t_m^a Fmi
        t1newA->gemm(true, false, FijA, t1A, -1.0, 1.0);

        // t_i^a <= \sum_{m,e} u_im^ae Fme
        U = std::make_shared<Tensor2d>("U2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        ccsd_u2_amps(U, t2);
        t1newA->gemv(false, U, FiaA, 1.0, 1.0);
        U.reset();

        // t_i^a <= \sum_{Q} t_Q b_ia^Q
        t1newA->gemv(true, bQiaA, T1c, 1.0, 1.0);

        // t_i^a <= -\sum_{Q,m} (T_ma^Q + t_ma^Q) b_mi^Q -= \sum_{Q,m} B(Qm,i) X(Qm,a)
        T = std::make_shared<Tensor2d>("T2 (Q|IA)", nQ, naoccA, navirA);
        T->read(psio_, PSIF_DFOCC_AMPS);
        T1 = std::make_shared<Tensor2d>("T1 (Q|IA)", nQ, naoccA, navirA);
        T1->read(psio_, PSIF_DFOCC_AMPS);
        U = std::make_shared<Tensor2d>("T1+T2 (Q|IA)", nQ, naoccA, navirA);
        U->copy(T);
        T.reset();
        U->add(T1);
        T1.reset();
        t1newA->contract(true, false, naoccA, navirA, nQ * naoccA, bQijA, U, -1.0, 1.0);
        U.reset();

        // t_i^a <= \sum_{Q,e} T_ie^Q b_ea^Q = \sum_{Qe} T^Q(i,e) B^Q(e,a)
        T = std::make_shared<Tensor2d>("T2 (Q|IA)", nQ, naoccA, navirA);
        T->read(psio_, PSIF_DFOCC_AMPS);
        Tau = std::make_shared<Tensor2d>("Temp (Q|AI)", nQ, navirA, naoccA);
        Tau->swap_3index_col(T);
        T.reset();
        t1newA->contract(true, false, naoccA, navirA, nQ * navirA, Tau, bQabA, 1.0, 1.0);
        Tau.reset();

        // Denom
        for (int i = 0; i < naoccA; ++i) {
            for (int a = 0; a < navirA; ++a) {
                double value = FockA->get(i + nfrzc, i + nfrzc) - FockA->get(a + noccA, a + noccA);
                t1newA->set(i, a, t1newA->get(i, a) / value);
            }
        }
        // t1newA->print();

    }  // if (reference_ == "RESTRICTED")

    // UHF
    else if (reference_ == "UNRESTRICTED") {
        SharedTensor2d T, X, Y, K;

        // Alpha Part
        // t(I,A) * D(I,A) = f(I,A)
        // D(I,A) = f(I,I) - f(A,A)
        for(int i=0; i<naoccA; ++i) {
            for(int a=0; a<navirA; ++a) {
               t1newA->set(i, a, FockA->get(i+nfrzc, a+noccA));
            }
        }

        // t(I,A) += \sum_(E) t(I,E) * F(A,E)
        t1newA->gemm(false, true, t1A, FabA, 1.0, 0.0);

        // t(I,A) -= \sum_(M) t(M,A) * F(M,I)
        t1newA->gemm(true, false, FijA, t1A, -1.0, 1.0);

        // t(I,A) += \sum_(M,E) T(IM,AE) * F(M,E)
        T = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        T->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        X = std::make_shared<Tensor2d>("Tx (IA|JB)", naoccA, navirA, naoccA, navirA);
        X->sort(1324, T, 1.0, 0.0);
        T.reset();
        t1newA->gemv(false, X, FiaA, 1.0, 1.0);
        X.reset();

        // t(I,A) += \sum_(m,e) T(Im,Ae) * F(m,e)
        T = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T->read(psio_, PSIF_DFOCC_AMPS);
        X = std::make_shared<Tensor2d>("Tx (IA|jb)", naoccA, navirA, naoccB, navirB);
        X->sort(1324, T, 1.0, 0.0);
        T.reset();
        t1newA->gemv(false, X, FiaB, 1.0, 1.0);
        X.reset();

        // t(I,A) += \sum_(Q) tQ * b(Q,AI)
        t1newA->gemv(true, bQiaA, T1c, 1.0, 1.0);

        // t(I,A) -= \sum_(Q,M) [T(Q,MA) + t(Q,MA)] b(Q,MI) -= \sum_(Q,M) b(QM,I) X(QM,A)
        T = std::make_shared<Tensor2d>("T2 (Q|IA)", nQ, naoccA, navirA);
        T->read(psio_, PSIF_DFOCC_AMPS);
        K = std::make_shared<Tensor2d>("T1 (Q|IA)", nQ, naoccA, navirA);
        K->read(psio_, PSIF_DFOCC_AMPS);
        X = std::make_shared<Tensor2d>("T1+T2 (Q|IA)", nQ, naoccA, navirA);
        X->copy(T);
        T.reset();
        X->add(K);
        K.reset();

        t1newA->contract(true, false, naoccA, navirA, nQ * naoccA, bQijA, X, -1.0, 1.0);
        X.reset();

        // t(I,A) += \sum_(Q,E) T(Q,IE) * b(Q,AE)
        T = std::make_shared<Tensor2d>("T2 (Q|IA)", nQ, naoccA, navirA);
        T->read(psio_, PSIF_DFOCC_AMPS);
        Y = std::make_shared<Tensor2d>("Temp (Q|AI)", nQ, navirA, naoccA);
        Y->swap_3index_col(T);
        T.reset();

        t1newA->contract(true, false, naoccA, navirA, nQ*navirA, Y, bQabA, 1.0, 1.0);
        Y.reset();

        // Denom
        for(int i=0; i<naoccA; ++i) {
             for(int a=0; a<navirA; ++a) {
                 double value = FockA->get(i+nfrzc, i+nfrzc) - FockA->get(a+noccA, a+noccA);
                 t1newA->set(i, a, t1newA->get(i,a)/value);
             }
        }
        //t1newA->print();

        // Beta Part
        // t(i,a) * D(i,a) = f(i,a)
        // D(i,a) = f(i,i) - f(a,a)
        for(int i=0; i<naoccB; ++i) {
            for(int a=0; a<navirB; ++a) {
               t1newB->set(i, a, FockB->get(i+nfrzc, a+noccB));
            }
        }

        // t(i,a) += \sum_(e) t(i,e) * F(a,e)
        t1newB->gemm(false, true, t1B, FabB, 1.0, 1.0);

        // t(i,a) -= \sum_(m) t(m,a) * F(m,i)
        t1newB->gemm(true, false, FijB, t1B, -1.0, 1.0);

        // t(i,a) += \sum_(m,e) T(im,ae) * F(m,e)
        T = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        T->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        X = std::make_shared<Tensor2d>("Tx (ia|jb)", naoccB, navirB, naoccB, navirB);
        X->sort(1324, T, 1.0, 0.0);
        T.reset();
        t1newB->gemv(false, X, FiaB, 1.0, 1.0);
        X.reset();

        // t(i,a) += \sum_(M,E) T(Mi,Ea) * F(M,E)
        T = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T->read(psio_, PSIF_DFOCC_AMPS);
        X = std::make_shared<Tensor2d>("Tx (ia|JB)", naoccB, navirB, naoccA, navirA);
        X->sort(2413, T, 1.0, 0.0);
        T.reset();
        t1newB->gemv(false, X, FiaA, 1.0, 1.0);
        X.reset();

        // t(i,a) += \sum_(Q) tQ * b(Q,ai)
        t1newB->gemv(true, bQiaB, T1c, 1.0, 1.0);

        // t(i,a) -= \sum_(Q,m) [T(Q,ma) + t(Q,ma)] b(Q,mi) -= \sum_(Q,m) b(Qm,i) X(Qm,a)
        T = std::make_shared<Tensor2d>("T2 (Q|ia)", nQ, naoccB, navirB);
        T->read(psio_, PSIF_DFOCC_AMPS);
        K = std::make_shared<Tensor2d>("T1 (Q|ia)", nQ, naoccB, navirB);
        K->read(psio_, PSIF_DFOCC_AMPS);
        X = std::make_shared<Tensor2d>("T1+T2 (Q|ia)", nQ, naoccB, navirB);
        X->copy(T);
        T.reset();
        X->add(K);
        K.reset();

        t1newB->contract(true, false, naoccB, navirB, nQ*naoccB, bQijB, X, -1.0, 1.0);
        X.reset();

        // t(i,a) += \sum_(Q,e) T(Q,ie) * b(Q,ae)
        T = std::make_shared<Tensor2d>("T2 (Q|ia)", nQ, naoccB, navirB);
        T->read(psio_, PSIF_DFOCC_AMPS);
        Y = std::make_shared<Tensor2d>("Temp (Q|ai)", nQ, navirB, naoccB);
        Y->swap_3index_col(T);
        T.reset();

        t1newB->contract(true, false, naoccB, navirB, nQ*navirB, Y, bQabB, 1.0, 1.0);
        Y.reset();

        // Denom
        for(int i=0; i<naoccB; ++i) {
             for(int a=0; a<navirB; ++a) {
                 double value = FockB->get(i+nfrzc, i+nfrzc) - FockB->get(a+noccB, a+noccB);
                 t1newB->set(i, a, t1newB->get(i,a)/value);
             }
        }


    }  // else if (reference_ == "UNRESTRICTED")

}  // end ccsd_t1_amps
}  // namespace dfoccwave
}  // namespace psi
