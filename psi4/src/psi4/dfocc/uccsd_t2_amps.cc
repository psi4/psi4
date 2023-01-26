/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

#include "psi4/libmints/matrix.h"
#include "psi4/libdiis/diismanager.h"
#include "dfocc.h"

namespace psi{
namespace dfoccwave {

void DFOCC::uccsd_t2AA_amps()
{
    SharedTensor2d X, Y, K, L, T, T2, T2new, G, J, W, U, Tau, R;

    // t(IJ,AB) = <IJ||AB>
    J = std::make_shared<Tensor2d>("J (IA|JB)", naoccA, navirA, naoccA, navirA);
    J->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    G = std::make_shared<Tensor2d>("G <IJ||AB>", naoccA, naoccA, navirA, navirA);
    G->sort(1324, J, 1.0, 0.0);
    G->sort(1342, J, -1.0, 1.0);
    J.reset();

    // New T2AA
    T2new = std::make_shared<Tensor2d>("New T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2new->copy(G);
    G.reset();

    // t(IJ,AB) += P_(AB) \sum_(E) t(IJ,AE) * Ft_(B,E)
    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X <IJ|AB>", naoccA, naoccA, navirA, navirA);
    X->contract(false, true, naoccA * naoccA * navirA, navirA, navirA, T2, FtabA, 1.0, 0.0);
    T2.reset();
    T2new->add(X);
    T2new->sort(1243, X, -1.0, 1.0);
    X.reset();

    // t(IJ,AB) -= P_(IJ) \sum_(M) t(IM,AB) * Ft_(M,J)
    // t(IJ,AB) -= P_(IJ) \sum_(M) t(MJ,AB) * Ft_(M,I)
    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Y = std::make_shared<Tensor2d>("T2 <MJ|AB>", naoccA, naoccA, navirA, navirA);
    Y->contract(true, false, naoccA, navirA * navirA * naoccA, naoccA, FtijA, T2, -1.0, 0.0);
    T2.reset();
    T2new->sort(2134, Y, -1.0, 1.0);
    T2new->axpy(Y, 1.0);
    Y.reset();

    // t(IJ,AB) += P_(IJ)P_(AB) \sum_(M,E) T(IM,AE) * W_MBEJ
    // X(IA,JB) = \sum(ME) T(IA,ME) W(ME,JB)
    // T(IJ,AB) <= P_(IJ) * P_(AB) X(IA,JB)
    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    U = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    U->sort(1324, T2, 1.0, 0.0);
    T2.reset();
    W = std::make_shared<Tensor2d>("W (ME|JB)", naoccA, navirA, naoccA, navirA);
    W->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X1 (IA|JB)", naoccA, navirA, naoccA, navirA);
    X->gemm(false, false, U, W, 1.0, 0.0);
    U.reset();
    W.reset();
    //T2new->P_ijab(X);
    T2new->sort(1324, X, 1.0, 1.0);
    T2new->sort(3124, X, -1.0, 1.0);
    T2new->sort(1342, X, -1.0, 1.0);
    T2new->sort(3142, X, 1.0, 1.0);
    X.reset();

    // t(IJ,AB) += P_(IJ)P_(AB) \sum_(m,e) T(Im,Ae) * W_mBeJ
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    U = std::make_shared<Tensor2d>("T2 (IA|jb)", naoccA, navirA, naoccB, navirB);
    U->sort(1324, T2, 1.0, 0.0);
    T2.reset();
    X = std::make_shared<Tensor2d>("X2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    W = std::make_shared<Tensor2d>("W (me|JB)", naoccB, navirB, naoccA, navirA);
    W->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(false, false, U, W, 1.0, 0.0);
    U.reset();
    W.reset();
    //T2new->P_ijab(X);
    T2new->sort(1324, X, 1.0, 1.0);
    T2new->sort(3124, X, -1.0, 1.0);
    T2new->sort(1342, X, -1.0, 1.0);
    T2new->sort(3142, X, 1.0, 1.0);
    X.reset();


    // t(IJ,AB) = P_(IJ)P_(AB) \sum_(Q) [t(Q,AJ) * t(Q,IB)] + [t(Qp,IA) * b(Q,JB)] 
    // X(IA,JB) = t(Qp,IA) * b(Q,JB)
    X = std::make_shared<Tensor2d>("X3 (IA|JB)", naoccA, navirA, naoccA, navirA);
    T = std::make_shared<Tensor2d>("T1p (Q|IA)", nQ, naoccA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(true, false, T, bQiaA, 1.0, 0.0);
    T.reset();

    // X(IA,JB) = t(Q,AI) * t(Q,JB)
    T = std::make_shared<Tensor2d>("T1 (Q|AI)", nQ, navirA, naoccA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("Temp T1 (Q|IA)", nQ, naoccA, navirA);
    K->swap_3index_col(T);
    T.reset();
    L = std::make_shared<Tensor2d>("T1 (Q|IA)", nQ, naoccA, navirA);
    L->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(true, false, K, L, -1.0, 1.0);   ///////////////////////////////////////////////////
    K.reset();
    L.reset();
    //T2new->P_ijab(X);
    T2new->sort(1324, X, 1.0, 1.0);
    T2new->sort(3124, X, -1.0, 1.0);
    T2new->sort(1342, X, -1.0, 1.0);
    T2new->sort(3142, X, 1.0, 1.0);
    X.reset();

    // Write and close
    T2new->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T2new.reset();

    // WmnijT2
    uccsd_W_MNIJT2AA();

    // WabefT2
    ccsd_Wabef2T2AA();

    // Denom
    T2new = std::make_shared<Tensor2d>("New T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2new->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T2new->apply_denom(nfrzc, noccA, FockA);
    T2new->write_anti_symm(psio_, PSIF_DFOCC_AMPS);

    // Error vector
    T = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    R = std::make_shared<Tensor2d>("RT2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    R->copy(T2new);
    T2new.reset();
    R->subtract(T);
    T.reset();
    R->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
}// end ccsd_t2AA_amps


void DFOCC::uccsd_t2BB_amps()
{
    SharedTensor2d X, Y, K, L, T, T2, T2new, G, J, W, U, Tau, R;

    // t(ij,ab) = <ij||ab>
    J = std::make_shared<Tensor2d>("J (ia|jb)", naoccB, navirB, naoccB, navirB);
    J->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    G = std::make_shared<Tensor2d>("G <ij||ab>", naoccB, naoccB, navirB, navirB);
    G->sort(1324, J, 1.0, 0.0);
    G->sort(1342, J, -1.0, 1.0);
    J.reset();

    // New T2BB
    T2new = std::make_shared<Tensor2d>("New T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2new->copy(G);
    G.reset();

    // t(ij,ab) += P_(ab) \sum_(e) t(ij,ae) * Ft_(b,e)
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X <ij|ab>", naoccB, naoccB, navirB, navirB);
    X->contract(false, true, naoccB * naoccB * navirB, navirB, navirB, T2, FtabB, 1.0, 0.0);
    T2.reset();
    T2new->axpy(X, 1.0);
    T2new->sort(1243, X, -1.0, 1.0);
    X.reset();

    // t(ij,ab) -= P_(ij) \sum_(m) t(im,ab) * Ft_(m,j)
    // t(ij,ab) -= P_(ij) \sum_(m) t(mj,ab) * Ft_(m,i)
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Y = std::make_shared<Tensor2d>("T2 <mj|ab>", naoccB, naoccB, navirB, navirB);
    Y->contract(true, false, naoccB, navirB * navirB * naoccB, naoccB, FtijB, T2, -1.0, 0.0);
    T2.reset();
    T2new->sort(2134, Y, -1.0, 1.0);
    T2new->axpy(Y, 1.0);
    Y.reset();

    // t(ij,ab) += P_(ij) * P_(ab) \sum_(m,e) T(im,ae) * W_mbej
    // X(ia,jb) = \sum(m,e) T(ia,me) W(me,jb)
    // T(ij,ab) <= P_(ij) * P_(ab) X(ia,jb)
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    U = std::make_shared<Tensor2d>("T2 (ia|jb)", naoccB, navirB, naoccB, navirB);
    U->sort(1324, T2, 1.0, 0.0);
    T2.reset();
    W = std::make_shared<Tensor2d>("W (me|jb)", naoccB, navirB, naoccB, navirB);
    W->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X1 (ia|jb)", naoccB, navirB, naoccB, navirB);
    X->gemm(false, false, U, W, 1.0, 0.0);
    U.reset();
    W.reset();
    //T2new->P_ijab(X);
    T2new->sort(1324, X, 1.0, 1.0);
    T2new->sort(3124, X, -1.0, 1.0);
    T2new->sort(1342, X, -1.0, 1.0);
    T2new->sort(3142, X, 1.0, 1.0);
    X.reset();

    // t(ij,ab) += P_(ij) * P_(ab) \sum_(m,e) T(Mi,Ea) * W_MbEj
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    U = std::make_shared<Tensor2d>("T2 (IA|jb)", naoccA, navirA, naoccB, navirB);
    U->sort(1324, T2, 1.0, 0.0);
    T2.reset();
    X = std::make_shared<Tensor2d>("X2 (ia|jb)", naoccB, navirB, naoccB, navirB);
    W = std::make_shared<Tensor2d>("W (ME|jb)", naoccA, navirA, naoccB, navirB);
    W->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(true, false, U, W, 1.0, 0.0);
    U.reset();
    W.reset();
    //T2new->P_ijab(X);
    T2new->sort(1324, X, 1.0, 1.0);
    T2new->sort(3124, X, -1.0, 1.0);
    T2new->sort(1342, X, -1.0, 1.0);
    T2new->sort(3142, X, 1.0, 1.0);
    X.reset();

    // t(ij,ab) = P_(ij) P_(ab) \sum_(Q) [t(Q,aj) * t(Q,ib)] + [t(Qp,ia) * b(Q,jb)] 
    X = std::make_shared<Tensor2d>("X3 (ia|jb)", naoccB, navirB, naoccB, navirB);

    // X(ia,jb) = t(Qp,ia) * b(Q,jb)
    T = std::make_shared<Tensor2d>("T1p (Q|ia)", nQ, naoccB, navirB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(true, false, T, bQiaB, 1.0, 0.0);
    T.reset();

    // X(ia,jb) = t(Q,ai) * t(Q,jb)
    T = std::make_shared<Tensor2d>("T1 (Q|ai)", nQ, navirB, naoccB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("Temp T1 (Q|ia)", nQ, naoccB, navirB);
    K->swap_3index_col(T);
    T.reset();
    L = std::make_shared<Tensor2d>("T1 (Q|ia)", nQ, naoccB, navirB);
    L->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(true, false, K, L, -1.0, 1.0);
    K.reset();
    T.reset();
    //T2new->P_ijab(X);
    T2new->sort(1324, X, 1.0, 1.0);
    T2new->sort(3124, X, -1.0, 1.0);
    T2new->sort(1342, X, -1.0, 1.0);
    T2new->sort(3142, X, 1.0, 1.0);
    X.reset();

    // Write and close
    T2new->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T2new.reset();

    // WmnijT2
    uccsd_W_mnijT2BB();

    // WabefT2
    ccsd_Wabef2T2BB();

    // Denom
    T2new = std::make_shared<Tensor2d>("New T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2new->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T2new->apply_denom(nfrzc, noccB, FockB);
    T2new->write_anti_symm(psio_, PSIF_DFOCC_AMPS);

    // Error vector
    T = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    R = std::make_shared<Tensor2d>("RT2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    R->copy(T2new);
    T2new.reset();
    R->subtract(T);
    T.reset();
    R->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
}// End ccsd_t2BB_amps


void DFOCC::uccsd_t2AB_amps()
{
    SharedTensor2d X, Y, K, L, T, T2, T2new, G, J, W, U, Tau, R;

    // t(Ij,Ab) = <Ij|Ab>
    J = std::make_shared<Tensor2d>("J (IA|jb)", naoccA, navirA, naoccB, navirB);
    J->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    T2new = std::make_shared<Tensor2d>("New T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2new->sort(1324, J, 1.0, 0.0);
    J.reset();

    // Read T2
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);

    // T(Ij,Ab) += \sum_(e) T(Ij,Ae) Ft_(b,e)
    T2new->contract424(4, 2, T2, FtabB, 1.0, 1.0);
    //T2new->contract(false, true, naoccA*naoccB*navirA, navirB, navirB, T2, FtabB, 1.0, 1.0);

    // T(Ij,Ab) += \sum_(E) T(Ij,Eb) Ft_(A,E)
    T2new->contract424(3, 2, T2, FtabA, 1.0, 1.0);

    // T(Ij,Ab) -= \sum_(m) T(Im,Ab) Ft_(m,j)
    T2new->contract424(2, 1, T2, FtijB, -1.0, 1.0);

    // T(Ij,Ab) -= \sum_(M) T(Mj,Ab) Ft_(M,I)
    T2new->contract424(1, 1, T2, FtijA, -1.0, 1.0);
    T2.reset();

    // t(Ij,Ab) += \sum_(Q) t(Qp,IA) * b(Q,jb) 
    X = std::make_shared<Tensor2d>("X (IA|jb)", naoccA, navirA, naoccB, navirB);
    T = std::make_shared<Tensor2d>("T1p (Q|IA)", nQ, naoccA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(true, false, T, bQiaB, 1.0, 0.0);
    T.reset();

    // t(Ij,Ab) -= \sum_(Q) t(Q,AI) * t(Q,jb) 
    T = std::make_shared<Tensor2d>("T1 (Q|AI)", nQ, navirA, naoccA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("Temp T1 (Q|IA)", nQ, naoccA, navirA);
    K->swap_3index_col(T);
    T.reset();
    L = std::make_shared<Tensor2d>("T1 (Q|ia)", nQ, navirB, naoccB);
    L->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(true, false, K, L, -1.0, 1.0);
    K.reset();
    L.reset();

    // t(Ij,Ab) -= \sum_(Q) t(Q,bj) * t(Q,IA) 
    T = std::make_shared<Tensor2d>("T1 (Q|ai)", nQ, navirB, naoccB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("Temp T1 (Q|ia)", nQ, naoccB, navirB);
    K->swap_3index_col(T);
    T.reset();
    L = std::make_shared<Tensor2d>("T1 (Q|IA)", nQ, naoccA, navirA);
    L->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(true, false, L, K, -1.0, 1.0);
    K.reset();
    L.reset();

    // t(Ij,Ab) += \sum_(Q) t(Qp,jb) * b(Q,IA) 
    T = std::make_shared<Tensor2d>("T1p (Q|ia)", nQ, naoccB, navirB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(true, false, bQiaA, T, 1.0, 1.0);
    T.reset();

    T2new->sort(1324, X, 1.0, 1.0);
    X.reset();

    //T2new->write(psio_, PSIF_DFOCC_AMPS);
    //T2new.reset();

///WmbejT2AB Terms
    // T(ME,jb) <= T(Mj,Eb)
    // X(IA,jb) += \sum(ME) T(Mj,Eb) W(MA,EI) = \sum(ME) T(ME,jb) W(ME,IA)
    // T(Ij,Ab) <= X(IA,jb)
    W = std::make_shared<Tensor2d>("W (ME|JB)", naoccA, navirA, naoccA, navirA);
    W->read(psio_, PSIF_DFOCC_AMPS);

    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("T2 (IA|jb)", naoccA, navirA, naoccB, navirB);
    K->sort(1324, T2, 1.0, 0.0);
    T2.reset();

    X = std::make_shared<Tensor2d>("X1 (IA|jb)", naoccA, navirA, naoccB, navirB);
    X->gemm(true, false, W, K, 1.0, 0.0);
    K.reset();
    W.reset();
    //T2new = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    //T2new->read(psio_, PSIF_DFOCC_AMPS);
    T2new->sort(1324, X, 1.0, 1.0);
    X.reset();
    //T2new->write(psio_, PSIF_DFOCC_AMPS);
    //T2new.reset();

    // X(IA,jb) = \sum(me) T(Im,Ae) W(mb,ej) = \sum(me) T(IA,me) W(me,jb)
    // T(Ij,Ab) <= X(IA,jb)
    W = std::make_shared<Tensor2d>("W (me|jb)", naoccB, navirB, naoccB, navirB);
    W->read(psio_, PSIF_DFOCC_AMPS);

    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("T2 (IA|jb)", naoccA, navirA, naoccB, navirB);
    K->sort(1324, T2, 1.0, 0.0);
    T2.reset();

    X = std::make_shared<Tensor2d>("X2 (IA|jb)", naoccA, navirA, naoccB, navirB);
    X->gemm(false, false, K, W, 1.0, 0.0);
    K.reset();
    W.reset();
    //T2new = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    //T2new->read(psio_, PSIF_DFOCC_AMPS);
    T2new->sort(1324, X, 1.0, 1.0);
    X.reset();
    //T2new->write(psio_, PSIF_DFOCC_AMPS);
    //T2new.reset();

    // X(IA,jb) = \sum(ME) T(IM,AE) W(Mb,Ej) = \sum(ME) T(IA,ME) W(ME,jb)
    // T(Ij,Ab) <= X(IA,jb)
    W = std::make_shared<Tensor2d>("W (ME|jb)", naoccA, navirA, naoccB, navirB);
    W->read(psio_, PSIF_DFOCC_AMPS);

    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    K->sort(1324, T2, 1.0, 0.0);
    T2.reset();

    X = std::make_shared<Tensor2d>("X3 (IA|jb)", naoccA, navirA, naoccB, navirB);
    X->gemm(false, false, K, W, 1.0, 0.0);
    K.reset();
    W.reset();
    //T2new = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    //T2new->read(psio_, PSIF_DFOCC_AMPS);
    T2new->sort(1324, X, 1.0, 1.0);
    X.reset();
    //T2new->write(psio_, PSIF_DFOCC_AMPS);
    //T2new.reset();

    // X(IA,jb) = \sum(me) T(jm,be) W(mA,eI) = \sum(me) T(me,jb) W(IA,me)
    // T(Ij,Ab) <= X(IA,jb)
    W = std::make_shared<Tensor2d>("W (me|JB)", naoccB, navirB, naoccA, navirA);
    W->read(psio_, PSIF_DFOCC_AMPS);
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("T2 (ia|jb)", naoccB, navirB, naoccB, navirB);
    K->sort(2413, T2, 1.0, 0.0);
    T2.reset();

    X = std::make_shared<Tensor2d>("X4 (IA|jb)", naoccA, navirA, naoccB, navirB);
    X->gemm(true, false, W, K, 1.0, 0.0);
    K.reset();
    W.reset();
    //T2new = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    //T2new->read(psio_, PSIF_DFOCC_AMPS);
    T2new->sort(1324, X, 1.0, 1.0);
    X.reset();
    //T2new->write(psio_, PSIF_DFOCC_AMPS);
    //T2new.reset();

    // X(Ib,jA) = \sum(mE) T(Im,Eb) W(mA,Ej) = \sum(mE) T'(Ib,mE) W(mE,jA)
    // T(Ij,Ab) <= X(Ib,jA)
    X = std::make_shared<Tensor2d>("X5 (Ib|jA)", naoccA, navirB, naoccB, navirA);
    W = std::make_shared<Tensor2d>("W (mE|jB)", naoccB, navirA, naoccB, navirA);
    W->read(psio_, PSIF_DFOCC_AMPS);

    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("T2 (Ib|jA)", naoccA, navirB, naoccB, navirA);
    K->sort(1423, T2, 1.0, 0.0);
    T2.reset(); 

    X->gemm(false, false, K, W, 1.0, 0.0);
    K.reset();
    W.reset();
    //T2new = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    //T2new->read(psio_, PSIF_DFOCC_AMPS);
    T2new->sort(1342, X, 1.0, 1.0);
    X.reset();
    //T2new->write(psio_, PSIF_DFOCC_AMPS);
    //T2new.reset();

    // X(Ib,jA) = \sum(Me) T(Mj,Ae) W(Mb,eI) = \sum(me) T'(Me,jA) W(Me,Ib)
    // T(Ij,Ab) <= X(Ib,jA)
    W = std::make_shared<Tensor2d>("W (Me|Jb)", naoccA, navirB, naoccA, navirB);
    W->read(psio_, PSIF_DFOCC_AMPS);

    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("T2 (Ib|jA)", naoccA, navirB, naoccB, navirA);
    K->sort(1423, T2, 1.0, 0.0);
    T2.reset();

    X = std::make_shared<Tensor2d>("X6 (Ib|jA)", naoccA, navirB, naoccB, navirA);
    X->gemm(true, false, W, K, 1.0, 0.0);
    K.reset();
    W.reset();
    //T2new = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    //T2new->read(psio_, PSIF_DFOCC_AMPS);
    T2new->sort(1342, X, 1.0, 1.0);
    X.reset();
    T2new->write(psio_, PSIF_DFOCC_AMPS);
    T2new.reset();

    // WmnijT2AB
    uccsd_W_MnIjT2AB();

    // WabefT2AB
    ccsd_Wabef2T2AB();

    // Denom
    T2new = std::make_shared<Tensor2d>("New T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2new->read(psio_, PSIF_DFOCC_AMPS);
    T2new->apply_denom_os(nfrzc, noccA, noccB, FockA, FockB);
    T2new->write(psio_, PSIF_DFOCC_AMPS);

    // Error vector
    T = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    R = std::make_shared<Tensor2d>("RT2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    R->copy(T2new);
    T2new.reset();
    R->subtract(T);
    T.reset();
    R->write(psio_, PSIF_DFOCC_AMPS);
    R.reset();


    // Error vectors of t1A and t1B
    SharedTensor2d RA = std::make_shared<Tensor2d>("RT1 <I|A>", naoccA, navirA);
    RA->copy(t1newA);
    RA->subtract(t1A);

    SharedTensor2d RB = std::make_shared<Tensor2d>("RT1 <i|a>", naoccB, navirB);
    RB->copy(t1newB);
    RB->subtract(t1B);

    // DIIS
    SharedTensor2d RAA, TAA, RBB, TBB, RAB, TAB, TA, TB;

        //if (orb_opt_ == "FALSE" || mo_optimized == 1) {
            // RAA
            RAA = std::make_shared<Tensor2d>("RT2 <IJ|AB>", ntri_anti_ijAA, ntri_anti_abAA);
            RAA->read(psio_, PSIF_DFOCC_AMPS);
            auto RT2AA = std::make_shared<Matrix>("RT2AA", ntri_anti_ijAA, ntri_anti_abAA);
            RAA->to_matrix(RT2AA);
            RAA.reset();
            // TAA
            TAA = std::make_shared<Tensor2d>("New T2 <IJ|AB>", ntri_anti_ijAA, ntri_anti_abAA);
            TAA->read(psio_, PSIF_DFOCC_AMPS);
            auto T2AA = std::make_shared<Matrix>("T2AA", ntri_anti_ijAA, ntri_anti_abAA);
            TAA->to_matrix(T2AA);

            // RBB
            RBB = std::make_shared<Tensor2d>("RT2 <ij|ab>", ntri_anti_ijBB, ntri_anti_abBB);
            RBB->read(psio_, PSIF_DFOCC_AMPS);
            auto RT2BB = std::make_shared<Matrix>("RT2BB", ntri_anti_ijBB, ntri_anti_abBB);
            RBB->to_matrix(RT2BB);
            RBB.reset();
            // TBB
            TBB = std::make_shared<Tensor2d>("New T2 <ij|ab>", ntri_anti_ijBB, ntri_anti_abBB);
            TBB->read(psio_, PSIF_DFOCC_AMPS);
            auto T2BB = std::make_shared<Matrix>("T2BB", ntri_anti_ijBB, ntri_anti_abBB);
            TBB->to_matrix(T2BB);

            // RAB
            RAB = std::make_shared<Tensor2d>("RT2 <Ij|Ab>", naoccA * naoccB, navirA * navirB);
            RAB->read(psio_, PSIF_DFOCC_AMPS);
            auto RT2AB = std::make_shared<Matrix>("RT2AB", naoccA * naoccB, navirA * navirB);
            RAB->to_matrix(RT2AB);
            RAB.reset();
            // TAB
            TAB = std::make_shared<Tensor2d>("New T2 <Ij|Ab>", naoccA * naoccB, navirA * navirB);
            TAB->read(psio_, PSIF_DFOCC_AMPS);
            auto T2AB = std::make_shared<Matrix>("T2AB", naoccA * naoccB, navirA * navirB);
            TAB->to_matrix(T2AB);

            // T1A
            auto RT1A = std::make_shared<Matrix>("RT1A", naoccA, navirA);
            RA->to_matrix(RT1A);
            RA.reset();
            auto T1A = std::make_shared<Matrix>("T1A", naoccA, navirA);
            t1A->to_matrix(T1A);

            // T1B
            auto RT1B = std::make_shared<Matrix>("RT1B", naoccB, navirB);
            RB->to_matrix(RT1B);
            RB.reset();
            auto T1B = std::make_shared<Matrix>("T1B", naoccB, navirB);
            t1B->to_matrix(T1B);


            // add entry
            if (do_diis_ == 1) 
                ccsdDiisManager->add_entry(RT2AA.get(), RT2BB.get(), RT2AB.get(), RT1A.get(),  RT1B.get(), T2AA.get(), T2BB.get(),
                                           T2AB.get(), T1A.get(), T1B.get());
            RT2AA.reset();
            RT2BB.reset();
            RT2BB.reset();
            RT1A.reset();
            RT1B.reset(); 

            // extrapolate
            if (do_diis_ == 1) {
                if (ccsdDiisManager->subspace_size() >= cc_mindiis_)
                    ccsdDiisManager->extrapolate(T2AA.get(), T2BB.get(), T2AB.get(), T1A.get(), T1B.get());
                TAA->set2(T2AA);
                TBB->set2(T2BB);
                TAB->set2(T2AB);
                t1A->set2(T1A);
                t1B->set2(T1B);
  
            }
            T2AA.reset();
            T2BB.reset();
            T2AB.reset();
            T1A.reset();
            T1B.reset();              

            TAA->write(psio_, PSIF_DFOCC_AMPS);
            TBB->write(psio_, PSIF_DFOCC_AMPS);
            TAB->write(psio_, PSIF_DFOCC_AMPS);
            TAA.reset();
            TBB.reset();
            TAB.reset();
       // }  // end if diis

}// end ccsd_t2AB_amps

}  // namespace dfoccwave
}  // namespace psi
