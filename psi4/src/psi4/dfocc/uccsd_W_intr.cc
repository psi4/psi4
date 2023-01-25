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

#include "dfocc.h"

using namespace psi;
using namespace std;

namespace psi {
namespace dfoccwave {

// ccsd_W_MNIJT2AA
void DFOCC::uccsd_W_MNIJT2AA()
{
    SharedTensor2d J, W, I, X, Y, T;
    SharedTensor2d T2, Tau, T2new;

    // W_MNIJ Alpha Block
    // W_MNIJ = <MN||IJ> + P_(MN) P_(IJ) \sum_(Q) t_IM^Q b_JN^Q + \sum_(EF) Tau(IJ,EF) * <MN|EF> 
    // W_MNIJ = <MN||IJ> 
    J = std::make_shared<Tensor2d>("J (IM|JN)", naoccA, naoccA, naoccA, naoccA);
    J->gemm(true, false, bQijA, bQijA, 1.0, 0.0);
    W = std::make_shared<Tensor2d>("W <MN|IJ>", naoccA, naoccA, naoccA, naoccA);
    W->sort(1324, J, 1.0, 0.0);
    W->sort(1342, J, -1.0, 1.0);
    J.reset();

    // W_MNIJ += P_(MN) P_(IJ) \sum_(Q) t_IM^Q b_JN^Q
    // X(IM,JN) = \sum(Q) t(Q,IM) b(Q,JN)
    // W_MNIJ += P_(MN) * P_(IJ) X(IM,JN)
    T = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (IM|JN)", naoccA, naoccA, naoccA, naoccA);
    X->gemm(true, false, T, bQijA, 1.0, 0.0);
    T.reset();
    //W->P_ijab(X);
    W->sort(2413, X, 1.0, 1.0);
    W->sort(4213, X, -1.0, 1.0);
    W->sort(2431, X, -1.0, 1.0);
    W->sort(4231, X, 1.0, 1.0);
    X.reset();

    //W_MNIJ += \sum_(EF) Tau(IJ,EF) * <MN|EF> 
    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA);
    uccsd_tau_amps(naoccA, naoccA, navirA, navirA, Tau, T2, t1A, t1A);
    T2.reset();
    J = std::make_shared<Tensor2d>("J (ME|NF)", naoccA, navirA, naoccA, navirA);
    J->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    I = std::make_shared<Tensor2d>("I <MN|EF>", naoccA, naoccA, navirA, navirA);
    I->sort(1324, J, 1.0, 0.0);
    J.reset();
    W->gemm(false, true, I, Tau, 1.0, 1.0);
    I.reset();
    // 0.5 * Tau * WMNIJ
    T2new = std::make_shared<Tensor2d>("New T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2new->read_anti_symm(psio_, PSIF_DFOCC_AMPS);  // ccsd_t2_amps() fonksiyonundan okuyor
    T2new->gemm(true, false, W, Tau, 0.5, 1.0);
    W.reset();
    Tau.reset();
    T2new->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
}// ccsd_W_MNIJT2AA

// ccsd_W_mnijT2BB
void DFOCC::uccsd_W_mnijT2BB()
{
    SharedTensor2d J, W, I, X, Y, T;
    SharedTensor2d T2, Tau, T2new;

    // W_mnij Beta Block
    // W_mnij = <mn||ij> + P_(mn) P_(ij) \sum_(Q) t_im^Q b_jn^Q + \sum_(ef) Tau(ij,ef) * <mn|ef> 
    // W_mnij = <mn||ij> 
    J = std::make_shared<Tensor2d>("J (im|jn)", naoccB, naoccB, naoccB, naoccB);
    J->gemm(true, false, bQijB, bQijB, 1.0, 0.0);
    W = std::make_shared<Tensor2d>("W <mn|ij>", naoccB, naoccB, naoccB, naoccB);
    W->sort(1324, J, 1.0, 0.0);
    W->sort(1342, J, -1.0, 1.0);
    J.reset();

    // W_mnij += P_(mn) P_(ij) \sum_(Q) t_im^Q b_jn^Q
    // X(im,jn) = \sum(Q) t(Q,im) b(Q,jn)
    // W_mnij += P_(mn) * P_(ij) X(im,jn)
    T = std::make_shared<Tensor2d>("T1 (Q|ij)", nQ, naoccB, naoccB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (im|jn)", naoccB, naoccB, naoccB, naoccB);
    X->gemm(true, false, T, bQijB, 1.0, 0.0);
    T.reset();
    //W->P_ijab(X);
    W->sort(2413, X, 1.0, 1.0);
    W->sort(4213, X, -1.0, 1.0);
    W->sort(2431, X, -1.0, 1.0);
    W->sort(4231, X, 1.0, 1.0);
    X.reset();

    //W_mnij += \sum_(ef) Tau(ij,ef) * <mn|ef> 
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <ij|ab>", naoccB, naoccB, navirB, navirB);
    uccsd_tau_amps(naoccB, naoccB, navirB, navirB, Tau, T2, t1B, t1B);
    T2.reset();
    J = std::make_shared<Tensor2d>("J (me|nf)", naoccB, navirB, naoccB, navirB);
    J->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    I = std::make_shared<Tensor2d>("I <mn|ef>", naoccB, naoccB, navirB, navirB);
    I->sort(1324, J, 1.0, 0.0);
    J.reset();
    W->gemm(false, true, I, Tau, 1.0, 1.0);
    I.reset();

    // 0.5 * Tau * Wmnij
    T2new = std::make_shared<Tensor2d>("New T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2new->read_anti_symm(psio_, PSIF_DFOCC_AMPS);  // ccsd_t2_amps() fonksiyonundan okuyor
    T2new->gemm(true, false, W, Tau, 0.5, 1.0);
    W.reset();
    Tau.reset();
    T2new->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
}// ccsd_W_mnijT2BB

// ccsd_W_MnIjT2AB
void DFOCC::uccsd_W_MnIjT2AB()
{
    SharedTensor2d J, W, I, X, Y, T;
    SharedTensor2d T2, Tau, T2new;

    // W_MnIj Alpha-Beta Block
    // W_MnIj = <Mn|Ij> + \sum_(Q) t_IM^Q b_jn^Q + \sum_(Ef) Tau(Ij,Ef) * <Mn|Ef> 
    // W_MnIj = <Mn|Ij> 
    W = std::make_shared<Tensor2d>("W <Mn|Ij>", naoccA, naoccB, naoccA, naoccB);
    J = std::make_shared<Tensor2d>("J (MI|nj)", naoccA, naoccA, naoccB, naoccB);
    J->gemm(true, false, bQijA, bQijB, 1.0, 0.0);
    W->sort(1324, J, 1.0, 0.0);
    J.reset();

    // X(IM,jn) = \sum(Q) t(Q,IM) b(Q,jn)
    // W_MnIj += X(IM,jn)
    X = std::make_shared<Tensor2d>("X (IM|jn)", naoccA, naoccA, naoccB, naoccB);
    T = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(true, false, T, bQijB, 1.0, 0.0);
    T.reset();
    W->sort(2413, X, 1.0, 1.0);
    X.reset();

    // X(IM,jn) = \sum(Q) t(Q,jn) b(Q,IM)
    // W_MnIj += X(IM,jn)
    X = std::make_shared<Tensor2d>("X (IM|jn)", naoccA, naoccA, naoccB, naoccB);
    T = std::make_shared<Tensor2d>("T1 (Q|ij)", nQ, naoccB, naoccB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(true, false, bQijA, T, 1.0, 0.0);
    T.reset();
    W->sort(2413, X, 1.0, 1.0);
    X.reset();

    //W_MnIj += \sum_(Ef) Tau(Ij,Ef) * <Mn|Ef> 
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    uccsd_tau_amps_OS(naoccA, naoccB, navirA, navirB, Tau, T2, t1A, t1B);
    T2.reset();
    J = std::make_shared<Tensor2d>("J (ME|nf)", naoccA, navirA, naoccB, navirB);
    J->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    I = std::make_shared<Tensor2d>("I <Mn|Ef>", naoccA, naoccB, navirA, navirB);
    I->sort(1324, J, 1.0, 0.0);
    J.reset();
    W->gemm(false, true, I, Tau, 1.0, 1.0);
    I.reset();

    // New T2AB += W_MnIj * Tau(Mn,Ab)
    T2new = std::make_shared<Tensor2d>("New T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2new->read(psio_, PSIF_DFOCC_AMPS);  // ccsd_t2_amps() fonksiyonundan okuyor
    T2new->gemm(true, false, W, Tau, 1.0, 1.0);
    W.reset();
    Tau.reset();
    T2new->write(psio_, PSIF_DFOCC_AMPS);
}// ccsd_W_MnIjT2AB

// ccsd_W_MBEJAAAA
void DFOCC::uccsd_W_MBEJAAAA()
{
    SharedTensor2d W, J, I, X, Y, Z, K, L, M, T;

    // W (MB,EJ) = <MB||EJ> + \sum_(Q)[t(Qp,JB)+0.5*T(Q,JB)]*b(Q,ME) + \sum_(Q)t(Q,BE)*[t(Q,JM)+b(Q,JM)] - \sum_(Q)t(Q,JM)*b(Q,BE) - 0.5*\sum_(N,F)t(JN,BF)*<EM|NF>
    // W_MBEJ = W(ME,JB)
    // W(ME,JB) = (ME|JB) - <ME|JB>
    W = std::make_shared<Tensor2d>("W (ME|JB)", naoccA, navirA, naoccA, navirA);
    W->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    L = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IJ|AB)", naoccA, naoccA, navirA, navirA);
    L->gemm(true, false, bQijA, bQabA, 1.0, 0.0);
    W->sort(1324, L, -1.0, 1.0);
    L.reset();

    // W (ME,JB) += \sum_(Q) [t(Qp,JB) + 0.5 * T(Q,JB)] * b(Q,ME)
    T = std::make_shared<Tensor2d>("T1p (Q|IA)", nQ, naoccA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("T1A + T2/2 (Q|IA)", nQ, naoccA, navirA);
    X->copy(T);
    T.reset();
    K = std::make_shared<Tensor2d>("T2 (Q|IA)", nQ, naoccA, navirA);
    K->read(psio_, PSIF_DFOCC_AMPS);
    X->axpy(K, 0.5);
    K.reset();
    W->gemm(true, false, bQiaA, X, 1.0, 1.0);
    X.reset();

    // W (ME,JB) += \sum_(Q) t(Q,BE) * [t(Q,JM) + b(Q,JM)]
    T = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("T1(Q|IJ) + b(Q|IJ)", nQ, naoccA, naoccA);
    X->copy(T);
    T.reset();
    X->axpy(bQijA, 1.0);
    T = std::make_shared<Tensor2d>("T1 (Q|AB)", nQ, navirA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    I = std::make_shared<Tensor2d>("I (JM|BE)", naoccA, naoccA, navirA, navirA);
    I->gemm(true, false, X, T, 1.0, 0.0);
    X.reset();
    T.reset();
    W->sort(2413, I, 1.0, 1.0);
    I.reset();

    // W (ME,JB) -= \sum_(Q) t(Q,JM) * b(Q,BE) 
    T = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    M = std::make_shared<Tensor2d>("M (JM|BE)", naoccA, naoccA, navirA, navirA);
    M->gemm(true, false, T, bQabA, 1.0, 0.0);
    T.reset();
    W->sort(2413, M, -1.0, 1.0);
    M.reset();

    // W (ME,JB) -= 0.5 * \sum_(NF) t(JN,BF) Y(ME,NF)
    // <me|fn> = (mf|ne) = Y(me,nf) (sort: 1432)
    // t <JN|BF> = t (NF|JB) (sort: 2413)
    K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA);
    K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    X = std::make_shared<Tensor2d>("X_ (IA|JB)", naoccA, navirA, naoccA, navirA);
    X->sort(1432, K, 1.0, 0.0);
    K.reset();
    T = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    L = std::make_shared<Tensor2d>("T2tilde (IA|JB)", naoccA, navirA, naoccA, navirA);
    L->sort(2413, T, 1.0, 0.0);
    T.reset();
    W->gemm(false, false, X, L, -0.5, 1.0);
    X.reset();
    L.reset();
    W->write(psio_, PSIF_DFOCC_AMPS);
}// End ccsd_W_MBEJAAAA

// ccsd_W_mbejBBBB
void DFOCC::uccsd_W_mbejBBBB()
{
    SharedTensor2d W, J, I, X, Y, Z, K, L, M, T;

    // W (mb,ej) = <mb||ej> + \sum_(Q)[t(Qp,jb)+0.5*T(Q,jb)]*b(Q,me) + \sum_(Q)t(Q,be)*[t(Q,jm)+b(Q,jm)] - \sum_(Q)t(Q,jm)*b(Q,be) - 0.5*\sum_(n,f)t(jn,bf)*<em|nf>
    // W_mbej = W(me,jb)
    // W(me,jb) = (me|jb) - <me|jb>
    W = std::make_shared<Tensor2d>("W (me|jb)", naoccB, navirB, naoccB, navirB);
    W->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    L = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (ij|ab)", naoccB, naoccB, navirB, navirB);
    L->gemm(true, false, bQijB, bQabB, 1.0, 0.0);
    W->sort(1324, L, -1.0, 1.0);
    L.reset();

    // W (me,jb) += \sum_(Q) [t(Qp,jb) + 0.5 * T(Q,jb)] * b(Q,me)
    T = std::make_shared<Tensor2d>("T1p (Q|ia)", nQ, naoccB, navirB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("T1B + T2/2 (Q|ia)", nQ, naoccB, navirB);
    X->copy(T);
    T.reset();
    K = std::make_shared<Tensor2d>("T2 (Q|ia)", nQ, naoccB, navirB);
    K->read(psio_, PSIF_DFOCC_AMPS);
    X->axpy(K, 0.5);
    K.reset();
    W->gemm(true, false, bQiaB, X, 1.0, 1.0);
    X.reset();

    // W (me,jb) += \sum_(Q) t(Q,be) * [t(Q,jm) + b(Q,jm)]
    T = std::make_shared<Tensor2d>("T1 (Q|ij)", nQ, naoccB, naoccB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("T1(Q|ij) + b(Q|ij)", nQ, naoccB, naoccB);
    X->copy(T);
    T.reset();
    X->axpy(bQijB, 1.0);
    T = std::make_shared<Tensor2d>("T1 (Q|ab)", nQ, navirB, navirB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    I = std::make_shared<Tensor2d>("I (jm|be)", naoccB, naoccB, navirB, navirB);
    I->gemm(true, false, X, T, 1.0, 0.0);
    X.reset();
    T.reset();
    W->sort(2413, I, 1.0, 1.0);
    I.reset();

    // W (me,jb) -= \sum_(Q) t(Q,jm) * b(Q,be) 
    T = std::make_shared<Tensor2d>("T1 (Q|ij)", nQ, naoccB, naoccB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    M = std::make_shared<Tensor2d>("M (jm|be)", naoccB, naoccB, navirB, navirB);
    M->gemm(true, false, T, bQabB, 1.0, 0.0);
    T.reset();
    W->sort(2413, M, -1.0, 1.0);
    M.reset();

    // W (me,jb) -= 0.5 * \sum_(nf) t(jn,bf) Y(me,nf)
    // <me|fn> = (mf|ne) = Y(me,nf) (sort: 1432)
    // t <jn|bf> = t (nf|jb) (sort: 2413)
    K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB);
    K->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    X = std::make_shared<Tensor2d>("X (ia|jb)", naoccB, navirB, naoccB, navirB);
    X->sort(1432, K, 1.0, 0.0);
    K.reset();
    T = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    L = std::make_shared<Tensor2d>("T2tilde (ia|jb)", naoccB, navirB, naoccB, navirB);
    L->sort(2413, T, 1.0, 0.0);
    T.reset();
    W->gemm(false, false, X, L, -0.5, 1.0);
    X.reset();
    L.reset();
    W->write(psio_, PSIF_DFOCC_AMPS);
}// ccsd_W_mbejBBBB

// ccsd_W_MbEjABAB
void DFOCC::uccsd_W_MbEjABAB()
{
    SharedTensor2d W, X, K, L, T;
    // W (Mb,Ej) = <Mb|Ej> + \sum_(Q) [t(Qp,jb) + 0.5 * T(Q,jb)] * b(Q,ME) - 0.5 * \sum_(N,F) t(Nj,Fb) * <EM|NF>
    // W_MbEj = W(ME,jb)
    // W(ME,jb) <= (ME|jb)
    W = std::make_shared<Tensor2d>("W (ME|jb)", naoccA, navirA, naoccB, navirB);
    W->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);

    // W (ME,jb) += \sum_(Q) [t(Qp,jb) + 0.5 * T(Q,jb)] * b(Q,ME)
    T = std::make_shared<Tensor2d>("T1p (Q|ia)", nQ, naoccB, navirB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("T1A + T2/2 (Q|ia)", nQ, naoccB, navirB);
    X->copy(T);
    T.reset();
    K = std::make_shared<Tensor2d>("T2 (Q|ia)", nQ, naoccB, navirB);
    K->read(psio_, PSIF_DFOCC_AMPS);
    X->axpy(K, 0.5);
    K.reset();
    W->gemm(true, false, bQiaA, X, 1.0, 1.0);
    X.reset();

    // W (ME,jb) -= 0.5 * \sum_(N,F) t(Nj,Fb) * <EM|NF>
    // <ME|FN> = (MF|NE) = Y(ME,NF) (sort: 1432)
    // t <Nj|Fb> = t (NF|jb) (sort: 1324)
    K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA);
    K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    X = std::make_shared<Tensor2d>("X (IA|JB)", naoccA, navirA, naoccA, navirA);
    X->sort(1432, K, 1.0, 0.0);
    K.reset();
    T = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    L = std::make_shared<Tensor2d>("T2tilde (IA|jb)", naoccA, navirA, naoccB, navirB);
    L->sort(1324, T, 1.0, 0.0);
    T.reset();
    W->gemm(false, false, X, L, -0.5, 1.0);
    X.reset();
    L.reset();
    W->write(psio_, PSIF_DFOCC_AMPS);
}// ccsd_W_MbEjABAB

// ccsd_W_mBeJBABA
void DFOCC::uccsd_W_mBeJBABA()
{
    SharedTensor2d W, X, K, L, T;

    // W (mB,eJ) = <Bm|Je> + \sum_(Q) [t(Qp,JB) + 0.5 * T(Q,JB)] * b(Q,me) - 0.5 * \sum_(n,f) t(Jn,Bf) * <em|nf>
    // W_mBeJ = W(me,JB)
    // W(me,JB) <= (me|JB) 
    // <Bm|Je> = <mB|eJ> = (me|BJ) = (me|JB) 
    W = std::make_shared<Tensor2d>("W (me|JB)", naoccB, navirB, naoccA, navirA);
    W->gemm(true, false, bQiaB, bQiaA, 1.0, 0.0);

    // W (me,JB) += \sum_(Q) [t(Qp,JB) + 0.5 * T(Q,JB)] * b(Q,me) 
    T = std::make_shared<Tensor2d>("T1p (Q|IA)", nQ, naoccA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("T1A + T2/2 (Q|IA)", nQ, naoccA, navirA);
    X->copy(T);
    T.reset();
    K = std::make_shared<Tensor2d>("T2 (Q|IA)", nQ, naoccA, navirA);
    K->read(psio_, PSIF_DFOCC_AMPS);
    X->axpy(K, 0.5);
    K.reset();
    W->gemm(true, false, bQiaB, X, 1.0, 1.0);
    X.reset();

    //  W (me,JB) -= 0.5 * \sum_(n,f) t(Jn,Bf) * <em|nf>
    // <em|nf> = <me|fn> = (mf|ne) = Y(me,nf) (sort: 1432)
    // t <Jn|Bf> = t (nf|JB) (sort: 2413)
    K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB);
    K->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    X = std::make_shared<Tensor2d>("X (ia|jb)", naoccB, navirB, naoccB, navirB);
    X->sort(1432, K, 1.0, 0.0);
    K.reset();
    T = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    L = std::make_shared<Tensor2d>("T2tilde (IA|jb)", naoccA, navirA, naoccB, navirB);
    L->sort(1324, T, 1.0, 0.0);
    T.reset();
    W->gemm(false, true, X, L, -0.5, 1.0);
    X.reset();
    L.reset();
    W->write(psio_, PSIF_DFOCC_AMPS);
}// ccsd_W_mBeJBABA

// ccsd_W_MbeJABBA
void DFOCC::uccsd_W_MbeJABBA()
{
    SharedTensor2d W, X, K, L, T, I, M;

    // W(Me,Jb) = - <Mb|Je> + \sum_(Q) t(Q,be) * [t(Q,JM) + b(Q,JM)] - \sum_(Q) t(Q,JM) * b(Q,be) + 0.5 * \sum_(n,F) t(Jn,Fb) * <Me|Fn> 
    // W_MbeJ = W (Me,Jb)
    // W(Me,Jb) = - <Me|Jb> = -(MJ|be) = Y (Me,Jb) (sort: 1423)
    W = std::make_shared<Tensor2d>("W (Me|Jb)", naoccA, navirB, naoccA, navirB);
    K = std::make_shared<Tensor2d>("Int (MJ|be)", naoccA, naoccA, navirB, navirB);
    K->gemm(true, false, bQijA, bQabB, 1.0, 0.0);
    W->sort(1423, K, -1.0, 0.0);
    K.reset();

    // W (Me,Jb) += \sum_(Q) t(Q,be) * [t(Q,JM) + b(Q,JM)]
    T = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("T1(Q|IJ) + b(Q|IJ)", nQ, naoccA, naoccA);
    X->copy(T);
    T.reset();
    X->axpy(bQijA, 1.0);
    T = std::make_shared<Tensor2d>("T1 (Q|ab)", nQ, navirB, navirB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    I = std::make_shared<Tensor2d>("I (JM|be)", naoccA, naoccA, navirB, navirB);
    I->gemm(true, false, X, T, 1.0, 0.0);
    X.reset();
    T.reset();
    W->sort(2413, I, 1.0, 1.0);
    I.reset();

    // W (Me,Jb) -= \sum_(Q) t(Q,JM) * b(Q,be) 
    T = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    M = std::make_shared<Tensor2d>("M (JM|be)", naoccA, naoccA, navirB, navirB);
    M->gemm(true, false, T, bQabB, 1.0, 0.0);
    T.reset();
    W->sort(2413, M, -1.0, 1.0);
    M.reset();

    // W (Me,Jb) += 0.5 * \sum_(n,F) t(Jn,Fb) * <Me|Fn>   
    // <Me|Fn> = (MF|en) = (MF|ne) = Y(Me,Fn) (sort: 1432)
    // t <Jn|Fb> = t (nF|Jb) (sort: 2314)
    K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (MF|ne)", naoccA, navirA, naoccB, navirB);
    K->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    X = std::make_shared<Tensor2d>("X (Me|Fn)", naoccA, navirB, navirA, naoccB);
    X->sort(1423, K, 1.0, 0.0);
    K.reset();
    T = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    L = std::make_shared<Tensor2d>("T2tilde (Aj|Ib)", navirA, naoccB, naoccA, navirB);
    L->sort(3214, T, 1.0, 0.0);
    T.reset();
    W->gemm(false, false, X, L, 0.5, 1.0);
    X.reset();
    L.reset();
    W->write(psio_, PSIF_DFOCC_AMPS);
}// ccsd_W_MbeJABBA

// ccsd_W_mBEjBAAB
void DFOCC::uccsd_W_mBEjBAAB()
{
    SharedTensor2d W, X, K, L, T, I, M;

    // W(mE,jB) = -<Bm|Ej> + \sum_(Q) t(Q,BE) * [t(Q,jm) + b(Q,jm)] - \sum_(Q) t(Q,jm) * b(Q,BE) - 0.5 * \sum_(N,f) t(Nj,Bf) * <Em|Nf>   
    // W_mBEj = W(mE,jB)
    // W(mE,jB) = - <mE|jB> = -(EB|mj)
    W = std::make_shared<Tensor2d>("W (mE|jB)", naoccB, navirA, naoccB, navirA);
    K = std::make_shared<Tensor2d>("Int (EB|mj)", navirA, navirA, naoccB, naoccB);
    K->gemm(true, false, bQabA, bQijB, 1.0, 0.0);
    W->sort(3142, K, -1.0, 0.0);
    K.reset();

    // W (mE,jB) += \sum_(Q) t(Q,BE) * [t(Q,jm) + b(Q,jm)]
    T = std::make_shared<Tensor2d>("T1 (Q|ij)", nQ, naoccB, naoccB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("T1(Q|ij) + b(Q|ij)", nQ, naoccB, naoccB);
    X->copy(T);
    T.reset();
    X->axpy(bQijB, 1.0);
    T = std::make_shared<Tensor2d>("T1 (Q|AB)", nQ, navirA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    I = std::make_shared<Tensor2d>("I (jm|BE)", naoccB, naoccB, navirA, navirA);
    I->gemm(true, false, X, T, 1.0, 0.0);
    X.reset();
    T.reset();
    W->sort(2413, I, 1.0, 1.0);
    I.reset();

    // W (mE,jB) -= \sum_(Q) t(Q,jm) * b(Q,BE) 
    T = std::make_shared<Tensor2d>("T1 (Q|ij)", nQ, naoccB, naoccB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    M = std::make_shared<Tensor2d>("M (jm|BE)", naoccB, naoccB, navirA, navirA);
    M->gemm(true, false, T, bQabA, 1.0, 0.0);
    T.reset();
    W->sort(2413, M, -1.0, 1.0);
    M.reset();

    // W (Me,Jb) += 0.5 * \sum_(n,F) t(Jn,Fb) * <Me|Fn>   
    // <Em|Nf> = <mE|fN> = (mf|NE) = Y(mE,Nf) (sort: 1432)
    // t <Nj|Bf> = t (Nf|jB) (sort: 1423)
    K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (ia|JB)", naoccB, navirB, naoccA, navirA);
    K->gemm(true, false, bQiaB, bQiaA, 1.0, 0.0);
    X = std::make_shared<Tensor2d>("X (iA|Jb)", naoccB, navirA, naoccA, navirB);
    X->sort(1432, K, 1.0, 0.0);
    K.reset();
    T = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    L = std::make_shared<Tensor2d>("T2tilde (Ia|jB)", naoccA, navirB, naoccB, navirA);
    L->sort(1423, T, 1.0, 0.0);
    T.reset();
    W->gemm(false, false, X, L, 0.5, 1.0);
    X.reset();
    L.reset();
    W->write(psio_, PSIF_DFOCC_AMPS);
}// ccsd_W_mBEjBAAB

}  // namespace dfoccwave
}  // namespace psi
