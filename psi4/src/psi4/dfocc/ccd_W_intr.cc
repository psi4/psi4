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

void DFOCC::ccd_WmnijT2() {
    // defs
    SharedTensor2d K, T, Tnew, U, Tau, W, X;
    SharedTensor2d M, L, I, Y, S, A;
    SharedTensor2d V, Vs, Va, Ts, Ta;

    timer_on("WmnijT2");

    // W_mnij = <mn|ij>
    W = std::make_shared<Tensor2d>("W <MN|IJ>", naoccA, naoccA, naoccA, naoccA);
    K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IJ|KL)", naoccA, naoccA, naoccA, naoccA);
    K->gemm(true, false, bQijA, bQijA, 1.0, 0.0);
    W->sort(1324, K, 1.0, 0.0);
    K.reset();

    // W_mnij = \sum_{ef} T_ij^ef <mn|ef>
    // (+)T(ij, ab) = 1/2 (T_ij^ab + T_ji^ab) * (2 - \delta_{ab})
    // (-)T(ij, ab) = 1/2 (T_ij^ab - T_ji^ab) * (2 - \delta_{ab})
    U = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    U->sort(1324, t2, 1.0, 0.0);
    Ts = std::make_shared<Tensor2d>("(+)tT [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
    Ta = std::make_shared<Tensor2d>("(-)tT [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
    Ts->symm_col_packed4(U);
    Ta->antisymm_col_packed4(U);
    U.reset();
    // Form <mn|ef>
    K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA);
    L = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA);
    L->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    Vs = std::make_shared<Tensor2d>("(+)V [M>=N|E>=F]", ntri_ijAA, ntri_abAA);
    Va = std::make_shared<Tensor2d>("(-)V [M>=N|E>=F]", ntri_ijAA, ntri_abAA);
    Vs->symm4(K);
    Va->antisymm4(K);
    K.reset();
    // Form S/A
    S = std::make_shared<Tensor2d>("S [M>=N|I>=J]", ntri_ijAA, ntri_ijAA);
    A = std::make_shared<Tensor2d>("A [M>=N|I>=J]", ntri_ijAA, ntri_ijAA);
    S->gemm(false, true, Vs, Ts, 1.0, 0.0);
    A->gemm(false, true, Va, Ta, 1.0, 0.0);
    Vs.reset();
    Va.reset();
// add to W(mn,ij)
#pragma omp parallel for
    for (int m = 0; m < naoccA; ++m) {
        for (int n = 0; n < naoccA; ++n) {
            int mn = index2(m, n);
            int mn2 = ij_idxAA->get(m, n);
            for (int i = 0; i < naoccA; ++i) {
                for (int j = 0; j < naoccA; ++j) {
                    int ij = index2(i, j);
                    int ij2 = ij_idxAA->get(i, j);
                    int perm1 = (i > j) ? 1 : -1;
                    int perm2 = (m > n) ? 1 : -1;
                    double value = S->get(mn, ij) + (perm1 * perm2 * A->get(mn, ij));
                    W->add(mn2, ij2, value);
                }
            }
        }
    }
    S.reset();
    A.reset();

    // t_ij^ab <= \sum_{m,n} T_mn^ab Wmnij
    // (+)T(ij, ab) = 1/2 (T_ij^ab + T_ji^ab) * (2 - \delta_{ij})
    // (-)T(ij, ab) = 1/2 (T_ij^ab - T_ji^ab) * (2 - \delta_{ij})
    U = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    U->sort(1324, t2, 1.0, 0.0);
    Ts->symm_row_packed4(U);
    Ta->antisymm_row_packed4(U);
    U.reset();

    // Form (+/-)W(m>=n, i>=j)
    Vs = std::make_shared<Tensor2d>("(+)W [M>=N|I>=J]", ntri_ijAA, ntri_ijAA);
    Va = std::make_shared<Tensor2d>("(-)W [M>=N|I>=J]", ntri_ijAA, ntri_ijAA);
    Vs->symm4(W);
    Va->antisymm4(W);
    W.reset();

    // Symmetric & Anti-symmetric contributions
    S = std::make_shared<Tensor2d>("S (I>=J, A>=B)", ntri_ijAA, ntri_abAA);
    A = std::make_shared<Tensor2d>("A (I>=J, A>=B)", ntri_ijAA, ntri_abAA);
    S->gemm(true, false, Vs, Ts, 1.0, 0.0);
    A->gemm(true, false, Va, Ta, 1.0, 0.0);
    Ts.reset();
    Ta.reset();
    Vs.reset();
    Va.reset();

    // T(ia,jb) <-- S(a>=b,i>=j) + A(a>=b,i>=j)
    Tnew = std::make_shared<Tensor2d>("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
#pragma omp parallel for
    for (int a = 0; a < navirA; ++a) {
        for (int b = 0; b < navirA; ++b) {
            int ab = index2(a, b);
            for (int i = 0; i < naoccA; ++i) {
                int ia = ia_idxAA->get(i, a);
                for (int j = 0; j < naoccA; ++j) {
                    int jb = ia_idxAA->get(j, b);
                    int ij = index2(i, j);
                    int perm1 = (i > j) ? 1 : -1;
                    int perm2 = (a > b) ? 1 : -1;
                    double value = S->get(ij, ab) + (perm1 * perm2 * A->get(ij, ab));
                    Tnew->add(ia, jb, value);
                }
            }
        }
    }
    S.reset();
    A.reset();
    Tnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    timer_off("WmnijT2");

}  // end ccd_WmnijT2

//======================================================================
//    WmnijT2AA
//======================================================================
void DFOCC::ccd_W_MNIJT2AA()
{
    SharedTensor2d J, W, I, X, Y, T;
    SharedTensor2d T2, Tau, T2new;

    // W_MNIJ Alpha Block
    // W_MNIJ = <MN||IJ> + \sum_(EF) T(IJ,EF) * <MN|EF>
    // W_MNIJ = <MN||IJ>
    J = std::make_shared<Tensor2d>("J (IM|JN)", naoccA, naoccA, naoccA, naoccA);
    J->gemm(true, false, bQijA, bQijA, 1.0, 0.0);
    W = std::make_shared<Tensor2d>("W <MN|IJ>", naoccA, naoccA, naoccA, naoccA);
    W->sort(1324, J, 1.0, 0.0);
    W->sort(1342, J, -1.0, 1.0);
    J.reset();

    //W_MNIJ += \sum_(EF) T(IJ,EF) * <MN|EF>
    Tau = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    Tau->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    J = std::make_shared<Tensor2d>("J (ME|NF)", naoccA, navirA, naoccA, navirA);
    J->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    I = std::make_shared<Tensor2d>("I <MN|EF>", naoccA, naoccA, navirA, navirA);
    I->sort(1324, J, 1.0, 0.0);
    J.reset();
    W->gemm(false, true, I, Tau, 1.0, 1.0);
    I.reset();

    // 0.5 * Tau * WMNIJ
    T2new = std::make_shared<Tensor2d>("New T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2new->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T2new->gemm(true, false, W, Tau, 0.5, 1.0);
    W.reset();
    Tau.reset();
    T2new->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
}// ccd_W_MNIJT2AA

//======================================================================
//    WmnijT2BB
//======================================================================
void DFOCC::ccd_W_mnijT2BB()
{
    SharedTensor2d J, W, I, X, Y, T;
    SharedTensor2d T2, Tau, T2new;

    // W_mnij Beta Block
    // W_mnij = <mn||ij> + \sum_(ef) Tau(ij,ef) * <mn|ef>
    // W_mnij = <mn||ij>
    J = std::make_shared<Tensor2d>("J (im|jn)", naoccB, naoccB, naoccB, naoccB);
    J->gemm(true, false, bQijB, bQijB, 1.0, 0.0);
    W = std::make_shared<Tensor2d>("W <mn|ij>", naoccB, naoccB, naoccB, naoccB);
    W->sort(1324, J, 1.0, 0.0);
    W->sort(1342, J, -1.0, 1.0);
    J.reset();

    //W_mnij += \sum_(ef) T(ij,ef) * <mn|ef>
    Tau = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    Tau->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    J = std::make_shared<Tensor2d>("J (me|nf)", naoccB, navirB, naoccB, navirB);
    J->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    I = std::make_shared<Tensor2d>("I <mn|ef>", naoccB, naoccB, navirB, navirB);
    I->sort(1324, J, 1.0, 0.0);
    J.reset();
    W->gemm(false, true, I, Tau, 1.0, 1.0);
    I.reset();

    // 0.5 * Tau * Wmnij
    T2new = std::make_shared<Tensor2d>("New T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2new->read_anti_symm(psio_, PSIF_DFOCC_AMPS);  // ccd_t2_amps() fonksiyonundan okuyor
    T2new->gemm(true, false, W, Tau, 0.5, 1.0);
    W.reset();
    Tau.reset();
    T2new->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T2new.reset();

}// ccd_W_mnijT2BB

//======================================================================
//    WmnijT2AB
//======================================================================
void DFOCC::ccd_W_MnIjT2AB()
{
    SharedTensor2d J, W, I, X, Y, T;
    SharedTensor2d T2, Tau, T2new;

    // W_MnIj Alpha-Beta Block
    // W_MnIj = <Mn|Ij> + \sum_(Ef) Tau(Ij,Ef) * <Mn|Ef>
    // W_MnIj = <Mn|Ij>
    W = std::make_shared<Tensor2d>("W <Mn|Ij>", naoccA, naoccB, naoccA, naoccB);
    J = std::make_shared<Tensor2d>("J (MI|nj)", naoccA, naoccA, naoccB, naoccB);
    J->gemm(true, false, bQijA, bQijB, 1.0, 0.0);
    W->sort(1324, J, 1.0, 0.0);
    J.reset();

    //W_MnIj += \sum_(Ef) Tau(Ij,Ef) * <Mn|Ef>
    Tau = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    Tau->read(psio_, PSIF_DFOCC_AMPS);
    J = std::make_shared<Tensor2d>("J (ME|nf)", naoccA, navirA, naoccB, navirB);
    J->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    I = std::make_shared<Tensor2d>("I <Mn|Ef>", naoccA, naoccB, navirA, navirB);
    I->sort(1324, J, 1.0, 0.0);
    J.reset();
    W->gemm(false, true, I, Tau, 1.0, 1.0);
    I.reset();

    // New T2AB += W_MnIj * Tau(Mn,Ab)
    T2new = std::make_shared<Tensor2d>("New T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2new->read(psio_, PSIF_DFOCC_AMPS);  // ccd_t2_amps() fonksiyonundan okuyor
    T2new->gemm(true, false, W, Tau, 1.0, 1.0);
    W.reset();
    Tau.reset();
    T2new->write(psio_, PSIF_DFOCC_AMPS);
    T2new.reset();

}// ccd_W_MnIjT2AB

//======================================================================
//    WmbejT2
//======================================================================
void DFOCC::ccd_WmbejT2() {
    // defs
    SharedTensor2d K, L, T, T1, Tnew, U, Tau, W, W2, X, Y;

    timer_on("WmbejT2");

    // W_mbej = W(me,jb)
    // W(me,jb) <= (me|jb)
    W = std::make_shared<Tensor2d>("W (ME|JB)", naoccA, navirA, naoccA, navirA);
    W->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);

    // W(me,jb) <= 1/2\sum_{Q} T_jb^Q b_me^Q
    T = std::make_shared<Tensor2d>("T2 (Q|IA)", nQ, naoccA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    W->gemm(true, false, bQiaA, T, 0.5, 1.0);
    T.reset();

    // W(me,jb) <= -1/2 \sum_{nf} t_jn^bf X(me,nf)
    // (mf|ne) = X(me,nf) (sort: 1432)
    K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA);
    K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    X = std::make_shared<Tensor2d>("X (IA|JB)", naoccA, navirA, naoccA, navirA);
    X->sort(1432, K, 1.0, 0.0);
    K.reset();
    W->gemm(false, false, X, t2, -0.5, 1.0);
    X.reset();
    W->write(psio_, PSIF_DFOCC_AMPS);
    W.reset();

    // W_mbje = W'(me,jb)
    // W'(me,jb) <= <me|jb>
    W = std::make_shared<Tensor2d>("Wp (ME|JB)", naoccA, navirA, naoccA, navirA);
    L = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IJ|AB)", naoccA, naoccA, navirA, navirA);
    L->gemm(true, false, bQijA, bQabA, 1.0, 0.0);
    W->sort(1324, L, 1.0, 0.0);
    L.reset();

    // W'(me,jb) <= -1/2 \sum_{nf} t_nj^bf X(me,nf)
    // (mf|ne) = X(me,nf) (sort: 1432)
    K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA);
    K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    X = std::make_shared<Tensor2d>("X (IA|JB)", naoccA, navirA, naoccA, navirA);
    X->sort(1432, K, 1.0, 0.0);
    K.reset();
    T = std::make_shared<Tensor2d>("T2p (IA|JB)", naoccA, navirA, naoccA, navirA);
    ccsd_t2_prime_amps(T, t2);
    W->gemm(false, false, X, T, -0.5, 1.0);
    X.reset();
    T.reset();

    // t_ij^ab <= 1/2*C(ia,jb) + 1/2*C(jb,ia) + C(ja,ib) + C(ib,ja)
    // t_ij^ab <= Ct(ia,jb) + 2*Ct(ib,ja)
    // C(ia,jb) = -\sum_{me} t_mi^ae W'(me,jb) = -\sum_{me} T'(ia,me) W'(me,jb)
    U = std::make_shared<Tensor2d>("T2p (IA|JB)", naoccA, navirA, naoccA, navirA);
    ccsd_t2_prime_amps(U, t2);
    Y = std::make_shared<Tensor2d>("C2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    Y->gemm(false, false, U, W, -1.0, 0.0);
    U.reset();
    X = std::make_shared<Tensor2d>("C2+D2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    X->sort(1432, Y, 1.0, 0.0);
    X->axpy(Y, 0.5);
    Y.reset();

    // t_ij^ab <= D(ia,jb) + D(jb,ia)
    // D_ij^ab = 1/2 \sum_{me} u_im^ae [2*W(me,jb) - W'(me,jb)]
    Y = std::make_shared<Tensor2d>("2*W-W' (ME|JB)", naoccA, navirA, naoccA, navirA);
    Y->axpy(W, -1.0);
    W.reset();
    W2 = std::make_shared<Tensor2d>("W (ME|JB)", naoccA, navirA, naoccA, navirA);
    W2->read(psio_, PSIF_DFOCC_AMPS);
    Y->axpy(W2, 2.0);
    W2.reset();
    U = std::make_shared<Tensor2d>("U2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    ccsd_u2_amps(U, t2);
    X->gemm(false, false, U, Y, 0.5, 1.0);
    U.reset();
    Y.reset();
    X->symmetrize();
    Tnew = std::make_shared<Tensor2d>("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->axpy(X, 2.0);
    X.reset();
    Tnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    timer_off("WmbejT2");

}  // end ccd_WmbejT2

//======================================================================
// ccd_W_MBEJAAAA
//======================================================================
void DFOCC::ccd_W_MBEJAAAA()
{
    SharedTensor2d W, J, I, X, Y, Z, K, L, M, T;

    // W (MB,EJ) = <MB||EJ> + 0.5*\sum_(Q)T(Q,JB)*b(Q,ME) - 0.5*\sum_(N,F)t(JN,BF)*<EM|NF>
    // W_MBEJ = W(ME,JB)
    // W(ME,JB) = (ME|JB) - <ME|JB>
    W = std::make_shared<Tensor2d>("W (ME|JB)", naoccA, navirA, naoccA, navirA);
    W->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    L = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IJ|AB)", naoccA, naoccA, navirA, navirA);
    L->gemm(true, false, bQijA, bQabA, 1.0, 0.0);
    W->sort(1324, L, -1.0, 1.0);
    L.reset();

    // W (ME,JB) += 0.5 * \sum_(Q) [T(Q,JB)] * b(Q,ME)
    X = std::make_shared<Tensor2d>("T2 (Q|IA)", nQ, naoccA, navirA);
    X->read(psio_, PSIF_DFOCC_AMPS);
    W->gemm(true, false, bQiaA, X, 0.5, 1.0);
    X.reset();

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
    W.reset();

}// End ccd_W_MBEJAAAA

//======================================================================
// ccd_W_mbejBBBB
//======================================================================
void DFOCC::ccd_W_mbejBBBB()
{
    SharedTensor2d W, J, I, X, Y, Z, K, L, M, T;

    // W (mb,ej) = <mb||ej> + 0.5*\sum_(Q)[T(Q,jb)]*b(Q,me) - 0.5*\sum_(n,f)t(jn,bf)*<em|nf>
    // W_mbej = W(me,jb)
    // W(me,jb) = (me|jb) - <me|jb>
    W = std::make_shared<Tensor2d>("W (me|jb)", naoccB, navirB, naoccB, navirB);
    W->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    L = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (ij|ab)", naoccB, naoccB, navirB, navirB);
    L->gemm(true, false, bQijB, bQabB, 1.0, 0.0);
    W->sort(1324, L, -1.0, 1.0);
    L.reset();

    // W (me,jb) += 0.5*\sum_(Q) [T(Q,jb)] * b(Q,me)
    X = std::make_shared<Tensor2d>("T2 (Q|ia)", nQ, naoccB, navirB);
    X->read(psio_, PSIF_DFOCC_AMPS);
    W->gemm(true, false, bQiaB, X, 0.5, 1.0);
    X.reset();

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
    W.reset();
}// ccd_W_mbejBBBB

//======================================================================
// ccd_W_MbEjABAB
//======================================================================
void DFOCC::ccd_W_MbEjABAB()
{
    SharedTensor2d W, X, K, L, T;
    // W (Mb,Ej) = <Mb|Ej> + 0.5*\sum_(Q) [T(Q,jb)] * b(Q,ME) - 0.5 * \sum_(N,F) t(Nj,Fb) * <EM|NF>
    // W_MbEj = W(ME,jb)
    // W(ME,jb) <= (ME|jb)
    W = std::make_shared<Tensor2d>("W (ME|jb)", naoccA, navirA, naoccB, navirB);
    W->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);

    // W (ME,jb) += 0.5*\sum_(Q) [T(Q,jb)] * b(Q,ME)
    X = std::make_shared<Tensor2d>("T2 (Q|ia)", nQ, naoccB, navirB);
    X->read(psio_, PSIF_DFOCC_AMPS);
    W->gemm(true, false, bQiaA, X, 0.5, 1.0);
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
    W.reset();
}// ccd_W_MbEjABAB

//======================================================================
// ccd_W_mBeJBABA
//======================================================================
void DFOCC::ccd_W_mBeJBABA()
{
    SharedTensor2d W, X, K, L, T;

    // W (mB,eJ) = <Bm|Je> + 1/2\sum_(Q) [T(Q,JB)] * b(Q,me) - 0.5 * \sum_(n,f) t(Jn,Bf) * <em|nf>
    // W_mBeJ = W(me,JB)
    // W(me,JB) <= (me|JB)
    // <Bm|Je> = <mB|eJ> = (me|BJ) = (me|JB)
    W = std::make_shared<Tensor2d>("W (me|JB)", naoccB, navirB, naoccA, navirA);
    W->gemm(true, false, bQiaB, bQiaA, 1.0, 0.0);

    // W (me,JB) += 1/2\sum_(Q) [T(Q,JB)] * b(Q,me)
    X = std::make_shared<Tensor2d>("T2 (Q|IA)", nQ, naoccA, navirA);
    X->read(psio_, PSIF_DFOCC_AMPS);
    K.reset();
    W->gemm(true, false, bQiaB, X, 0.5, 1.0);
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
    W.reset();
}// ccd_W_mBeJBABA

//======================================================================
// ccd_W_MbeJABBA
//======================================================================
void DFOCC::ccd_W_MbeJABBA()
{
    SharedTensor2d W, X, K, L, T, I, M;

    // W(Me,Jb) = - <Mb|Je> + 0.5 * \sum_(n,F) t(Jn,Fb) * <Me|Fn>
    // W_MbeJ = W (Me,Jb)
    // W(Me,Jb) = - <Me|Jb> = -(MJ|be) = Y (Me,Jb) (sort: 1423)
    W = std::make_shared<Tensor2d>("W (Me|Jb)", naoccA, navirB, naoccA, navirB);
    K = std::make_shared<Tensor2d>("Int (MJ|be)", naoccA, naoccA, navirB, navirB);
    K->gemm(true, false, bQijA, bQabB, 1.0, 0.0);
    W->sort(1423, K, -1.0, 0.0);
    K.reset();

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
    W.reset();
}// ccd_W_MbeJABBA

//======================================================================
// ccd_W_mBEjBAAB
//======================================================================
void DFOCC::ccd_W_mBEjBAAB()
{
    SharedTensor2d W, X, K, L, T, I, M;

    // W(mE,jB) = -<Bm|Ej> - 0.5 * \sum_(N,f) t(Nj,Bf) * <Em|Nf>
    // W_mBEj = W(mE,jB)
    // W(mE,jB) = - <mE|jB> = -(EB|mj)
    W = std::make_shared<Tensor2d>("W (mE|jB)", naoccB, navirA, naoccB, navirA);
    K = std::make_shared<Tensor2d>("Int (EB|mj)", navirA, navirA, naoccB, naoccB);
    K->gemm(true, false, bQabA, bQijB, 1.0, 0.0);
    W->sort(3142, K, -1.0, 0.0);
    K.reset();

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
    W.reset();
}// ccd_W_mBEjBAAB

//======================================================================
//    WabefT2
//======================================================================
void DFOCC::ccd_WabefT2() {
    // defs
    SharedTensor2d K, M, L, I, T, Tnew, U, Tau, W, X, Y, S, A;
    SharedTensor2d V, Vs, Ts, Va, Ta;

    timer_on("WabefT2");

    // t_ij^ab <= \sum_{ef} T_ij^ef <ab|ef>
    // (+)T(ij, ab) = 1/2 (T_ij^ab + T_ji^ab) * (2 - \delta_{ab})
    // (-)T(ij, ab) = 1/2 (T_ij^ab - T_ji^ab) * (2 - \delta_{ab})
    U = std::make_shared<Tensor2d>("(+)T [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
    T = std::make_shared<Tensor2d>("(-)T [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
#pragma omp parallel for
    for (int i = 0; i < naoccA; ++i) {
        for (int j = 0; j <= i; ++j) {
            int ij = index2(i, j);
            for (int a = 0; a < navirA; ++a) {
                int ia = ia_idxAA->get(i, a);
                int ja = ia_idxAA->get(j, a);
                for (int b = 0; b <= a; ++b) {
                    double perm = (a == b ? 1.0 : 2.0);
                    int ab = index2(a, b);
                    int jb = ia_idxAA->get(j, b);
                    int ib = ia_idxAA->get(i, b);
                    double value1 = 0.5 * perm * (t2->get(ia, jb) + t2->get(ja, ib));
                    double value2 = 0.5 * perm * (t2->get(ia, jb) - t2->get(ja, ib));
                    U->set(ij, ab, value1);
                    T->set(ij, ab, value2);
                }
            }
        }
    }

    // Read B(Q,ab)
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (AB|Q)", navirA * navirA, nQ);
    K = bQabA->transpose();

    // malloc
    I = std::make_shared<Tensor2d>("I[A] <BF|E>", navirA * navirA, navirA);
    Vs = std::make_shared<Tensor2d>("(+)V[A] (B, E>=F)", navirA, ntri_abAA);
    Va = std::make_shared<Tensor2d>("(-)V[A] (B, E>=F)", navirA, ntri_abAA);
    Ts = std::make_shared<Tensor2d>("(+)T[A] (B, I>=J)", navirA, ntri_ijAA);
    Ta = std::make_shared<Tensor2d>("(-)T[B] (B, I>=J)", navirA, ntri_ijAA);

    // Symmetric & Anti-symmetric contributions
    S = std::make_shared<Tensor2d>("S (A>=B, I>=J)", ntri_abAA, ntri_ijAA);
    A = std::make_shared<Tensor2d>("A (A>=B, I>=J)", ntri_abAA, ntri_ijAA);
    // Main loop
    for (int a = 0; a < navirA; ++a) {
        int nb = a + 1;

        // Form V[a](bf,e) = \sum_{Q} B(bfQ)*B(aeQ) cost = V^4N/2
        I->contract(false, true, navirA * nb, navirA, nQ, K, K, 0, a * navirA * nQ, 1.0, 0.0);

// Form (+)V[a](b, e>=f)
#pragma omp parallel for
        for (int b = 0; b <= a; ++b) {
            for (int e = 0; e < navirA; ++e) {
                int be = e + (b * navirA);
                for (int f = 0; f <= e; ++f) {
                    int ef = index2(e, f);
                    int bf = f + (b * navirA);
                    double value1 = 0.5 * (I->get(bf, e) + I->get(be, f));
                    double value2 = 0.5 * (I->get(bf, e) - I->get(be, f));
                    Vs->set(b, ef, value1);
                    Va->set(b, ef, value2);
                }
            }
        }

        // Form T[a](b, i>=j) = \sum_{e>=f} Tau(i>=j,e>=f) V[a](b, e>=f)
        Ts->contract(false, true, nb, ntri_ijAA, ntri_abAA, Vs, U, 1.0, 0.0);
        Ta->contract(false, true, nb, ntri_ijAA, ntri_abAA, Va, T, 1.0, 0.0);

// Form S(ij,ab) & A(ij,ab)
#pragma omp parallel for
        for (int b = 0; b <= a; ++b) {
            int ab = index2(a, b);
            for (int i = 0; i < naoccA; ++i) {
                for (int j = 0; j <= i; ++j) {
                    int ij = index2(i, j);
                    S->add(ab, ij, Ts->get(b, ij));
                    A->add(ab, ij, Ta->get(b, ij));
                }
            }
        }
    }
    K.reset();
    I.reset();
    Vs.reset();
    Va.reset();
    Ts.reset();
    Ta.reset();
    U.reset();
    T.reset();

    // T(ia,jb) <-- S(a>=b,i>=j) + A(a>=b,i>=j)
    Tnew = std::make_shared<Tensor2d>("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
#pragma omp parallel for
    for (int a = 0; a < navirA; ++a) {
        for (int b = 0; b < navirA; ++b) {
            int ab = index2(a, b);
            for (int i = 0; i < naoccA; ++i) {
                int ia = ia_idxAA->get(i, a);
                for (int j = 0; j < naoccA; ++j) {
                    int jb = ia_idxAA->get(j, b);
                    int ij = index2(i, j);
                    int perm1 = (i > j) ? 1 : -1;
                    int perm2 = (a > b) ? 1 : -1;
                    double value = S->get(ab, ij) + (perm1 * perm2 * A->get(ab, ij));
                    Tnew->add(ia, jb, value);
                }
            }
        }
    }
    S.reset();
    A.reset();
    Tnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    timer_off("WabefT2");

}  // end ccd_WabefT2

//======================================================================
//    WabefT2: High memory version
//======================================================================
void DFOCC::ccd_WabefT2_high_mem() {
    // defs
    SharedTensor2d K, M, L, I, T, Tnew, U, Tau, W, X, Y, S, A;
    SharedTensor2d V, Vs, Ts, Va, Ta, J, T1;

    timer_on("WabefT2");

    // t_ij^ab <= \sum_{ef} T_ij^ef <ab|ef>
    // (+)T(ij, ab) = 1/2 (T_ij^ab + T_ji^ab) * (2 - \delta_{ab})
    // (-)T(ij, ab) = 1/2 (T_ij^ab - T_ji^ab) * (2 - \delta_{ab})
    U = std::make_shared<Tensor2d>("(+)T [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
    T = std::make_shared<Tensor2d>("(-)T [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
#pragma omp parallel for
    for (int i = 0; i < naoccA; ++i) {
        for (int j = 0; j <= i; ++j) {
            int ij = index2(i, j);
            for (int a = 0; a < navirA; ++a) {
                int ia = ia_idxAA->get(i, a);
                int ja = ia_idxAA->get(j, a);
                for (int b = 0; b <= a; ++b) {
                    double perm = (a == b ? 1.0 : 2.0);
                    int ab = index2(a, b);
                    int jb = ia_idxAA->get(j, b);
                    int ib = ia_idxAA->get(i, b);
                    double value1 = 0.5 * perm * (t2->get(ia, jb) + t2->get(ja, ib));
                    double value2 = 0.5 * perm * (t2->get(ia, jb) - t2->get(ja, ib));
                    U->set(ij, ab, value1);
                    T->set(ij, ab, value2);
                }
            }
        }
    }

    // Read B(Q,a>=b)
    bQabA.reset();
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, ntri_abAA);
    K->read(psio_, PSIF_DFOCC_INTS);

    // Form (A>=E|B>=F) : cost = V4N/4
    J = std::make_shared<Tensor2d>("J (A>=E|B>=F)", ntri_abAA, ntri_abAA);
    J->gemm(true, false, K, K, 1.0, 0.0);
    K.reset();

    // malloc
    Vs = std::make_shared<Tensor2d>("(+)V[A] (B, E>=F)", navirA, ntri_abAA);
    Va = std::make_shared<Tensor2d>("(-)V[A] (B, E>=F)", navirA, ntri_abAA);
    Ts = std::make_shared<Tensor2d>("(+)T[A] (B, I>=J)", navirA, ntri_ijAA);
    Ta = std::make_shared<Tensor2d>("(-)T[B] (B, I>=J)", navirA, ntri_ijAA);

    // Symmetric & Anti-symmetric contributions
    S = std::make_shared<Tensor2d>("S (A>=B, I>=J)", ntri_abAA, ntri_ijAA);
    A = std::make_shared<Tensor2d>("A (A>=B, I>=J)", ntri_abAA, ntri_ijAA);
    // Main loop
    for (int a = 0; a < navirA; ++a) {
        int nb = a + 1;

// Form (+)V[a](b, e>=f)
#pragma omp parallel for
        for (int b = 0; b <= a; ++b) {
            for (int e = 0; e < navirA; ++e) {
                int ae = index2(a, e);
                int be = index2(b, e);
                for (int f = 0; f <= e; ++f) {
                    int ef = index2(e, f);
                    int bf = index2(b, f);
                    int af = index2(a, f);
                    double value1 = 0.5 * (J->get(ae, bf) + J->get(af, be));
                    double value2 = 0.5 * (J->get(ae, bf) - J->get(af, be));
                    Vs->set(b, ef, value1);
                    Va->set(b, ef, value2);
                }
            }
        }

        // Form T[a](b, i>=j) = \sum_{e>=f} Tau(i>=j,e>=f) V[a](b, e>=f)
        Ts->contract(false, true, nb, ntri_ijAA, ntri_abAA, Vs, U, 1.0, 0.0);
        Ta->contract(false, true, nb, ntri_ijAA, ntri_abAA, Va, T, 1.0, 0.0);

// Form S(ij,ab) & A(ij,ab)
#pragma omp parallel for
        for (int b = 0; b <= a; ++b) {
            int ab = index2(a, b);
            for (int i = 0; i < naoccA; ++i) {
                for (int j = 0; j <= i; ++j) {
                    int ij = index2(i, j);
                    S->add(ab, ij, Ts->get(b, ij));
                    A->add(ab, ij, Ta->get(b, ij));
                }
            }
        }
    }
    J.reset();
    Vs.reset();
    Va.reset();
    Ts.reset();
    Ta.reset();
    U.reset();
    T.reset();

    // T(ia,jb) <-- S(a>=b,i>=j) + A(a>=b,i>=j)
    Tnew = std::make_shared<Tensor2d>("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
#pragma omp parallel for
    for (int a = 0; a < navirA; ++a) {
        for (int b = 0; b < navirA; ++b) {
            int ab = index2(a, b);
            for (int i = 0; i < naoccA; ++i) {
                int ia = ia_idxAA->get(i, a);
                for (int j = 0; j < naoccA; ++j) {
                    int jb = ia_idxAA->get(j, b);
                    int ij = index2(i, j);
                    int perm1 = (i > j) ? 1 : -1;
                    int perm2 = (a > b) ? 1 : -1;
                    double value = S->get(ab, ij) + (perm1 * perm2 * A->get(ab, ij));
                    Tnew->add(ia, jb, value);
                }
            }
        }
    }
    S.reset();
    A.reset();
    Tnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    // Read B(Q,ab)
    bQabA = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA);
    bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);

    timer_off("WabefT2");

}  // end ccd_WabefT2_high_mem

}  // namespace dfoccwave
}  // namespace psi
