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

void DFOCC::ccsd_WmnijT2_low() {
    // defs
    SharedTensor2d K, T, Tnew, U, Tau, W, X;
    SharedTensor2d M, L, I, Y, S, A;
    SharedTensor2d V, Vs, Va, Ts, Ta;

    timer_on("WmnijT2");

    // W_mnij = <mn|ij>
    W = std::make_shared<Tensor2d>("W <MN|IJ>", naoccA, naoccA, naoccA, naoccA);
    tei_ijkl_phys_directAA(W);

    // W_mnij += X(im,jn) + X(jn,im) += 2Xt(im,jn)
    // X_imjn = \sum_{Q} t_im^Q b_jn^Q
    X = std::make_shared<Tensor2d>("X <MN|IJ>", naoccA, naoccA, naoccA, naoccA);
    T = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA);
    K->read(psio_, PSIF_DFOCC_INTS);
    X->gemm(true, false, T, K, 1.0, 0.0);
    T.reset();
    K.reset();
    X->symmetrize();
    W->sort(2413, X, 2.0, 1.0);
    X.reset();

    // W_mnij = \sum_{ef} Tau_ij^ef <mn|ef>
    // (+)Tau(ij, ab) = 1/2 (Tau_ij^ab + Tau_ji^ab) * (2 - \delta_{ab})
    // (-)Tau(ij, ab) = 1/2 (Tau_ij^ab - Tau_ji^ab) * (2 - \delta_{ab})
    Tau = std::make_shared<Tensor2d>("Tau (IA|JB)", naoccA, navirA, naoccA, navirA);
    Tau->read_symm(psio_, PSIF_DFOCC_AMPS);
    U = std::make_shared<Tensor2d>("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA);
    U->sort(1324, Tau, 1.0, 0.0);
    Tau.reset();
    Ts = std::make_shared<Tensor2d>("(+)tTau [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
    Ta = std::make_shared<Tensor2d>("(-)tTau [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
    Ts->symm_col_packed4(U);
    Ta->antisymm_col_packed4(U);
    // Form <mn|ef>
    K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA);
    tei_ijab_phys_directAA(K);
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

    // t_ij^ab <= \sum_{m,n} Tau_mn^ab Wmnij
    // (+)Tau(ij, ab) = 1/2 (Tau_ij^ab + Tau_ji^ab) * (2 - \delta_{ij})
    // (-)Tau(ij, ab) = 1/2 (Tau_ij^ab - Tau_ji^ab) * (2 - \delta_{ij})
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

}  // end ccsd_WmnijT2_low

//======================================================================
//    WmnijT2
//======================================================================
/*
void DFOCC::ccsd_WmnijT2_low()
{
    // defs
    SharedTensor2d K, T, Tnew, U, Tau, W, X;

    timer_on("WmnijT2");

    // W_mnij = <mn|ij>
    W = std::make_shared<Tensor2d>("W <MN|IJ>", naoccA, naoccA, naoccA, naoccA);
    tei_ijkl_phys_directAA(W);

    // W_mnij += X(im,jn) + X(jn,im) += 2Xt(im,jn)
    // X_imjn = \sum_{Q} t_im^Q b_jn^Q
    X = std::make_shared<Tensor2d>("X <MN|IJ>", naoccA, naoccA, naoccA, naoccA);
    T = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA);
    K->read(psio_, PSIF_DFOCC_INTS);
    X->gemm(true, false, T, K, 1.0, 0.0);
    T.reset();
    K.reset();
    X->symmetrize();
    W->sort(2413, X, 2.0, 1.0);
    X.reset();

    // W_mnij = \sum_{ef} Tau_ij^ef <mn|ef>
    Tau = std::make_shared<Tensor2d>("Tau (IA|JB)", naoccA, navirA, naoccA, navirA);
    Tau->read_symm(psio_, PSIF_DFOCC_AMPS);
    U = std::make_shared<Tensor2d>("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA);
    U->sort(1324, Tau, 1.0, 0.0);
    Tau.reset();
    K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA);
    tei_ijab_phys_directAA(K);
    W->gemm(false, true, K, U, 1.0, 1.0);
    K.reset();

    // t_ij^ab <= \sum_{m,n} Tau_mn^ab Wmnij
    T = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T->gemm(true, false, W, U, 1.0, 0.0);
    U.reset();
    W.reset();
    Tnew = std::make_shared<Tensor2d>("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->sort(1324, T, 1.0, 1.0);
    T.reset();
    Tnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    timer_off("WmnijT2");

}// end ccsd_WmnijT2_low
*/

//======================================================================
//    WmbejT2
//======================================================================
void DFOCC::ccsd_WmbejT2_low() {
    // defs
    SharedTensor2d K, T, T1, Tnew, U, Tau, W, W2, X, Y;

    timer_on("WmbejT2");

    // W_mbej = W(me,jb)
    // W(me,jb) <= (me|jb)
    W = std::make_shared<Tensor2d>("W (ME|JB)", naoccA, navirA, naoccA, navirA);
    tei_iajb_chem_directAA(W);

    // W(me,jb) <= \sum_{Q} (t_jb^Q' + 1/2 T_jb^Q) b_me^Q
    T1 = std::make_shared<Tensor2d>("T1p (Q|IA)", nQ, naoccA, navirA);
    T1->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T2 (Q|IA)", nQ, naoccA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    U = std::make_shared<Tensor2d>("T1 + T2/2 (Q|IA)", nQ, naoccA, navirA);
    U->copy(T);
    T.reset();
    U->scale(0.5);
    U->add(T1);
    T1.reset();
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
    K->read(psio_, PSIF_DFOCC_INTS);
    W->gemm(true, false, K, U, 1.0, 1.0);
    K.reset();
    U.reset();

    // W(me,jb) <= -1/2 \sum_{nf} t_jn^bf X(me,nf)
    // (mf|ne) = X(me,nf) (sort: 1432)
    K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA);
    tei_iajb_chem_directAA(K);
    X = std::make_shared<Tensor2d>("X (IA|JB)", naoccA, navirA, naoccA, navirA);
    X->sort(1432, K, 1.0, 0.0);
    K.reset();
    T = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    T->read_symm(psio_, PSIF_DFOCC_AMPS);
    W->gemm(false, false, X, T, -0.5, 1.0);
    T.reset();
    X.reset();
    W->write(psio_, PSIF_DFOCC_AMPS);
    W.reset();

    // W_mbje = W'(me,jb)
    // W'(me,jb) <= <me|jb>
    W = std::make_shared<Tensor2d>("Wp (ME|JB)", naoccA, navirA, naoccA, navirA);
    tei_iajb_phys_directAA(W);

    // X(jm,be) <= -\sum_{Q} t_be^Q ( t_jm^Q + b_jm^Q )
    T = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA);
    K->read(psio_, PSIF_DFOCC_INTS);
    K->add(T);
    T.reset();
    T = std::make_shared<Tensor2d>("T1 (Q|AB)", nQ, navirA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (IJ|AB)", naoccA, naoccA, navirA, navirA);
    X->gemm(true, false, K, T, -1.0, 0.0);
    T.reset();
    K.reset();
    // X(jm,be) <= \sum_{Q} t_jm^Q b_be^Q
    T = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    X->gemm(true, false, T, K, 1.0, 1.0);
    // W'(me,jb) <= X(jm,be)
    W->sort(2413, X, 1.0, 1.0);
    X.reset();

    // W'(me,jb) <= -1/2 \sum_{nf} t_nj^bf X(me,nf)
    // (mf|ne) = X(me,nf) (sort: 1432)
    K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA);
    tei_iajb_chem_directAA(K);
    X = std::make_shared<Tensor2d>("X (IA|JB)", naoccA, navirA, naoccA, navirA);
    X->sort(1432, K, 1.0, 0.0);
    K.reset();
    T = std::make_shared<Tensor2d>("T2p (IA|JB)", naoccA, navirA, naoccA, navirA);
    T->read_symm(psio_, PSIF_DFOCC_AMPS);
    W->gemm(false, false, X, T, -0.5, 1.0);
    X.reset();
    T.reset();

    // t_ij^ab <= 1/2*C(ia,jb) + 1/2*C(jb,ia) + C(ja,ib) + C(ib,ja)
    // t_ij^ab <= Ct(ia,jb) + 2*Ct(ib,ja)
    // C(ia,jb) = -\sum_{me} t_mi^ae W'(me,jb) = -\sum_{me} T'(ia,me) W'(me,jb)
    U = std::make_shared<Tensor2d>("T2p (IA|JB)", naoccA, navirA, naoccA, navirA);
    U->read_symm(psio_, PSIF_DFOCC_AMPS);
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
    U->read_symm(psio_, PSIF_DFOCC_AMPS);
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

}  // end ccsd_WmbejT2_low

//======================================================================
//    WijamT2
//======================================================================
void DFOCC::ccsd_WijamT2_low() {
    // defs
    SharedTensor2d K, M, L, I, T, Tnew, U, Tau, W, X, Y, S, A;
    SharedTensor2d V, Vs, Ts, Va, Ta;

    timer_on("WijamT2");

    // W_ijam = W(ij,am)
    W = std::make_shared<Tensor2d>("W (IJ|AM)", naoccA, naoccA, navirA, naoccA);

    // W_ijam = \sum_{ef} Tau_ij^ef <am|ef>
    Tau = std::make_shared<Tensor2d>("Tau (IA|JB)", naoccA, navirA, naoccA, navirA);
    Tau->read_symm(psio_, PSIF_DFOCC_AMPS);
    // (+)Tau(ij, ab) = 1/2 (Tau_ij^ab + Tau_ji^ab) * (2 - \delta_{ab})
    // (-)Tau(ij, ab) = 1/2 (Tau_ij^ab - Tau_ji^ab) * (2 - \delta_{ab})
    U = std::make_shared<Tensor2d>("(+)Tau [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
    T = std::make_shared<Tensor2d>("(-)Tau [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
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
                    double value1 = 0.5 * perm * (Tau->get(ia, jb) + Tau->get(ja, ib));
                    double value2 = 0.5 * perm * (Tau->get(ia, jb) - Tau->get(ja, ib));
                    U->set(ij, ab, value1);
                    T->set(ij, ab, value2);
                }
            }
        }
    }
    Tau.reset();

    // Read B(Q,ab)
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    X = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
    X->read(psio_, PSIF_DFOCC_INTS);
    L = std::make_shared<Tensor2d>("DF_BASIS_CC B (IA|Q)", naoccA * navirA, nQ);
    L = X->transpose();
    X.reset();
    I = std::make_shared<Tensor2d>("I[M] <AE|F>", navirA * navirA, navirA);

    // Symmetric & Anti-symmetric contributions
    Vs = std::make_shared<Tensor2d>("(+)V[M] (A, E>=F)", navirA, ntri_abAA);
    Va = std::make_shared<Tensor2d>("(-)V[M] (A, E>=F)", navirA, ntri_abAA);
    S = std::make_shared<Tensor2d>("S (AM, I>=J)", navirA, ntri_ijAA);
    A = std::make_shared<Tensor2d>("A (AM, I>=J)", navirA, ntri_ijAA);
    // Main loop
    for (int m = 0; m < naoccA; ++m) {
        // Form V[m](ae,f) = \sum_{Q} b(Q,ae) B(mfQ)
        I->contract(true, true, navirA * navirA, navirA, nQ, K, L, 0, m * navirA * nQ, 1.0, 0.0);

// Form (+)V[m](a, e>=f)
#pragma omp parallel for
        for (int a = 0; a < navirA; ++a) {
            for (int e = 0; e < navirA; ++e) {
                int ae = ab_idxAA->get(a, e);
                for (int f = 0; f <= e; ++f) {
                    int af = ab_idxAA->get(a, f);
                    int ef = index2(e, f);
                    double value1 = 0.5 * (I->get(ae, f) + I->get(af, e));
                    double value2 = 0.5 * (I->get(ae, f) - I->get(af, e));
                    Vs->set(a, ef, value1);
                    Va->set(a, ef, value2);
                }
            }
        }

        // Form S[m](a, i>=j) = \sum_{e>=f} Tau(i>=j,e>=f) V[m](a, e>=f)
        S->gemm(false, true, Vs, U, 1.0, 0.0);
        A->gemm(false, true, Va, T, 1.0, 0.0);

// Form S(am,ij) & A(ab,ij)-->W(ijam)
#pragma omp parallel for
        for (int a = 0; a < navirA; ++a) {
            int am = ai_idxAA->get(a, m);
            for (int i = 0; i < naoccA; ++i) {
                for (int j = 0; j < naoccA; ++j) {
                    int ij = ij_idxAA->get(i, j);
                    int ij2 = index2(i, j);
                    int perm = (i > j) ? 1 : -1;
                    double value = S->get(a, ij2) + (perm * A->get(a, ij2));
                    W->set(ij, am, value);
                }
            }
        }
    }
    K.reset();
    I.reset();
    Vs.reset();
    Va.reset();
    U.reset();
    T.reset();
    S.reset();
    A.reset();
    L.reset();

    // Y(ij,ab) <= -\sum_{m} t_m^b Wijam
    // X(ia,jb) = Y(ij,ab)
    // t_ij^ab <= X(ia,jb) + X(jb,ia)
    Y = std::make_shared<Tensor2d>("Y <IJ|AB>", naoccA, naoccA, navirA, navirA);
    Y->contract(false, false, naoccA * naoccA * navirA, navirA, naoccA, W, t1A, -1.0, 0.0);
    W.reset();
    X = std::make_shared<Tensor2d>("X (IA|JB)", naoccA, navirA, naoccA, navirA);
    X->sort(1324, Y, 1.0, 0.0);
    Y.reset();
    X->symmetrize();
    Tnew = std::make_shared<Tensor2d>("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->axpy(X, 2.0);
    X.reset();
    Tnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    timer_off("WijamT2");

}  // end ccsd_WijamT2_low

//======================================================================
//    WabefT2
//======================================================================
void DFOCC::ccsd_WabefT2_low() {
    // defs
    SharedTensor2d K, M, L, I, T, Tnew, U, Tau, W, X, Y, S, A;
    SharedTensor2d V, Vs, Ts, Va, Ta;

    timer_on("WabefT2");

    // t_ij^ab <= \sum_{ef} Tau_ij^ef <ab|ef>
    Tau = std::make_shared<Tensor2d>("Tau (IA|JB)", naoccA, navirA, naoccA, navirA);
    Tau->read_symm(psio_, PSIF_DFOCC_AMPS);
    // (+)Tau(ij, ab) = 1/2 (Tau_ij^ab + Tau_ji^ab) * (2 - \delta_{ab})
    // (-)Tau(ij, ab) = 1/2 (Tau_ij^ab - Tau_ji^ab) * (2 - \delta_{ab})
    U = std::make_shared<Tensor2d>("(+)Tau [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
    T = std::make_shared<Tensor2d>("(-)Tau [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
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
                    double value1 = 0.5 * perm * (Tau->get(ia, jb) + Tau->get(ja, ib));
                    double value2 = 0.5 * perm * (Tau->get(ia, jb) - Tau->get(ja, ib));
                    U->set(ij, ab, value1);
                    T->set(ij, ab, value2);
                }
            }
        }
    }
    Tau.reset();

    // Read B(Q,ab)
    L = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA);
    L->read(psio_, PSIF_DFOCC_INTS, true, true);
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (AB|Q)", navirA * navirA, nQ);
    K = L->transpose();
    L.reset();

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

}  // end ccsd_WabefT2_low

//======================================================================
//    Wabef2T2: 2nd version
//======================================================================
void DFOCC::ccsd_Wabef2T2_low() {
    // defs
    SharedTensor2d K, M, L, I, T, Tnew, U, Tau, W, X, Y, S, A;
    SharedTensor2d V, Vs, Ts, Va, Ta, J, T1;

    timer_on("WabefT2");

    // t_ij^ab <= \sum_{ef} Tau_ij^ef <ab|ef>
    Tau = std::make_shared<Tensor2d>("Tau (IA|JB)", naoccA, navirA, naoccA, navirA);
    Tau->read_symm(psio_, PSIF_DFOCC_AMPS);
    // (+)Tau(ij, ab) = 1/2 (Tau_ij^ab + Tau_ji^ab) * (2 - \delta_{ab})
    // (-)Tau(ij, ab) = 1/2 (Tau_ij^ab - Tau_ji^ab) * (2 - \delta_{ab})
    U = std::make_shared<Tensor2d>("(+)Tau [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
    T = std::make_shared<Tensor2d>("(-)Tau [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
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
                    double value1 = 0.5 * perm * (Tau->get(ia, jb) + Tau->get(ja, ib));
                    double value2 = 0.5 * perm * (Tau->get(ia, jb) - Tau->get(ja, ib));
                    U->set(ij, ab, value1);
                    T->set(ij, ab, value2);
                }
            }
        }
    }
    Tau.reset();

    // Read B(Q,ab)
    L = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA);
    L->read(psio_, PSIF_DFOCC_INTS, true, true);
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (AB|Q)", navirA * navirA, nQ);
    K = L->transpose();
    L.reset();
    // T1Q
    X = std::make_shared<Tensor2d>("T1 (Q|AB)", nQ, navirA, navirA);
    X->read(psio_, PSIF_DFOCC_AMPS);
    Y = std::make_shared<Tensor2d>("T1 (AB|Q)", navirA * navirA, nQ);
    Y = X->transpose();
    X.reset();
    X = std::make_shared<Tensor2d>("B-T1 (AB|Q)", navirA * navirA, nQ);
    X->copy(K);
    X->subtract(Y);
    Y.reset();
    // B(iaQ)
    M = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
    M->read(psio_, PSIF_DFOCC_INTS);
    L = std::make_shared<Tensor2d>("DF_BASIS_CC B (IA|Q)", naoccA * navirA, nQ);
    L = M->transpose();
    M.reset();
    // T(a,i)
    T1 = std::make_shared<Tensor2d>("T1 (A|I)", navirA, naoccA);
    T1 = t1A->transpose();

    // malloc
    I = std::make_shared<Tensor2d>("I[A] <BF|E>", navirA * navirA, navirA);
    J = std::make_shared<Tensor2d>("J[A] <MF|E>", naoccA * navirA, navirA);
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

        // Form J[a](bf,e) = \sum_{Q} B(bfQ)*[B(aeQ)-T(aeQ)] cost = V^4N/2
        I->contract(false, true, navirA * nb, navirA, nQ, K, X, 0, a * navirA * nQ, 1.0, 0.0);

        // Form J[a](mf,e) = \sum_{Q} B(mfQ)*B(aeQ) cost = OV^3N
        J->contract(false, true, navirA * naoccA, navirA, nQ, L, K, 0, a * navirA * nQ, 1.0, 0.0);

        // J[a](b,fe) -= \sum_{m} t(m,b) * J[a](mf,e)
        I->contract(false, false, nb, navirA * navirA, naoccA, T1, J, -1.0, 1.0);

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
    X.reset();
    Vs.reset();
    Va.reset();
    Ts.reset();
    Ta.reset();
    U.reset();
    T.reset();
    J.reset();
    L.reset();
    T1.reset();

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

}  // end ccsd_Wabef2T2_low

}  // namespace dfoccwave
}  // namespace psi
