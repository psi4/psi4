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

void DFOCC::ccsd_WmnijT2()
{
    // defs
    SharedTensor2d K, T, Tnew, U, Tau, W, X;
    SharedTensor2d M, L, I, Y, S, A;
    SharedTensor2d V, Vs, Va, Ts, Ta;

    timer_on("WmnijT2");

    // W_mnij = <mn|ij>
    W = SharedTensor2d(new Tensor2d("W <MN|IJ>", naoccA, naoccA, naoccA, naoccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|KL)", naoccA, naoccA, naoccA, naoccA));
    K->gemm(true, false, bQijA, bQijA, 1.0, 0.0);
    W->sort(1324, K, 1.0, 0.0);
    K.reset();

    // W_mnij += X(im,jn) + X(jn,im) += 2Xt(im,jn)
    // X_imjn = \sum_{Q} t_im^Q b_jn^Q
    X = SharedTensor2d(new Tensor2d("X <MN|IJ>", naoccA, naoccA, naoccA, naoccA));
    T = SharedTensor2d(new Tensor2d("T1 (Q|IJ)", nQ, naoccA, naoccA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(true, false, T, bQijA, 1.0, 0.0);
    T.reset();
    X->symmetrize();
    W->sort(2413, X, 2.0, 1.0);
    X.reset();

    // W_mnij = \sum_{ef} Tau_ij^ef <mn|ef>
    // (+)Tau(ij, ab) = 1/2 (Tau_ij^ab + Tau_ji^ab) * (2 - \delta_{ab})
    // (-)Tau(ij, ab) = 1/2 (Tau_ij^ab - Tau_ji^ab) * (2 - \delta_{ab})
    Tau = SharedTensor2d(new Tensor2d("Tau (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_tau_amps(Tau,t2);
    U = SharedTensor2d(new Tensor2d("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA));
    U->sort(1324, Tau, 1.0, 0.0);
    Tau.reset();
    Ts = SharedTensor2d(new Tensor2d("(+)tTau [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    Ta = SharedTensor2d(new Tensor2d("(-)tTau [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    Ts->symm_col_packed4(U);
    Ta->antisymm_col_packed4(U);
    U.reset();
    // Form <mn|ef>
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA));
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    L->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    Vs = SharedTensor2d(new Tensor2d("(+)V [M>=N|E>=F]", ntri_ijAA, ntri_abAA));
    Va = SharedTensor2d(new Tensor2d("(-)V [M>=N|E>=F]", ntri_ijAA, ntri_abAA));
    Vs->symm4(K);
    Va->antisymm4(K);
    K.reset();
    // Form S/A
    S = SharedTensor2d(new Tensor2d("S [M>=N|I>=J]", ntri_ijAA, ntri_ijAA));
    A = SharedTensor2d(new Tensor2d("A [M>=N|I>=J]", ntri_ijAA, ntri_ijAA));
    S->gemm(false, true, Vs, Ts, 1.0, 0.0);
    A->gemm(false, true, Va, Ta, 1.0, 0.0);
    Vs.reset();
    Va.reset();
    // add to W(mn,ij)
    #pragma omp parallel for
    for(int m = 0 ; m < naoccA; ++m){
        for(int n = 0 ; n < naoccA; ++n){
            int mn = index2(m,n);
            int mn2 = ij_idxAA->get(m,n);
            for(int i = 0 ; i < naoccA; ++i){
                for(int j = 0 ; j < naoccA; ++j){
                    int ij = index2(i,j);
                    int ij2 = ij_idxAA->get(i,j);
                    int perm1 = ( i > j ) ? 1 : -1;
                    int perm2 = ( m > n ) ? 1 : -1;
                    double value = S->get(mn,ij) + (perm1 * perm2 * A->get(mn,ij));
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
    Tau = SharedTensor2d(new Tensor2d("Tau (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_tau_amps(Tau,t2);
    U = SharedTensor2d(new Tensor2d("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA));
    U->sort(1324, Tau, 1.0, 0.0);
    Tau.reset();
    Ts->symm_row_packed4(U);
    Ta->antisymm_row_packed4(U);
    U.reset();

    // Form (+/-)W(m>=n, i>=j)
    Vs = SharedTensor2d(new Tensor2d("(+)W [M>=N|I>=J]", ntri_ijAA, ntri_ijAA));
    Va = SharedTensor2d(new Tensor2d("(-)W [M>=N|I>=J]", ntri_ijAA, ntri_ijAA));
    Vs->symm4(W);
    Va->antisymm4(W);
    W.reset();

    // Symmetric & Anti-symmetric contributions
    S = SharedTensor2d(new Tensor2d("S (I>=J, A>=B)", ntri_ijAA, ntri_abAA));
    A = SharedTensor2d(new Tensor2d("A (I>=J, A>=B)", ntri_ijAA, ntri_abAA));
    S->gemm(true, false, Vs, Ts, 1.0, 0.0);
    A->gemm(true, false, Va, Ta, 1.0, 0.0);
    Ts.reset();
    Ta.reset();
    Vs.reset();
    Va.reset();

    // T(ia,jb) <-- S(a>=b,i>=j) + A(a>=b,i>=j)
    Tnew = SharedTensor2d(new Tensor2d("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    #pragma omp parallel for
    for(int a = 0 ; a < navirA; ++a){
        for(int b = 0 ; b < navirA; ++b){
            int ab = index2(a,b);
            for(int i = 0 ; i < naoccA; ++i){
                int ia = ia_idxAA->get(i,a);
                for(int j = 0 ; j < naoccA; ++j){
                    int jb = ia_idxAA->get(j,b);
                    int ij = index2(i,j);
                    int perm1 = ( i > j ) ? 1 : -1;
                    int perm2 = ( a > b ) ? 1 : -1;
                    double value = S->get(ij,ab) + (perm1 * perm2 * A->get(ij,ab));
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

}// end ccsd_WmnijT2

//======================================================================
//    WmbejT2
//======================================================================
void DFOCC::ccsd_WmbejT2()
{
    // defs
    SharedTensor2d K, L, T, T1, Tnew, U, Tau, W, W2, X, Y;

    timer_on("WmbejT2");

    // W_mbej = W(me,jb)
    // W(me,jb) <= (me|jb)
    W = SharedTensor2d(new Tensor2d("W (ME|JB)", naoccA, navirA, naoccA, navirA));
    W->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);

    // W(me,jb) <= \sum_{Q} (t_jb^Q' + 1/2 T_jb^Q) b_me^Q
    T1 = SharedTensor2d(new Tensor2d("T1p (Q|IA)", nQ, naoccA, navirA));
    T1->read(psio_, PSIF_DFOCC_AMPS);
    T = SharedTensor2d(new Tensor2d("T2 (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    U = SharedTensor2d(new Tensor2d("T1 + T2/2 (Q|IA)", nQ, naoccA, navirA));
    U->copy(T);
    T.reset();
    U->scale(0.5);
    U->add(T1);
    T1.reset();
    W->gemm(true, false, bQiaA, U, 1.0, 1.0);
    U.reset();

    // W(me,jb) <= -1/2 \sum_{nf} t_jn^bf X(me,nf)
    // (mf|ne) = X(me,nf) (sort: 1432)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    X = SharedTensor2d(new Tensor2d("X (IA|JB)", naoccA, navirA, naoccA, navirA));
    X->sort(1432, K, 1.0, 0.0);
    K.reset();
    W->gemm(false, false, X, t2, -0.5, 1.0);
    X.reset();
    W->write(psio_, PSIF_DFOCC_AMPS);
    W.reset();

    // W_mbje = W'(me,jb)
    // W'(me,jb) <= <me|jb>
    W = SharedTensor2d(new Tensor2d("Wp (ME|JB)", naoccA, navirA, naoccA, navirA));
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|AB)", naoccA, naoccA, navirA, navirA));
    L->gemm(true, false, bQijA, bQabA, 1.0, 0.0);
    W->sort(1324, L, 1.0, 0.0);
    L.reset();

    // X(jm,be) <= -\sum_{Q} t_be^Q ( t_jm^Q + b_jm^Q )
    T = SharedTensor2d(new Tensor2d("T1 (Q|IJ)", nQ, naoccA, naoccA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA));
    K->copy(bQijA);
    K->add(T);
    T.reset();
    T = SharedTensor2d(new Tensor2d("T1 (Q|AB)", nQ, navirA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    X = SharedTensor2d(new Tensor2d("X (IJ|AB)", naoccA, naoccA, navirA, navirA));
    X->gemm(true, false, K, T, -1.0, 0.0);
    T.reset();
    K.reset();
    // X(jm,be) <= \sum_{Q} t_jm^Q b_be^Q
    T = SharedTensor2d(new Tensor2d("T1 (Q|IJ)", nQ, naoccA, naoccA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(true, false, T, bQabA, 1.0, 1.0);
    T.reset();
    // W'(me,jb) <= X(jm,be)
    W->sort(2413, X, 1.0, 1.0);
    X.reset();

    // W'(me,jb) <= -1/2 \sum_{nf} t_nj^bf X(me,nf)
    // (mf|ne) = X(me,nf) (sort: 1432)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    X = SharedTensor2d(new Tensor2d("X (IA|JB)", naoccA, navirA, naoccA, navirA));
    X->sort(1432, K, 1.0, 0.0);
    K.reset();
    T = SharedTensor2d(new Tensor2d("T2p (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_t2_prime_amps(T,t2);
    W->gemm(false, false, X, T, -0.5, 1.0);
    X.reset();
    T.reset();

    // t_ij^ab <= 1/2*C(ia,jb) + 1/2*C(jb,ia) + C(ja,ib) + C(ib,ja)
    // t_ij^ab <= Ct(ia,jb) + 2*Ct(ib,ja)
    // C(ia,jb) = -\sum_{me} t_mi^ae W'(me,jb) = -\sum_{me} T'(ia,me) W'(me,jb)
    U = SharedTensor2d(new Tensor2d("T2p (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_t2_prime_amps(U,t2);
    Y = SharedTensor2d(new Tensor2d("C2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Y->gemm(false, false, U, W, -1.0, 0.0);
    U.reset();
    X = SharedTensor2d(new Tensor2d("C2+D2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    X->sort(1432, Y, 1.0, 0.0);
    X->axpy(Y, 0.5);
    Y.reset();

    // t_ij^ab <= D(ia,jb) + D(jb,ia)
    // D_ij^ab = 1/2 \sum_{me} u_im^ae [2*W(me,jb) - W'(me,jb)]
    Y = SharedTensor2d(new Tensor2d("2*W-W' (ME|JB)", naoccA, navirA, naoccA, navirA));
    Y->axpy(W, -1.0);
    W.reset();
    W2 = SharedTensor2d(new Tensor2d("W (ME|JB)", naoccA, navirA, naoccA, navirA));
    W2->read(psio_, PSIF_DFOCC_AMPS);
    Y->axpy(W2, 2.0);
    W2.reset();
    U = SharedTensor2d(new Tensor2d("U2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_u2_amps(U,t2);
    X->gemm(false, false, U, Y, 0.5, 1.0);
    U.reset();
    Y.reset();
    X->symmetrize();
    Tnew = SharedTensor2d(new Tensor2d("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->axpy(X, 2.0);
    X.reset();
    Tnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    timer_off("WmbejT2");

}// end ccsd_WmbejT2

//======================================================================
//    WijamT2
//======================================================================
void DFOCC::ccsd_WijamT2()
{
    // defs
    SharedTensor2d K, M, L, I, T, Tnew, U, Tau, W, X, Y, S, A;
    SharedTensor2d V, Vs, Ts, Va, Ta;

    timer_on("WijamT2");

    // W_ijam = W(ij,am)
    W = SharedTensor2d(new Tensor2d("W (IJ|AM)", naoccA, naoccA, navirA, naoccA));

    // W_ijam = \sum_{ef} Tau_ij^ef <am|ef>
    Tau = SharedTensor2d(new Tensor2d("Tau (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_tau_amps(Tau,t2);
    // (+)Tau(ij, ab) = 1/2 (Tau_ij^ab + Tau_ji^ab) * (2 - \delta_{ab})
    // (-)Tau(ij, ab) = 1/2 (Tau_ij^ab - Tau_ji^ab) * (2 - \delta_{ab})
    U = SharedTensor2d(new Tensor2d("(+)Tau [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    T = SharedTensor2d(new Tensor2d("(-)Tau [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    #pragma omp parallel for
    for(int i = 0 ; i < naoccA; ++i){
        for(int j = 0 ; j <= i; ++j){
            int ij = index2(i,j);
            for(int a = 0 ; a < navirA; ++a){
                int ia = ia_idxAA->get(i,a);
                int ja = ia_idxAA->get(j,a);
                for(int b = 0 ; b <= a; ++b){
                    double perm = (a == b ? 1.0 : 2.0);
                    int ab = index2(a,b);
                    int jb = ia_idxAA->get(j,b);
                    int ib = ia_idxAA->get(i,b);
                    double value1 = 0.5 * perm * ( Tau->get(ia,jb) + Tau->get(ja,ib) );
                    double value2 = 0.5 * perm * ( Tau->get(ia,jb) - Tau->get(ja,ib) );
                    U->set(ij,ab,value1);
                    T->set(ij,ab,value2);
                }
            }
        }
    }
    Tau.reset();

    // Read B(Q,ab)
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (IA|Q)", naoccA * navirA, nQ));
    L = bQiaA->transpose();
    I = SharedTensor2d(new Tensor2d("I[M] <AE|F>", navirA * navirA, navirA));

    // Symmetric & Anti-symmetric contributions
    Vs = SharedTensor2d(new Tensor2d("(+)V[M] (A, E>=F)", navirA, ntri_abAA));
    Va = SharedTensor2d(new Tensor2d("(-)V[M] (A, E>=F)", navirA, ntri_abAA));
    S = SharedTensor2d(new Tensor2d("S (AM, I>=J)", navirA, ntri_ijAA));
    A = SharedTensor2d(new Tensor2d("A (AM, I>=J)", navirA, ntri_ijAA));
    // Main loop
    for(int m = 0 ; m < naoccA; ++m){

            // Form V[m](ae,f) = \sum_{Q} b(Q,ae) B(mfQ)
            I->contract(true, true, navirA * navirA, navirA, nQ, bQabA, L, 0, m*navirA*nQ, 1.0, 0.0);

            // Form (+)V[m](a, e>=f)
            #pragma omp parallel for
            for(int a = 0 ; a < navirA; ++a){
                for(int e = 0 ; e < navirA; ++e){
                    int ae = ab_idxAA->get(a,e);
                    for(int f = 0 ; f <= e; ++f){
                        int af = ab_idxAA->get(a,f);
                        int ef = index2(e,f);
                        double value1 = 0.5 * ( I->get(ae, f) + I->get(af, e) );
                        double value2 = 0.5 * ( I->get(ae, f) - I->get(af, e) );
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
            for(int a = 0 ; a < navirA; ++a){
                int am = ai_idxAA->get(a,m);
                for(int i = 0 ; i < naoccA; ++i){
                    for(int j = 0 ; j < naoccA; ++j){
                        int ij = ij_idxAA->get(i,j);
                        int ij2 = index2(i,j);
                        int perm = ( i > j ) ? 1 : -1;
                        double value = S->get(a,ij2) + (perm * A->get(a,ij2));
                        W->set(ij, am, value);
                    }
                }
            }

    }
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
    Y = SharedTensor2d(new Tensor2d("Y <IJ|AB>", naoccA, naoccA, navirA, navirA));
    Y->contract(false, false, naoccA * naoccA * navirA, navirA, naoccA, W, t1A, -1.0, 0.0);
    W.reset();
    X = SharedTensor2d(new Tensor2d("X (IA|JB)", naoccA, navirA, naoccA, navirA));
    X->sort(1324, Y, 1.0, 0.0);
    Y.reset();
    X->symmetrize();
    Tnew = SharedTensor2d(new Tensor2d("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->axpy(X, 2.0);
    X.reset();
    Tnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    timer_off("WijamT2");

}// end ccsd_WijamT2

//======================================================================
//    WijamT2: HIGH MEMORY
//======================================================================
void DFOCC::ccsd_WijamT2_high_mem()
{
    // defs
    SharedTensor2d K, M, L, I, T, Tnew, U, Tau, W, X, Y, S, A;
    SharedTensor2d V, Vs, Ts, Va, Ta, J;

    timer_on("WijamT2");

    // W_ijam = \sum_{ef} Tau_ij^ef <am|ef>
    Tau = SharedTensor2d(new Tensor2d("Tau (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_tau_amps(Tau,t2);
    // (+)Tau(ij, ab) = 1/2 (Tau_ij^ab + Tau_ji^ab) * (2 - \delta_{ab})
    // (-)Tau(ij, ab) = 1/2 (Tau_ij^ab - Tau_ji^ab) * (2 - \delta_{ab})
    U = SharedTensor2d(new Tensor2d("(+)Tau [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    T = SharedTensor2d(new Tensor2d("(-)Tau [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    #pragma omp parallel for
    for(int i = 0 ; i < naoccA; ++i){
        for(int j = 0 ; j <= i; ++j){
            int ij = index2(i,j);
            for(int a = 0 ; a < navirA; ++a){
                int ia = ia_idxAA->get(i,a);
                int ja = ia_idxAA->get(j,a);
                for(int b = 0 ; b <= a; ++b){
                    double perm = (a == b ? 1.0 : 2.0);
                    int ab = index2(a,b);
                    int jb = ia_idxAA->get(j,b);
                    int ib = ia_idxAA->get(i,b);
                    double value1 = 0.5 * perm * ( Tau->get(ia,jb) + Tau->get(ja,ib) );
                    double value2 = 0.5 * perm * ( Tau->get(ia,jb) - Tau->get(ja,ib) );
                    U->set(ij,ab,value1);
                    T->set(ij,ab,value2);
                }
            }
        }
    }
    Tau.reset();

    // Read in B(Q,a>=b)
    bQabA.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, ntri_abAA));
    K->read(psio_, PSIF_DFOCC_INTS);

    // J(mf,a>=e)
    J = SharedTensor2d(new Tensor2d("J (MF, A>=E)", naoccA*navirA, ntri_abAA));
    J->gemm(true, false, bQiaA, K, 1.0, 0.0);
    K.reset();

    // Build <am|ef>
    I = SharedTensor2d(new Tensor2d("I (AM,EF)", navirA, naoccA, navirA, navirA));
    #pragma omp parallel for
    for(int a = 0 ; a < navirA; ++a){
         for(int m = 0 ; m < naoccA; ++m){
             int am = ai_idxAA->get(a,m);
             for(int e = 0 ; e < navirA; ++e){
                 int ae = index2(a,e);
                 for(int f = 0 ; f < navirA; ++f){
                     int ef = ab_idxAA->get(e,f);
                     int mf = ia_idxAA->get(m,f);
		     I->set(am, ef, J->get(mf,ae));
		 }
	     }
	 }
    }
    J.reset();

    // Symmetric & Anti-symmetric contributions
    Vs = SharedTensor2d(new Tensor2d("(+)V (AM, E>=F)", navirA*naoccA, ntri_abAA));
    Va = SharedTensor2d(new Tensor2d("(-)V (AM, E>=F)", navirA*naoccA, ntri_abAA));
    #pragma omp parallel for
    for(int a = 0 ; a < navirA; ++a){
         for(int m = 0 ; m < naoccA; ++m){
             int am = ai_idxAA->get(a,m);
             for(int e = 0 ; e < navirA; ++e){
                 int ae = index2(a,e);
                 for(int f = 0 ; f <= e; ++f){
                     int ef = ab_idxAA->get(e,f);
                     int fe = ab_idxAA->get(f,e);
                     int ef2 = index2(e,f);
                     int mf = ia_idxAA->get(m,f);
                     double value1 = 0.5 * ( I->get(am, ef) + I->get(am, fe) );
                     double value2 = 0.5 * ( I->get(am, ef) - I->get(am, fe) );
		     Vs->set(am, ef2, value1);
		     Va->set(am, ef2, value2);
		 }
	     }
	 }
    }
    I.reset();

    // Form S(am, i>=j) = \sum_{e>=f} Tau(i>=j,e>=f) V(am, e>=f)
    S = SharedTensor2d(new Tensor2d("S (AM, I>=J)", navirA*naoccA, ntri_ijAA));
    A = SharedTensor2d(new Tensor2d("A (AM, I>=J)", navirA*naoccA, ntri_ijAA));
    S->gemm(false, true, Vs, U, 1.0, 0.0);
    A->gemm(false, true, Va, T, 1.0, 0.0);
    Vs.reset();
    Va.reset();
    U.reset();
    T.reset();

    // Form S(am,ij) & A(am,ij)-->W(ijam)
    // W_ijam = W(ij,am)
    W = SharedTensor2d(new Tensor2d("W (IJ|AM)", naoccA, naoccA, navirA, naoccA));
    #pragma omp parallel for
    for(int a = 0 ; a < navirA; ++a){
        for(int m = 0 ; m < naoccA; ++m){
            int am = ai_idxAA->get(a,m);
            for(int i = 0 ; i < naoccA; ++i){
                for(int j = 0 ; j < naoccA; ++j){
                    int ij = ij_idxAA->get(i,j);
                    int ij2 = index2(i,j);
                    int perm = ( i > j ) ? 1 : -1;
                    double value = S->get(am,ij2) + (perm * A->get(am,ij2));
                    W->set(ij, am, value);
		}
	    }
	}
    }
    S.reset();
    A.reset();

    // Y(ij,ab) <= -\sum_{m} t_m^b Wijam
    // X(ia,jb) = Y(ij,ab)
    // t_ij^ab <= X(ia,jb) + X(jb,ia)
    Y = SharedTensor2d(new Tensor2d("Y <IJ|AB>", naoccA, naoccA, navirA, navirA));
    Y->contract(false, false, naoccA * naoccA * navirA, navirA, naoccA, W, t1A, -1.0, 0.0);
    W.reset();
    X = SharedTensor2d(new Tensor2d("X (IA|JB)", naoccA, navirA, naoccA, navirA));
    X->sort(1324, Y, 1.0, 0.0);
    Y.reset();
    X->symmetrize();
    Tnew = SharedTensor2d(new Tensor2d("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->axpy(X, 2.0);
    X.reset();
    Tnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    // Read in B(Q,ab)
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);

    timer_off("WijamT2");

}// end ccsd_WijamT2_high_mem

//======================================================================
//    WabefT2
//======================================================================
void DFOCC::ccsd_WabefT2()
{
    // defs
    SharedTensor2d K, M, L, I, T, Tnew, U, Tau, W, X, Y, S, A;
    SharedTensor2d V, Vs, Ts, Va, Ta;

    timer_on("WabefT2");

    // t_ij^ab <= \sum_{ef} Tau_ij^ef <ab|ef>
    Tau = SharedTensor2d(new Tensor2d("Tau (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_tau_amps(Tau,t2);
    // (+)Tau(ij, ab) = 1/2 (Tau_ij^ab + Tau_ji^ab) * (2 - \delta_{ab})
    // (-)Tau(ij, ab) = 1/2 (Tau_ij^ab - Tau_ji^ab) * (2 - \delta_{ab})
    U = SharedTensor2d(new Tensor2d("(+)Tau [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    T = SharedTensor2d(new Tensor2d("(-)Tau [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    #pragma omp parallel for
    for(int i = 0 ; i < naoccA; ++i){
        for(int j = 0 ; j <= i; ++j){
            int ij = index2(i,j);
            for(int a = 0 ; a < navirA; ++a){
                int ia = ia_idxAA->get(i,a);
                int ja = ia_idxAA->get(j,a);
                for(int b = 0 ; b <= a; ++b){
                    double perm = (a == b ? 1.0 : 2.0);
                    int ab = index2(a,b);
                    int jb = ia_idxAA->get(j,b);
                    int ib = ia_idxAA->get(i,b);
                    double value1 = 0.5 * perm * ( Tau->get(ia,jb) + Tau->get(ja,ib) );
                    double value2 = 0.5 * perm * ( Tau->get(ia,jb) - Tau->get(ja,ib) );
                    U->set(ij,ab,value1);
                    T->set(ij,ab,value2);
                }
            }
        }
    }
    Tau.reset();

    // Read B(Q,ab)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (AB|Q)", navirA * navirA, nQ));
    K = bQabA->transpose();
    bQabA.reset();

    // malloc
    I = SharedTensor2d(new Tensor2d("I[A] <BF|E>", navirA * navirA, navirA));
    Vs = SharedTensor2d(new Tensor2d("(+)V[A] (B, E>=F)", navirA, ntri_abAA));
    Va = SharedTensor2d(new Tensor2d("(-)V[A] (B, E>=F)", navirA, ntri_abAA));
    Ts = SharedTensor2d(new Tensor2d("(+)T[A] (B, I>=J)", navirA, ntri_ijAA));
    Ta = SharedTensor2d(new Tensor2d("(-)T[B] (B, I>=J)", navirA, ntri_ijAA));

    // Symmetric & Anti-symmetric contributions
    S = SharedTensor2d(new Tensor2d("S (A>=B, I>=J)", ntri_abAA, ntri_ijAA));
    A = SharedTensor2d(new Tensor2d("A (A>=B, I>=J)", ntri_abAA, ntri_ijAA));
    // Main loop
    for(int a = 0 ; a < navirA; ++a){
            int nb = a+1;

            // Form V[a](bf,e) = \sum_{Q} B(bfQ)*B(aeQ) cost = V^4N/2
            I->contract(false, true, navirA*nb, navirA, nQ, K, K, 0, a*navirA*nQ, 1.0, 0.0);

            // Form (+)V[a](b, e>=f)
            #pragma omp parallel for
            for(int b = 0 ; b <= a; ++b){
                for(int e = 0 ; e < navirA; ++e){
                    int be = e + (b * navirA);
                    for(int f = 0 ; f <= e; ++f){
                        int ef = index2(e,f);
                        int bf = f + (b * navirA);
                        double value1 = 0.5 * ( I->get(bf, e) + I->get(be, f) );
                        double value2 = 0.5 * ( I->get(bf, e) - I->get(be, f) );
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
            for(int b = 0 ; b <=a; ++b){
                int ab = index2(a,b);
                for(int i = 0 ; i < naoccA; ++i){
                    for(int j = 0 ; j <= i; ++j){
                        int ij = index2(i,j);
                        S->add(ab, ij, Ts->get(b,ij));
                        A->add(ab, ij, Ta->get(b,ij));
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
    Tnew = SharedTensor2d(new Tensor2d("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    #pragma omp parallel for
    for(int a = 0 ; a < navirA; ++a){
        for(int b = 0 ; b < navirA; ++b){
            int ab = index2(a,b);
            for(int i = 0 ; i < naoccA; ++i){
                int ia = ia_idxAA->get(i,a);
                for(int j = 0 ; j < naoccA; ++j){
                    int jb = ia_idxAA->get(j,b);
                    int ij = index2(i,j);
                    int perm1 = ( i > j ) ? 1 : -1;
                    int perm2 = ( a > b ) ? 1 : -1;
                    double value = S->get(ab,ij) + (perm1 * perm2 * A->get(ab,ij));
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
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);

    timer_off("WabefT2");

}// end ccsd_WabefT2

//======================================================================
//    WabefT2: 2nd version, consistent with Scuseria & Schaefer
//======================================================================
void DFOCC::ccsd_Wabef2T2()
{
    // defs
    SharedTensor2d K, M, L, I, T, Tnew, U, Tau, W, X, Y, S, A;
    SharedTensor2d V, Vs, Ts, Va, Ta, J, T1;

    timer_on("WabefT2");

    // t_ij^ab <= \sum_{ef} Tau_ij^ef W_abef
    Tau = SharedTensor2d(new Tensor2d("Tau (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_tau_amps(Tau,t2);
    // (+)Tau(ij, ab) = 1/2 (Tau_ij^ab + Tau_ji^ab) * (2 - \delta_{ab})
    // (-)Tau(ij, ab) = 1/2 (Tau_ij^ab - Tau_ji^ab) * (2 - \delta_{ab})
    U = SharedTensor2d(new Tensor2d("(+)Tau [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    T = SharedTensor2d(new Tensor2d("(-)Tau [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    #pragma omp parallel for
    for(int i = 0 ; i < naoccA; ++i){
        for(int j = 0 ; j <= i; ++j){
            int ij = index2(i,j);
            for(int a = 0 ; a < navirA; ++a){
                int ia = ia_idxAA->get(i,a);
                int ja = ia_idxAA->get(j,a);
                for(int b = 0 ; b <= a; ++b){
                    double perm = (a == b ? 1.0 : 2.0);
                    int ab = index2(a,b);
                    int jb = ia_idxAA->get(j,b);
                    int ib = ia_idxAA->get(i,b);
                    double value1 = 0.5 * perm * ( Tau->get(ia,jb) + Tau->get(ja,ib) );
                    double value2 = 0.5 * perm * ( Tau->get(ia,jb) - Tau->get(ja,ib) );
                    U->set(ij,ab,value1);
                    T->set(ij,ab,value2);
                }
            }
        }
    }
    Tau.reset();

    // Read B(Q,ab) and B(Q,ia)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (AB|Q)", navirA * navirA, nQ));
    K = bQabA->transpose();
    //bQabA.reset();
    // T1Q
    X = SharedTensor2d(new Tensor2d("T1 (Q|AB)", nQ, navirA, navirA));
    X->read(psio_, PSIF_DFOCC_AMPS);
    Y = SharedTensor2d(new Tensor2d("T1 (AB|Q)", navirA * navirA, nQ));
    Y = X->transpose();
    X.reset();
    X = SharedTensor2d(new Tensor2d("B-T1 (AB|Q)", navirA * navirA, nQ));
    X->copy(K);
    X->subtract(Y);
    Y.reset();
    // B(iaQ)
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA));
    M->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (IA|Q)", naoccA * navirA, nQ));
    L = M->transpose();
    M.reset();
    // T(a,i)
    T1 = SharedTensor2d(new Tensor2d("T1 (A|I)", navirA, naoccA));
    T1 = t1A->transpose();

    // malloc
    I = SharedTensor2d(new Tensor2d("I[A] <BF|E>", navirA * navirA, navirA));
    J = SharedTensor2d(new Tensor2d("J[A] <MF|E>", naoccA * navirA, navirA));
    Vs = SharedTensor2d(new Tensor2d("(+)V[A] (B, E>=F)", navirA, ntri_abAA));
    Va = SharedTensor2d(new Tensor2d("(-)V[A] (B, E>=F)", navirA, ntri_abAA));
    Ts = SharedTensor2d(new Tensor2d("(+)T[A] (B, I>=J)", navirA, ntri_ijAA));
    Ta = SharedTensor2d(new Tensor2d("(-)T[B] (B, I>=J)", navirA, ntri_ijAA));

    // Symmetric & Anti-symmetric contributions
    S = SharedTensor2d(new Tensor2d("S (A>=B, I>=J)", ntri_abAA, ntri_ijAA));
    A = SharedTensor2d(new Tensor2d("A (A>=B, I>=J)", ntri_abAA, ntri_ijAA));
    // Main loop
    for(int a = 0 ; a < navirA; ++a){
        int nb = a+1;

        // Form J[a](bf,e) = \sum_{Q} B(bfQ)*[B(aeQ)-T(aeQ)] cost = V^4N/2
        I->contract(false, true, navirA*nb, navirA, nQ, K, X, 0, a*navirA*nQ, 1.0, 0.0);

        // Form J[a](mf,e) = \sum_{Q} B(mfQ)*B(aeQ) cost = OV^3N
        J->contract(false, true, navirA*naoccA, navirA, nQ, L, K, 0, a*navirA*nQ, 1.0, 0.0);

        // J[a](b,fe) -= \sum_{m} t(b,m) * J[a](m,fe)
        I->contract(false, false, nb, navirA*navirA, naoccA, T1, J, -1.0, 1.0);

        // Form (+)V[a](b, e>=f)
        #pragma omp parallel for
        for(int b = 0 ; b <= a; ++b){
            for(int e = 0 ; e < navirA; ++e){
                int be = e + (b * navirA);
                for(int f = 0 ; f <= e; ++f){
                    int ef = index2(e,f);
                    int bf = f + (b * navirA);
                    double value1 = 0.5 * ( I->get(bf, e) + I->get(be, f) );
                    double value2 = 0.5 * ( I->get(bf, e) - I->get(be, f) );
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
        for(int b = 0 ; b <=a; ++b){
            int ab = index2(a,b);
            for(int i = 0 ; i < naoccA; ++i){
                for(int j = 0 ; j <= i; ++j){
                    int ij = index2(i,j);
                     S->add(ab, ij, Ts->get(b,ij));
                     A->add(ab, ij, Ta->get(b,ij));
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
    Tnew = SharedTensor2d(new Tensor2d("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    #pragma omp parallel for
    for(int a = 0 ; a < navirA; ++a){
        for(int b = 0 ; b < navirA; ++b){
            int ab = index2(a,b);
            for(int i = 0 ; i < naoccA; ++i){
                int ia = ia_idxAA->get(i,a);
                for(int j = 0 ; j < naoccA; ++j){
                    int jb = ia_idxAA->get(j,b);
                    int ij = index2(i,j);
                    int perm1 = ( i > j ) ? 1 : -1;
                    int perm2 = ( a > b ) ? 1 : -1;
                    double value = S->get(ab,ij) + (perm1 * perm2 * A->get(ab,ij));
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
    //bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    //bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);

    timer_off("WabefT2");

}// end ccsd_Wabef2T2

//======================================================================
//    WabefT2: High memory version
//======================================================================
void DFOCC::ccsd_WabefT2_high_mem()
{

    // defs
    SharedTensor2d K, M, L, I, T, Tnew, U, Tau, W, X, Y, S, A;
    SharedTensor2d V, Vs, Ts, Va, Ta, J, T1;

    timer_on("WabefT2");

    // t_ij^ab <= \sum_{ef} Tau_ij^ef W_abef
    Tau = SharedTensor2d(new Tensor2d("Tau (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_tau_amps(Tau,t2);
    // (+)Tau(ij, ab) = 1/2 (Tau_ij^ab + Tau_ji^ab) * (2 - \delta_{ab})
    // (-)Tau(ij, ab) = 1/2 (Tau_ij^ab - Tau_ji^ab) * (2 - \delta_{ab})
    U = SharedTensor2d(new Tensor2d("(+)Tau [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    T = SharedTensor2d(new Tensor2d("(-)Tau [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    #pragma omp parallel for
    for(int i = 0 ; i < naoccA; ++i){
        for(int j = 0 ; j <= i; ++j){
            int ij = index2(i,j);
            for(int a = 0 ; a < navirA; ++a){
                int ia = ia_idxAA->get(i,a);
                int ja = ia_idxAA->get(j,a);
                for(int b = 0 ; b <= a; ++b){
                    double perm = (a == b ? 1.0 : 2.0);
                    int ab = index2(a,b);
                    int jb = ia_idxAA->get(j,b);
                    int ib = ia_idxAA->get(i,b);
                    double value1 = 0.5 * perm * ( Tau->get(ia,jb) + Tau->get(ja,ib) );
                    double value2 = 0.5 * perm * ( Tau->get(ia,jb) - Tau->get(ja,ib) );
                    U->set(ij,ab,value1);
                    T->set(ij,ab,value2);
                }
            }
        }
    }
    Tau.reset();

    // Read B(Q,a>=b)
    bQabA.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, ntri_abAA));
    K->read(psio_, PSIF_DFOCC_INTS);

    // Form (A>=E|B>=F) : cost = V4N/4
    J = SharedTensor2d(new Tensor2d("J (A>=E|B>=F)", ntri_abAA, ntri_abAA));
    J->gemm(true, false, K, K, 1.0, 0.0);
    K.reset();

    // malloc
    Vs = SharedTensor2d(new Tensor2d("(+)V[A] (B, E>=F)", navirA, ntri_abAA));
    Va = SharedTensor2d(new Tensor2d("(-)V[A] (B, E>=F)", navirA, ntri_abAA));
    Ts = SharedTensor2d(new Tensor2d("(+)T[A] (B, I>=J)", navirA, ntri_ijAA));
    Ta = SharedTensor2d(new Tensor2d("(-)T[B] (B, I>=J)", navirA, ntri_ijAA));

    // Symmetric & Anti-symmetric contributions
    S = SharedTensor2d(new Tensor2d("S (A>=B, I>=J)", ntri_abAA, ntri_ijAA));
    A = SharedTensor2d(new Tensor2d("A (A>=B, I>=J)", ntri_abAA, ntri_ijAA));
    // Main loop
    for(int a = 0 ; a < navirA; ++a){
        int nb = a+1;

        // Form (+)V[a](b, e>=f)
        #pragma omp parallel for
        for(int b = 0 ; b <= a; ++b){
            for(int e = 0 ; e < navirA; ++e){
                int ae = index2(a,e);
                int be = index2(b,e);
                for(int f = 0 ; f <= e; ++f){
                    int ef = index2(e,f);
                    int bf = index2(b,f);
                    int af = index2(a,f);
                    double value1 = 0.5 * ( J->get(ae, bf) + J->get(af, be) );
                    double value2 = 0.5 * ( J->get(ae, bf) - J->get(af, be) );
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
        for(int b = 0 ; b <=a; ++b){
            int ab = index2(a,b);
            for(int i = 0 ; i < naoccA; ++i){
                for(int j = 0 ; j <= i; ++j){
                    int ij = index2(i,j);
                     S->add(ab, ij, Ts->get(b,ij));
                     A->add(ab, ij, Ta->get(b,ij));
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
    Tnew = SharedTensor2d(new Tensor2d("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    #pragma omp parallel for
    for(int a = 0 ; a < navirA; ++a){
        for(int b = 0 ; b < navirA; ++b){
            int ab = index2(a,b);
            for(int i = 0 ; i < naoccA; ++i){
                int ia = ia_idxAA->get(i,a);
                for(int j = 0 ; j < naoccA; ++j){
                    int jb = ia_idxAA->get(j,b);
                    int ij = index2(i,j);
                    int perm1 = ( i > j ) ? 1 : -1;
                    int perm2 = ( a > b ) ? 1 : -1;
                    double value = S->get(ab,ij) + (perm1 * perm2 * A->get(ab,ij));
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
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);

    timer_off("WabefT2");

}// end ccsd_WabefT2_high_mem

//======================================================================
//    CD-WabefT2:
//======================================================================
void DFOCC::ccsd_WabefT2_cd()
{

    timer_on("WabefT2");

    // defs
    SharedTensor2d K, M, L, I, T, Tnew, U, Tau, W, X, Y, S, A;
    SharedTensor2d V, Vs, Ts, Va, Ta, bQ;


    // t_ij^ab <= \sum_{ef} Tau_ij^ef <ab|ef>
    Tau = SharedTensor2d(new Tensor2d("Tau (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_tau_amps(Tau,t2);
    // (+)Tau(ij, ab) = 1/2 (Tau_ij^ab + Tau_ji^ab) * (2 - \delta_{ab})
    // (-)Tau(ij, ab) = 1/2 (Tau_ij^ab - Tau_ji^ab) * (2 - \delta_{ab})
    U = SharedTensor2d(new Tensor2d("(+)Tau [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    T = SharedTensor2d(new Tensor2d("(-)Tau [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    #pragma omp parallel for
    for(int i = 0 ; i < naoccA; ++i){
        for(int j = 0 ; j <= i; ++j){
            int ij = index2(i,j);
            for(int a = 0 ; a < navirA; ++a){
                int ia = ia_idxAA->get(i,a);
                int ja = ia_idxAA->get(j,a);
                for(int b = 0 ; b <= a; ++b){
                    double perm = (a == b ? 1.0 : 2.0);
                    int ab = index2(a,b);
                    int jb = ia_idxAA->get(j,b);
                    int ib = ia_idxAA->get(i,b);
                    double value1 = 0.5 * perm * ( Tau->get(ia,jb) + Tau->get(ja,ib) );
                    double value2 = 0.5 * perm * ( Tau->get(ia,jb) - Tau->get(ja,ib) );
                    U->set(ij,ab,value1);
                    T->set(ij,ab,value2);
                }
            }
        }
    }
    Tau.reset();

    // Read B(Q,ab)
    bQ = SharedTensor2d(new Tensor2d("L <Q|AB>", nQ_cd, navirA, navirA));
    bQ->read(psio_, PSIF_DFOCC_INTS, true, true);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (AB|Q)", navirA * navirA, nQ_cd));
    K = bQ->transpose();
    bQ.reset();

    // malloc
    I = SharedTensor2d(new Tensor2d("I[A] <BF|E>", navirA * navirA, navirA));
    Vs = SharedTensor2d(new Tensor2d("(+)V[A] (B, E>=F)", navirA, ntri_abAA));
    Va = SharedTensor2d(new Tensor2d("(-)V[A] (B, E>=F)", navirA, ntri_abAA));
    Ts = SharedTensor2d(new Tensor2d("(+)T[A] (B, I>=J)", navirA, ntri_ijAA));
    Ta = SharedTensor2d(new Tensor2d("(-)T[B] (B, I>=J)", navirA, ntri_ijAA));

    // Symmetric & Anti-symmetric contributions
    S = SharedTensor2d(new Tensor2d("S (A>=B, I>=J)", ntri_abAA, ntri_ijAA));
    A = SharedTensor2d(new Tensor2d("A (A>=B, I>=J)", ntri_abAA, ntri_ijAA));
    // Main loop
    for(int a = 0 ; a < navirA; ++a){
            int nb = a+1;

            // Form V[a](bf,e) = \sum_{Q} B(bfQ)*B(aeQ) cost = V^4N/2
            I->contract(false, true, navirA*nb, navirA, nQ_cd, K, K, 0, a*navirA*nQ_cd, 1.0, 0.0);

            // Form (+)V[a](b, e>=f)
            #pragma omp parallel for
            for(int b = 0 ; b <= a; ++b){
                for(int e = 0 ; e < navirA; ++e){
                    int be = e + (b * navirA);
                    for(int f = 0 ; f <= e; ++f){
                        int ef = index2(e,f);
                        int bf = f + (b * navirA);
                        double value1 = 0.5 * ( I->get(bf, e) + I->get(be, f) );
                        double value2 = 0.5 * ( I->get(bf, e) - I->get(be, f) );
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
            for(int b = 0 ; b <=a; ++b){
                int ab = index2(a,b);
                for(int i = 0 ; i < naoccA; ++i){
                    for(int j = 0 ; j <= i; ++j){
                        int ij = index2(i,j);
                        S->add(ab, ij, Ts->get(b,ij));
                        A->add(ab, ij, Ta->get(b,ij));
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
    Tnew = SharedTensor2d(new Tensor2d("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    #pragma omp parallel for
    for(int a = 0 ; a < navirA; ++a){
        for(int b = 0 ; b < navirA; ++b){
            int ab = index2(a,b);
            for(int i = 0 ; i < naoccA; ++i){
                int ia = ia_idxAA->get(i,a);
                for(int j = 0 ; j < naoccA; ++j){
                    int jb = ia_idxAA->get(j,b);
                    int ij = index2(i,j);
                    int perm1 = ( i > j ) ? 1 : -1;
                    int perm2 = ( a > b ) ? 1 : -1;
                    double value = S->get(ab,ij) + (perm1 * perm2 * A->get(ab,ij));
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

}// end ccsd_WabefT2_cd

//======================================================================
//    LDL-WabefT2
//======================================================================
void DFOCC::ccsd_WabefT2_ldl()
{

    // defs
    SharedTensor2d L, T, Tnew, U, Tau, X;

    timer_on("WabefT2");

    // t_ij^ab <= \sum_{ef} Tau_ij^ef K_abef
    T = SharedTensor2d(new Tensor2d("Tau (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_tau_amps(T,t2);
    Tau = SharedTensor2d(new Tensor2d("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA));
    Tau->sort(1324, T, 1.0, 1.0);
    T.reset();

    // X(ij,Q) = \sum(ef) Tau(ij,ef) * U(Q,ef)
    U = SharedTensor2d(new Tensor2d("U <Q|CD>", nQ_cd, navirA*navirA));
    U->read(psio_, PSIF_DFOCC_INTS);
    X = SharedTensor2d(new Tensor2d("X <IJ|Q>", naoccA*naoccA, nQ_cd));
    X->gemm(false, true, Tau, U, 1.0, 0.0);
    Tau.reset();
    U.reset();

    // t_ij^ab <= \sum_{Q} X(ij,Q) L(ab,Q)
    L = SharedTensor2d(new Tensor2d("L <AB|Q>", navirA*navirA, nQ_cd));
    L->read(psio_, PSIF_DFOCC_INTS);
    T = SharedTensor2d(new Tensor2d("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    T->gemm(false, true, X, L, 1.0, 0.0);
    L.reset();
    X.reset();

    // T(IJ,AB) --> T(IA,JB)
    Tnew = SharedTensor2d(new Tensor2d("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->sort(1324, T, 1.0, 1.0);
    T.reset();
    Tnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    timer_off("WabefT2");

}// end ccsd_WabefT2_ldl

//======================================================================
//    WabefT2: AO Basis: EXPERIMENTAL & NOT FINISHED
//======================================================================
/*
void DFOCC::ccsd_WabefT2_ao_basis()
{

    // defs
    SharedTensor2d K, M, L, I, T, Tnew, U, Tau, W, X, Y, S, A;
    SharedTensor2d V, Vs, Ts, Va, Ta, G, D, Tt, J;

    timer_on("WabefT2");

    // Read SO integrals
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mn)", nQ, nso_, nso_));
    bQso->read(psio_, PSIF_DFOCC_INTS, true, true);
    G = SharedTensor2d(new Tensor2d("DF_BASIS_CC J (mn|mn)", nso_, nso_));

    #pragma omp parallel for
    for(int m = 0 ; m < nso_; ++m){
        for(int n = 0 ; n < nso_; ++n){
            int mn = n + (m * nso_);
	    double sum = 0.0;
            for(int Q = 0 ; Q < nQ; ++Q){
                sum += bQso->get(Q,mn) * bQso->get(Q,mn);
            }
	    double value = sqrt(sum);
	    G->set(m, n, value);
        }
    }
    //bQso.reset();

    // Transpose B
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (mn|Q)", nso_*nso_, nQ));
    K->trans(bQso);
    bQso.reset();

    // Tau_ij^es = \sum(f) Tau_ij^ef C_sf
    Tau = SharedTensor2d(new Tensor2d("Tau (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_tau_amps(Tau,t2);
    U = SharedTensor2d(new Tensor2d("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA));
    U->sort(1324, Tau, 1.0, 0.0);
    Tau.reset();
    Tau = SharedTensor2d(new Tensor2d("Tau <IJ|ES>", naoccA, naoccA, navirA, nso_));
    Tau->contract(false, true, naoccA*naoccA*navirA, nso_, navirA, U, CavirA, 1.0, 0.0);
    U.reset();

    // Tau_ij^ls = \sum(e) Tau_ij^es C_le
    U = SharedTensor2d(new Tensor2d("Tau <IJ|SE>", naoccA, naoccA, nso_, navirA));
    U->sort(1243, Tau, 1.0, 0.0);
    Tau.reset();
    X = SharedTensor2d(new Tensor2d("Tau <IJ|SL>", naoccA, naoccA, nso_, nso_));
    X->contract(false, true, naoccA*naoccA*nso_, nso_, navirA, U, CavirA, 1.0, 0.0);
    U.reset();
    Tau = SharedTensor2d(new Tensor2d("Tau <IJ|LS>", naoccA, naoccA, nso_, nso_));
    Tau->sort(1243, X, 1.0, 0.0);
    X.reset();
    //U->cont424("IJLS", "IJES", "LE", false, Tau, CavirA, 1.0, 0.0);// do not work!

    // tTau_ij^mn <= \sum_{ls} Tau_ij^ls <mn|ls>
    // (+)Tau(ij, ls) = 1/2 (Tau_ij^ls + Tau_ij^sl) * (2 - \delta_{ls})
    // (-)Tau(ij, ls) = 1/2 (Tau_ij^ls - Tau_ij^sl) * (2 - \delta_{ls})
    U = SharedTensor2d(new Tensor2d("(+)Tau [I>=J|L>=S]", ntri_ijAA, ntri_so));
    T = SharedTensor2d(new Tensor2d("(-)Tau [I>=J|L>=S]", ntri_ijAA, ntri_so));
    U->symm_col_packed4(Tau);
    T->antisymm_col_packed4(Tau);
    D = SharedTensor2d(new Tensor2d("Tau[LS] (IJ)", naoccA, naoccA));
    //Tau.reset();
    //Tau->print();

    // malloc
    I = SharedTensor2d(new Tensor2d("I[m] <ns|l>", nso_ * nso_, nso_));
    Vs = SharedTensor2d(new Tensor2d("(+)V[m] (n, l>=s)", nso_, ntri_so));
    Va = SharedTensor2d(new Tensor2d("(-)V[m] (n, l>=s)", nso_, ntri_so));
    Ts = SharedTensor2d(new Tensor2d("(+)T[m] (n, I>=J)", nso_, ntri_ijAA));
    Ta = SharedTensor2d(new Tensor2d("(-)T[m] (n, I>=J)", nso_, ntri_ijAA));

    // Symmetric & Anti-symmetric contributions
    S = SharedTensor2d(new Tensor2d("S (m>=n, I>=J)", ntri_so, ntri_ijAA));
    A = SharedTensor2d(new Tensor2d("A (m>=n, I>=J)", ntri_so, ntri_ijAA));
    J = SharedTensor2d(new Tensor2d("J <ML|NS>", 1, 1));
    // Main loop
    // Here: a = mu, b = nu, e = lambda, f = sigma
    ndf_nz = 0;
    double tau_max = Tau->get_max_element();
    #pragma omp parallel for
    for(int a = 0 ; a < nso_; ++a){
            int nb = a+1;

	    // Form V[a](bf,e)
            #pragma omp parallel for
            for(int b = 0 ; b <= a; ++b){
                for(int e = 0 ; e < nso_; ++e){
                    int ae = e + (a * nso_);
                    for(int f = 0 ; f < nso_; ++f){
                        int bf = f + (b * nso_);

  	                // prescreening
			double Gaebf = tau_max * G->get(a,e) * G->get(b,f);
			if (fabs(Gaebf) > int_cutoff_) {
			      ndf_nz++;
                              J->contract(false, true, 1, 1, nQ, K, K, ae*nQ, bf*nQ, 1.0, 0.0);
			      I->set(bf, e, J->get(0,0));
			}// end if
			else I->set(bf, e, 0.0);

                    }
                }
            }

            // Form (+)V[a](b, e>=f)
            #pragma omp parallel for
            for(int b = 0 ; b <= a; ++b){
                for(int e = 0 ; e < nso_; ++e){
                    int be = e + (b * nso_);
                    for(int f = 0 ; f <= e; ++f){
                        int ef = index2(e,f);
                        int bf = f + (b * nso_);
                        double value1 = 0.5 * ( I->get(bf, e) + I->get(be, f) );
                        double value2 = 0.5 * ( I->get(bf, e) - I->get(be, f) );
                        Vs->set(b, ef, value1);
                        Va->set(b, ef, value2);
                    }
                }
            }

            // Form T[a](b, i>=j) = \sum_{e>=f} Tau(i>=j,e>=f) V[a](b, e>=f)
            //#pragma omp single
            Ts->contract(false, true, nb, ntri_ijAA, ntri_so, Vs, U, 1.0, 0.0);
            Ta->contract(false, true, nb, ntri_ijAA, ntri_so, Va, T, 1.0, 0.0);

            // Form S(ij,ab) & A(ij,ab)
            #pragma omp parallel for
            for(int b = 0 ; b <=a; ++b){
                int ab = index2(a,b);
                for(int i = 0 ; i < naoccA; ++i){
                    for(int j = 0 ; j <= i; ++j){
                        int ij = index2(i,j);
                        S->add(ab, ij, Ts->get(b,ij));
                        A->add(ab, ij, Ta->get(b,ij));
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
    Tau.reset();
    D.reset();
    J.reset();
    outfile->Printf("\tNumber of prescreened integrals in Wabef    : %3d\n", ndf_nz);

    // Tt(ij,mn) <-- S(m>=n,i>=j) + A(m>=n,i>=j)
    Tau = SharedTensor2d(new Tensor2d("Tau <IJ|MN>", naoccA, naoccA, nso_, nso_));
    // Here: a = mu, b = nu, e = lambda, f = sigma
    #pragma omp parallel for
    for(int a = 0 ; a < nso_; ++a){
        for(int b = 0 ; b < nso_; ++b){
            int ab = index2(a,b);
	    int ab2 = b + (a*nso_);
            for(int i = 0 ; i < naoccA; ++i){
                for(int j = 0 ; j < naoccA; ++j){
                    int ij2 = ij_idxAA->get(i,j);
                    int ij = index2(i,j);
                    int perm1 = ( i > j ) ? 1 : -1;
                    int perm2 = ( a > b ) ? 1 : -1;
                    double value = S->get(ab,ij) + (perm1 * perm2 * A->get(ab,ij));
                    Tau->set(ij2, ab2, value);
                }
            }
        }
    }
    S.reset();
    A.reset();

    // Taut_ij^mb = \sum(n) Taut_ij^mn C(n,b)
    T = SharedTensor2d(new Tensor2d("Tau <IJ|MB>", naoccA, naoccA, nso_, navirA));
    T->contract(false, false, naoccA*naoccA*nso_, navirA, nso_, Tau, CavirA, 1.0, 0.0);
    Tau.reset();


    // T_ij^ab = \sum(m) Taut_ij^mb C(m,a)
    U = SharedTensor2d(new Tensor2d("Tau <IJ|BM>", naoccA, naoccA, navirA, nso_));
    U->sort(1243, T, 1.0, 0.0);
    T.reset();
    T = SharedTensor2d(new Tensor2d("New T2 <IJ|BA>", naoccA, naoccA, navirA, navirA));
    T->contract(false, false, naoccA*naoccA*navirA, navirA, nso_, U, CavirA, 1.0, 0.0);
    U.reset();
    Tnew = SharedTensor2d(new Tensor2d("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->sort(1423, T, 1.0, 1.0);
    T.reset();
    Tnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    // Read B(Q,ab)
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);

    // 2nd version
    // Read SO integrals
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mn)", nQ, nso_, nso_));
    bQso->read(psio_, PSIF_DFOCC_INTS, true, true);
    G = SharedTensor2d(new Tensor2d("DF_BASIS_CC J (mn|mn)", nso_, nso_));

    #pragma omp parallel for
    for(int m = 0 ; m < nso_; ++m){
        for(int n = 0 ; n < nso_; ++n){
            int mn = n + (m * nso_);
	    double sum = 0.0;
            for(int Q = 0 ; Q < nQ; ++Q){
                sum += bQso->get(Q,mn) * bQso->get(Q,mn);
            }
	    double value = sqrt(sum);
	    G->set(m, n, value);
        }
    }
    bQso.reset();

    // Transpose B
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mn)", nQ, ntri_so));
    bQso->read(psio_, PSIF_DFOCC_INTS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (m>=n|Q)", ntri_so, nQ));
    K->trans(bQso);
    bQso.reset();

    // Tau_ij^es = \sum(f) Tau_ij^ef C_sf
    Tau = SharedTensor2d(new Tensor2d("Tau (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_tau_amps(Tau,t2);
    U = SharedTensor2d(new Tensor2d("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA));
    U->sort(1324, Tau, 1.0, 0.0);
    Tau.reset();
    Tau = SharedTensor2d(new Tensor2d("Tau <IJ|ES>", naoccA, naoccA, navirA, nso_));
    Tau->contract(false, true, naoccA*naoccA*navirA, nso_, navirA, U, CavirA, 1.0, 0.0);
    U.reset();

    // Tau_ij^ls = \sum(e) Tau_ij^es C_le
    U = SharedTensor2d(new Tensor2d("Tau <IJ|SE>", naoccA, naoccA, nso_, navirA));
    U->sort(1243, Tau, 1.0, 0.0);
    Tau.reset();
    X = SharedTensor2d(new Tensor2d("Tau <IJ|SL>", naoccA, naoccA, nso_, nso_));
    X->contract(false, true, naoccA*naoccA*nso_, nso_, navirA, U, CavirA, 1.0, 0.0);
    U.reset();
    Tau = SharedTensor2d(new Tensor2d("Tau <IJ|LS>", naoccA, naoccA, nso_, nso_));
    Tau->sort(1243, X, 1.0, 0.0);
    X.reset();

    // tTau_ij^mn <= \sum_{ls} Tau_ij^ls <mn|ls>
    // tTau_ij^mn <= \sum_{ls} Tau_ij^ls (ml|ns)
    Tt = SharedTensor2d(new Tensor2d("Tau <IJ|MN>", naoccA, naoccA, nso_, nso_));
    SharedTensor2i nz_ints = SharedTensor2i(new Tensor2i("nz_ints", naoccA, naoccA));
    double Jmlns = 0.0;
    J = SharedTensor2d(new Tensor2d("J <ML|NS>", 1, 1));
    ndf_nz = 0;
    double tau_max = Tau->get_max_element();
    outfile->Printf("\tTau_max: %12.10f \n", tau_max);
    #pragma omp parallel for
    for(int m = 0 ; m < nso_; ++m){
        for(int l = 0 ; l <=m; ++l){
            int ml = index2(m,l);
            for(int n = 0 ; n < nso_; ++n){
                for(int s = 0 ; s <=n; ++s){
                    int ns = index2(n,s);

		    if (ml >= ns) {

		        // prescreening
			double Gmlns = tau_max * G->get(m,l) * G->get(n,s);
			//double Gmlns = Sso->get(m,l) * Sso->get(n,s) / tau_max;
			if (fabs(Gmlns) > int_cutoff_) {
			    ndf_nz++;
                            J->contract(false, true, 1, 1, nQ, K, K, ml*nQ, ns*nQ, 1.0, 0.0);
		            Jmlns = J->get(0,0);

			    // adressing
                            int ls = s + (l * nso_);
                            int ms = s + (m * nso_);
                            int mn = n + (m * nso_);
                            int ln = n + (l * nso_);
                            //int ls = s + (l * nso_);
                            //int ms = s + (m * nso_);
                            //int sl = l + (s * nso_);
                            //int sm = m + (s * nso_);
                            //int nl = l + (n * nso_);
                            //int nm = m + (n * nso_);

                            for(int i = 0 ; i < naoccA; ++i){
                                for(int j = 0 ; j < naoccA; ++j){
                                    int ij = j + (i * naoccA);

				    // tTau_ij^mn <= Tau_ij^ls * (ml|ns) (1)
				    Tt->add(ij, mn, Tau->get(ij,ls) * Jmlns);

				    // tTau_ij^ln <= Tau_ij^ms * (ml|ns) (2)
				    Tt->add(ij, ln, Tau->get(ij,ms) * Jmlns);

				    // we have 6 more terms!

			        }//j
	                    }// i

			}// end if cutoff


		    }// end if ml >= ns



		}// s
            }// n
        }// l
    }// m
    Tau.reset();
    K.reset();
    J.reset();
    outfile->Printf("\tNumber of prescreened integrals in Wabef    : %3d\n", ndf_nz);

    // Taut_ij^mb = \sum(n) Taut_ij^mn C(n,b)
    T = SharedTensor2d(new Tensor2d("Tau <IJ|MB>", naoccA, naoccA, nso_, navirA));
    T->contract(false, false, naoccA*naoccA*nso_, navirA, nso_, Tt, CavirA, 1.0, 0.0);
    Tt.reset();

    // T_ij^ab = \sum(m) Taut_ij^mb C(m,a)
    U = SharedTensor2d(new Tensor2d("Tau <IJ|BM>", naoccA, naoccA, navirA, nso_));
    U->sort(1243, T, 1.0, 0.0);
    T.reset();
    T = SharedTensor2d(new Tensor2d("New T2 <IJ|BA>", naoccA, naoccA, navirA, navirA));
    T->contract(false, false, naoccA*naoccA*navirA, navirA, nso_, U, CavirA, 1.0, 0.0);
    U.reset();
    Tnew = SharedTensor2d(new Tensor2d("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->sort(1423, T, 1.0, 1.0);
    T.reset();
    Tnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    // Read B(Q,ab)
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);

    timer_off("WabefT2");

}// end ccsd_WabefT2_basis
*/


}} // End Namespaces



