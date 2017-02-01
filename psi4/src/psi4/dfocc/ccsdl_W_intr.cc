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

void DFOCC::ccsdl_VmnijL2()
{
    // defs
    SharedTensor2d K, T, Lnew, U, Tau, W, X;
    SharedTensor2d M, L, I, Y, S, A;
    SharedTensor2d V, Vs, Va, Ts, Ta;

    timer_on("VmnijL2");

    // Read Vmnij
    V = SharedTensor2d(new Tensor2d("V <IJ|KL>", naoccA, naoccA, naoccA, naoccA));
    V->read(psio_, PSIF_DFOCC_AMPS);

    // Form (ma|nb)
    K = SharedTensor2d(new Tensor2d("Int (MA|NB)", naoccA, navirA, naoccA, navirA));
    K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);

    // Form <mn|ab>
    U = SharedTensor2d(new Tensor2d("Int <MN|AB>", naoccA, naoccA, navirA, navirA));
    U->sort(1324, K, 1.0, 0.0);
    K.reset();

    // l_ij^ab <= \sum_{m,n} I_mn^ab Vmnij
    // I_mn^ab = <mn|ab>
    // (+)I(ij, ab) = 1/2 (I_ij^ab + I_ji^ab) * (2 - \delta_{ij})
    // (-)I(ij, ab) = 1/2 (I_ij^ab - I_ji^ab) * (2 - \delta_{ij})
    Ts = SharedTensor2d(new Tensor2d("(+)tI [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    Ta = SharedTensor2d(new Tensor2d("(-)tI [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    Ts->symm_row_packed4(U);
    Ta->antisymm_row_packed4(U);
    U.reset();

    // Form (+/-)V(m>=n, i>=j)
    Vs = SharedTensor2d(new Tensor2d("(+)V [M>=N|I>=J]", ntri_ijAA, ntri_ijAA));
    Va = SharedTensor2d(new Tensor2d("(-)V [M>=N|I>=J]", ntri_ijAA, ntri_ijAA));
    Vs->symm4(V);
    Va->antisymm4(V);
    V.reset();

    // Symmetric & Anti-symmetric contributions
    S = SharedTensor2d(new Tensor2d("S (I>=J, A>=B)", ntri_ijAA, ntri_abAA));
    A = SharedTensor2d(new Tensor2d("A (I>=J, A>=B)", ntri_ijAA, ntri_abAA));
    S->gemm(true, false, Vs, Ts, 1.0, 0.0);
    A->gemm(true, false, Va, Ta, 1.0, 0.0);
    Ts.reset();
    Ta.reset();
    Vs.reset();
    Va.reset();

    // L(ia,jb) <-- S(a>=b,i>=j) + A(a>=b,i>=j)
    Lnew = SharedTensor2d(new Tensor2d("New L2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Lnew->read_symm(psio_, PSIF_DFOCC_AMPS);
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
                    Lnew->add(ia, jb, value);
                }
            }
        }
    }
    S.reset();
    A.reset();
    Lnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Lnew.reset();

    /*
    // Read Vmnij
    V = SharedTensor2d(new Tensor2d("V <IJ|KL>", naoccA, naoccA, naoccA, naoccA));
    V->read(psio_, PSIF_DFOCC_AMPS);

    // Form (ma|nb)
    K = SharedTensor2d(new Tensor2d("Int (MA|NB)", naoccA, navirA, naoccA, navirA));
    K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);

    // Form <mn|ab>
    M = SharedTensor2d(new Tensor2d("Int <MN|AB>", naoccA, naoccA, navirA, navirA));
    M->sort(1324, K, 1.0, 0.0);
    K.reset();

    // L_ij^ab
    L = SharedTensor2d(new Tensor2d("New L2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    L->gemm(true, false, V, M, 1.0, 0.0);
    V.reset();
    M.reset();

    // New L2
    Lnew = SharedTensor2d(new Tensor2d("New L2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Lnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Lnew->sort(1324, L, 1.0, 1.0);
    L.reset();
    Lnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Lnew.reset();
    */


    timer_off("VmnijL2");

}// end ccsdl_VmnijL2

//======================================================================
//    Wmnij
//======================================================================
void DFOCC::ccsdl_Wmnij()
{
    // defs
    SharedTensor2d K, T, Lnew, U, Tau, W, X;
    SharedTensor2d M, L, I, Y, S, A;
    SharedTensor2d V, Vs, Va, Ts, Ta;

    timer_on("Wmnij");

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
    //Tau->read_symm(psio_, PSIF_DFOCC_AMPS);
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
    W->write(psio_, PSIF_DFOCC_AMPS);
    W.reset();

    timer_off("Wmnij");

}// end ccsdl_Wmnij

//======================================================================
//    WijmnL2
//======================================================================
void DFOCC::ccsdl_WijmnL2()
{
    // defs
    SharedTensor2d K, T, Lnew, U, Tau, W, X;
    SharedTensor2d M, L, I, Y, S, A;
    SharedTensor2d V, Vs, Va, Ts, Ta;

    timer_on("WijmnL2");

    // W_mnij = <mn|ij>
    W = SharedTensor2d(new Tensor2d("W <MN|IJ>", naoccA, naoccA, naoccA, naoccA));
    W->read(psio_, PSIF_DFOCC_AMPS);
    /*
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
    Tau->read_symm(psio_, PSIF_DFOCC_AMPS);
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
    */

    // l_ij^ab <= \sum_{m,n} L_mn^ab W_ijmn
    // (+)L(ij, ab) = 1/2 (L_ij^ab + L_ji^ab) * (2 - \delta_{ij})
    // (-)L(ij, ab) = 1/2 (L_ij^ab - L_ji^ab) * (2 - \delta_{ij})
    U = SharedTensor2d(new Tensor2d("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    U->sort(1324, l2, 1.0, 0.0);
    Ts = SharedTensor2d(new Tensor2d("(+)tTau [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    Ta = SharedTensor2d(new Tensor2d("(-)tTau [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    Ts->symm_row_packed4(U);
    Ta->antisymm_row_packed4(U);
    U.reset();

    // Form (+/-)W(i>=j, m>=n)
    Vs = SharedTensor2d(new Tensor2d("(+)W [I>=J|M>=N]", ntri_ijAA, ntri_ijAA));
    Va = SharedTensor2d(new Tensor2d("(-)W [I>=J|M>=N]", ntri_ijAA, ntri_ijAA));
    Vs->symm4(W);
    Va->antisymm4(W);
    W.reset();

    // Symmetric & Anti-symmetric contributions
    S = SharedTensor2d(new Tensor2d("S (I>=J, A>=B)", ntri_ijAA, ntri_abAA));
    A = SharedTensor2d(new Tensor2d("A (I>=J, A>=B)", ntri_ijAA, ntri_abAA));
    S->gemm(false, false, Vs, Ts, 1.0, 0.0);
    A->gemm(false, false, Va, Ta, 1.0, 0.0);
    Ts.reset();
    Ta.reset();
    Vs.reset();
    Va.reset();

    // L(ia,jb) <-- S(a>=b,i>=j) + A(a>=b,i>=j)
    Lnew = SharedTensor2d(new Tensor2d("New L2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Lnew->read_symm(psio_, PSIF_DFOCC_AMPS);
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
                    Lnew->add(ia, jb, value);
                }
            }
        }
    }
    S.reset();
    A.reset();
    Lnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Lnew.reset();

    timer_off("WijmnL2");

}// end ccsdl_WijmnL2

//======================================================================
//    Wmbej
//======================================================================
void DFOCC::ccsdl_Wmbej()
{
    // defs
    SharedTensor2d K, L, T, T1, Tnew, U, Tau, W, X, Y, Z;

    timer_on("Wmbej");

    // Z_mbej = Z(me,jb)
    // Z(me,jb) <= (me|jb)
    Z = SharedTensor2d(new Tensor2d("Z (ME|JB)", naoccA, navirA, naoccA, navirA));
    Z->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);

    // Z(me,jb) <= \sum_{Q} T_jb^Q b_me^Q
    T = SharedTensor2d(new Tensor2d("T2 (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    Z->gemm(true, false, bQiaA, T, 1.0, 1.0);
    T.reset();

    // Z(me,jb) <= -\sum_{nf} t_jn^bf X(me,nf)
    // (mf|ne) = X(me,nf) (sort: 1432)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    X = SharedTensor2d(new Tensor2d("X (IA|JB)", naoccA, navirA, naoccA, navirA));
    X->sort(1432, K, 1.0, 0.0);
    K.reset();
    Z->gemm(false, false, X, t2, -1.0, 1.0);
    X.reset();
    Z->write(psio_, PSIF_DFOCC_AMPS);

    // Wl_mbej = Wl(me,jb)
    // Wl(me,jb) = Z(me,jb)
    W = SharedTensor2d(new Tensor2d("WL (ME|JB)", naoccA, navirA, naoccA, navirA));
    W->copy(Z);
    Z.reset();

    // Wl(me,jb) <= \sum_{Q} t_jb^Q' b_me^Q
    T1 = SharedTensor2d(new Tensor2d("T1p (Q|IA)", nQ, naoccA, navirA));
    T1->read(psio_, PSIF_DFOCC_AMPS);
    W->gemm(true, false, bQiaA, T1, 1.0, 1.0);
    T1.reset();
    W->write(psio_, PSIF_DFOCC_AMPS);
    //W->print();
    W.reset();

    timer_off("Wmbej");

}// end ccsdl_Wmbej

//======================================================================
//    Wmbje
//======================================================================
void DFOCC::ccsdl_Wmbje()
{
    // defs
    SharedTensor2d K, L, T, T1, Tnew, U, Tau, W, X, Y, Z;

    timer_on("Wmbje");

    // Z_mbje = Z'(me,jb)
    // Z'(me,jb) <= <me|jb>
    Z = SharedTensor2d(new Tensor2d("Zp (ME|JB)", naoccA, navirA, naoccA, navirA));
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|AB)", naoccA, naoccA, navirA, navirA));
    L->gemm(true, false, bQijA, bQabA, 1.0, 0.0);
    Z->sort(1324, L, 1.0, 0.0);
    L.reset();

    // Z'(me,jb) <= -\sum_{nf} t_nj^bf X(me,nf)
    // (mf|ne) = X(me,nf) (sort: 1432)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    X = SharedTensor2d(new Tensor2d("X (IA|JB)", naoccA, navirA, naoccA, navirA));
    X->sort(1432, K, 1.0, 0.0);
    K.reset();
    T = SharedTensor2d(new Tensor2d("T2p (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_t2_prime_amps(T,t2);
    Z->gemm(false, false, X, T, -1.0, 1.0);
    X.reset();
    T.reset();
    Z->write(psio_, PSIF_DFOCC_AMPS);

    // W_mbje = W'(me,jb)
    // W'(me,jb) <= Z'(me,jb)
    W = SharedTensor2d(new Tensor2d("WLp (ME|JB)", naoccA, navirA, naoccA, navirA));
    W->copy(Z);
    Z.reset();

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
    W->write(psio_, PSIF_DFOCC_AMPS);
    //W->print();
    W.reset();

    timer_off("Wmbje");

}// end ccsdl_Wmbje

//======================================================================
//    Wmnie
//======================================================================
void DFOCC::ccsdl_Wmnie()
{
    // defs
    SharedTensor2d T, W, X;

    timer_on("Wmnie");

    // X(im,ne) <= \sum_{Q} b_ne^Q ( t_im^Q + b_im^Q )
    T = SharedTensor2d(new Tensor2d("T1 (Q|IJ)", nQ, naoccA, naoccA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    T->add(bQijA);
    X = SharedTensor2d(new Tensor2d("X (IM|NE)", naoccA, naoccA, naoccA, navirA));
    X->gemm(true, false, T, bQiaA, 1.0, 0.0);
    T.reset();

    // W_mnie = W(mn,ie)
    W = SharedTensor2d(new Tensor2d("WL (MN|IE)", naoccA, naoccA, naoccA, navirA));
    // W(mn,ie) <= X(im,ne)
    W->sort(2314, X, 1.0, 0.0);
    X.reset();
    W->write(psio_, PSIF_DFOCC_AMPS);
    //W->print();
    W.reset();

    timer_off("Wmnie");

}// end ccsd_Wmnie

//======================================================================
//    Wmnie DIRECT
//======================================================================
void DFOCC::ccsdl_Wmnie_direct(SharedTensor2d &W)
{
    // defs
    SharedTensor2d T, X;

    timer_on("Wmnie");

    // X(im,ne) <= \sum_{Q} b_ne^Q ( t_im^Q + b_im^Q )
    T = SharedTensor2d(new Tensor2d("T1 (Q|IJ)", nQ, naoccA, naoccA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    T->add(bQijA);
    X = SharedTensor2d(new Tensor2d("X (IM|NE)", naoccA, naoccA, naoccA, navirA));
    X->gemm(true, false, T, bQiaA, 1.0, 0.0);
    T.reset();

    // W_mnie = W(mn,ie)
    // W(mn,ie) <= X(im,ne)
    W->sort(2314, X, 1.0, 0.0);
    X.reset();

    timer_off("Wmnie");

}// end ccsdl_Wmnie_direct

//======================================================================
//    Wmbij
//======================================================================
void DFOCC::ccsdl_Wmbij()
{
    // defs
    SharedTensor2d K, M, L, I, T, Tnew, U, Tau, W, WL, X, Y, Z, S, A;
    SharedTensor2d V, Vs, Ts, Va, Ta;

    timer_on("Wmbij");

    // WL_mbij = WL(mb,ij)
    // WL(mb,ij) = <mb|ij>
    WL = SharedTensor2d(new Tensor2d("WL (MB|IJ)", naoccA, navirA, naoccA, naoccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC (MI|JB)", naoccA, naoccA, naoccA, navirA));
    K->gemm(true, false, bQijA, bQiaA, 1.0, 0.0);
    WL->sort(1423, K, 1.0, 0.0);
    K.reset();

    // WL(mb,ij) = \sum_{e} t_ij^eb Ft_m = \sum_{e} Ft(m,e) T'(e,bij)
    //t2 = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    //t2->read_symm(psio_, PSIF_DFOCC_AMPS);
    T = SharedTensor2d(new Tensor2d("T2 <AB|IJ>)", navirA, navirA, naoccA, naoccA));
    T->sort(2413, t2, 1.0, 0.0);
    WL->contract(false, false, naoccA, navirA * naoccA * naoccA, navirA, FiaA, T, 1.0, 1.0);
    T.reset();

    // W_mnij = <mn|ij>
    W = SharedTensor2d(new Tensor2d("W <MN|IJ>", naoccA, naoccA, naoccA, naoccA));
    W->read(psio_, PSIF_DFOCC_AMPS);
    /*
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
    */

    // WL(mb,ij) = -\sum_{n} t_n^b Wmnij = -\sum_{n} W(mn,ij) T(n,b)
    WL->contract424(2, 1, W, t1A, -1.0, 1.0);
    W.reset();

    /*
    // WL(mb,ij) = \sum_{ef} tau_ij^ef <mb|ef>
    Tau = SharedTensor2d(new Tensor2d("Tau (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_tau_amps(Tau,t2);
    U = SharedTensor2d(new Tensor2d("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA));
    U->sort(1324, Tau, 1.0, 0.0);
    Tau.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ME|BF)", naoccA, navirA, navirA, navirA));
    K->gemm(true, false, bQiaA, bQabA, 1.0, 0.0);
    I = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <MB|EF>", naoccA, navirA, navirA, navirA));
    I->sort(1324, K, 1.0, 0.0);
    K.reset();
    WL->gemm(false, true, I, U, 1.0, 1.0);
    I.reset();
    */

    // WL(mb,ij) = \sum_{ef} tau_ij^ef <mb|ef>
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
    I = SharedTensor2d(new Tensor2d("I[M] <AF|E>", navirA * navirA, navirA));

    // Symmetric & Anti-symmetric contributions
    Vs = SharedTensor2d(new Tensor2d("(+)V[M] (A, E>=F)", navirA, ntri_abAA));
    Va = SharedTensor2d(new Tensor2d("(-)V[M] (A, E>=F)", navirA, ntri_abAA));
    S = SharedTensor2d(new Tensor2d("S[M] (A, I>=J)", navirA, ntri_ijAA));
    A = SharedTensor2d(new Tensor2d("A[M] (A, I>=J)", navirA, ntri_ijAA));
    // Main loop
    for(int m = 0 ; m < naoccA; ++m){

            // Form V[m](af,e) = \sum_{Q} b(Q,af) B(meQ)
            I->contract(true, true, navirA * navirA, navirA, nQ, bQabA, L, 0, m*navirA*nQ, 1.0, 0.0);

            // Form (+)V[m](a, e>=f)
            #pragma omp parallel for
            for(int a = 0 ; a < navirA; ++a){
                for(int e = 0 ; e < navirA; ++e){
                    int ae = ab_idxAA->get(a,e);
                    for(int f = 0 ; f <= e; ++f){
                        int af = ab_idxAA->get(a,f);
                        int ef = index2(e,f);
                        double value1 = 0.5 * ( I->get(af, e) + I->get(ae, f) );
                        double value2 = 0.5 * ( I->get(af, e) - I->get(ae, f) );
                        Vs->set(a, ef, value1);
                        Va->set(a, ef, value2);
                    }
                }
            }

            // Form S[m](a, i>=j) = \sum_{e>=f} Tau(i>=j,e>=f) V[m](a, e>=f)
            S->gemm(false, true, Vs, U, 1.0, 0.0);
            A->gemm(false, true, Va, T, 1.0, 0.0);

            // Form S(am,ij) & A(am,ij)-->W(maij)
            #pragma omp parallel for
            for(int a = 0 ; a < navirA; ++a){
                int ma = ia_idxAA->get(m,a);
                for(int i = 0 ; i < naoccA; ++i){
                    for(int j = 0 ; j < naoccA; ++j){
                        int ij = ij_idxAA->get(i,j);
                        int ij2 = index2(i,j);
                        int perm = ( i > j ) ? 1 : -1;
                        double value = S->get(a,ij2) + (perm * A->get(a,ij2));
                        WL->add(ma, ij, value);
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

    //outfile->Printf("\tI am here.\n");

    // WL(mb,ij) = \sum_{e} t_i^e Zmbej
    Z = SharedTensor2d(new Tensor2d("Z (ME|JB)", naoccA, navirA, naoccA, navirA));
    Z->read(psio_, PSIF_DFOCC_AMPS);
    X = SharedTensor2d(new Tensor2d("Z (MB|EJ)", naoccA, navirA, navirA, naoccA));
    X->sort(1423, Z, 1.0, 0.0);
    Z.reset();
    WL->contract424(3, 2, X, t1A, 1.0, 1.0);
    X.reset();

    // WL(mb,ij) = \sum_{e} t_j^e Zmbie = \sum_{m} Z'(me,ib) T(j,e)
    Z = SharedTensor2d(new Tensor2d("Zp (ME|JB)", naoccA, navirA, naoccA, navirA));
    Z->read(psio_, PSIF_DFOCC_AMPS);
    X = SharedTensor2d(new Tensor2d("Zp (MB|JE)", naoccA, navirA, naoccA, navirA));
    X->sort(1432, Z, 1.0, 0.0);
    Z.reset();
    //WL->contract424(4, 2, X, t1A, 1.0, 1.0);
    WL->contract(false, true, naoccA*navirA*naoccA, naoccA, navirA, X, t1A, 1.0, 1.0);
    X.reset();

    // WL(mb,ij) = 1/2\sum_{ne} u_jb^ne [ 2<mn|ie> - <nm|ie>]
    // WL(mb,ij) = 1/2\sum_{ne} U(jb,ne) [ 2(mi|ne) - (ni|me)]
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (MI|NE)", naoccA, naoccA, naoccA, navirA));
    K->gemm(true, false, bQijA, bQiaA, 1.0, 0.0);
    L = SharedTensor2d(new Tensor2d("2(MI|NE) - (NI|ME)", naoccA, naoccA, naoccA, navirA));
    L->tei_cs4_anti_symm(K,K);
    U = SharedTensor2d(new Tensor2d("U2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_u2_amps(U,t2);
    X = SharedTensor2d(new Tensor2d("X (MI|JB)", naoccA, naoccA, naoccA, navirA));
    X->gemm(false, false, L, U, 0.5, 1.0);
    L.reset();
    U.reset();
    WL->sort(1423, X, 1.0, 1.0);
    X.reset();

    // WL(mb,ij) = X_jmib + 1/2 X_imjb
    // X_jmib = -\sum_{ne} <jm|ne> T'(ne,ib)
    // <jm|ne> = (jn|me)
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <JM|NE>", naoccA, naoccA, naoccA, navirA));
    L->sort(1324, K, 1.0, 0.0);
    K.reset();
    X = SharedTensor2d(new Tensor2d("X (JM|IB)", naoccA, naoccA, naoccA, navirA));
    U = SharedTensor2d(new Tensor2d("Tp (IB|JA)", naoccA, navirA, naoccA, navirA));
    ccsd_t2_prime_amps(U,t2);
    X->gemm(false, false, L, U, -1.0, 0.0);
    L.reset();
    U.reset();
    WL->sort(2413, X, 0.5, 1.0);
    WL->sort(2431, X, 1.0, 1.0);
    X.reset();

    // Write and delete
    WL->write(psio_, PSIF_DFOCC_AMPS);
    //WL->print();
    WL.reset();

    timer_off("Wmbij");

}// end ccsd_Wmbij

//======================================================================
//    WmbejL2
//======================================================================
void DFOCC::ccsdl_WmbejL2()
{
    // defs
    SharedTensor2d K, L, T, Lnew, U, Tau, W, W2, X, Y;

    timer_on("WmbejL2");

    // W_mbje = W'(me,jb)
    W = SharedTensor2d(new Tensor2d("WLp (ME|JB)", naoccA, navirA, naoccA, navirA));
    W->read(psio_, PSIF_DFOCC_AMPS);

    // l_ij^ab <= 1/2*C(ia,jb) + 1/2*C(jb,ia) + C(ja,ib) + C(ib,ja)
    // l_ij^ab <= Ct(ia,jb) + 2*Ct(ib,ja)
    // C(ia,jb) = -\sum_{me} l_mi^ae W'(jb,me) = -\sum_{me} L'(ia,me) W'(jb,me)
    U = SharedTensor2d(new Tensor2d("L2p (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_t2_prime_amps(U,l2);
    Y = SharedTensor2d(new Tensor2d("C2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Y->gemm(false, true, U, W, -1.0, 0.0);
    U.reset();
    X = SharedTensor2d(new Tensor2d("C2+D2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    X->sort(1432, Y, 1.0, 0.0);
    X->axpy(Y, 0.5);
    Y.reset();

    // l_ij^ab <= D(ia,jb) + D(jb,ia)
    // D_ij^ab = 1/2 \sum_{me} Ut_im^ae [2*W(jb,me) - W'(jb,me)]
    Y = SharedTensor2d(new Tensor2d("2*W-W' (ME|JB)", naoccA, navirA, naoccA, navirA));
    Y->axpy(W, -1.0);
    W.reset();
    // W_mbej = W(me,jb)
    W = SharedTensor2d(new Tensor2d("WL (ME|JB)", naoccA, navirA, naoccA, navirA));
    W->read(psio_, PSIF_DFOCC_AMPS);
    Y->axpy(W, 2.0);
    W.reset();
    U = SharedTensor2d(new Tensor2d("Ut2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_u2_amps(U,l2);
    X->gemm(false, true, U, Y, 0.5, 1.0);
    U.reset();
    Y.reset();
    X->symmetrize();
    Lnew = SharedTensor2d(new Tensor2d("New L2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Lnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Lnew->axpy(X, 2.0);
    X.reset();
    Lnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Lnew.reset();

    timer_off("WmbejL2");

}// end ccsdl_WmbejL2

//======================================================================
//    LijmeL2: HIGH MEMORY
//======================================================================
void DFOCC::ccsdl_LijmeL2_high_mem()
{
    // defs
    SharedTensor2d K, M, L, I, T, Lnew, U, Tau, W, X, Y, S, A;
    SharedTensor2d V, Vs, Ts, Va, Ta, J;

    timer_on("LijmeL2");

    L = SharedTensor2d(new Tensor2d("L <IJ|KA>", naoccA, naoccA, naoccA, navirA));
    L->read(psio_, PSIF_DFOCC_AMPS);

    // Read in B(Q,a>=b)
    bQabA.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, ntri_abAA));
    K->read(psio_, PSIF_DFOCC_INTS);

    // J(ma,e>=b)
    J = SharedTensor2d(new Tensor2d("J (MF, A>=E)", naoccA*navirA, ntri_abAA));
    J->gemm(true, false, bQiaA, K, 1.0, 0.0);
    K.reset();

    // Build <me|ab>
    I = SharedTensor2d(new Tensor2d("I (ME,AB)", naoccA, navirA, navirA, navirA));
    #pragma omp parallel for
    for(int m = 0 ; m < naoccA; ++m){
         for(int e = 0 ; e < navirA; ++e){
             int me = ia_idxAA->get(m,e);
             for(int a = 0 ; a < navirA; ++a){
                 int ma = ia_idxAA->get(m,a);
                 for(int b = 0 ; b < navirA; ++b){
                     int eb = index2(e,b);
                     int ab = ab_idxAA->get(a,b);
		     I->set(me, ab, J->get(ma,eb));
		 }
	     }
	 }
    }
    J.reset();

    // Y(ij,ab) <= -\sum_{me} L(ij,me) <me|ab>
    // X(ia,jb) = Y(ij,ab)
    // l_ij^ab <= X(ia,jb) + X(jb,ia)
    Y = SharedTensor2d(new Tensor2d("Y <IJ|AB>", naoccA, naoccA, navirA, navirA));
    Y->gemm(false, false, L, I, -1.0, 0.0);
    I.reset();
    L.reset();
    X = SharedTensor2d(new Tensor2d("X (IA|JB)", naoccA, navirA, naoccA, navirA));
    X->sort(1324, Y, 1.0, 0.0);
    Y.reset();
    X->symmetrize();
    Lnew = SharedTensor2d(new Tensor2d("New L2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Lnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Lnew->axpy(X, 2.0);
    X.reset();
    Lnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Lnew.reset();

    // Read in B(Q,ab)
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);

    timer_off("LijmeL2");

}// end ccsd_LijmeL2_high_mem

//======================================================================
//    WabefL2
//======================================================================
void DFOCC::ccsdl_WabefL2()
{
    // defs
    SharedTensor2d K, M, L, I, I2, T, Lnew, U, Tau, W, X, Y, Z, S, A;
    SharedTensor2d V, Vs, Ts, Va, Ta, J, J2, T1;

    timer_on("WabefL2");

    // l_ij^ab <= \sum_{ef} l_ij^ef W_efab
    // (+)l(ij, ab) = 1/2 (l_ij^ab + l_ji^ab) * (2 - \delta_{ab})
    // (-)l(ij, ab) = 1/2 (l_ij^ab - l_ji^ab) * (2 - \delta_{ab})
    U = SharedTensor2d(new Tensor2d("(+)L [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    T = SharedTensor2d(new Tensor2d("(-)L [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
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
                    double value1 = 0.5 * perm * ( l2->get(ia,jb) + l2->get(ja,ib) );
                    double value2 = 0.5 * perm * ( l2->get(ia,jb) - l2->get(ja,ib) );
                    U->set(ij,ab,value1);
                    T->set(ij,ab,value2);
                }
            }
        }
    }

    // Read B(Q,ab) and B(Q,ia)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (AB|Q)", navirA * navirA, nQ));
    K = bQabA->transpose();
    // T1Q
    X = SharedTensor2d(new Tensor2d("T1 (Q|AB)", nQ, navirA, navirA));
    X->read(psio_, PSIF_DFOCC_AMPS);
    Y = SharedTensor2d(new Tensor2d("T1 (Q|BA)", nQ, navirA, navirA));
    Y->swap_3index_col(X);
    X.reset();
    Z = SharedTensor2d(new Tensor2d("T1 (BA|Q)", navirA * navirA, nQ));
    Z = Y->transpose();
    Y.reset();
    X = SharedTensor2d(new Tensor2d("B(AB|Q) - T1(BA|Q)", navirA * navirA, nQ));
    X->copy(Z);
    Z.reset();
    X->scale(-1.0);
    X->add(K);
    // B(aiQ)
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AI)", nQ, navirA, naoccA));
    M->swap_3index_col(bQiaA);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (AI|Q)", naoccA * navirA, nQ));
    L = M->transpose();
    M.reset();

    // malloc
    I = SharedTensor2d(new Tensor2d("I[A] <BF|E>", navirA * navirA, navirA));
    I2 = SharedTensor2d(new Tensor2d("I[A] <BE|F>", navirA * navirA, navirA));
    J = SharedTensor2d(new Tensor2d("J[A] <BM|E>", navirA * naoccA, navirA));
    J2 = SharedTensor2d(new Tensor2d("J[A] <BE|M>", navirA * navirA, naoccA));
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

        // Form J[a](bf,e) = \sum_{Q} B(bfQ)*[B(aeQ)-T(eaQ)] cost = V^4N/2
        I->contract(false, true, navirA*nb, navirA, nQ, K, X, 0, a*navirA*nQ, 1.0, 0.0);

        // Form J[a](bm,e) = \sum_{Q} B(bmQ)*B(aeQ) cost = OV^3N
        J->contract(false, true, nb*naoccA, navirA, nQ, L, K, 0, a*navirA*nQ, 1.0, 0.0);

        // J[a](be,m) = J[a](bm,e)
        J2->sort3b(132, navirA, naoccA, navirA, J, 1.0, 0.0);

        // J[a](bef) -= \sum_{m} J[a](be,m) * t(m,f)
        I2->contract(false, false, navirA*nb, navirA, naoccA, J2, t1A, -1.0, 0.0);

        // J[a](bf,e) += J[a](be,f)
        I->sort3b(132, navirA, navirA, navirA, I2, 1.0, 1.0);

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

        // Form L[a](b, i>=j) = \sum_{e>=f} L(i>=j,e>=f) V[a](b, e>=f)
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
    I2.reset();
    X.reset();
    Vs.reset();
    Va.reset();
    Ts.reset();
    Ta.reset();
    U.reset();
    T.reset();
    J.reset();
    J2.reset();
    L.reset();

    // L(ia,jb) <-- S(a>=b,i>=j) + A(a>=b,i>=j)
    Lnew = SharedTensor2d(new Tensor2d("New L2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Lnew->read_symm(psio_, PSIF_DFOCC_AMPS);
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
                    Lnew->add(ia, jb, value);
                }
            }
        }
    }
    S.reset();
    A.reset();
    Lnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Lnew.reset();

    /*
    // X(ae,bf) = W(ab,ef)
    X = SharedTensor2d(new Tensor2d("W (AE|BF)", navirA, navirA, navirA, navirA));
    // X(ae,bf) = -\sum(Q) B(Q,ae) * T(Q,bf)
    T = SharedTensor2d(new Tensor2d("T1 (Q|AB)", nQ, navirA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(true, false, bQabA, T, -1.0, 0.0);
    // X(ae,bf) += \sum(Q) [B(Q,ae)-T(Q,ae)] * B(Q,bf)
    U = SharedTensor2d(new Tensor2d("B-T1 (Q|AB)", nQ, navirA, navirA));
    U->copy(bQabA);
    U->axpy(T,-1.0);
    T.reset();
    X->gemm(true, false, U, bQabA, 1.0, 1.0);
    U.reset();
    // Wabef
    W = SharedTensor2d(new Tensor2d("W <AB|EF>", navirA, navirA, navirA, navirA));
    W->sort(1324, X, 1.0, 0.0);
    X.reset();

    // L_ij^ab
    U = SharedTensor2d(new Tensor2d("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    U->sort(1324, l2, 1.0, 0.0);
    L = SharedTensor2d(new Tensor2d("New L2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    L->gemm(false, false, U, W, 1.0, 0.0);
    W.reset();
    U.reset();

    // New L2
    Lnew = SharedTensor2d(new Tensor2d("New L2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Lnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Lnew->sort(1324, L, 1.0, 1.0);
    L.reset();
    Lnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Lnew.reset();
    */

    timer_off("WabefL2");

}// end ccsdl_WabefL2

//======================================================================
//    WabefL2: High memory version
//======================================================================
void DFOCC::ccsdl_WabefL2_high_mem()
{

    // defs
    SharedTensor2d K, M, L, I, T, Lnew, U, Tau, W, X, Y, S, A;
    SharedTensor2d V, Vs, Ts, Va, Ta, J, T1;

    timer_on("WabefL2");

    // l_ij^ab <= \sum_{ef} l_ij^ef W_efab
    // (+)l(ij, ab) = 1/2 (l_ij^ab + l_ji^ab) * (2 - \delta_{ab})
    // (-)l(ij, ab) = 1/2 (l_ij^ab - l_ji^ab) * (2 - \delta_{ab})
    U = SharedTensor2d(new Tensor2d("(+)L [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    T = SharedTensor2d(new Tensor2d("(-)L [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
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
                    double value1 = 0.5 * perm * ( l2->get(ia,jb) + l2->get(ja,ib) );
                    double value2 = 0.5 * perm * ( l2->get(ia,jb) - l2->get(ja,ib) );
                    U->set(ij,ab,value1);
                    T->set(ij,ab,value2);
                }
            }
        }
    }

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

    // L(ia,jb) <-- S(a>=b,i>=j) + A(a>=b,i>=j)
    Lnew = SharedTensor2d(new Tensor2d("New L2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Lnew->read_symm(psio_, PSIF_DFOCC_AMPS);
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
                    Lnew->add(ia, jb, value);
                }
            }
        }
    }
    S.reset();
    A.reset();
    Lnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Lnew.reset();

    // Read B(Q,ab)
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);

    timer_off("WabefL2");

}// end ccsdl_WabefL2_high_mem

}} // End Namespaces

