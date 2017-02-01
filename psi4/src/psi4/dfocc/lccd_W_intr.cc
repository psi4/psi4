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

void DFOCC::lccd_WmnijT2()
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

    // t_ij^ab <= \sum_{m,n} T_mn^ab Wmnij
    // (+)T(ij, ab) = 1/2 (T_ij^ab + T_ji^ab) * (2 - \delta_{ij})
    // (-)T(ij, ab) = 1/2 (T_ij^ab - T_ji^ab) * (2 - \delta_{ij})
    U = SharedTensor2d(new Tensor2d("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    U->sort(1324, t2, 1.0, 0.0);
    Ts = SharedTensor2d(new Tensor2d("(+)tT [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    Ta = SharedTensor2d(new Tensor2d("(-)tT [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
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

}// end lccd_WmnijT2

//======================================================================
//    WmnijT2AA
//======================================================================
void DFOCC::lccd_WmnijT2AA()
{
    // defs
    SharedTensor2d K, T, Tnew, U, Tau, W, X;
    SharedTensor2d M, L, I, Y, S, A;
    SharedTensor2d V, Vs, Va, Ts, Ta;

    timer_on("WmnijT2");

    // W_mnij = <mn||ij>
    W = SharedTensor2d(new Tensor2d("W <MN|IJ>", naoccA, naoccA, naoccA, naoccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|KL)", naoccA, naoccA, naoccA, naoccA));
    K->gemm(true, false, bQijA, bQijA, 1.0, 0.0);
    W->sort(1324, K, 1.0, 0.0);
    W->sort(1342, K, -1.0, 1.0);
    K.reset();

    // t_ij^ab <= 1/2 \sum_{m,n} T_mn^ab Wmnij
    // (-)T(ij, ab) = 1/2 (T_ij^ab - T_ji^ab) * (2 - \delta_{ij})
    U = SharedTensor2d(new Tensor2d("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    U->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Ta = SharedTensor2d(new Tensor2d("(-)tT [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    Ta->antisymm_row_packed4(U);
    U.reset();

    // Form (-)W(m>=n, i>=j)
    Va = SharedTensor2d(new Tensor2d("(-)W [M>=N|I>=J]", ntri_ijAA, ntri_ijAA));
    Va->antisymm4(W);
    W.reset();

    // Anti-symmetric contributions
    A = SharedTensor2d(new Tensor2d("A (I>=J, A>=B)", ntri_ijAA, ntri_abAA));
    A->gemm(true, false, Va, Ta, 0.5, 0.0);
    Ta.reset();
    Va.reset();

    // T(ia,jb) <-- A(a>=b,i>=j)
    Tnew = SharedTensor2d(new Tensor2d("New T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    Tnew->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    #pragma omp parallel for
    for(int a = 0 ; a < navirA; ++a){
        for(int b = 0 ; b < navirA; ++b){
            int ab = index2(a,b);
            int ab2 = ab_idxAA->get(a,b);
            for(int i = 0 ; i < naoccA; ++i){
                for(int j = 0 ; j < naoccA; ++j){
                    int ij2 = ij_idxAA->get(i,j);
                    int ij = index2(i,j);
                    int perm1 = ( i > j ) ? 1 : -1;
                    int perm2 = ( a > b ) ? 1 : -1;
                    double value = perm1 * perm2 * A->get(ij,ab);
                    Tnew->add(ij2, ab2, value);
                }
            }
        }
    }
    A.reset();
    Tnew->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    timer_off("WmnijT2");

}// end lccd_WmnijT2AA

//======================================================================
//    WmnijT2BB
//======================================================================
void DFOCC::lccd_WmnijT2BB()
{
    // defs
    SharedTensor2d K, T, Tnew, U, Tau, W, X;
    SharedTensor2d M, L, I, Y, S, A;
    SharedTensor2d V, Vs, Va, Ts, Ta;

    timer_on("WmnijT2");

    // W_mnij = <mn||ij>
    W = SharedTensor2d(new Tensor2d("W <mn|ij>", naoccB, naoccB, naoccB, naoccB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ij|kl)", naoccB, naoccB, naoccB, naoccB));
    K->gemm(true, false, bQijB, bQijB, 1.0, 0.0);
    W->sort(1324, K, 1.0, 0.0);
    W->sort(1342, K, -1.0, 1.0);
    K.reset();

    // t_ij^ab <= 1/2 \sum_{m,n} T_mn^ab Wmnij
    // (-)T(ij, ab) = 1/2 (T_ij^ab - T_ji^ab) * (2 - \delta_{ij})
    U = SharedTensor2d(new Tensor2d("T2 <ij|ab>", naoccB, naoccB, navirB, navirB));
    U->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Ta = SharedTensor2d(new Tensor2d("(-)tT [I>=J|A>=B]", ntri_ijBB, ntri_abBB));
    Ta->antisymm_row_packed4(U);
    U.reset();

    // Form (-)W(m>=n, i>=j)
    Va = SharedTensor2d(new Tensor2d("(-)W [M>=N|I>=J]", ntri_ijBB, ntri_ijBB));
    Va->antisymm4(W);
    W.reset();

    // Anti-symmetric contributions
    A = SharedTensor2d(new Tensor2d("A (I>=J, A>=B)", ntri_ijBB, ntri_abBB));
    A->gemm(true, false, Va, Ta, 1.0, 0.0);
    Ta.reset();
    Va.reset();

    // T(ia,jb) <-- A(a>=b,i>=j)
    Tnew = SharedTensor2d(new Tensor2d("New T2 <ij|ab>", naoccB, naoccB, navirB, navirB));
    Tnew->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    #pragma omp parallel for
    for(int a = 0 ; a < navirB; ++a){
        for(int b = 0 ; b < navirB; ++b){
            int ab = index2(a,b);
            int ab2 = ab_idxBB->get(a,b);
            for(int i = 0 ; i < naoccB; ++i){
                for(int j = 0 ; j < naoccB; ++j){
                    int ij2 = ij_idxBB->get(i,j);
                    int ij = index2(i,j);
                    int perm1 = ( i > j ) ? 1 : -1;
                    int perm2 = ( a > b ) ? 1 : -1;
                    double value = perm1 * perm2 * A->get(ij,ab);
                    Tnew->add(ij2, ab2, 0.5*value);
                }
            }
        }
    }
    A.reset();
    Tnew->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    timer_off("WmnijT2");

}// end lccd_WmnijT2BB

//======================================================================
//    WmnijT2AB
//======================================================================
void DFOCC::lccd_WmnijT2AB()
{
    // defs
    SharedTensor2d K, T, Tnew, U, Tau, W, X;

    timer_on("WmnijT2");

    // W_mnij = <mn|ij>
    W = SharedTensor2d(new Tensor2d("W <Mn|Ij>", naoccA, naoccB, naoccA, naoccB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|kl)", naoccA, naoccA, naoccB, naoccB));
    K->gemm(true, false, bQijA, bQijB, 1.0, 0.0);
    W->sort(1324, K, 1.0, 0.0);
    K.reset();

    // t_Ij^Ab <= \sum_{M,n} T_Mn^Ab W_MnIj
    T = SharedTensor2d(new Tensor2d("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    T->read(psio_, PSIF_DFOCC_AMPS);
    Tnew = SharedTensor2d(new Tensor2d("New T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    Tnew->read(psio_, PSIF_DFOCC_AMPS);
    Tnew->gemm(true, false, W, T, 1.0, 1.0);
    T.reset();
    W.reset();
    Tnew->write(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    timer_off("WmnijT2");

}// end lccd_WmnijT2AB

//======================================================================
//    WmbejT2
//======================================================================
void DFOCC::lccd_WmbejT2()
{
    // defs
    SharedTensor2d K, L, T, T1, Tnew, U, Tau, W, W2, X, Y;

    timer_on("WmbejT2");

    // W_mbej = W(me,jb)
    // W(me,jb) <= (me|jb)
    W = SharedTensor2d(new Tensor2d("W (ME|JB)", naoccA, navirA, naoccA, navirA));
    W->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    W->write(psio_, PSIF_DFOCC_AMPS);
    W.reset();

    // W_mbje = W'(me,jb)
    // W'(me,jb) <= <me|jb>
    W = SharedTensor2d(new Tensor2d("Wp (ME|JB)", naoccA, navirA, naoccA, navirA));
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|AB)", naoccA, naoccA, navirA, navirA));
    L->gemm(true, false, bQijA, bQabA, 1.0, 0.0);
    W->sort(1324, L, 1.0, 0.0);
    L.reset();

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

}// end lccd_WmbejT2

//======================================================================
//    WmbejT2AA
//======================================================================
void DFOCC::lccd_WmbejT2AA()
{
    // defs
    SharedTensor2d K, L, T, T1, Tnew, U, Tau, W, W2, X, Y;

    timer_on("WmbejT2");

    // W_MBEJ = W(ME,JB)
    // W(ME,JB) = (ME|JB) - <ME|JB>
    W = SharedTensor2d(new Tensor2d("W (ME|JB)", naoccA, navirA, naoccA, navirA));
    W->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|AB)", naoccA, naoccA, navirA, navirA));
    L->gemm(true, false, bQijA, bQabA, 1.0, 0.0);
    W->sort(1324, L, -1.0, 1.0);
    L.reset();

    // X(IA,JB) = \sum(ME) T(IA,ME) W(ME,JB)
    // T(IJ,AB) <= P_(IJ) * P_(AB) X(IA,JB)
    T = SharedTensor2d(new Tensor2d("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    T->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    U = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    U->sort(1324, T, 1.0, 0.0);
    T.reset();
    X = SharedTensor2d(new Tensor2d("X (IA|JB)", naoccA, navirA, naoccA, navirA));
    X->gemm(false, false, U, W, 1.0, 0.0);
    U.reset();
    W.reset();
    Tnew = SharedTensor2d(new Tensor2d("New T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    Tnew->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->P_ijab(X);
    X.reset();
    Tnew->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    // W_MbEJ = W(ME,jb)
    // W(ME,jb) = (ME|jb)
    W = SharedTensor2d(new Tensor2d("W (ME|jb)", naoccA, navirA, naoccB, navirB));
    W->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);

    // X(IA,JB) = \sum(me) T(IA,me) W(JB,me)
    // T(IJ,AB) <= P_(IJ) * P_(AB) X(IA,JB)
    T = SharedTensor2d(new Tensor2d("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    T->read(psio_, PSIF_DFOCC_AMPS);
    U = SharedTensor2d(new Tensor2d("T2 (IA|jb)", naoccA, navirA, naoccB, navirB));
    U->sort(1324, T, 1.0, 0.0);
    T.reset();
    X = SharedTensor2d(new Tensor2d("X (IA|JB)", naoccA, navirA, naoccA, navirA));
    X->gemm(false, true, U, W, 1.0, 0.0);
    U.reset();
    W.reset();
    Tnew = SharedTensor2d(new Tensor2d("New T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    Tnew->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->P_ijab(X);
    X.reset();
    Tnew->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    timer_off("WmbejT2");

}// end lccd_WmbejT2AA

//======================================================================
//    WmbejT2BB
//======================================================================
void DFOCC::lccd_WmbejT2BB()
{
    // defs
    SharedTensor2d K, L, T, T1, Tnew, U, Tau, W, W2, X, Y;

    timer_on("WmbejT2");

    // W_mbej = W(me,jb)
    // W(me,jb) = (me|jb) - <me|jb>
    W = SharedTensor2d(new Tensor2d("W (me|jb)", naoccB, navirB, naoccB, navirB));
    W->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ij|ab)", naoccB, naoccB, navirB, navirB));
    L->gemm(true, false, bQijB, bQabB, 1.0, 0.0);
    W->sort(1324, L, -1.0, 1.0);
    L.reset();

    // X(ia,jb) = \sum(me) T(ia,me) W(me,jb)
    // T(ij,ab) <= P_(ij) * P_(ab) X(ia,jb)
    T = SharedTensor2d(new Tensor2d("T2 <ij|ab>", naoccB, naoccB, navirB, navirB));
    T->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    U = SharedTensor2d(new Tensor2d("T2 (ia|jb)", naoccB, navirB, naoccB, navirB));
    U->sort(1324, T, 1.0, 0.0);
    T.reset();
    X = SharedTensor2d(new Tensor2d("X (ia|jb)", naoccB, navirB, naoccB, navirB));
    X->gemm(false, false, U, W, 1.0, 0.0);
    U.reset();
    W.reset();
    Tnew = SharedTensor2d(new Tensor2d("New T2 <ij|ab>", naoccB, naoccB, navirB, navirB));
    Tnew->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->P_ijab(X);
    X.reset();
    Tnew->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    // W_MbEJ = W(ME,jb)
    // W(ME,jb) = (ME|jb)
    W = SharedTensor2d(new Tensor2d("W (ME|jb)", naoccA, navirA, naoccB, navirB));
    W->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);

    // X(ia,jb) = \sum(ME) T(ME,ia) W(ME,jb)
    // T(ij,ab) <= P_(ij) * P_(ab) X(ia,jb)
    T = SharedTensor2d(new Tensor2d("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    T->read(psio_, PSIF_DFOCC_AMPS);
    U = SharedTensor2d(new Tensor2d("T2 (IA|jb)", naoccA, navirA, naoccB, navirB));
    U->sort(1324, T, 1.0, 0.0);
    T.reset();
    X = SharedTensor2d(new Tensor2d("X (ia|jb)", naoccB, navirB, naoccB, navirB));
    X->gemm(true, false, U, W, 1.0, 0.0);
    U.reset();
    W.reset();
    Tnew = SharedTensor2d(new Tensor2d("New T2 <ij|ab>", naoccB, naoccB, navirB, navirB));
    Tnew->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->P_ijab(X);
    X.reset();
    Tnew->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    timer_off("WmbejT2");

}// end lccd_WmbejT2BB

//======================================================================
//    WmbejT2AB
//======================================================================
void DFOCC::lccd_WmbejT2AB()
{
    // defs
    SharedTensor2d K, L, T, T1, Tnew, U, Tau, W, W2, X, Y;

    timer_on("WmbejT2");

    // Read
    Tnew = SharedTensor2d(new Tensor2d("New T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    Tnew->read(psio_, PSIF_DFOCC_AMPS);

    // W_MBEJ = W(ME,JB)
    // W(ME,JB) = (ME|JB) - <ME|JB>
    W = SharedTensor2d(new Tensor2d("W (ME|JB)", naoccA, navirA, naoccA, navirA));
    W->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|AB)", naoccA, naoccA, navirA, navirA));
    L->gemm(true, false, bQijA, bQabA, 1.0, 0.0);
    W->sort(1324, L, -1.0, 1.0);
    L.reset();

    // X(IA,jb) = \sum(ME) T(ME,jb) W(ME,IA)
    // T(Ij,Ab) <= X(IA,jb)
    T = SharedTensor2d(new Tensor2d("T2 (IA|jb)", naoccA, navirA, naoccB, navirB));
    T->read(psio_, PSIF_DFOCC_AMPS);
    X = SharedTensor2d(new Tensor2d("X (IA|jb)", naoccA, navirA, naoccB, navirB));
    X->gemm(true, false, W, T, 1.0, 0.0);
    T.reset();
    W.reset();
    Tnew->sort(1324, X, 1.0, 1.0);
    X.reset();


    // W_mbej = W(me,jb)
    // W(me,jb) = (me|jb) - <me|jb>
    W = SharedTensor2d(new Tensor2d("W (me|jb)", naoccB, navirB, naoccB, navirB));
    W->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ij|ab)", naoccB, naoccB, navirB, navirB));
    L->gemm(true, false, bQijB, bQabB, 1.0, 0.0);
    W->sort(1324, L, -1.0, 1.0);
    L.reset();

    // X(IA,jb) = \sum(me) T(IA,me) W(me,jb)
    // T(Ij,Ab) <= X(IA,jb)
    T = SharedTensor2d(new Tensor2d("T2 (IA|jb)", naoccA, navirA, naoccB, navirB));
    T->read(psio_, PSIF_DFOCC_AMPS);
    X = SharedTensor2d(new Tensor2d("X (IA|jb)", naoccA, navirA, naoccB, navirB));
    X->gemm(false, false, T, W, 1.0, 0.0);
    T.reset();
    W.reset();
    Tnew->sort(1324, X, 1.0, 1.0);
    X.reset();

    // W_MbEJ = W(ME,jb)
    // W(ME,jb) = (ME|jb)
    W = SharedTensor2d(new Tensor2d("W (ME|jb)", naoccA, navirA, naoccB, navirB));
    W->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);

    // X(IA,jb) = \sum(ME) T(IA,ME) W(ME,jb)
    // T(Ij,Ab) <= X(IA,jb)
    T = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    T->read_symm(psio_, PSIF_DFOCC_AMPS);
    X = SharedTensor2d(new Tensor2d("X (IA|jb)", naoccA, navirA, naoccB, navirB));
    X->gemm(false, false, T, W, 1.0, 0.0);
    T.reset();
    Tnew->sort(1324, X, 1.0, 1.0);
    X.reset();

    // X(IA,jb) = \sum(me) T(me,jb) W(IA,me)
    // T(Ij,Ab) <= X(IA,jb)
    T = SharedTensor2d(new Tensor2d("T2 (ia|jb)", naoccB, navirB, naoccB, navirB));
    T->read_symm(psio_, PSIF_DFOCC_AMPS);
    X = SharedTensor2d(new Tensor2d("X (IA|jb)", naoccA, navirA, naoccB, navirB));
    X->gemm(false, false, W, T, 1.0, 0.0);
    T.reset();
    W.reset();
    Tnew->sort(1324, X, 1.0, 1.0);
    X.reset();

    // W_mBEj = W(mE,jB)
    // W(mE,jB) = - <mE|jB> = -(EB|mj)
    W = SharedTensor2d(new Tensor2d("W (mE|jB)", naoccB, navirA, naoccB, navirA));
    K = SharedTensor2d(new Tensor2d("Int (EB|mj)", navirA, navirA, naoccB, naoccB));
    K->gemm(true, false, bQabA, bQijB, 1.0, 0.0);
    W->sort(3142, K, -1.0, 0.0);
    K.reset();

    // X(Ib,jA) = \sum(mE) T(Im,Eb) W(mE,jA) = \sum(me) T'(Ib,mE) W(mE,jA)
    // T(Ij,Ab) <= X(Ib,jA)
    U = SharedTensor2d(new Tensor2d("T2 (IA|jb)", naoccA, navirA, naoccB, navirB));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T = SharedTensor2d(new Tensor2d("T2 (Ib|jA)", naoccA, navirB, naoccB, navirA));
    T->sort(1432, U, 1.0, 0.0);
    U.reset();
    X = SharedTensor2d(new Tensor2d("X (Ib|jA)", naoccA, navirB, naoccB, navirA));
    X->gemm(false, false, T, W, 1.0, 0.0);
    T.reset();
    W.reset();
    Tnew->sort(1342, X, 1.0, 1.0);
    X.reset();

    // W_MbeJ = W(Me,Jb)
    // W(Me,Jb) = - <Me|Jb> = -(MJ|eb)
    W = SharedTensor2d(new Tensor2d("W (Me|Jb)", naoccA, navirB, naoccA, navirB));
    K = SharedTensor2d(new Tensor2d("Int (MJ|eb)", naoccA, naoccA, navirB, navirB));
    K->gemm(true, false, bQijA, bQabB, 1.0, 0.0);
    W->sort(1324, K, -1.0, 0.0);
    K.reset();

    // X(Ib,jA) = \sum(Me) T(Mj,Ae) W(Me,Ib) = \sum(me) T'(Me,jA) W(Me,Ib)
    // T(Ij,Ab) <= X(Ib,jA)
    U = SharedTensor2d(new Tensor2d("T2 (IA|jb)", naoccA, navirA, naoccB, navirB));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T = SharedTensor2d(new Tensor2d("T2 (Ib|jA)", naoccA, navirB, naoccB, navirA));
    T->sort(1432, U, 1.0, 0.0);
    U.reset();
    X = SharedTensor2d(new Tensor2d("X (Ib|jA)", naoccA, navirB, naoccB, navirA));
    X->gemm(true, false, W, T, 1.0, 0.0);
    T.reset();
    W.reset();
    Tnew->sort(1342, X, 1.0, 1.0);
    X.reset();

    // Write
    Tnew->write(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    timer_off("WmbejT2");

}// end lccd_WmbejT2AB

//======================================================================
//    WabefT2
//======================================================================
void DFOCC::lccd_WabefT2()
{
    // defs
    SharedTensor2d K, M, L, I, T, Tnew, U, Tau, W, X, Y, S, A;
    SharedTensor2d V, Vs, Ts, Va, Ta;

    timer_on("WabefT2");

    // t_ij^ab <= \sum_{ef} T_ij^ef <ab|ef>
    // (+)T(ij, ab) = 1/2 (T_ij^ab + T_ji^ab) * (2 - \delta_{ab})
    // (-)T(ij, ab) = 1/2 (T_ij^ab - T_ji^ab) * (2 - \delta_{ab})
    U = SharedTensor2d(new Tensor2d("(+)T [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    T = SharedTensor2d(new Tensor2d("(-)T [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
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
                    double value1 = 0.5 * perm * ( t2->get(ia,jb) + t2->get(ja,ib) );
                    double value2 = 0.5 * perm * ( t2->get(ia,jb) - t2->get(ja,ib) );
                    U->set(ij,ab,value1);
                    T->set(ij,ab,value2);
                }
            }
        }
    }

    // Read B(Q,ab)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (AB|Q)", navirA * navirA, nQ));
    K = bQabA->transpose();

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

    timer_off("WabefT2");

}// end lccd_WabefT2

//======================================================================
//    WabefT2AA
//======================================================================
void DFOCC::lccd_WabefT2AA()
{
    // defs
    SharedTensor2d K, M, L, I, T, Tnew, U, Tau, W, X, Y, S, A;
    SharedTensor2d V, Vs, Ts, Va, Ta;

    timer_on("WabefT2");

    // t_ij^ab <= 1/2 \sum_{ef} T_ij^ef <ab||ef> = \sum_{ef} T_ij^ef <ab|ef>
    // (-)T(ij, ab) = 1/2 (T_ij^ab - T_ji^ab) * (2 - \delta_{ab})
    X = SharedTensor2d(new Tensor2d("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    X->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T = SharedTensor2d(new Tensor2d("(-)T [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    T->antisymm_col_packed4(X);
    X.reset();

    // Read B(Q,ab)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (AB|Q)", navirA * navirA, nQ));
    K = bQabA->transpose();

    // malloc
    I = SharedTensor2d(new Tensor2d("I[A] <BF|E>", navirA * navirA, navirA));
    Va = SharedTensor2d(new Tensor2d("(-)V[A] (B, E>=F)", navirA, ntri_abAA));
    Ta = SharedTensor2d(new Tensor2d("(-)T[B] (B, I>=J)", navirA, ntri_ijAA));

    // Anti-symmetric contributions
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
                        double value2 = 0.5 * ( I->get(bf, e) - I->get(be, f) );
                        Va->set(b, ef, value2);
                    }
                }
            }

            // Form T[a](b, i>=j) = 1/2\sum_{e>=f} Tau(i>=j,e>=f) V[a](b, e>=f)
            Ta->contract(false, true, nb, ntri_ijAA, ntri_abAA, Va, T, 1.0, 0.0);

            // Form A(ij,ab)
            #pragma omp parallel for
            for(int b = 0 ; b <=a; ++b){
                int ab = index2(a,b);
                for(int i = 0 ; i < naoccA; ++i){
                    for(int j = 0 ; j <= i; ++j){
                        int ij = index2(i,j);
                        A->add(ab, ij, Ta->get(b,ij));
                    }
                }
            }

    }
    K.reset();
    I.reset();
    Va.reset();
    Ta.reset();
    T.reset();

    // T(ia,jb) <-- A(a>=b,i>=j)
    Tnew = SharedTensor2d(new Tensor2d("New T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    Tnew->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    #pragma omp parallel for
    for(int a = 0 ; a < navirA; ++a){
        for(int b = 0 ; b < navirA; ++b){
            int ab = index2(a,b);
            int ab2 = ab_idxAA->get(a,b);
            for(int i = 0 ; i < naoccA; ++i){
                for(int j = 0 ; j < naoccA; ++j){
                    int ij2 = ij_idxAA->get(i,j);
                    int ij = index2(i,j);
                    int perm1 = ( i > j ) ? 1 : -1;
                    int perm2 = ( a > b ) ? 1 : -1;
                    double value = perm1 * perm2 * A->get(ab,ij);
                    Tnew->add(ij2, ab2, value);
                }
            }
        }
    }
    A.reset();
    Tnew->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    timer_off("WabefT2");

}// end lccd_WabefT2AA

//======================================================================
//    WabefT2BB
//======================================================================
void DFOCC::lccd_WabefT2BB()
{
    // defs
    SharedTensor2d K, M, L, I, T, Tnew, U, Tau, W, X, Y, S, A;
    SharedTensor2d V, Vs, Ts, Va, Ta;

    timer_on("WabefT2");

    // t_ij^ab <= 1/2 \sum_{ef} T_ij^ef <ab||ef> = \sum_{ef} T_ij^ef <ab|ef>
    // (-)T(ij, ab) = 1/2 (T_ij^ab - T_ji^ab) * (2 - \delta_{ab})
    X = SharedTensor2d(new Tensor2d("T2 <ij|ab>", naoccB, naoccB, navirB, navirB));
    X->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T = SharedTensor2d(new Tensor2d("(-)T [I>=J|A>=B]", ntri_ijBB, ntri_abBB));
    T->antisymm_col_packed4(X);
    X.reset();

    // Read B(Q,ab)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (ab|Q)", navirB * navirB, nQ));
    K = bQabB->transpose();

    // malloc
    I = SharedTensor2d(new Tensor2d("I[A] <BF|E>", navirB * navirB, navirB));
    Va = SharedTensor2d(new Tensor2d("(-)V[A] (B, E>=F)", navirB, ntri_abBB));
    Ta = SharedTensor2d(new Tensor2d("(-)T[B] (B, I>=J)", navirB, ntri_ijBB));

    // Anti-symmetric contributions
    A = SharedTensor2d(new Tensor2d("A (A>=B, I>=J)", ntri_abBB, ntri_ijBB));
    // Main loop
    for(int a = 0 ; a < navirB; ++a){
            int nb = a+1;

            // Form V[a](bf,e) = \sum_{Q} B(bfQ)*B(aeQ) cost = V^4N/2
            I->contract(false, true, navirB*nb, navirB, nQ, K, K, 0, a*navirB*nQ, 1.0, 0.0);

            // Form (+)V[a](b, e>=f)
            #pragma omp parallel for
            for(int b = 0 ; b <= a; ++b){
                for(int e = 0 ; e < navirB; ++e){
                    int be = e + (b * navirB);
                    for(int f = 0 ; f <= e; ++f){
                        int ef = index2(e,f);
                        int bf = f + (b * navirB);
                        double value2 = 0.5 * ( I->get(bf, e) - I->get(be, f) );
                        Va->set(b, ef, value2);
                    }
                }
            }

            // Form T[a](b, i>=j) = \sum_{e>=f} Tau(i>=j,e>=f) V[a](b, e>=f)
            Ta->contract(false, true, nb, ntri_ijBB, ntri_abBB, Va, T, 1.0, 0.0);

            // Form A(ij,ab)
            #pragma omp parallel for
            for(int b = 0 ; b <=a; ++b){
                int ab = index2(a,b);
                for(int i = 0 ; i < naoccB; ++i){
                    for(int j = 0 ; j <= i; ++j){
                        int ij = index2(i,j);
                        A->add(ab, ij, Ta->get(b,ij));
                    }
                }
            }

    }
    K.reset();
    I.reset();
    Va.reset();
    Ta.reset();
    T.reset();

    // T(ia,jb) <-- A(a>=b,i>=j)
    Tnew = SharedTensor2d(new Tensor2d("New T2 <ij|ab>", naoccB, naoccB, navirB, navirB));
    Tnew->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    #pragma omp parallel for
    for(int a = 0 ; a < navirB; ++a){
        for(int b = 0 ; b < navirB; ++b){
            int ab = index2(a,b);
            int ab2 = ab_idxBB->get(a,b);
            for(int i = 0 ; i < naoccB; ++i){
                for(int j = 0 ; j < naoccB; ++j){
                    int ij2 = ij_idxBB->get(i,j);
                    int ij = index2(i,j);
                    int perm1 = ( i > j ) ? 1 : -1;
                    int perm2 = ( a > b ) ? 1 : -1;
                    double value = perm1 * perm2 * A->get(ab,ij);
                    Tnew->add(ij2, ab2, value);
                }
            }
        }
    }
    A.reset();
    Tnew->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    timer_off("WabefT2");

}// end lccd_WabefT2BB

//======================================================================
//    WabefT2AB
//======================================================================
void DFOCC::lccd_WabefT2AB()
{
    // defs
    SharedTensor2d K, M, L, I, J, T, Tnew, U, Tau, W, X, Y, S, A;

    timer_on("WabefT2");

    // Read
    Tnew = SharedTensor2d(new Tensor2d("New T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    Tnew->read(psio_, PSIF_DFOCC_AMPS);

    // t_Ij^Ab <= \sum_{Ef} T_Ij^Ef <Ab|Ef>
    T = SharedTensor2d(new Tensor2d("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    T->read(psio_, PSIF_DFOCC_AMPS);

    // malloc
    J = SharedTensor2d(new Tensor2d("J[A] <E|bf>", navirA, navirB*navirB));
    I = SharedTensor2d(new Tensor2d("I[A] <b|Ef>", navirB, navirA*navirB));
    X = SharedTensor2d(new Tensor2d("T[A] <b|Ij>", navirB, naoccA*naoccB));
    K = SharedTensor2d(new Tensor2d("B[A] <E|Q>", navirA, nQ));

    // Main loop
    for(int a = 0 ; a < navirA; ++a){

	    // Form B[A](e,Q)
            #pragma omp parallel for
            for(int Q = 0 ; Q < nQ; ++Q){
                for(int e = 0 ; e < navirA; ++e){
                    int ae = ab_idxAA->get(a,e);
                    K->set(e, Q, bQabA->get(Q,ae));
                }
            }

            // Form J[A](E,bf) = \sum_{Q} B[A](e,Q) * B(Q,bf)
            J->gemm(false, false, K, bQabB, 1.0, 0.0);

            // Form I[A](b,Ef)
            #pragma omp parallel for
            for(int b = 0 ; b < navirB; ++b){
                for(int e = 0 ; e < navirA; ++e){
                    for(int f = 0 ; f < navirB; ++f){
                        int bf = f + (b * navirB);
                        int ef = ab_idxAB->get(e,f);
                        I->set(b, ef, J->get(e,bf));
                    }
                }
            }

            // Form T[A](b,Ij) = \sum_{Ef} I[A](b, Ef) T(Ij,Ef)
            X->gemm(false, true, I, T, 1.0, 0.0);

            // T[A](b,Ij) --> T(Ij,Ab)
            #pragma omp parallel for
            for(int b = 0 ; b < navirB; ++b){
                int ab = ab_idxAB->get(a,b);
                for(int i = 0 ; i < naoccA; ++i){
                    for(int j = 0 ; j < naoccB; ++j){
                        int ij = ij_idxAB->get(i,j);
                        Tnew->add(ij, ab, X->get(b,ij));
                    }
                }
            }

    }
    K.reset();
    J.reset();
    I.reset();
    T.reset();
    X.reset();

    // Write
    Tnew->write(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    timer_off("WabefT2");

}// end lccd_WabefT2AB


}} // End Namespaces


