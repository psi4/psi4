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

void DFOCC::ccdl_VmnijL2() {
    // defs
    SharedTensor2d K, T, Lnew, U, Tau, W, X;
    SharedTensor2d M, L, I, Y, S, A;
    SharedTensor2d V, Vs, Va, Ts, Ta;

    timer_on("VmnijL2");

    // Read Vmnij
    V = std::make_shared<Tensor2d>("V <IJ|KL>", naoccA, naoccA, naoccA, naoccA);
    V->read(psio_, PSIF_DFOCC_AMPS);

    // Form (ma|nb)
    K = std::make_shared<Tensor2d>("Int (MA|NB)", naoccA, navirA, naoccA, navirA);
    K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);

    // Form <mn|ab>
    U = std::make_shared<Tensor2d>("Int <MN|AB>", naoccA, naoccA, navirA, navirA);
    U->sort(1324, K, 1.0, 0.0);
    K.reset();

    // l_ij^ab <= \sum_{m,n} I_mn^ab Vmnij
    // I_mn^ab = <mn|ab>
    // (+)I(ij, ab) = 1/2 (I_ij^ab + I_ji^ab) * (2 - \delta_{ij})
    // (-)I(ij, ab) = 1/2 (I_ij^ab - I_ji^ab) * (2 - \delta_{ij})
    Ts = std::make_shared<Tensor2d>("(+)tI [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
    Ta = std::make_shared<Tensor2d>("(-)tI [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
    Ts->symm_row_packed4(U);
    Ta->antisymm_row_packed4(U);
    U.reset();

    // Form (+/-)V(m>=n, i>=j)
    Vs = std::make_shared<Tensor2d>("(+)V [M>=N|I>=J]", ntri_ijAA, ntri_ijAA);
    Va = std::make_shared<Tensor2d>("(-)V [M>=N|I>=J]", ntri_ijAA, ntri_ijAA);
    Vs->symm4(V);
    Va->antisymm4(V);
    V.reset();

    // Symmetric & Anti-symmetric contributions
    S = std::make_shared<Tensor2d>("S (I>=J, A>=B)", ntri_ijAA, ntri_abAA);
    A = std::make_shared<Tensor2d>("A (I>=J, A>=B)", ntri_ijAA, ntri_abAA);
    S->gemm(true, false, Vs, Ts, 1.0, 0.0);
    A->gemm(true, false, Va, Ta, 1.0, 0.0);
    Ts.reset();
    Ta.reset();
    Vs.reset();
    Va.reset();

    // L(ia,jb) <-- S(a>=b,i>=j) + A(a>=b,i>=j)
    Lnew = std::make_shared<Tensor2d>("New L2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    Lnew->read_symm(psio_, PSIF_DFOCC_AMPS);
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
    V = std::make_shared<Tensor2d>("V <IJ|KL>", naoccA, naoccA, naoccA, naoccA);
    V->read(psio_, PSIF_DFOCC_AMPS);

    // Form (ma|nb)
    K = std::make_shared<Tensor2d>("Int (MA|NB)", naoccA, navirA, naoccA, navirA);
    K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);

    // Form <mn|ab>
    M = std::make_shared<Tensor2d>("Int <MN|AB>", naoccA, naoccA, navirA, navirA);
    M->sort(1324, K, 1.0, 0.0);
    K.reset();

    // L_ij^ab
    L = std::make_shared<Tensor2d>("New L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    L->gemm(true, false, V, M, 1.0, 0.0);
    V.reset();
    M.reset();

    // New L2
    Lnew = std::make_shared<Tensor2d>("New L2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    Lnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Lnew->sort(1324, L, 1.0, 1.0);
    L.reset();
    Lnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Lnew.reset();
    */

    timer_off("VmnijL2");

}  // end ccdl_VmnijL2

//======================================================================
//    Wmnij
//======================================================================
void DFOCC::ccdl_Wmnij() {
    // defs
    SharedTensor2d K, T, Lnew, U, Tau, W, X;
    SharedTensor2d M, L, I, Y, S, A;
    SharedTensor2d V, Vs, Va, Ts, Ta;

    timer_on("Wmnij");

    // W_mnij = <mn|ij>
    W = std::make_shared<Tensor2d>("W <MN|IJ>", naoccA, naoccA, naoccA, naoccA);
    K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IJ|KL)", naoccA, naoccA, naoccA, naoccA);
    K->gemm(true, false, bQijA, bQijA, 1.0, 0.0);
    W->sort(1324, K, 1.0, 0.0);
    K.reset();

    // W_mnij = \sum_{ef} Tau_ij^ef <mn|ef>
    // (+)Tau(ij, ab) = 1/2 (Tau_ij^ab + Tau_ji^ab) * (2 - \delta_{ab})
    // (-)Tau(ij, ab) = 1/2 (Tau_ij^ab - Tau_ji^ab) * (2 - \delta_{ab})
    // Tau = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    // Tau->read_symm(psio_, PSIF_DFOCC_AMPS);
    U = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    U->sort(1324, t2, 1.0, 0.0);
    // Tau.reset();
    Ts = std::make_shared<Tensor2d>("(+)tTau [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
    Ta = std::make_shared<Tensor2d>("(-)tTau [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
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
    W->write(psio_, PSIF_DFOCC_AMPS);
    W.reset();

    timer_off("Wmnij");

}  // end ccdl_Wmnij

//======================================================================
//    WijmnL2
//======================================================================
void DFOCC::ccdl_WijmnL2() {
    // defs
    SharedTensor2d K, T, Lnew, U, Tau, W, X;
    SharedTensor2d M, L, I, Y, S, A;
    SharedTensor2d V, Vs, Va, Ts, Ta;

    timer_on("WijmnL2");

    // W_mnij = <mn|ij>
    W = std::make_shared<Tensor2d>("W <MN|IJ>", naoccA, naoccA, naoccA, naoccA);
    W->read(psio_, PSIF_DFOCC_AMPS);
    /*
    K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IJ|KL)", naoccA, naoccA, naoccA, naoccA);
    K->gemm(true, false, bQijA, bQijA, 1.0, 0.0);
    W->sort(1324, K, 1.0, 0.0);
    K.reset();

    // W_mnij = \sum_{ef} Tau_ij^ef <mn|ef>
    // (+)Tau(ij, ab) = 1/2 (Tau_ij^ab + Tau_ji^ab) * (2 - \delta_{ab})
    // (-)Tau(ij, ab) = 1/2 (Tau_ij^ab - Tau_ji^ab) * (2 - \delta_{ab})
    Tau = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    Tau->read_symm(psio_, PSIF_DFOCC_AMPS);
    U = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    U->sort(1324, Tau, 1.0, 0.0);
    Tau.reset();
    Ts = std::make_shared<Tensor2d>("(+)tTau [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
    Ta = std::make_shared<Tensor2d>("(-)tTau [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
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
    U = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    U->sort(1324, l2, 1.0, 0.0);
    Ts = std::make_shared<Tensor2d>("(+)tTau [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
    Ta = std::make_shared<Tensor2d>("(-)tTau [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
    Ts->symm_row_packed4(U);
    Ta->antisymm_row_packed4(U);
    U.reset();

    // Form (+/-)W(i>=j, m>=n)
    Vs = std::make_shared<Tensor2d>("(+)W [I>=J|M>=N]", ntri_ijAA, ntri_ijAA);
    Va = std::make_shared<Tensor2d>("(-)W [I>=J|M>=N]", ntri_ijAA, ntri_ijAA);
    Vs->symm4(W);
    Va->antisymm4(W);
    W.reset();

    // Symmetric & Anti-symmetric contributions
    S = std::make_shared<Tensor2d>("S (I>=J, A>=B)", ntri_ijAA, ntri_abAA);
    A = std::make_shared<Tensor2d>("A (I>=J, A>=B)", ntri_ijAA, ntri_abAA);
    S->gemm(false, false, Vs, Ts, 1.0, 0.0);
    A->gemm(false, false, Va, Ta, 1.0, 0.0);
    Ts.reset();
    Ta.reset();
    Vs.reset();
    Va.reset();

    // L(ia,jb) <-- S(a>=b,i>=j) + A(a>=b,i>=j)
    Lnew = std::make_shared<Tensor2d>("New L2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    Lnew->read_symm(psio_, PSIF_DFOCC_AMPS);
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

}  // end ccdl_WijmnL2

//======================================================================
//    WmnijAA
//======================================================================
void DFOCC::ccdl_WmnijAA()
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
    W->write(psio_, PSIF_DFOCC_AMPS);
}// ccd_WmnijAA

//======================================================================
//    WmnijBB
//======================================================================
void DFOCC::ccdl_WmnijBB()
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
    W->write(psio_, PSIF_DFOCC_AMPS);
}// ccd_WmnijBB

//======================================================================
//    WmnijAB
//======================================================================
void DFOCC::ccdl_WmnijAB()
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
    W->write(psio_, PSIF_DFOCC_AMPS);

}// ccd_WmnijAB

//======================================================================
//    Wmbej
//======================================================================
void DFOCC::ccdl_Wmbej() {
    // defs
    SharedTensor2d K, L, T, T1, Tnew, U, Tau, W, X, Y, Z;

    timer_on("Wmbej");

    // Z_mbej = Z(me,jb)
    // Z(me,jb) <= (me|jb)
    // Wl_mbej = Wl(me,jb)
    // Wl(me,jb) = Z(me,jb)
    Z = std::make_shared<Tensor2d>("WL (ME|JB)", naoccA, navirA, naoccA, navirA);
    Z->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);

    // Z(me,jb) <= \sum_{Q} T_jb^Q b_me^Q
    T = std::make_shared<Tensor2d>("T2 (Q|IA)", nQ, naoccA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    Z->gemm(true, false, bQiaA, T, 1.0, 1.0);
    T.reset();

    // Z(me,jb) <= -\sum_{nf} t_jn^bf X(me,nf)
    // (mf|ne) = X(me,nf) (sort: 1432)
    K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA);
    K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    X = std::make_shared<Tensor2d>("X (IA|JB)", naoccA, navirA, naoccA, navirA);
    X->sort(1432, K, 1.0, 0.0);
    K.reset();
    Z->gemm(false, false, X, t2, -1.0, 1.0);
    X.reset();
    Z->write(psio_, PSIF_DFOCC_AMPS);
    Z.reset();

    timer_off("Wmbej");

}  // end ccdl_Wmbej

//======================================================================
//    Wmbje
//======================================================================
void DFOCC::ccdl_Wmbje() {
    // defs
    SharedTensor2d K, L, T, T1, Tnew, U, Tau, W, X, Y, Z;

    timer_on("Wmbje");

    // W_mbje = W'(me,jb)
    // W'(me,jb) = Z'(me,jb)
    // Z_mbje = Z'(me,jb)
    // Z'(me,jb) <= <me|jb>
    Z = std::make_shared<Tensor2d>("WLp (ME|JB)", naoccA, navirA, naoccA, navirA);
    L = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IJ|AB)", naoccA, naoccA, navirA, navirA);
    L->gemm(true, false, bQijA, bQabA, 1.0, 0.0);
    Z->sort(1324, L, 1.0, 0.0);
    L.reset();

    // Z'(me,jb) <= -\sum_{nf} t_nj^bf X(me,nf)
    // (mf|ne) = X(me,nf) (sort: 1432)
    K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA);
    K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    X = std::make_shared<Tensor2d>("X (IA|JB)", naoccA, navirA, naoccA, navirA);
    X->sort(1432, K, 1.0, 0.0);
    K.reset();
    T = std::make_shared<Tensor2d>("T2p (IA|JB)", naoccA, navirA, naoccA, navirA);
    ccsd_t2_prime_amps(T, t2);
    Z->gemm(false, false, X, T, -1.0, 1.0);
    X.reset();
    T.reset();
    Z->write(psio_, PSIF_DFOCC_AMPS);
    Z.reset();

    timer_off("Wmbje");

}  // end ccdl_Wmbje

//======================================================================
//    WmbejL2
//======================================================================
void DFOCC::ccdl_WmbejL2() {
    // defs
    SharedTensor2d K, L, T, Lnew, U, Tau, W, W2, X, Y;

    timer_on("WmbejL2");

    // W_mbje = W'(me,jb)
    W = std::make_shared<Tensor2d>("WLp (ME|JB)", naoccA, navirA, naoccA, navirA);
    W->read(psio_, PSIF_DFOCC_AMPS);

    // l_ij^ab <= 1/2*C(ia,jb) + 1/2*C(jb,ia) + C(ja,ib) + C(ib,ja)
    // l_ij^ab <= Ct(ia,jb) + 2*Ct(ib,ja)
    // C(ia,jb) = -\sum_{me} l_mi^ae W'(jb,me) = -\sum_{me} L'(ia,me) W'(jb,me)
    U = std::make_shared<Tensor2d>("L2p (IA|JB)", naoccA, navirA, naoccA, navirA);
    ccsd_t2_prime_amps(U, l2);
    Y = std::make_shared<Tensor2d>("C2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    Y->gemm(false, true, U, W, -1.0, 0.0);
    U.reset();
    X = std::make_shared<Tensor2d>("C2+D2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    X->sort(1432, Y, 1.0, 0.0);
    X->axpy(Y, 0.5);
    Y.reset();

    // l_ij^ab <= D(ia,jb) + D(jb,ia)
    // D_ij^ab = 1/2 \sum_{me} Ut_im^ae [2*W(jb,me) - W'(jb,me)]
    Y = std::make_shared<Tensor2d>("2*W-W' (ME|JB)", naoccA, navirA, naoccA, navirA);
    Y->axpy(W, -1.0);
    W.reset();
    // W_mbej = W(me,jb)
    W = std::make_shared<Tensor2d>("WL (ME|JB)", naoccA, navirA, naoccA, navirA);
    W->read(psio_, PSIF_DFOCC_AMPS);
    Y->axpy(W, 2.0);
    W.reset();
    U = std::make_shared<Tensor2d>("Ut2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    ccsd_u2_amps(U, l2);
    X->gemm(false, true, U, Y, 0.5, 1.0);
    U.reset();
    Y.reset();
    X->symmetrize();
    Lnew = std::make_shared<Tensor2d>("New L2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    Lnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Lnew->axpy(X, 2.0);
    X.reset();
    Lnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Lnew.reset();

    timer_off("WmbejL2");

}  // end ccdl_WmbejL2

//======================================================================
//    UHF Wmbej = Zmbej for CCD
//======================================================================
void DFOCC::ccdl_WMBEJ_AAAA()
{
    SharedTensor2d J, W, I, K, X, Y, T, Z, L;
    SharedTensor2d T2, Tau, T2new;
    // AAAA Block
    // Z_MBEJ =  <MB||EJ> - \sum_{NF} t_NJ^BF <MN||EF> + \sum_{nf} t_Jn^Bf <Mn|Ef>     (91)
    // Z_MBEJ =  <MB||EJ>
    // Z(ME,JB) = (ME|JB) - <ME|JB>
    Z = std::make_shared<Tensor2d>("WL (ME|JB)", naoccA, navirA, naoccA, navirA);
    Z->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    L = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IJ|AB)", naoccA, naoccA, navirA, navirA);
    L->gemm(true, false, bQijA, bQabA, 1.0, 0.0);
    Z->sort(1324, L, -1.0, 1.0);
    L.reset();
    // Z_MBEJ -= \sum_{NF} t_NJ^BF <MN||EF>
    // <MN||EF> = <MN|EF> - <MN|FE> = (ME|FN) - <ME|FN>
    // t <JN|BF> = t (NF|JB) (sort: 2413)
    J = std::make_shared<Tensor2d>("J (ME|NF)", naoccA, navirA, naoccA, navirA);
    J->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    K = std::make_shared<Tensor2d>("K (ME|NF)", naoccA, navirA, naoccA, navirA);
    K->sort(1432, J, -1.0, 0.0);
    K->axpy(J,1.0);
    J.reset();
    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T_ <NF|JB>", naoccA, navirA, naoccA, navirA);
    T->sort(1423, T2, 1.0, 0.0);
    T2.reset();
    Z->gemm(false, false, K, T, -1.0, 1.0);
    K.reset();
    T.reset();
    // Z_MBEJ += \sum_{nf} t_Jn^Bf <Mn|Ef>
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T <nf|JB>", naoccB, navirB, naoccA, navirA);
    T->sort(2413, T2, 1.0, 0.0);
    T2.reset();
    J = std::make_shared<Tensor2d>("J (ME|nf)", naoccA, navirA, naoccB, navirB);
    J->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    Z->gemm(false, false, J, T, 1.0, 1.0);
    Z->write(psio_, PSIF_DFOCC_AMPS);
} //ccsdl_ZMBEJ_AAAA()

void DFOCC::ccdl_Wmbej_BBBB()
{
    SharedTensor2d J, W, I, K, X, Y, T, Z, L;
    SharedTensor2d T2, Tau, T2new;
    // BBBB Block
    // Z_mbej =  <mb||ej> - \sum_{nf} t_nj^bf <mn||ef> + \sum_{NF} t_Nj^Fb <Nm|Fe>     (92)
    // Z_mbej =  <mb||ej>
    // Z(me,jb) = (me|jb) - <me|jb>
    Z = std::make_shared<Tensor2d>("WL (me|jb)", naoccB, navirB, naoccB, navirB);
    Z->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    L = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (ij|ab)", naoccB, naoccB, navirB, navirB);
    L->gemm(true, false, bQijB, bQabB, 1.0, 0.0);
    Z->sort(1324, L, -1.0, 1.0);
    L.reset();
    // Z_mbej -=  \sum_{nf} t_nj^bf <mn||ef>
    J = std::make_shared<Tensor2d>("J (me|nf)", naoccB, navirB, naoccB, navirB);
    J->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    K = std::make_shared<Tensor2d>("K (me|nf)", naoccB, navirB, naoccB, navirB);
    K->sort(1432, J, -1.0, 0.0);
    K->axpy(J,1.0);
    J.reset();
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T_ <nf|jb>", naoccB, navirB, naoccB, navirB);
    T->sort(1423, T2, 1.0, 0.0);
    T2.reset();
    Z->gemm(false, false, K, T, -1.0, 1.0);
    K.reset();
    T.reset();
    // Z_mbej +=  \sum_{NF} t_Nj^Fb <Nm|Fe>
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T <NF|jb>", naoccA, navirA, naoccB, navirB);
    T->sort(1324, T2, 1.0, 0.0);
    T2.reset();
    J = std::make_shared<Tensor2d>("J (me|NF)", naoccB, navirB, naoccA, navirA);
    J->gemm(true, false, bQiaB, bQiaA, 1.0, 0.0);
    Z->gemm(false, false, J, T, 1.0, 1.0);
    Z->write(psio_, PSIF_DFOCC_AMPS);
} // ccsdl_Zmbej_BBBB()

void DFOCC::ccdl_WMbEj_ABAB()
{
    SharedTensor2d J, W, I, K, X, Y, T, Z, L;
    SharedTensor2d T2, Tau, T2new;
    // ABAB Block
    // Z_MbEj =  <Mb|Ej> - \sum_{nf} t_nj^bf <Mn|Ef> + \sum_{NF} t_Nj^Fb <MN||EF>     (93)
    // Z_MbEj =  <Mb|Ej>
    // Z(ME,jb) = (ME|jb)
    Z = std::make_shared<Tensor2d>("WL (ME|jb)", naoccA, navirA, naoccB, navirB);
    Z->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    // Z_MbEj -= \sum_{nf} t_nj^bf <Mn|Ef>
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T_ <nf|jb>", naoccB, navirB, naoccB, navirB);
    T->sort(1423, T2, 1.0, 0.0);
    T2.reset();
    J = std::make_shared<Tensor2d>("J (ME|nf)", naoccA, navirA, naoccB, navirB);
    J->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    Z->gemm(false, false, J, T, -1.0, 1.0);
    J.reset();
    T.reset();
    // Z_MbEj += \sum_{NF} t_Nj^Fb <MN||EF>
    J = std::make_shared<Tensor2d>("J (ME|NF)", naoccA, navirA, naoccA, navirA);
    J->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    K = std::make_shared<Tensor2d>("K (ME|NF)", naoccA, navirA, naoccA, navirA);
    K->sort(1432, J, -1.0, 0.0);
    K->axpy(J,1.0);
    J.reset();
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T <NF|jb>", naoccA, navirA, naoccB, navirB);
    T->sort(1324, T2, 1.0, 0.0);
    T2.reset();
    Z->gemm(false, false, K, T, 1.0, 1.0);
    Z->write(psio_, PSIF_DFOCC_AMPS);
} // ccsdl_ZMbEj_ABAB()

void DFOCC::ccdl_WmBeJ_BABA()
{
    SharedTensor2d J, W, I, K, X, Y, T, Z, L;
    SharedTensor2d T2, Tau, T2new;
    // BABA Block
    // Z_mBeJ =  <Bm|Je> - \sum_{NF} t_NJ^BF <Nm|Fe> + \sum_{nf} t_Jn^Bf <mn||ef>      (94)
    // Z_mBeJ =  <Bm|Je>
    // <Bm|Je> = <mB|eJ> = (me|BJ) = (me|JB)
    Z = std::make_shared<Tensor2d>("WL (me|JB)", naoccB, navirB, naoccA, navirA);
    Z->gemm(true, false, bQiaB, bQiaA, 1.0, 0.0);
    // Z_mBeJ -= \sum_{NF} t_NJ^BF <Nm|Fe>
    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T_ <NF|JB>", naoccA, navirA, naoccA, navirA);
    T->sort(1423, T2, 1.0, 0.0);
    T2.reset();
    J = std::make_shared<Tensor2d>("J (me|NF)", naoccB, navirB, naoccA, navirA);
    J->gemm(true, false, bQiaB, bQiaA, 1.0, 0.0);//UB
    Z->gemm(false, false, J, T, -1.0, 1.0);
    J.reset();
    T.reset();
    // Z_mBeJ += \sum_{nf} t_Jn^Bf <mn||ef>
    J = std::make_shared<Tensor2d>("J (me|nf)", naoccB, navirB, naoccB, navirB);
    J->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    K = std::make_shared<Tensor2d>("K (me|nf)", naoccB, navirB, naoccB, navirB);
    K->sort(1432, J, -1.0, 0.0);
    K->axpy(J,1.0);
    J.reset();
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T <nf|JB>", naoccB, navirB, naoccA, navirA);
    T->sort(2413, T2, 1.0, 0.0);
    T2.reset();
    Z->gemm(false, false, K, T, 1.0, 1.0);
    Z->write(psio_, PSIF_DFOCC_AMPS);
} // ccsdl_ZmBeJ_BABA()

void DFOCC::ccdl_WMbeJ_ABBA()
{
    SharedTensor2d J, W, I, K, X, Y, T, Z, L;
    SharedTensor2d T2, Tau, T2new;
    // ABBA Block
    // Z_MbeJ =  - <Mb|Je> + \sum_{nF} t_Jn^Fb <Mn|Fe>     (95)
    // Z_MbeJ =  - <Mb|Je>
    // W(Me,Jb) = - <Me|Jb> = -(MJ|be) = Y (Me,Jb) (sort: 1423)
    J = std::make_shared<Tensor2d>("Int (MJ|be)", naoccA, naoccA, navirB, navirB);
    J->gemm(true, false, bQijA, bQabB, 1.0, 0.0);
    Z = std::make_shared<Tensor2d>("WL (Me|Jb)", naoccA, navirB, naoccA, navirB);
    Z->sort(1423, J, -1.0, 0.0);
    J.reset();
    // Z_MbeJ +=  \sum_{nF} t_Jn^Fb <Mn|Fe>
    J = std::make_shared<Tensor2d>("J (MF|ne)", naoccA, navirA, naoccB, navirB);
    J->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    X = std::make_shared<Tensor2d>("X (Me|nF)", naoccA, navirB, naoccB, navirA);
    X->sort(1432, J, 1.0, 0.0);
    J.reset();
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T <nF|Jb>", naoccB, navirA, naoccA, navirB);
    T->sort(2314, T2, 1.0, 0.0);
    T2.reset();
    Z->gemm(false, false, X, T, 1.0, 1.0);
    Z->write(psio_, PSIF_DFOCC_AMPS);
}  // ccsdl_ZMbeJ_ABBA()

void DFOCC::ccdl_WmBEj_BAAB()
{
    SharedTensor2d J, W, I, K, X, Y, T, Z, L;
    SharedTensor2d T2, Tau, T2new;
    // BAAB Block
    // Z_mBEj =  - <Bm|Ej> - \sum_{Nf} t_Nj^Bf <Nm|Ef>     (96)
    // Z_mBEj =  - <Bm|Ej>
    // W(mE,jB) = - <mE|jB> = -(EB|mj)
    J = std::make_shared<Tensor2d>("Int (EB|mj)", navirA, navirA, naoccB, naoccB);
    J->gemm(true, false, bQabA, bQijB, 1.0, 0.0);
    Z = std::make_shared<Tensor2d>("WL (mE|jB)", naoccB, navirA, naoccB, navirA);
    Z->sort(3142, J, -1.0, 0.0);
    J.reset();
    // Z_mBEj -= \sum_{Nf} t_Nj^Bf <Nm|Ef>
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T <Nf|jB>", naoccA, navirB, naoccB, navirA);
    T->sort(1423, T2, 1.0, 0.0);
    T2.reset();
    J = std::make_shared<Tensor2d>("J (NE|mf)", naoccA, navirA, naoccB, navirB);
    J->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    X = std::make_shared<Tensor2d>("X (mE|Nf)", naoccB, navirA, naoccA, navirB);
    X->sort(3214, J, -1.0, 0.0);
    J.reset();
    Z->gemm(false, false, X, T, -1.0, 1.0);
    Z->write(psio_, PSIF_DFOCC_AMPS);
}// ccsdl_ZmBEj_BAAB()

//======================================================================
//    WabefL2
//======================================================================
void DFOCC::ccdl_WabefL2() {
    // defs
    SharedTensor2d K, M, L, I, T, Lnew, U, Tau, W, X, Y, Z, S, A;
    SharedTensor2d V, Vs, Ts, Va, Ta, J, T1;

    timer_on("WabefL2");

    // l_ij^ab <= \sum_{ef} l_ij^ef W_efab
    // (+)l(ij, ab) = 1/2 (l_ij^ab + l_ji^ab) * (2 - \delta_{ab})
    // (-)l(ij, ab) = 1/2 (l_ij^ab - l_ji^ab) * (2 - \delta_{ab})
    U = std::make_shared<Tensor2d>("(+)L [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
    T = std::make_shared<Tensor2d>("(-)L [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
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
                    double value1 = 0.5 * perm * (l2->get(ia, jb) + l2->get(ja, ib));
                    double value2 = 0.5 * perm * (l2->get(ia, jb) - l2->get(ja, ib));
                    U->set(ij, ab, value1);
                    T->set(ij, ab, value2);
                }
            }
        }
    }

    // Read B(Q,ab) and B(Q,ia)
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (AB|Q)", navirA * navirA, nQ);
    K = bQabA->transpose();
    // B(iaQ)
    M = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AI)", nQ, navirA, naoccA);
    M->swap_3index_col(bQiaA);
    L = std::make_shared<Tensor2d>("DF_BASIS_CC B (AI|Q)", naoccA * navirA, nQ);
    L = M->transpose();
    M.reset();

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

        // Form J[a](bf,e) = \sum_{Q} B(bfQ)*B(aeQ) cost = V^4N/2
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

        // Form L[a](b, i>=j) = \sum_{e>=f} L(i>=j,e>=f) V[a](b, e>=f)
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
    L.reset();

    // L(ia,jb) <-- S(a>=b,i>=j) + A(a>=b,i>=j)
    Lnew = std::make_shared<Tensor2d>("New L2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    Lnew->read_symm(psio_, PSIF_DFOCC_AMPS);
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
                    Lnew->add(ia, jb, value);
                }
            }
        }
    }
    S.reset();
    A.reset();
    Lnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Lnew.reset();

    timer_off("WabefL2");

}  // end ccdl_WabefL2

}  // namespace dfoccwave
}  // namespace psi
