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

/** Standard library includes */
#include "psi4/libqt/qt.h"
#include "defines.h"
#include "dfocc.h"

using namespace psi;

namespace psi {
namespace dfoccwave {

void DFOCC::ccd_tpdm() {
    timer_on("tpdm");
    SharedTensor2d T, U, Tau, G, G2, V, X, Y, Z;
    SharedTensor2d Ts, Ta, Vs, Va, S, A;

    // if (reference_ == "RESTRICTED") {

    //============================
    // OO-Block Correlation TPDM
    //============================

    // G_ij^Q += P+(ij) V_ij^Q
    G = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|IJ)", nQ, naoccA, naoccA);
    V = std::make_shared<Tensor2d>("V (Q|IJ)", nQ, naoccA, naoccA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, 1.0);
    V.reset();
    // G_ij^Q -= P+(ij) 2*V'_ij^Q
    V = std::make_shared<Tensor2d>("Vp (Q|IJ)", nQ, naoccA, naoccA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, -2.0);
    V.reset();

    // symmetrize
    G->symmetrize3(G);
    G->scale(2.0);
    G2 = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|OO)", nQ, noccA, noccA);
    G2->set3_act_oo(nfrzc, G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS);
    if (print_ > 3) G2->print();
    G2.reset();

    //============================
    // OV-Block Correlation TPDM
    //============================

    // G_ia^Q = T_ia^Q
    G = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|IA)", nQ, naoccA, navirA);
    T = std::make_shared<Tensor2d>("T2 (Q|IA)", nQ, naoccA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, 1.0);
    T.reset();
    // G_ia^Q += L_ia^Q
    T = std::make_shared<Tensor2d>("L2 (Q|IA)", nQ, naoccA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, 1.0);
    T.reset();
    // G_ia^Q += 2*y_ia^Q
    T = std::make_shared<Tensor2d>("Y (Q|IA)", nQ, naoccA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, 2.0);
    T.reset();

    // G_ia^Q -= \sum(m) T_ma^Q G_im
    T = std::make_shared<Tensor2d>("T2 (Q|IA)", nQ, naoccA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->contract233(false, false, naoccA, navirA, GijA, T, -1.0, 1.0);

    // G_ia^Q += \sum(e) T_ie^Q G_ea
    G->contract(false, false, nQ * naoccA, navirA, navirA, T, GabA, 1.0, 1.0);
    // G->cont323("IA", "IE", "EA", false, T, GabA, 1.0, 1.0); // it works
    T.reset();

    // G_ia^Q += \sum(me) U(ia,me) (G_em^Q - G_me^Q)
    U = std::make_shared<Tensor2d>("G (Q|AI)", nQ, navirA, naoccA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("Temp (Q|IA)", nQ, naoccA, navirA);
    T->swap_3index_col(U);
    U.reset();
    U = std::make_shared<Tensor2d>("G (Q|IA)", nQ, naoccA, navirA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, -1.0);
    U.reset();
    U = std::make_shared<Tensor2d>("U2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    U->read_symm(psio_, PSIF_DFOCC_AMPS);
    G->gemm(false, true, T, U, 1.0, 1.0);
    U.reset();
    T.reset();

    // Form overall OV Block
    G2 = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|OV)", nQ, noccA, nvirA);
    G2->set3_act_ov(nfrzc, naoccA, navirA, nvirA, G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS);
    if (print_ > 3) G2->print();

    // Form G_vo^Q
    G = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|VO)", nQ, nvirA, noccA);
    G->swap_3index_col(G2);
    G2.reset();
    G->write(psio_, PSIF_DFOCC_DENS);
    if (print_ > 3) G->print();
    G.reset();

    //============================
    // VV-Block Correlation TPDM
    //============================

    // G_ab^Q -= P+(ab) 2*V_ab^Q
    G = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|AB)", nQ, navirA, navirA);
    V = std::make_shared<Tensor2d>("V (Q|AB)", nQ, navirA, navirA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, -2.0);
    V.reset();

    // G_ab^Q += P+(ab) \sum(ef) Vt_aebf b_ef^Q
    // Read b_ef^Q
    bQabA = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA);
    bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);

    // (+)Ut(mn, ae) = 1/2 (Ut_mn^ae + Ut_nm^ae) * (2 - \delta_{mn})
    // (-)Ut(mn, ae) = 1/2 (Ut_mn^ae - Ut_nm^ae) * (2 - \delta_{mn})
    T = std::make_shared<Tensor2d>("Ut2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    T->read_symm(psio_, PSIF_DFOCC_AMPS);
    U = std::make_shared<Tensor2d>("Ut2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    U->sort(1324, T, 1.0, 0.0);
    T.reset();
    Ts = std::make_shared<Tensor2d>("(+)Ut [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
    Ta = std::make_shared<Tensor2d>("(-)Ut [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
    Ts->symm_row_packed4(U);
    Ta->antisymm_row_packed4(U);
    U.reset();

    // Symmetric & Anti-symmetric contributions
    U = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    U->read_symm(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    Tau->sort(1324, U, 1.0, 0.0);
    U.reset();
    Vs = std::make_shared<Tensor2d>("(+)T[B] (F,M>=N)", navirA, ntri_ijAA);
    Va = std::make_shared<Tensor2d>("(-)T[B] (F,M>=N)", navirA, ntri_ijAA);
    S = std::make_shared<Tensor2d>("S[B] (F,A>=E)", navirA, ntri_abAA);
    A = std::make_shared<Tensor2d>("A[B] (F,A>=E)", navirA, ntri_abAA);
    V = std::make_shared<Tensor2d>("V[B] (A,EF)", navirA, navirA, navirA);
    X = std::make_shared<Tensor2d>("TPDM[B] (Q|A)", nQ, navirA);
    // Main loop
    for (int b = 0; b < navirA; ++b) {
// Form (+)Tau[b](f, m>=n)
#pragma omp parallel for
        for (int m = 0; m < naoccA; ++m) {
            for (int n = 0; n <= m; ++n) {
                int mn2 = index2(m, n);
                int mn = n + (m * naoccA);
                int nm = m + (n * naoccA);
                for (int f = 0; f < navirA; ++f) {
                    int bf = f + (b * navirA);
                    double value1 = 0.5 * (Tau->get(mn, bf) + Tau->get(nm, bf));
                    double value2 = 0.5 * (Tau->get(mn, bf) - Tau->get(nm, bf));
                    Vs->set(f, mn2, value1);
                    Va->set(f, mn2, value2);
                }
            }
        }

        // Form S[b](f,a>=e) = \sum(m>=n) (+)Ut(m>=n,a>=e) * Tau[b](f,m>=n)
        S->contract(false, false, navirA, ntri_abAA, ntri_ijAA, Vs, Ts, 1.0, 0.0);
        A->contract(false, false, navirA, ntri_abAA, ntri_ijAA, Va, Ta, 1.0, 0.0);

// Form V[b](a,ef)
#pragma omp parallel for
        for (int a = 0; a < navirA; ++a) {
            for (int e = 0; e < navirA; ++e) {
                int ae = index2(a, e);
                for (int f = 0; f < navirA; ++f) {
                    int ef = ab_idxAA->get(e, f);
                    int perm1 = (a > e) ? 1 : -1;
                    double value = S->get(f, ae) + (perm1 * A->get(f, ae));
                    V->set(a, ef, value);
                }
            }
        }

        // G2[b](Q,a) = \sum(ef) b(Q,ef) * V[b](a,ef)
        X->gemm(false, true, bQabA, V, 1.0, 0.0);

// Form G2
#pragma omp parallel for
        for (int Q = 0; Q < nQ; ++Q) {
            for (int a = 0; a < navirA; ++a) {
                int ab = ab_idxAA->get(a, b);
                G->add(Q, ab, X->get(Q, a));
            }
        }
    }
    Tau.reset();
    Ts.reset();
    Ta.reset();
    Vs.reset();
    Va.reset();
    S.reset();
    A.reset();
    X.reset();
    V.reset();
    bQabA.reset();

    // symmetrize
    G->symmetrize3(G);
    G->scale(2.0);
    G2 = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|VV)", nQ, nvirA, nvirA);
    G2->set3_act_vv(G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS, true, true);
    if (print_ > 3) G2->print();
    G2.reset();

    //}// end if (reference_ == "RESTRICTED")
    // else if (reference_ == "UNRESTRICTED") {
    //}// else if (reference_ == "UNRESTRICTED")
    timer_off("tpdm");
}  // end ccd_tpdm

//======================================================================
//    PPL Alpha
//======================================================================
void DFOCC::ccd_tpdm_pplA(SharedTensor2d& G, std::string amps) {

    SharedTensor2d T, T2, L2, Tau, X, Y, Z, V, U, L, G2;
    SharedTensor2d T1, Ts, Ta, Vs, Va, S, A;

        // G_AB^Q += P+(ab) 1/2\sum(EF) V_AEBF b_EF^Q
        // V_AEBF = 1/2 \sum(MN) L_MN^AE T_MN^BF
        // Read b_EF^Q
        bQabA = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA);
        bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);

        // (-)Tt(mn, ae) = 1/2 (L_mn^ae - L_nm^ae) * (2 - \delta_{mn})
        T1 = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        T1->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        Ta = std::make_shared<Tensor2d>("(-)Ut [I>=J|A>=B]", ntri_ijAA, ntri_abAA);
        Ta->antisymm_row_packed4(T1);

        // Anti-symmetric contributions
        if (amps == "T2") {
            T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
            T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        }
        else if (amps == "Tau") {
            Tau = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
            Tau->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
            T2 = std::make_shared<Tensor2d>("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA);
            uccsd_tau_amps(naoccA, naoccA, navirA, navirA, T2, Tau, t1A, t1A);
            Tau.reset();
        }
        else {
             std::cout << "ccd_tpdm_pplA ---> Unrecognized amps, it should be T2 or Tau ! \n";
             exit(1);
        }
        Va = std::make_shared<Tensor2d>("(-)T[B] (F,M>=N)", navirA, ntri_ijAA);
        A = std::make_shared<Tensor2d>("A[B] (F,A>=E)", navirA, ntri_abAA);
        V = std::make_shared<Tensor2d>("V[B] (A,EF)", navirA, navirA, navirA);
        X = std::make_shared<Tensor2d>("TPDM[B] (Q|A)", nQ, navirA);
        // Main loop
        for (int b = 0; b < navirA; ++b) {
// Form (+)T[b](f, m>=n)
#pragma omp parallel for
            for (int m = 0; m < naoccA; ++m) {
                for (int n = 0; n <= m; ++n) {
                    int mn2 = index2(m, n);
                    int mn = n + (m * naoccA);
                    int nm = m + (n * naoccA);
                    for (int f = 0; f < navirA; ++f) {
                        int bf = f + (b * navirA);
                        double value2 = 0.5 * (T2->get(mn, bf) - T2->get(nm, bf));
                        Va->set(f, mn2, value2);
                    }
                }
            }

            // Form A[b](f,a>=e) = 1/2 \sum(m>=n) (+)Tt(m>=n,a>=e) * T[b](f,m>=n)
            A->contract(false, false, navirA, ntri_abAA, ntri_ijAA, Va, Ta, 0.5, 0.0);

// Form V[b](a,ef)
#pragma omp parallel for
            for (int a = 0; a < navirA; ++a) {
                for (int e = 0; e < navirA; ++e) {
                    int ae = index2(a, e);
                    for (int f = 0; f < navirA; ++f) {
                        int ef = ab_idxAA->get(e, f);
                        int perm1 = (a > e) ? 1 : -1;
                        double value = perm1 * A->get(f, ae);
                        V->set(a, ef, value);
                    }
                }
            }

            // G2[b](Q,a) = 1/2\sum(ef) b(Q,ef) * V[b](a,ef)
            X->gemm(false, true, bQabA, V, 0.5, 0.0);

// Form G2
#pragma omp parallel for
            for (int Q = 0; Q < nQ; ++Q) {
                for (int a = 0; a < navirA; ++a) {
                    int ab = ab_idxAA->get(a, b);
                    G->add(Q, ab, X->get(Q, a));
                }
            }
        }
        T1.reset();
        T2.reset();
        Ta.reset();
        Va.reset();
        A.reset();
        X.reset();
        V.reset();
        bQabA.reset();

        // G_AB^Q += P+(AB) 1/2 \sum(ef) V_AeBf b_ef^Q
        // V_AeBf = \sum(Mn) L_Mn^Ae T_Mn^Bf
        // Read b_ef^Q
        bQabB = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|ab)", nQ, navirB, navirB);
        bQabB->read(psio_, PSIF_DFOCC_INTS, true, true);

        T1 = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T1->read(psio_, PSIF_DFOCC_AMPS);
        if (amps == "T2") {
            T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
            T2->read(psio_, PSIF_DFOCC_AMPS);
        }
        else if (amps == "Tau") {
            Tau = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
            Tau->read(psio_, PSIF_DFOCC_AMPS);
            T2 = std::make_shared<Tensor2d>("Tau <Ij|Ab>", naoccA, naoccB, navirA, navirB);
            uccsd_tau_amps_OS(naoccA, naoccB, navirA, navirB, T2, Tau, t1A, t1B);
            Tau.reset();
        }
        else {
             std::cout << "ccd_tpdm_pplA ---> Unrecognized amps, it should be T2 or Tau ! \n";
             exit(1);
        }
        Ts = std::make_shared<Tensor2d>("T[B] (Mn,f)", naoccA * naoccB, navirB);
        V = std::make_shared<Tensor2d>("V[B] (A,ef)", navirA, navirB * navirB);
        X = std::make_shared<Tensor2d>("TPDM[B] (Q|A)", nQ, navirA);
        // Main loop
        for (int b = 0; b < navirA; ++b) {
// Form (+)T[B](Mn,f)
#pragma omp parallel for
            for (int m = 0; m < naoccA; ++m) {
                for (int n = 0; n < naoccB; ++n) {
                    int mn = n + (m * naoccB);
                    for (int f = 0; f < navirB; ++f) {
                        int bf = f + (b * navirB);
                        Ts->set(mn, f, T2->get(mn, bf));
                    }
                }
            }

            // Form V[B](A,ef) = \sum(Mn) L(Mn,Ae) T[B](Mn,f)
            V->contract(true, false, navirA * navirB, navirB, naoccA * naoccB, T1, Ts, 1.0, 0.0);

            // G2[B](Q,A) = 1/2\sum(ef) b(Q,ef) * V[B](A,ef)
            X->gemm(false, true, bQabB, V, 0.5, 0.0);

// Form G2
#pragma omp parallel for
            for (int Q = 0; Q < nQ; ++Q) {
                for (int a = 0; a < navirA; ++a) {
                    int ab = ab_idxAA->get(a, b);
                    G->add(Q, ab, X->get(Q, a));
                }
            }
        }// main b-loop
        T1.reset();
        T2.reset();
        Ts.reset();
        X.reset();
        V.reset();
        bQabB.reset();

}  // end ccd_tpdm_pplA

//======================================================================
//    PPL Beta
//======================================================================
void DFOCC::ccd_tpdm_pplB(SharedTensor2d& G, std::string amps) {

    SharedTensor2d T, T2, L2, Tau, X, Y, Z, V, U, L, G2;
    SharedTensor2d T1, Ts, Ta, Vs, Va, S, A;

        // G_ab^Q += P+(ab) 1/2 \sum(ef) V_aebf b_ef^Q
        // V_aebf = 1/2 \sum(MN) L_mn^ae T_mn^bf
        // Read b_ef^Q
        bQabB = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|ab)", nQ, navirB, navirB);
        bQabB->read(psio_, PSIF_DFOCC_INTS, true, true);

        // (-)Tt(mn, ae) = 1/2 (T_mn^ae - T_nm^ae) * (2 - \delta_{mn})
        T1 = std::make_shared<Tensor2d>("L2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        T1->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        Ta = std::make_shared<Tensor2d>("(-)Ut [i>=j|a>=b]", ntri_ijBB, ntri_abBB);
        Ta->antisymm_row_packed4(T1);

        // Anti-symmetric contributions
        if (amps == "T2") {
            T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
            T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        }
        else if (amps == "Tau") {
            Tau = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
            Tau->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
            T2 = std::make_shared<Tensor2d>("Tau <ij|ab>", naoccB, naoccB, navirB, navirB);
            uccsd_tau_amps(naoccB, naoccB, navirB, navirB, T2, Tau, t1B, t1B);
            Tau.reset();
        }
        else {
             std::cout << "ccd_tpdm_pplB --> Unrecognized amps, it should be T2 or Tau ! \n";
             exit(1);
        }
        Va = std::make_shared<Tensor2d>("(-)T[b] (f,m>=n)", navirB, ntri_ijBB);
        A = std::make_shared<Tensor2d>("A[b] (f,a>=e)", navirB, ntri_abBB);
        V = std::make_shared<Tensor2d>("V[b] (a,ef)", navirB, navirB, navirB);
        X = std::make_shared<Tensor2d>("TPDM[b] (Q|a)", nQ, navirB);
        // Main loop
        for (int b = 0; b < navirB; ++b) {
// Form (+)T[b](f, m>=n)
#pragma omp parallel for
            for (int m = 0; m < naoccB; ++m) {
                for (int n = 0; n <= m; ++n) {
                    int mn2 = index2(m, n);
                    int mn = n + (m * naoccB);
                    int nm = m + (n * naoccB);
                    for (int f = 0; f < navirB; ++f) {
                        int bf = f + (b * navirB);
                        double value2 = 0.5 * (T2->get(mn, bf) - T2->get(nm, bf));
                        Va->set(f, mn2, value2);
                    }
                }
            }

            // Form A[b](f,a>=e) = 1/2 \sum(m>=n) (+)Tt(m>=n,a>=e) * T[b](f,m>=n)
            A->contract(false, false, navirB, ntri_abBB, ntri_ijBB, Va, Ta, 0.5, 0.0);

// Form V[b](a,ef)
#pragma omp parallel for
            for (int a = 0; a < navirB; ++a) {
                for (int e = 0; e < navirB; ++e) {
                    int ae = index2(a, e);
                    for (int f = 0; f < navirB; ++f) {
                        int ef = ab_idxBB->get(e, f);
                        int perm1 = (a > e) ? 1 : -1;
                        double value = perm1 * A->get(f, ae);
                        V->set(a, ef, value);
                    }
                }
            }

            // G2[b](Q,a) = 1/2\sum(ef) b(Q,ef) * V[b](a,ef)
            X->gemm(false, true, bQabB, V, 0.5, 0.0);

// Form G2
#pragma omp parallel for
            for (int Q = 0; Q < nQ; ++Q) {
                for (int a = 0; a < navirB; ++a) {
                    int ab = ab_idxBB->get(a, b);
                    G->add(Q, ab, X->get(Q, a));
                }
            }
        }
        T1.reset();
        T2.reset();
        Ta.reset();
        Va.reset();
        A.reset();
        X.reset();
        V.reset();
        bQabB.reset();

        // G_ab^Q += P+(ab) 1/2 \sum(ef) V_EaFb b_EF^Q
        // V_EaFb = \sum(Mn) L_Mn^Ea T_Mn^Fb
        // Read b_EF^Q
        bQabA = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA);
        bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);

        T1 = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T1->read(psio_, PSIF_DFOCC_AMPS);
        if (amps == "T2") {
            T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
            T2->read(psio_, PSIF_DFOCC_AMPS);
        }
        else if (amps == "Tau") {
            Tau = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
            Tau->read(psio_, PSIF_DFOCC_AMPS);
            T2 = std::make_shared<Tensor2d>("Tau <Ij|Ab>", naoccA, naoccB, navirA, navirB);
            uccsd_tau_amps_OS(naoccA, naoccB, navirA, navirB, T2, Tau, t1A, t1B);
            Tau.reset();
        }
        else {
             std::cout << "ccd_tpdm_pplB ---> Unrecognized amps, it should be T2 or Tau ! \n";
             exit(1);
        }
        Ts = std::make_shared<Tensor2d>("T[b] (Mn,F)", naoccA * naoccB, navirA);
        V = std::make_shared<Tensor2d>("V[b] (Ea,F)", navirA * navirB, navirA);
        Vs = std::make_shared<Tensor2d>("V[b] (EF,a)", navirA * navirA, navirB);
        X = std::make_shared<Tensor2d>("TPDM[b] (Q|a)", nQ, navirB);
        // Main loop
        for (int b = 0; b < navirB; ++b) {
// Form (+)T[b](Mn,F)
#pragma omp parallel for
            for (int m = 0; m < naoccA; ++m) {
                for (int n = 0; n < naoccB; ++n) {
                    int mn = n + (m * naoccB);
                    for (int f = 0; f < navirA; ++f) {
                        int fb = (f * navirB) + b;
                        Ts->set(mn, f, T2->get(mn, fb));
                    }
                }
            }

            // Form V[b](Ea,F) = \sum(Mn) L(Mn,Ea) T[b](Mn,F)
            V->gemm(true, false, T1, Ts, 1.0, 0.0);

// Form Vs[b](EF,a)
#pragma omp parallel for
            for (int e = 0; e < navirA; ++e) {
                for (int f = 0; f < navirA; ++f) {
                    int ef = ab_idxAA->get(e, f);
                    for (int a = 0; a < navirB; ++a) {
                        int ea = (e * navirB) + a;
                        Vs->set(ef, a, V->get(ea, f));
                    }
                }
            }

            // G2[b](Q,a) = 1/2 \sum(EF) b(Q,EF) * Vs[b](EF,a)
            X->gemm(false, false, bQabA, Vs, 0.5, 0.0);

// Form G2
#pragma omp parallel for
            for (int Q = 0; Q < nQ; ++Q) {
                for (int a = 0; a < navirB; ++a) {
                    int ab = ab_idxBB->get(a, b);
                    G->add(Q, ab, X->get(Q, a));
                }
            }
        }
        T1.reset();
        T2.reset();
        Ts.reset();
        X.reset();
        V.reset();
        Vs.reset();
        bQabA.reset();

}  // end ccd_tpdm_pplB

}  // namespace dfoccwave
}  // namespace psi
