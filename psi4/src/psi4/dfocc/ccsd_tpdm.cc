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

/** Standard library includes */
#include "psi4/libqt/qt.h"
#include "defines.h"
#include "dfocc.h"


using namespace psi;
using namespace std;


namespace psi{ namespace dfoccwave{

void DFOCC::ccsd_tpdm()
{

    timer_on("tpdm");
    SharedTensor2d T, U, Tau, G, G2, V, X, Y, Z;
    SharedTensor2d Ts, Ta, Vs, Va, S, A;

//if (reference_ == "RESTRICTED") {

    //============================
    // OO-Block Correlation TPDM
    //============================

    // G_ij^Q = P+(ij) 2*\cal(V)_ij^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|IJ)", nQ, naoccA, naoccA));
    V = SharedTensor2d(new Tensor2d("calV (Q|IJ)", nQ, naoccA, naoccA));
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, 2.0);
    V.reset();
    // G_ij^Q += P+(ij) V_ij^Q
    V = SharedTensor2d(new Tensor2d("V (Q|IJ)", nQ, naoccA, naoccA));
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, 1.0);
    V.reset();
    // G_ij^Q += P+(ij) Vt_ij^Q
    V = SharedTensor2d(new Tensor2d("Vt (Q|IJ)", nQ, naoccA, naoccA));
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, 1.0);
    V.reset();
    // G_ij^Q -= P+(ij) 2*V'_ij^Q
    V = SharedTensor2d(new Tensor2d("Vp (Q|IJ)", nQ, naoccA, naoccA));
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, -2.0);
    V.reset();
    // G_ij^Q -= P+(ij) Z_ij^Q
    V = SharedTensor2d(new Tensor2d("Zeta (Q|IJ)", nQ, naoccA, naoccA));
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, -1.0);
    V.reset();
    // G_ij^Q -= P+(ij) G_ij t_Q
    G->dirprd123(T1c, GijA, -1.0, 1.0);

    // G_ij^Q += P+(ij) \sum(m) G_mj t_im^Q
    T = SharedTensor2d(new Tensor2d("T1 (Q|IJ)", nQ, naoccA, naoccA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->contract(false, false, nQ*naoccA, naoccA, naoccA, T, GijA, 1.0, 1.0);
    T.reset();

    // G_ij^Q -= P+(ij) \sum(e) l_j^e (t_ie^Q + Tau_ie^Q)
    T = SharedTensor2d(new Tensor2d("T1 (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    U = SharedTensor2d(new Tensor2d("Tau (Q|IA)", nQ, naoccA, navirA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();
    G->contract(false, true, nQ*naoccA, naoccA, navirA, T, l1A, -1.0, 1.0);
    T.reset();

    // G_ij^Q -= P+(ij) \sum(e) t_j^e (L_ie^Q + 2*V_ei^Q)
    V = SharedTensor2d(new Tensor2d("V (Q|AI)", nQ, navirA, naoccA));
    V->read(psio_, PSIF_DFOCC_AMPS);
    T = SharedTensor2d(new Tensor2d("Temp (Q|IA)", nQ, naoccA, navirA));
    T->swap_3index_col(V);
    V.reset();
    T->scale(2.0);
    U = SharedTensor2d(new Tensor2d("L2 (Q|IA)", nQ, naoccA, navirA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->add(U);
    U.reset();
    G->contract(false, true, nQ*naoccA, naoccA, navirA, T, t1A, -1.0, 1.0);
    T.reset();

    // symmetrize
    G->symmetrize3(G);
    G->scale(2.0);
    G2 = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OO)", nQ, noccA, noccA));
    G2->set3_act_oo(nfrzc, G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS);
    if(print_ > 3) G2->print();
    G2.reset();

    //============================
    // OV-Block Correlation TPDM
    //============================

    // G_ia^Q = Tau_ia^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|IA)", nQ, naoccA, navirA));
    T = SharedTensor2d(new Tensor2d("Tau (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, 1.0);
    T.reset();
    // G_ia^Q += L_ia^Q
    T = SharedTensor2d(new Tensor2d("L2 (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, 1.0);
    T.reset();
    // G_ia^Q += Z_ia^Q
    T = SharedTensor2d(new Tensor2d("Zeta (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, 1.0);
    T.reset();
    // G_ia^Q += 2*y_ia^Q
    T = SharedTensor2d(new Tensor2d("Y (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, 2.0);
    T.reset();
    // G_ia^Q -= tEta_ia^Q
    T = SharedTensor2d(new Tensor2d("Eta2 (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, -1.0);
    T.reset();
    // G_ia^Q += l_i^a t_Q
    G->dirprd123(T1c, l1A, 1.0, 1.0);
    // G_ia^Q += t_i^a (l_Q - G_Q - Gt_Q)
    SharedTensor1d Tt = SharedTensor1d(new Tensor1d("TEMP", nQ));
    Tt->copy(L1c);
    Tt->axpy(gQ, -1.0);
    Tt->axpy(gQt, -1.0);
    G->dirprd123(Tt, t1A, 1.0, 1.0);
    Tt.reset();

    // G_ia^Q += \sum(m) t_m^a (G_im^Q + V_im^Q - 2*Vp_im^Q - n_mi^Q)
    T = SharedTensor2d(new Tensor2d("Temp (Q|IJ)", nQ, naoccA, naoccA));
    U = SharedTensor2d(new Tensor2d("Eta (Q|IJ)", nQ, naoccA, naoccA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->swap_3index_col(U);
    U.reset();
    T->scale(-1.0);
    U = SharedTensor2d(new Tensor2d("G (Q|IJ)", nQ, naoccA, naoccA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();
    U = SharedTensor2d(new Tensor2d("V (Q|IJ)", nQ, naoccA, naoccA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();
    U = SharedTensor2d(new Tensor2d("Vp (Q|IJ)", nQ, naoccA, naoccA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, -2.0);
    U.reset();
    G->contract(false, false, nQ*naoccA, navirA, naoccA, T, t1A, 1.0, 1.0);
    T.reset();

    // G_ia^Q -= \sum(m) Tau_ma^Q Gt_im
    T = SharedTensor2d(new Tensor2d("Tau (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->contract233(false, false, naoccA, navirA, GtijA, T, -1.0, 1.0);

    // G_ia^Q += \sum(e) Tau_ie^Q Gt_ea
    G->contract(false, false, nQ*naoccA, navirA, navirA, T, GtabA, 1.0, 1.0);
    T.reset();

    // G_ia^Q += 2\sum(e) t_i^e V_ea^Q
    T = SharedTensor2d(new Tensor2d("V (Q|AB)", nQ, navirA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->contract233(false, false, naoccA, navirA, t1A, T, 2.0, 1.0);
    T.reset();

    // G_ia^Q += \sum(e) t_ie^Q G_ea
    T = SharedTensor2d(new Tensor2d("T1 (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->contract(false, false, nQ*naoccA, navirA, navirA, T, GabA, 1.0, 1.0);
    T.reset();

    // G_ia^Q += \sum(me) { Ut(ia,me) - 2*X(ia,me) } (t_me^Q - t_em^Q)
    // X(ia,me) = Y(ie,am)
    // Y(ie,am) = 2*V_ieam - V_iema
    T = SharedTensor2d(new Tensor2d("Temp (Q|IA)", nQ, naoccA, navirA));
    U = SharedTensor2d(new Tensor2d("T1 (Q|AI)", nQ, navirA, naoccA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->swap_3index_col(U);
    U.reset();
    T->scale(-1.0);
    U = SharedTensor2d(new Tensor2d("T1 (Q|IA)", nQ, naoccA, navirA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();
    Y = SharedTensor2d(new Tensor2d("Ytemp (IA|BJ)", naoccA, navirA, navirA, naoccA));
    V = SharedTensor2d(new Tensor2d("V (IA|BJ)", naoccA, navirA, navirA, naoccA));
    V->read(psio_, PSIF_DFOCC_AMPS);
    Y->axpy(V, 2.0);
    V.reset();
    V = SharedTensor2d(new Tensor2d("V (IA|JB)", naoccA, navirA, naoccA, navirA));
    V->read(psio_, PSIF_DFOCC_AMPS);
    Y->sort(1243, V, -1.0, 1.0);
    V.reset();
    X = SharedTensor2d(new Tensor2d("X (IA|ME)", naoccA, navirA, naoccA, navirA));
    X->sort(1342, Y, 1.0, 0.0);
    Y.reset();
    U = SharedTensor2d(new Tensor2d("Ut2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    U->read_symm(psio_, PSIF_DFOCC_AMPS);
    U->axpy(X, -2.0);
    X.reset();
    G->gemm(false, true, T, U, 1.0, 1.0);
    U.reset();
    T.reset();

    // G_ia^Q += \sum(me) Ubar(ia,me) (Gt_em^Q - Gt_me^Q + l_me^Q - l_em^Q)
    U = SharedTensor2d(new Tensor2d("Gt (Q|AI)", nQ, navirA, naoccA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    X = SharedTensor2d(new Tensor2d("L1 (Q|AI)", nQ, navirA, naoccA));
    X->read(psio_, PSIF_DFOCC_AMPS);
    U->axpy(X, -1.0);
    X.reset();
    T = SharedTensor2d(new Tensor2d("Temp (Q|IA)", nQ, naoccA, navirA));
    T->swap_3index_col(U);
    U.reset();
    U = SharedTensor2d(new Tensor2d("Gt (Q|IA)", nQ, naoccA, navirA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, -1.0);
    U.reset();
    U = SharedTensor2d(new Tensor2d("L1 (Q|IA)", nQ, naoccA, navirA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();
    U = SharedTensor2d(new Tensor2d("2*Tau(ia,jb) - Tau(ib,ja)", naoccA, navirA, naoccA, navirA));
    U->read_symm(psio_, PSIF_DFOCC_AMPS);
    G->gemm(false, true, T, U, 1.0, 1.0);
    U.reset();
    T.reset();

    // Form overall OV Block
    G2 = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA, nvirA));
    G2->set3_act_ov(nfrzc, naoccA, navirA, nvirA, G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS);
    if(print_ > 3) G2->print();

    // Form G_vo^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VO)", nQ, nvirA, noccA));
    G->swap_3index_col(G2);
    G2.reset();
    G->write(psio_, PSIF_DFOCC_DENS);
    if(print_ > 3) G->print();
    G.reset();

    //============================
    // VV-Block Correlation TPDM
    //============================

    // G_ab^Q -= P+(ab) 2*V_ab^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|AB)", nQ, navirA, navirA));
    V = SharedTensor2d(new Tensor2d("V (Q|AB)", nQ, navirA, navirA));
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, -2.0);
    V.reset();
    // G_ab^Q -= P+(ab) 2*Vt_ab^Q
    V = SharedTensor2d(new Tensor2d("Vt (Q|AB)", nQ, navirA, navirA));
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, -2.0);
    V.reset();
    // G_ab^Q -= P+(ab) tEta_ab^Q
    V = SharedTensor2d(new Tensor2d("Eta2 (Q|AB)", nQ, navirA, navirA));
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, -1.0);
    V.reset();
    // G_ab^Q -= P+(ab) G_ab t_Q
    G->dirprd123(T1c, GabA, -1.0, 1.0);

    // G_ab^Q -= P+(ab) \sum(m) (t_am^Q - Tau_ma^Q) l_m^b
    T = SharedTensor2d(new Tensor2d("Temp (Q|AI)", nQ, navirA, naoccA));
    U = SharedTensor2d(new Tensor2d("Tau (Q|IA)", nQ, naoccA, navirA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->swap_3index_col(U);
    U.reset();
    T->scale(-1.0);
    U = SharedTensor2d(new Tensor2d("T1 (Q|AI)", nQ, navirA, naoccA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();
    G->contract(false, false, nQ*navirA, navirA, naoccA, T, l1A, -1.0, 1.0);
    T.reset();

    // G_ab^Q += P+(ab) \sum(m) (Eta_ma^Q + G_am^Q + L_ma^Q + 2*V_am^Q) t_m^b
    U = SharedTensor2d(new Tensor2d("Eta (Q|IA)", nQ, naoccA, navirA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T = SharedTensor2d(new Tensor2d("L2 (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    U->add(T);
    T.reset();
    T = SharedTensor2d(new Tensor2d("Temp (Q|AI)", nQ, navirA, naoccA));
    T->swap_3index_col(U);
    U.reset();
    U = SharedTensor2d(new Tensor2d("G (Q|AI)", nQ, navirA, naoccA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();
    U = SharedTensor2d(new Tensor2d("V (Q|AI)", nQ, navirA, naoccA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 2.0);
    U.reset();
    G->contract(false, false, nQ*navirA, navirA, naoccA, T, t1A, 1.0, 1.0);
    T.reset();


    // G_ab^Q += P+(ab) \sum(ef) Vt_aebf b_ef^Q
    // Read b_ef^Q
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);

    // (+)Ut(mn, ae) = 1/2 (Ut_mn^ae + Ut_nm^ae) * (2 - \delta_{mn})
    // (-)Ut(mn, ae) = 1/2 (Ut_mn^ae - Ut_nm^ae) * (2 - \delta_{mn})
    T = SharedTensor2d(new Tensor2d("Ut2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    T->read_symm(psio_, PSIF_DFOCC_AMPS);
    U = SharedTensor2d(new Tensor2d("Ut2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    U->sort(1324, T, 1.0, 0.0);
    T.reset();
    Ts = SharedTensor2d(new Tensor2d("(+)Ut [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    Ta = SharedTensor2d(new Tensor2d("(-)Ut [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    Ts->symm_row_packed4(U);
    Ta->antisymm_row_packed4(U);
    U.reset();

    // Symmetric & Anti-symmetric contributions
    U = SharedTensor2d(new Tensor2d("Tau (IA|JB)", naoccA, navirA, naoccA, navirA));
    U->read_symm(psio_, PSIF_DFOCC_AMPS);
    Tau = SharedTensor2d(new Tensor2d("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA));
    Tau->sort(1324, U, 1.0, 0.0);
    U.reset();
    Vs = SharedTensor2d(new Tensor2d("(+)Tau[B] (F,M>=N)", navirA, ntri_ijAA));
    Va = SharedTensor2d(new Tensor2d("(-)Tau[B] (F,M>=N)", navirA, ntri_ijAA));
    S = SharedTensor2d(new Tensor2d("S[B] (F,A>=E)", navirA, ntri_abAA));
    A = SharedTensor2d(new Tensor2d("A[B] (F,A>=E)", navirA, ntri_abAA));
    V = SharedTensor2d(new Tensor2d("V[B] (A,EF)", navirA, navirA, navirA));
    X = SharedTensor2d(new Tensor2d("TPDM[B] (Q|A)", nQ, navirA));
    // Main loop
    for(int b = 0 ; b < navirA; ++b){

	    // does not work!
	    /*
	    int nf = b+1;
            // Form (+)Tau[b](f, m>=n)
            #pragma omp parallel for
            for(int m = 0 ; m < naoccA; ++m){
                for(int n = 0 ; n <= m; ++n){
                    int mn2 = index2(m,n);
                    int mn = n + (m * naoccA);
                    int nm = m + (n * naoccA);
                    for(int f = 0 ; f <= b; ++f){
                        int bf = f + (b * navirA);
                        double value1 = 0.5 * ( Tau->get(mn, bf) + Tau->get(nm, bf) );
                        double value2 = 0.5 * ( Tau->get(mn, bf) - Tau->get(nm, bf) );
                        Vs->set(f, mn2, value1);
                        Va->set(f, mn2, value2);
                    }
                }
            }

            // Form S[b](f,a>=e) = \sum(m>=n) (+)Ut(m>=n,a>=e) * Tau[b](f,m>=n)
	    S->contract(false, false, nf, ntri_abAA, ntri_ijAA, Vs, Ts, 1.0, 0.0);
	    A->contract(false, false, nf, ntri_abAA, ntri_ijAA, Va, Ta, 1.0, 0.0);

	    // Form V[b](a,ef)
            #pragma omp parallel for
            for(int a = 0 ; a < navirA; ++a){
                for(int e = 0 ; e < navirA; ++e){
                    int ae = index2(a,e);
                    for(int f = 0 ; f <= b; ++f){
                        int ef = ab_idxAA->get(e,f);
			int perm1 = ( a > e ) ? 1 : -1;
			int perm2 = ( b > f ) ? 1 : -1;
			int perm3 = ( b != f) ? 2 : 1;
			double value = S->get(f,ae) + (perm1 * perm2 * A->get(f,ae));
			value *= perm3;
                        V->set(a, ef, value);
                    }
                }
            }
	    */
	    //Does not work!

            // Form (+)Tau[b](f, m>=n)
            #pragma omp parallel for
            for(int m = 0 ; m < naoccA; ++m){
                for(int n = 0 ; n <= m; ++n){
                    int mn2 = index2(m,n);
                    int mn = n + (m * naoccA);
                    int nm = m + (n * naoccA);
                    for(int f = 0 ; f < navirA; ++f){
                        int bf = f + (b * navirA);
                        double value1 = 0.5 * ( Tau->get(mn, bf) + Tau->get(nm, bf) );
                        double value2 = 0.5 * ( Tau->get(mn, bf) - Tau->get(nm, bf) );
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
            for(int a = 0 ; a < navirA; ++a){
                for(int e = 0 ; e < navirA; ++e){
                    int ae = index2(a,e);
                    for(int f = 0 ; f < navirA; ++f){
                        int ef = ab_idxAA->get(e,f);
			int perm1 = ( a > e ) ? 1 : -1;
			double value = S->get(f,ae) + (perm1 * A->get(f,ae));
                        V->set(a, ef, value);
                    }
                }
            }

	    // G2[b](Q,a) = \sum(ef) b(Q,ef) * V[b](a,ef)
            X->gemm(false, true, bQabA, V, 1.0, 0.0);

	    // Form G2
            #pragma omp parallel for
            for(int Q = 0 ; Q < nQ; ++Q){
                for(int a = 0 ; a < navirA; ++a){
                    int ab = ab_idxAA->get(a,b);
		    G->add(Q, ab, X->get(Q,a));
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

    /*
    // Read Tau
    U = SharedTensor2d(new Tensor2d("Tau (IA|JB)", naoccA, navirA, naoccA, navirA));
    U->read_symm(psio_, PSIF_DFOCC_AMPS);
    Tau = SharedTensor2d(new Tensor2d("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA));
    Tau->sort(1324, U, 1.0, 0.0);
    U.reset();
    // Read Ut
    T = SharedTensor2d(new Tensor2d("Ut2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    T->read_symm(psio_, PSIF_DFOCC_AMPS);
    U = SharedTensor2d(new Tensor2d("Ut2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    U->sort(1324, T, 1.0, 0.0);
    T.reset();
    // Vt_(ab,cd) = \sum(mn) Ut(mn,ab) Tau(mn,cd)
    V = SharedTensor2d(new Tensor2d("V <AB|CD>", navirA, navirA, navirA, navirA));
    V->gemm(true, false, U, Tau, 1.0, 0.0);
    Tau.reset();
    U.reset();
    X = SharedTensor2d(new Tensor2d("X <BD|AC>", navirA, navirA, navirA, navirA));
    X->sort(2413, V, 1.0, 0.0);
    V.reset();
    // Read b_ef^Q
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);
    // G_ab^Q += P+(ab) \sum(ef) Vt_aebf b_ef^Q
    G->gemm(false, false, bQabA, X, 1.0, 1.0);
    X.reset();
    bQabA.reset();
    */

    // symmetrize
    //G->print();
    G->symmetrize3(G);
    G->scale(2.0);
    G2 = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VV)", nQ, nvirA, nvirA));
    G2->set3_act_vv(G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS, true, true);
    if(print_ > 3) G2->print();
    G2.reset();

//}// end if (reference_ == "RESTRICTED")
//else if (reference_ == "UNRESTRICTED") {
//}// else if (reference_ == "UNRESTRICTED")
    timer_off("tpdm");
} // end ccsd_tpdm

}} // End Namespaces
