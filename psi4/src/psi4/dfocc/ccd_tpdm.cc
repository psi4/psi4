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

void DFOCC::ccd_tpdm()
{

    timer_on("tpdm");
    SharedTensor2d T, U, Tau, G, G2, V, X, Y, Z;
    SharedTensor2d Ts, Ta, Vs, Va, S, A;

//if (reference_ == "RESTRICTED") {

    //============================
    // OO-Block Correlation TPDM
    //============================

    // G_ij^Q += P+(ij) V_ij^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|IJ)", nQ, naoccA, naoccA));
    V = SharedTensor2d(new Tensor2d("V (Q|IJ)", nQ, naoccA, naoccA));
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, 1.0);
    V.reset();
    // G_ij^Q -= P+(ij) 2*V'_ij^Q
    V = SharedTensor2d(new Tensor2d("Vp (Q|IJ)", nQ, naoccA, naoccA));
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, -2.0);
    V.reset();

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

    // G_ia^Q = T_ia^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|IA)", nQ, naoccA, navirA));
    T = SharedTensor2d(new Tensor2d("T2 (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, 1.0);
    T.reset();
    // G_ia^Q += L_ia^Q
    T = SharedTensor2d(new Tensor2d("L2 (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, 1.0);
    T.reset();
    // G_ia^Q += 2*y_ia^Q
    T = SharedTensor2d(new Tensor2d("Y (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, 2.0);
    T.reset();

    // G_ia^Q -= \sum(m) T_ma^Q G_im
    T = SharedTensor2d(new Tensor2d("T2 (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->contract233(false, false, naoccA, navirA, GijA, T, -1.0, 1.0);

    // G_ia^Q += \sum(e) T_ie^Q G_ea
    G->contract(false, false, nQ*naoccA, navirA, navirA, T, GabA, 1.0, 1.0);
    //G->cont323("IA", "IE", "EA", false, T, GabA, 1.0, 1.0); // it works
    T.reset();

    // G_ia^Q += \sum(me) U(ia,me) (G_em^Q - G_me^Q)
    U = SharedTensor2d(new Tensor2d("G (Q|AI)", nQ, navirA, naoccA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T = SharedTensor2d(new Tensor2d("Temp (Q|IA)", nQ, naoccA, navirA));
    T->swap_3index_col(U);
    U.reset();
    U = SharedTensor2d(new Tensor2d("G (Q|IA)", nQ, naoccA, navirA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, -1.0);
    U.reset();
    U = SharedTensor2d(new Tensor2d("U2 (IA|JB)", naoccA, navirA, naoccA, navirA));
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
    U = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    U->read_symm(psio_, PSIF_DFOCC_AMPS);
    Tau = SharedTensor2d(new Tensor2d("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    Tau->sort(1324, U, 1.0, 0.0);
    U.reset();
    Vs = SharedTensor2d(new Tensor2d("(+)T[B] (F,M>=N)", navirA, ntri_ijAA));
    Va = SharedTensor2d(new Tensor2d("(-)T[B] (F,M>=N)", navirA, ntri_ijAA));
    S = SharedTensor2d(new Tensor2d("S[B] (F,A>=E)", navirA, ntri_abAA));
    A = SharedTensor2d(new Tensor2d("A[B] (F,A>=E)", navirA, ntri_abAA));
    V = SharedTensor2d(new Tensor2d("V[B] (A,EF)", navirA, navirA, navirA));
    X = SharedTensor2d(new Tensor2d("TPDM[B] (Q|A)", nQ, navirA));
    // Main loop
    for(int b = 0 ; b < navirA; ++b){

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

    // symmetrize
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
} // end ccd_tpdm

}} // End Namespaces
