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

void DFOCC::olccd_tpdm()
{

    timer_on("tpdm");
    SharedTensor2d T, T1, U, Tau, G, G2, V, X, Y, Z;
    SharedTensor2d Ts, Ta, Vs, Va, S, A;

if (reference_ == "RESTRICTED") {

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

    // G_ia^Q = 2*T_ia^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|IA)", nQ, naoccA, navirA));
    T = SharedTensor2d(new Tensor2d("T2 (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, 2.0);
    T.reset();
    // G_ia^Q += 2*y_ia^Q
    T = SharedTensor2d(new Tensor2d("Y (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, 2.0);
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

    // (+)Ut(mn, ae) = 1/2 (U_mn^ae + U_nm^ae) * (2 - \delta_{mn})
    // (-)Ut(mn, ae) = 1/2 (U_mn^ae - U_nm^ae) * (2 - \delta_{mn})
    T1 = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    T1->read_symm(psio_, PSIF_DFOCC_AMPS);
    T = SharedTensor2d(new Tensor2d("U2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_u2_amps(T,T1);
    U = SharedTensor2d(new Tensor2d("U2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    U->sort(1324, T, 1.0, 0.0);
    T.reset();
    Ts = SharedTensor2d(new Tensor2d("(+)Ut [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    Ta = SharedTensor2d(new Tensor2d("(-)Ut [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    Ts->symm_row_packed4(U);
    Ta->antisymm_row_packed4(U);
    U.reset();

    // Symmetric & Anti-symmetric contributions
    Tau = SharedTensor2d(new Tensor2d("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    Tau->sort(1324, T1, 1.0, 0.0);
    T1.reset();
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

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    //============================
    // OO-Block Correlation TPDM
    //============================

    // Build G_IJ^Q
    // G_IJ^Q += P+(IJ) V_IJ^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|IJ)", nQ, naoccA, naoccA));
    V = SharedTensor2d(new Tensor2d("V (Q|IJ)", nQ, naoccA, naoccA));
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, 1.0);
    V.reset();
    // G_IJ^Q -= P+(IJ) 2*V'_IJ^Q
    V = SharedTensor2d(new Tensor2d("Vp (Q|IJ)", nQ, naoccA, naoccA));
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, -2.0);
    V.reset();

    // Active to full
    //G->symmetrize3(G);
    //G->scale(2.0);
    G2 = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OO)", nQ, noccA, noccA));
    G2->set3_act_oo(nfrzc, G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS);
    if(print_ > 3) G2->print();
    G2.reset();

    // Build G_ij^Q
    // G_ij^Q += P+(ij) V_ij^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|ij)", nQ, naoccB, naoccB));
    V = SharedTensor2d(new Tensor2d("V (Q|ij)", nQ, naoccB, naoccB));
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, 1.0);
    V.reset();
    // G_ij^Q -= P+(ij) 2*V'_ij^Q
    V = SharedTensor2d(new Tensor2d("Vp (Q|ij)", nQ, naoccB, naoccB));
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, -2.0);
    V.reset();

    // Active to full
    //G->symmetrize3(G);
    //G->scale(2.0);
    G2 = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|oo)", nQ, noccB, noccB));
    G2->set3_act_oo(nfrzc, G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS);
    if(print_ > 3) G2->print();
    G2.reset();

    //outfile->Printf("\tOO block done.\n");

    //============================
    // OV-Block Correlation TPDM
    //============================

    // Build G_IA^Q
    // G_IA^Q = T_IA^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|IA)", nQ, naoccA, navirA));
    T = SharedTensor2d(new Tensor2d("T2 (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, 1.0);
    T.reset();
    // G_IA^Q += y_IA^Q
    T = SharedTensor2d(new Tensor2d("Y (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, 1.0);
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

    // Build G_ia^Q
    // G_ia^Q = T_ia^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|ia)", nQ, naoccB, navirB));
    T = SharedTensor2d(new Tensor2d("T2 (Q|ia)", nQ, naoccB, navirB));
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, 1.0);
    T.reset();
    // G_ia^Q += y_ia^Q
    T = SharedTensor2d(new Tensor2d("Y (Q|ia)", nQ, naoccB, navirB));
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, 1.0);
    T.reset();

    // Form overall OV Block
    G2 = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|ov)", nQ, noccB, nvirB));
    G2->set3_act_ov(nfrzc, naoccB, navirB, nvirB, G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS);
    if(print_ > 3) G2->print();

    // Form G_vo^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|vo)", nQ, nvirB, noccB));
    G->swap_3index_col(G2);
    G2.reset();
    G->write(psio_, PSIF_DFOCC_DENS);
    if(print_ > 3) G->print();
    G.reset();

    //outfile->Printf("\tOV block done.\n");

    //============================
    // VV-Block Correlation TPDM
    //============================

    // Build G_AB^Q
    // G_AB^Q -= 2*V_AB^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|AB)", nQ, navirA, navirA));
    V = SharedTensor2d(new Tensor2d("V (Q|AB)", nQ, navirA, navirA));
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, -2.0);
    V.reset();

    // G_AB^Q += \sum(EF) V_AEBF b_EF^Q
    // V_AEBF = 1/2 \sum(MN) L_MN^AE T_MN^BF
    // Read b_EF^Q
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);

    // (-)Tt(mn, ae) = 1/2 (T_mn^ae - T_nm^ae) * (2 - \delta_{mn})
    T1 = SharedTensor2d(new Tensor2d("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    T1->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Ta = SharedTensor2d(new Tensor2d("(-)Ut [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    Ta->antisymm_row_packed4(T1);

    // Anti-symmetric contributions
    Va = SharedTensor2d(new Tensor2d("(-)T[B] (F,M>=N)", navirA, ntri_ijAA));
    A = SharedTensor2d(new Tensor2d("A[B] (F,A>=E)", navirA, ntri_abAA));
    V = SharedTensor2d(new Tensor2d("V[B] (A,EF)", navirA, navirA, navirA));
    X = SharedTensor2d(new Tensor2d("TPDM[B] (Q|A)", nQ, navirA));
    // Main loop
    for(int b = 0 ; b < navirA; ++b){

            // Form (+)T[b](f, m>=n)
            #pragma omp parallel for
            for(int m = 0 ; m < naoccA; ++m){
                for(int n = 0 ; n <= m; ++n){
                    int mn2 = index2(m,n);
                    int mn = n + (m * naoccA);
                    int nm = m + (n * naoccA);
                    for(int f = 0 ; f < navirA; ++f){
                        int bf = f + (b * navirA);
                        double value2 = 0.5 * ( T1->get(mn, bf) - T1->get(nm, bf) );
                        Va->set(f, mn2, value2);
                    }
                }
            }

            // Form A[b](f,a>=e) = 1/2 \sum(m>=n) (+)Tt(m>=n,a>=e) * T[b](f,m>=n)
	    A->contract(false, false, navirA, ntri_abAA, ntri_ijAA, Va, Ta, 0.5, 0.0);

	    // Form V[b](a,ef)
            #pragma omp parallel for
            for(int a = 0 ; a < navirA; ++a){
                for(int e = 0 ; e < navirA; ++e){
                    int ae = index2(a,e);
                    for(int f = 0 ; f < navirA; ++f){
                        int ef = ab_idxAA->get(e,f);
			int perm1 = ( a > e ) ? 1 : -1;
			double value = perm1 * A->get(f,ae);
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
    T1.reset();
    Ta.reset();
    Va.reset();
    A.reset();
    X.reset();
    V.reset();
    bQabA.reset();

    // G_AB^Q += \sum(ef) V_AeBf b_ef^Q
    // V_AeBf = \sum(Mn) L_Mn^Ae T_Mn^Bf
    // Read b_ef^Q
    bQabB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ab)", nQ, navirB, navirB));
    bQabB->read(psio_, PSIF_DFOCC_INTS, true, true);

    T1 = SharedTensor2d(new Tensor2d("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    T1->read(psio_, PSIF_DFOCC_AMPS);
    Ts = SharedTensor2d(new Tensor2d("T[B] (Mn,f)", naoccA*naoccB, navirB));
    V = SharedTensor2d(new Tensor2d("V[B] (A,ef)", navirA, navirB*navirB));
    X = SharedTensor2d(new Tensor2d("TPDM[B] (Q|A)", nQ, navirA));
    // Main loop
    for(int b = 0 ; b < navirA; ++b){

            // Form (+)T[B](Mn,f)
            #pragma omp parallel for
            for(int m = 0 ; m < naoccA; ++m){
                for(int n = 0 ; n < naoccB; ++n){
                    int mn = n + (m * naoccB);
                    for(int f = 0 ; f < navirB; ++f){
                        int bf = f + (b * navirB);
                        Ts->set(mn, f, T1->get(mn,bf));
                    }
                }
            }

	    // Form V[B](A,ef) = \sum(Mn) L(Mn,Ae) T[B](Mn,f)
	    V->contract(true, false, navirA*navirB, navirB, naoccA*naoccB, T1, Ts, 1.0, 0.0);

	    // G2[B](Q,A) = \sum(ef) b(Q,ef) * V[B](A,ef)
            X->gemm(false, true, bQabB, V, 1.0, 0.0);

	    // Form G2
            #pragma omp parallel for
            for(int Q = 0 ; Q < nQ; ++Q){
                for(int a = 0 ; a < navirA; ++a){
                    int ab = ab_idxAA->get(a,b);
		    G->add(Q, ab, X->get(Q,a));
                }
            }

    }
    T1.reset();
    Ts.reset();
    X.reset();
    V.reset();
    bQabB.reset();

    // Active to full
    //G->symmetrize3(G);
    //G->scale(2.0);
    G2 = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VV)", nQ, nvirA, nvirA));
    G2->set3_act_vv(G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS, true, true);
    if(print_ > 3) G2->print();
    G2.reset();

    //outfile->Printf("\tVV block done.\n");

    // Build G_ab^Q
    // G_ab^Q -= P+(ab) 2*V_ab^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|ab)", nQ, navirB, navirB));
    V = SharedTensor2d(new Tensor2d("V (Q|ab)", nQ, navirB, navirB));
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, -2.0);
    V.reset();

    // G_ab^Q += \sum(ef) V_aebf b_ef^Q
    // V_aebf = 1/2 \sum(MN) L_mn^ae T_mn^bf
    // Read b_ef^Q
    bQabB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ab)", nQ, navirB, navirB));
    bQabB->read(psio_, PSIF_DFOCC_INTS, true, true);

    // (-)Tt(mn, ae) = 1/2 (T_mn^ae - T_nm^ae) * (2 - \delta_{mn})
    T1 = SharedTensor2d(new Tensor2d("T2 <ij|ab>", naoccB, naoccB, navirB, navirB));
    T1->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Ta = SharedTensor2d(new Tensor2d("(-)Ut [i>=j|a>=b]", ntri_ijBB, ntri_abBB));
    Ta->antisymm_row_packed4(T1);

    // Anti-symmetric contributions
    Va = SharedTensor2d(new Tensor2d("(-)T[b] (f,m>=n)", navirB, ntri_ijBB));
    A = SharedTensor2d(new Tensor2d("A[b] (f,a>=e)", navirB, ntri_abBB));
    V = SharedTensor2d(new Tensor2d("V[b] (a,ef)", navirB, navirB, navirB));
    X = SharedTensor2d(new Tensor2d("TPDM[b] (Q|a)", nQ, navirB));
    // Main loop
    for(int b = 0 ; b < navirB; ++b){

            // Form (+)T[b](f, m>=n)
            #pragma omp parallel for
            for(int m = 0 ; m < naoccB; ++m){
                for(int n = 0 ; n <= m; ++n){
                    int mn2 = index2(m,n);
                    int mn = n + (m * naoccB);
                    int nm = m + (n * naoccB);
                    for(int f = 0 ; f < navirB; ++f){
                        int bf = f + (b * navirB);
                        double value2 = 0.5 * ( T1->get(mn, bf) - T1->get(nm, bf) );
                        Va->set(f, mn2, value2);
                    }
                }
            }

            // Form A[b](f,a>=e) = 1/2 \sum(m>=n) (+)Tt(m>=n,a>=e) * T[b](f,m>=n)
	    A->contract(false, false, navirB, ntri_abBB, ntri_ijBB, Va, Ta, 0.5, 0.0);

	    // Form V[b](a,ef)
            #pragma omp parallel for
            for(int a = 0 ; a < navirB; ++a){
                for(int e = 0 ; e < navirB; ++e){
                    int ae = index2(a,e);
                    for(int f = 0 ; f < navirB; ++f){
                        int ef = ab_idxBB->get(e,f);
			int perm1 = ( a > e ) ? 1 : -1;
			double value = perm1 * A->get(f,ae);
                        V->set(a, ef, value);
                    }
                }
            }

	    // G2[b](Q,a) = \sum(ef) b(Q,ef) * V[b](a,ef)
            X->gemm(false, true, bQabB, V, 1.0, 0.0);

	    // Form G2
            #pragma omp parallel for
            for(int Q = 0 ; Q < nQ; ++Q){
                for(int a = 0 ; a < navirB; ++a){
                    int ab = ab_idxBB->get(a,b);
		    G->add(Q, ab, X->get(Q,a));
                }
            }

    }
    T1.reset();
    Ta.reset();
    Va.reset();
    A.reset();
    X.reset();
    V.reset();
    bQabB.reset();

    // G_ab^Q += \sum(ef) V_EaFb b_EF^Q
    // V_EaFb = \sum(Mn) L_Mn^Ea T_Mn^Fb
    // Read b_EF^Q
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);

    T1 = SharedTensor2d(new Tensor2d("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    T1->read(psio_, PSIF_DFOCC_AMPS);
    Ts = SharedTensor2d(new Tensor2d("T[b] (Mn,F)", naoccA*naoccB, navirA));
    V = SharedTensor2d(new Tensor2d("V[b] (Ea,F)", navirA*navirB, navirA));
    Vs = SharedTensor2d(new Tensor2d("V[b] (EF,a)", navirA*navirA, navirB));
    X = SharedTensor2d(new Tensor2d("TPDM[b] (Q|a)", nQ, navirB));
    // Main loop
    for(int b = 0 ; b < navirB; ++b){

            // Form (+)T[b](Mn,F)
            #pragma omp parallel for
            for(int m = 0 ; m < naoccA; ++m){
                for(int n = 0 ; n < naoccB; ++n){
                    int mn = n + (m * naoccB);
                    for(int f = 0 ; f < navirA; ++f){
                        int fb = (f * navirB) + b;
                        Ts->set(mn, f, T1->get(mn,fb));
                    }
                }
            }

	    // Form V[b](Ea,F) = \sum(Mn) L(Mn,Ea) T[b](Mn,F)
	    V->gemm(true, false, T1, Ts, 1.0, 0.0);

            // Form Vs[b](EF,a)
            #pragma omp parallel for
            for(int e = 0 ; e < navirA; ++e){
                for(int f = 0 ; f < navirA; ++f){
                    int ef = ab_idxAA->get(e,f);
                    for(int a = 0 ; a < navirB; ++a){
                        int ea = (e * navirB) + a;
                        Vs->set(ef, a, V->get(ea,f));
                    }
                }
            }

	    // G2[b](Q,a) = \sum(EF) b(Q,EF) * Vs[b](EF,a)
            X->gemm(false, false, bQabA, Vs, 1.0, 0.0);

	    // Form G2
            #pragma omp parallel for
            for(int Q = 0 ; Q < nQ; ++Q){
                for(int a = 0 ; a < navirB; ++a){
                    int ab = ab_idxBB->get(a,b);
		    G->add(Q, ab, X->get(Q,a));
                }
            }

    }
    T1.reset();
    Ts.reset();
    X.reset();
    V.reset();
    Vs.reset();
    bQabA.reset();

    // Active to full
    //G->symmetrize3(G);
    //G->scale(2.0);
    G2 = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|vv)", nQ, nvirB, nvirB));
    G2->set3_act_vv(G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS, true, true);
    if(print_ > 3) G2->print();
    G2.reset();

    //outfile->Printf("\tvv block done.\n");

}// else if (reference_ == "UNRESTRICTED")
    timer_off("tpdm");
} // end olccd_tpdm

}} // End Namespaces
