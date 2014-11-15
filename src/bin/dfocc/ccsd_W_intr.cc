/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#include <libqt/qt.h>
#include "defines.h"
#include "dfocc.h"

using namespace psi;
using namespace std;


namespace psi{ namespace dfoccwave{
  
void DFOCC::ccsd_WmnijT2()
{
    // defs
    SharedTensor2d K, T, Tnew, U, Tau, W, X;

    timer_on("WmnijT2");

    // W_mnij = <mn|ij>
    W = SharedTensor2d(new Tensor2d("W <MN|IJ>", naoccA, naoccA, naoccA, naoccA));
    tei_ijkl_phys_directAA(W);

    // W_mnij += X(im,jn) + X(jn,im) += 2Xt(im,jn) 
    // X_imjn = \sum_{Q} t_im^Q b_jn^Q
    X = SharedTensor2d(new Tensor2d("X <MN|IJ>", naoccA, naoccA, naoccA, naoccA));
    T = SharedTensor2d(new Tensor2d("T1 (Q|IJ)", nQ, naoccA, naoccA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    X->gemm(true, false, T, K, 1.0, 0.0);
    T.reset();
    K.reset();
    X->symmetrize();
    W->sort(2413, X, 2.0, 1.0);
    X.reset();

    // W_mnij = \sum_{ef} Tau_ij^ef <mn|ef>
    Tau = SharedTensor2d(new Tensor2d("Tau (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tau->read_symm(psio_, PSIF_DFOCC_AMPS);
    U = SharedTensor2d(new Tensor2d("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA));
    U->sort(1324, Tau, 1.0, 0.0);
    Tau.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA));
    tei_ijab_phys_directAA(K);
    W->gemm(false, true, K, U, 1.0, 1.0);
    K.reset();

    // t_ij^ab <= \sum_{m,n} Tau_mn^ab Wmnij
    T = SharedTensor2d(new Tensor2d("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    T->gemm(true, false, W, U, 1.0, 0.0);
    U.reset();
    W.reset();
    Tnew = SharedTensor2d(new Tensor2d("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->sort(1324, T, 1.0, 1.0);
    T.reset();
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
    SharedTensor2d K, T, T1, Tnew, U, Tau, W, X, Y;

    timer_on("WmbejT2");

    // W_mbej = W(me,jb)
    // W(me,jb) <= (me|jb)
    W = SharedTensor2d(new Tensor2d("W (ME|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(W);

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
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA));
    K->read(psio_, PSIF_DFOCC_INTS);
    W->gemm(true, false, K, U, 1.0, 1.0);
    K.reset();
    U.reset();

    // W(me,jb) <= -1/2 \sum_{nf} t_jn^bf X(me,nf)
    // (mf|ne) = X(me,nf) (sort: 1432)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(K);
    X = SharedTensor2d(new Tensor2d("X (IA|JB)", naoccA, navirA, naoccA, navirA));
    X->sort(1432, K, 1.0, 0.0);
    K.reset();
    T = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    T->read_symm(psio_, PSIF_DFOCC_AMPS);
    W->gemm(false, false, X, T, -0.5, 1.0);
    T.reset();
    X.reset();

    // t_ij^ab <= X(ia,jb) + X(jb,ia)
    // X(ia,jb) = \sum_{me} u_im^ae W(me,jb) = \sum_{me} U(ia,me) W(me,jb)
    U = SharedTensor2d(new Tensor2d("U2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    X = SharedTensor2d(new Tensor2d("X (IA|JB)", naoccA, navirA, naoccA, navirA));
    U->read_symm(psio_, PSIF_DFOCC_AMPS);
    X->gemm(false, false, U, W, 1.0, 0.0);
    U.reset();
    W.reset();
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
//    WmbjeT2
//======================================================================             
void DFOCC::ccsd_WmbjeT2()
{
    // defs
    SharedTensor2d K, T, T1, Tnew, U, Tau, W, X, Y;

    timer_on("WmbjeT2");

    // W_mbje = W'(me,jb)
    // W'(me,jb) <= <me|jb>
    W = SharedTensor2d(new Tensor2d("Wp (ME|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_phys_directAA(W);

    // X(jm,be) <= -\sum_{Q} t_be^Q ( t_jm^Q + b_jm^Q )
    T = SharedTensor2d(new Tensor2d("T1 (Q|IJ)", nQ, naoccA, naoccA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA));
    K->read(psio_, PSIF_DFOCC_INTS);
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
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    X->gemm(true, false, T, K, 1.0, 1.0);
    // W'(me,jb) <= X(jm,be)
    W->sort(2413, X, 1.0, 1.0);
    X.reset();

    // W'(me,jb) <= -1/2 \sum_{nf} t_nj^bf X(me,nf)
    // (mf|ne) = X(me,nf) (sort: 1432)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(K);
    X = SharedTensor2d(new Tensor2d("X (IA|JB)", naoccA, navirA, naoccA, navirA));
    X->sort(1432, K, 1.0, 0.0);
    K.reset();
    T = SharedTensor2d(new Tensor2d("T2p (IA|JB)", naoccA, navirA, naoccA, navirA));
    T->read_symm(psio_, PSIF_DFOCC_AMPS);
    W->gemm(false, false, X, T, -0.5, 1.0);
    X.reset();

    // t_ij^ab <= Y(ib,ja) + Y(ja,ib)
    // Y(ib,ja) = -\sum_{me} t_mi^be W'(me,ja) = -\sum_{me} T'(ib,me) W'(me,ja)
    Y = SharedTensor2d(new Tensor2d("Y (IA|JB)", naoccA, navirA, naoccA, navirA));
    Y->gemm(false, false, T, W, -1.0, 0.0);
    T.reset();
    X = SharedTensor2d(new Tensor2d("X (IA|JB)", naoccA, navirA, naoccA, navirA));
    X->sort(1432, Y, 1.0, 0.0);
    Y.reset();
    // t_ij^ab <= X(ia,jb) + X(jb,ia)
    // X(ia,jb) = -\sum_{me} t_im^ae W'(me,jb) = -\sum_{me} T(ia,me) W'(me,jb)
    U = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    U->read_symm(psio_, PSIF_DFOCC_AMPS);
    X->gemm(false, false, U, W, -1.0, 1.0);
    U.reset();
    W.reset();
    X->symmetrize();
    Tnew = SharedTensor2d(new Tensor2d("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->axpy(X, 2.0);
    X.reset();
    Tnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();
    
    timer_off("WmbjeT2");

}// end ccsd_WmjeT2

//======================================================================
//    WijamT2
//======================================================================             
void DFOCC::ccsd_WijamT2()
{
    // defs
    SharedTensor2d K, L, I, T, T1, Tnew, U, Tau, W, Wrv, X, Y;

    timer_on("WijamT2");

    // W_ijam = W(ij,am)
    W = SharedTensor2d(new Tensor2d("W (IJ|AM)", naoccA, naoccA, navirA, naoccA));

    // W_ijam = \sum_{ef} Tau_ij^ef <am|ef>
    Tau = SharedTensor2d(new Tensor2d("Tau (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tau->read_symm(psio_, PSIF_DFOCC_AMPS);
    U = SharedTensor2d(new Tensor2d("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA));
    U->sort(1324, Tau, 1.0, 0.0);
    Tau.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA));
    K->read(psio_, PSIF_DFOCC_INTS);
    //L = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    //L->read(psio_, PSIF_DFOCC_INTS, true, true);
    // Read B(Q, a>=b)
    int ntri_vv = 0.5 * navirA * (navirA + 1);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, ntri_vv));
    L->read(psio_, PSIF_DFOCC_INTS);
    I = SharedTensor2d(new Tensor2d("I[AM] <E|F>", navirA, navirA));
    Wrv = SharedTensor2d(new Tensor2d("W[AM] <I|J>", naoccA, naoccA));

    // Main loop
    for(int a = 0 ; a < navirA; ++a){
        for(int m = 0 ; m < naoccA; ++m){
            int am = ai_idxAA->get(a,m);
            // Form I[am]
            #pragma omp parallel for
            for(int e = 0 ; e < navirA; ++e){
                int ae = index2(a,e); 
                for(int f = 0 ; f < navirA; ++f){
                    int mf = ia_idxAA->get(m,f);
                    double value = 0.0;
                    for(int Q = 0 ; Q < nQ; ++Q){
                        value += L->get(Q,ae) * K->get(Q,mf);
                    }
                    I->set(e,f,value);
                }
            }
            // Form W[am]
            Wrv->gemv(false, U, I, 1.0, 0.0);
            // Form W(ij,am)
            W->set_column(Wrv,am);
        }
    }
    L.reset();
    U.reset();
    K.reset();
    I.reset();
    Wrv.reset();

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
//    WabefT2
//======================================================================             
void DFOCC::ccsd_WabefT2()
{
    // defs
    SharedTensor2d K, I, T, Tnew, U, Tau, W, X, Y, S, A;
    SharedTensor1d Vs, Ts, Va, Ta;

    timer_on("WabefT2");

    // t_ij^ab <= \sum_{ef} Tau_ij^ef <ab|ef>
    Tau = SharedTensor2d(new Tensor2d("Tau (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tau->read_symm(psio_, PSIF_DFOCC_AMPS);
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

    // Read B(Q, a>=b)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, ntri_abAA));
    K->read(psio_, PSIF_DFOCC_INTS);
    Vs = SharedTensor1d(new Tensor1d("(+)V[AB] (E>=F)", ntri_abAA));
    Ts = SharedTensor1d(new Tensor1d("(+)T[AB] (I>=J)", ntri_ijAA));
    Va = SharedTensor1d(new Tensor1d("(-)V[AB] (E>=F)", ntri_abAA));
    Ta = SharedTensor1d(new Tensor1d("(-)T[AB] (I>=J)", ntri_ijAA));

    // Symmetric & Anti-symmetric contributions
    S = SharedTensor2d(new Tensor2d("S [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    A = SharedTensor2d(new Tensor2d("A [I>=J|A>=B]", ntri_ijAA, ntri_abAA));
    for(int a = 0 ; a < navirA; ++a){
        for(int b = 0 ; b <=a; ++b){
            int ab = index2(a,b); 

            // Form V[ab]
            #pragma omp parallel for
            for(int e = 0 ; e < navirA; ++e){
                int ae = index2(a,e); 
                int be = index2(b,e); 
                for(int f = 0 ; f <= e; ++f){
                    int ef = index2(e,f); 
                    int af = index2(a,f); 
                    int bf = index2(b,f); 
                    double value1 = 0.0;
                    double value2 = 0.0;
                    for(int Q = 0 ; Q < nQ; ++Q){
                        value1 += 0.5 * K->get(Q,ae) * K->get(Q,bf);
                        value2 += 0.5 * K->get(Q,af) * K->get(Q,be);
                    }
                    Vs->set(ef,value1+value2);
                    Va->set(ef,value1-value2);
                }
            }

            // Form T[ab]
            Ts->gemv(false, U, Vs, 1.0, 0.0);
            Ta->gemv(false, T, Va, 1.0, 0.0);

            // Form S(ij,ab) & A(ij,ab)
            #pragma omp parallel for
            for(int i = 0 ; i < naoccA; ++i){
                for(int j = 0 ; j <= i; ++j){
                    int ij = index2(i,j); 
                    S->add(ij,ab,Ts->get(ij));
                    A->add(ij,ab,Ta->get(ij));
                }
            }
        }
    }
    K.reset();
    U.reset();
    T.reset();
    Vs.reset();
    Ts.reset();
    Va.reset();
    Ta.reset();

    // T(ia,jb) <-- S(i>=j,a>=b) + A(i>=j,a>=b)
    Tnew = SharedTensor2d(new Tensor2d("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    #pragma omp parallel for
    for(int i = 0 ; i < naoccA; ++i){
        for(int j = 0 ; j < naoccA; ++j){
            int ij = index2(i,j); 
            for(int a = 0 ; a < navirA; ++a){
                int ia = ia_idxAA->get(i,a);
                for(int b = 0 ; b < navirA; ++b){
                    int perm = 0;
                    int ab = index2(a,b); 
                    int jb = ia_idxAA->get(j,b);
                    if (i >= j && a >= b) perm = 1;
                    else if (j >= i && a >= b) perm = -1;
                    else if (i >= j && b >= a) perm = -1;
                    else if (j >= i && b >= a) perm = 1;
                    double value = S->get(ij,ab) + (perm * A->get(ij,ab));
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

}// end ccsd_WabefT2


}} // End Namespaces



