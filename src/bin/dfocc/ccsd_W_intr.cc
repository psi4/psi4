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
    SharedTensor2d K, L, I, M, T, T1, Tnew, U, Tau, W, Wrv, X, Y;

    timer_on("WijamT2");

    // W_ijam = W(ij,am)
    W = SharedTensor2d(new Tensor2d("W (IJ|AM)", naoccA, naoccA, navirA, naoccA));

    // W_ijam = \sum_{ef} Tau_ij^ef <am|ef> = \sum_{ef} X_ij^ef <ma|ef>
    Tau = SharedTensor2d(new Tensor2d("Tau (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tau->read_symm(psio_, PSIF_DFOCC_AMPS);
    U = SharedTensor2d(new Tensor2d("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA));
    U->sort(1342, Tau, 1.0, 0.0);
    Tau.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    L->read(psio_, PSIF_DFOCC_INTS, true, true);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (E|Q)", navirA, nQ));
    I = SharedTensor2d(new Tensor2d("I[M] <E|FA>", navirA, navirA * navirA));
    Wrv = SharedTensor2d(new Tensor2d("W[M] <IJ|A>", naoccA * naoccA, navirA));

    // Main loop
    for(int m = 0 ; m < naoccA; ++m){

            // Form b[m](e,Q)
            #pragma omp parallel for
            for(int e = 0 ; e < navirA; ++e){
                int me = ia_idxAA->get(m,e);
                for(int Q = 0 ; Q < nQ; ++Q){
                    M->set(e, Q, K->get(Q,me));
                }
            }

            // Form J[m](e,fa) = \sum_{Q} b[m](e,Q) b(Q,fa)
            I->gemm(false, false, M, L, 1.0, 0.0);

            // Form W[m](ij,a) = \sum_{ef} X(ij,ef) J[m](ef,a)
            Wrv->contract(false, false, naoccA * naoccA, navirA, navirA * navirA, U, I, 1.0, 0.0);

            // Form W(ij,am) = W[m](ij,a)
            #pragma omp parallel for
            for(int i = 0 ; i < naoccA; ++i){
                for(int j = 0 ; j < naoccA; ++j){
                    int ij = ij_idxAA->get(i,j);
                    for(int a = 0 ; a < navirA; ++a){
                        int am = ai_idxAA->get(a,m);
                        W->set(ij, am, Wrv->get(ij,a));
                    }
                }
            }
    }
    L.reset();
    I.reset();
    U.reset();
    K.reset();
    M.reset();
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
//    WabefT2 High Mem
//======================================================================             
void DFOCC::ccsd_WabefT2_high()
{
    // defs
    SharedTensor2d K, M, L, I, T, Tnew, U, Tau, W, X, Y, S, A;
    SharedTensor2d V, Vs, Ts, Va, Ta;

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

    // Read B(Q,ab)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|E)", nQ, navirA));

    // Symmetric & Anti-symmetric contributions
    S = SharedTensor2d(new Tensor2d("S (A>=B, I>=J)", ntri_abAA, ntri_ijAA));
    A = SharedTensor2d(new Tensor2d("A (A>=B, I>=J)", ntri_abAA, ntri_ijAA));
    // Main loop
    for(int a = 0 ; a < navirA; ++a){
            int nb = a+1;
            I = SharedTensor2d(new Tensor2d("I[A] <E|FB>", navirA, navirA * nb));
            L = SharedTensor2d(new Tensor2d("B (Q|FB)", nQ, navirA * nb));

            // Form b[a](Qe) = b_ea^Q, cost = V^2N
            M->pcopy(K,1,navirA-1,a);

            // Form V[a](e,fb) = \sum_{Q} b[a](e,Q) b(Q,fb), cost = V^4N/2
            L->pcopy(K,nb,navirA-nb);
            I->gemm(true, false, M, L, 1.0, 0.0);
            L.reset();

            // memalloc
            Vs = SharedTensor2d(new Tensor2d("(+)V[A] (B, E>=F)", nb, ntri_abAA));
            Va = SharedTensor2d(new Tensor2d("(-)V[A] (B, E>=F)", nb, ntri_abAA));
            Ts = SharedTensor2d(new Tensor2d("(+)T[A] (B, I>=J)", nb, ntri_ijAA));
            Ta = SharedTensor2d(new Tensor2d("(-)T[B] (B, I>=J)", nb, ntri_ijAA));

            // Form (+)V[a](b, e>=f), cost = V^4/4
            #pragma omp parallel for
            for(int b = 0 ; b <= a; ++b){
                for(int e = 0 ; e < navirA; ++e){
                    int eb = b + (e * nb);
                    for(int f = 0 ; f <= e; ++f){
                        int ef = index2(e,f); 
                        int fb = b + (f * nb);
                        double value1 = 0.5 * ( I->get(e, fb) + I->get(f, eb) );
                        double value2 = 0.5 * ( I->get(e, fb) - I->get(f, eb) );
                        Vs->set(b, ef, value1);
                        Va->set(b, ef, value2);
                    }
                }
            }
            I.reset();

            // Form T[a](b, i>=j) = \sum_{e>=f} Tau(i>=j,e>=f) V[a](b, e>=f), cost = O^2V^4/2 
            Ts->gemm(false, true, Vs, U, 1.0, 0.0);
            Ta->gemm(false, true, Va, T, 1.0, 0.0);
            Vs.reset();
            Va.reset();

            // Form S(ij,ab) & A(ij,ab), cost = O^2V^2/2
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
            Ts.reset();
            Ta.reset();

    }
    K.reset();
    I.reset();
    U.reset();
    T.reset();
    M.reset();

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

}// end ccsd_WabefT2_high

//======================================================================
//    WabefT2 Low Mem
//======================================================================             
void DFOCC::ccsd_WabefT2_low()
{
    // defs
    SharedTensor2d K, M, I, T, Tnew, U, Tau, W, X, Y, S, A;
    SharedTensor2d V, Vs, Ts, Va, Ta;

    timer_on("WabefT2");

    // t_ij^ab <= \sum_{ef} Tau_ij^ef <ab|ef>
    Tau = SharedTensor2d(new Tensor2d("Tau (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tau->read_symm(psio_, PSIF_DFOCC_AMPS);
    U = SharedTensor2d(new Tensor2d("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA));
    U->sort(1324, Tau, 1.0, 0.0);
    Tau.reset();

    // Read B(Q,ab)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (E|Q)", navirA, nQ));
    I = SharedTensor2d(new Tensor2d("I[A] <E|FB>", navirA, navirA * navirA));
    T = SharedTensor2d(new Tensor2d("T[A] <IJ|B>", naoccA * naoccA, navirA));
    X = SharedTensor2d(new Tensor2d("X (IA|JB)", naoccA, navirA, naoccA, navirA));

    // Main loop
    for(int a = 0 ; a < navirA; ++a){

            // Form b[a](e,Q), cost = V^2N
            #pragma omp parallel for
            for(int e = 0 ; e < navirA; ++e){
                int ae = ab_idxAA->get(a,e);
                for(int Q = 0 ; Q < nQ; ++Q){
                    M->set(e, Q, K->get(Q,ae));
                }
            }

            // Form V[a](e,fb) = \sum_{Q} b[a](e,Q) b(Q,fb), cost = V^4N
            I->gemm(false, false, M, K, 1.0, 0.0);

            // T[a](ij,b) = \sum_{ef} Tau(ij,ef) V[a](ef,b)
            T->contract(false, false, naoccA * naoccA, navirA, navirA * navirA, U, I, 1.0, 0.0);

            // Form X(ij,ab) = T[a](ij,b)
            #pragma omp parallel for
            for(int i = 0 ; i < naoccA; ++i){
                int ia = ia_idxAA->get(i,a);
                for(int j = 0 ; j < naoccA; ++j){
                    int ij = ij_idxAA->get(i,j);
                    for(int b = 0 ; b < navirA; ++b){
                        int jb = ia_idxAA->get(j,b);
                        X->set(ia, jb, T->get(ij,b));
                    }
                }
            }

    }
    K.reset();
    I.reset();
    U.reset();
    T.reset();
    M.reset();

    // T(ia,jb) += X(ia,jb)
    Tnew = SharedTensor2d(new Tensor2d("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->add(X);
    X.reset();
    Tnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    timer_off("WabefT2");

}// end ccsd_WabefT2_low

}} // End Namespaces



