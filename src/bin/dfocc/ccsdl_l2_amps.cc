/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

#include <libqt/qt.h>
#include "defines.h"
#include "dfocc.h"

using namespace psi;
using namespace std;


namespace psi{ namespace dfoccwave{
  
void DFOCC::ccsdl_l2_amps()
{

    // defs
    SharedTensor2d K, I, M, L, Lnew, T, U, Tau, W, X, Y, Z;

    // l_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = l_i^a Ft_jb
    X = SharedTensor2d(new Tensor2d("X (IA|JB)", naoccA, navirA, naoccA, navirA));
    X->dirprd224(l1A, FiaA);

    // l_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = \sum_{e} l_ij^ae Ft_eb = \sum_{e} L(ia,je) Ft_eb
    X->contract(false, false, naoccA * navirA * naoccA, navirA, navirA, l2, FtabA, 1.0, 1.0);

    // l_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = -\sum_{m} l_mj^ab Ft_im = -\sum_{m} Ft(i,m) L(ma,jb)
    X->contract(false, false, naoccA, naoccA * navirA * navirA, naoccA, FtijA, l2, -1.0, 1.0);

    // l_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = -\sum_{m} l_m^a W_ijmb 
    W = SharedTensor2d(new Tensor2d("WL (MN|IE)", naoccA, naoccA, naoccA, navirA));
    //W->read(psio_, PSIF_DFOCC_AMPS);
    ccsdl_Wmnie_direct(W);
    Y = SharedTensor2d(new Tensor2d("Y (MN|EI)", naoccA, naoccA, navirA, naoccA));
    Y->sort(1243, W, 1.0, 0.0);
    W.reset();
    T = SharedTensor2d(new Tensor2d("Temp <IJ|BA>", naoccA, naoccA, navirA, navirA));
    T->contract(false, false, naoccA*naoccA*navirA, navirA, naoccA, Y, l1A, -1.0, 0.0);
    Y.reset();
    X->sort(1423, T, 1.0, 1.0);
    T.reset();

    // l_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = \sum_{Q} (G_ai^Q - G_ia^Q + l_ia^Q - lt_ia^Q) b_jb^Q
    U = SharedTensor2d(new Tensor2d("G (Q|AI)", nQ, navirA, naoccA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T = SharedTensor2d(new Tensor2d("Temp (Q|IA)", nQ, naoccA, navirA));
    T->swap_3index_col(U);
    U.reset();
    U = SharedTensor2d(new Tensor2d("G (Q|IA)", nQ, naoccA, navirA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, -1.0);
    U.reset();
    U = SharedTensor2d(new Tensor2d("L1 (Q|IA)", nQ, naoccA, navirA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();
    U = SharedTensor2d(new Tensor2d("L1t (Q|IA)", nQ, naoccA, navirA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, -1.0);
    U.reset();
    X->gemm(true, false, T, bQiaA, 1.0, 1.0);
    T.reset();
    X->symmetrize();

    // l_ij^ab <= <ij|ab>
    Lnew = SharedTensor2d(new Tensor2d("New L2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Lnew->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0); 

    // Contributions of X
    Lnew->axpy(X, 2.0);
    X.reset();

    // Write and close
    Lnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Lnew.reset();

    // VmnijL2
    ccsdl_VmnijL2();

    // WijmnL2
    ccsdl_WijmnL2();

    // WmbejTL2
    ccsdl_WmbejL2();

    // WabefL2
    if (Wabef_type_ == "AUTO") {
	if (!do_ppl_hm) ccsdl_WabefL2();
	else {
	    ccsdl_LijmeL2_high_mem();
	    ccsdl_WabefL2_high_mem();
	}
    }
    else if (Wabef_type_ == "LOW_MEM") ccsdl_WabefL2();
    else if (Wabef_type_ == "HIGH_MEM") {
	ccsdl_LijmeL2_high_mem();
	ccsdl_WabefL2_high_mem();
    }

    // Denom
    Lnew = SharedTensor2d(new Tensor2d("New L2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Lnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Lnew->apply_denom_chem(nfrzc, noccA, FockA);

    // Reset T1
    rms_t1 = l1newA->rms(l1A); 
    SharedTensor2d Rl1A = SharedTensor2d(new Tensor2d("RL1 <I|A>", naoccA, navirA));
    Rl1A->copy(l1newA);
    Rl1A->subtract(l1A);
    l1A->copy(l1newA);

    // Reset T2
    rms_t2 = Lnew->rms(l2); 
    // Error vector
    Tau = SharedTensor2d(new Tensor2d("RL2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tau->copy(Lnew);
    Tau->subtract(l2);
    l2->copy(Lnew);
    Lnew.reset();

    // DIIS
    boost::shared_ptr<Matrix> RL2(new Matrix("RL2", naoccA*navirA, naoccA*navirA));
    Tau->to_matrix(RL2);
    Tau.reset();
    boost::shared_ptr<Matrix> L2(new Matrix("L2", naoccA*navirA, naoccA*navirA));
    l2->to_matrix(L2);
    boost::shared_ptr<Matrix> RL1(new Matrix("RL1", naoccA, navirA));
    Rl1A->to_matrix(RL1);
    Rl1A.reset();
    boost::shared_ptr<Matrix> L1(new Matrix("L1", naoccA, navirA));
    l1A->to_matrix(L1);

    // add entry
    if (do_diis_ == 1) ccsdlDiisManager->add_entry(4, RL2.get(), RL1.get(), L2.get(), L1.get());
    RL2.reset();
    RL1.reset();

    // extrapolate
    if (do_diis_ == 1) {
        if (ccsdlDiisManager->subspace_size() >= cc_mindiis_) ccsdlDiisManager->extrapolate(2, L2.get(), L1.get());
        l2->set2(L2);
        l1A->set2(L1);
    }
    L1.reset();
    L2.reset();

    // Form Tau(ia,jb) = L(ia,jb) + l(ia) * l(jb)
    Tau = SharedTensor2d(new Tensor2d("Tau (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsdl_tau_amps(Tau,l2);
 
    // Energy
    U = SharedTensor2d(new Tensor2d("2*Tau(ia,jb) - Tau(ib,ja)", naoccA, navirA, naoccA, navirA));
    U->sort(1432, Tau, 1.0, 0.0);
    U->scale(-1.0);
    U->axpy(Tau, 2.0);
    Tau.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    EcorrL = U->vector_dot(K);
    U.reset();
    K.reset();
    EccsdL = Escf + EcorrL;

    // print
    //l2->print();

}// end ccsdl_l2_amps

//======================================================================
//    CCSD: Tau_ij^ab = l_ij^ab + l_i^a l_j^b
//======================================================================             
void DFOCC::ccsdl_tau_amps(SharedTensor2d &U, SharedTensor2d &T)
{
    U->dirprd224(l1A, l1A);
    U->add(T);
}// end ccsdl_tau_amps

}} // End Namespaces
