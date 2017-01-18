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
#include "psi4/libmints/matrix.h"
#include "psi4/libdiis/diismanager.h"
using namespace std;


namespace psi{ namespace dfoccwave{

void DFOCC::ccdl_l2_amps()
{

    // defs
    SharedTensor2d K, I, L, Lnew, T, U, Tau, W, X, Y, Z;


    // l_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = \sum_{e} l_ij^ae F_eb = \sum_{e} L(ia,je) F_eb
    X = SharedTensor2d(new Tensor2d("X (IA|JB)", naoccA, navirA, naoccA, navirA));
    X->contract(false, false, naoccA * navirA * naoccA, navirA, navirA, l2, FabA, 1.0, 0.0);

    // l_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = -\sum_{m} l_mj^ab F_im = -\sum_{m} F(i,m) L(ma,jb)
    X->contract(false, false, naoccA, naoccA * navirA * navirA, naoccA, FijA, l2, -1.0, 1.0);

    // l_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = \sum_{Q} (G_ai^Q - G_ia^Q) b_jb^Q
    U = SharedTensor2d(new Tensor2d("G (Q|AI)", nQ, navirA, naoccA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    T = SharedTensor2d(new Tensor2d("Temp (Q|IA)", nQ, naoccA, navirA));
    T->swap_3index_col(U);
    U.reset();
    U = SharedTensor2d(new Tensor2d("G (Q|IA)", nQ, naoccA, navirA));
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
    ccdl_VmnijL2();

    // WijmnL2
    ccdl_WijmnL2();

    // WmbejTL2
    ccdl_WmbejL2();

    // WabefL2
    if (Wabef_type_ == "AUTO") {
	if (!do_ppl_hm) ccdl_WabefL2();
	else ccsdl_WabefL2_high_mem();
    }
    else if (Wabef_type_ == "LOW_MEM") ccdl_WabefL2();
    else if (Wabef_type_ == "HIGH_MEM") ccsdl_WabefL2_high_mem();

    // Denom
    Lnew = SharedTensor2d(new Tensor2d("New L2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Lnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Lnew->apply_denom_chem(nfrzc, noccA, FockA);

    // Reset T2
    rms_t2 = Lnew->rms(l2);
    // Error vector
    Tau = SharedTensor2d(new Tensor2d("RL2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tau->copy(Lnew);
    Tau->subtract(l2);
    l2->copy(Lnew);
    Lnew.reset();

    // DIIS
    std::shared_ptr<Matrix> RL2(new Matrix("RL2", naoccA*navirA, naoccA*navirA));
    Tau->to_matrix(RL2);
    Tau.reset();
    std::shared_ptr<Matrix> L2(new Matrix("L2", naoccA*navirA, naoccA*navirA));
    l2->to_matrix(L2);

    // add entry
    if (do_diis_ == 1) ccsdlDiisManager->add_entry(2, RL2.get(), L2.get());
    RL2.reset();

    // extrapolate
    if (do_diis_ == 1) {
        if (ccsdlDiisManager->subspace_size() >= cc_mindiis_) ccsdlDiisManager->extrapolate(1, L2.get());
        l2->set2(L2);
    }
    L2.reset();

    // Energy
    U = SharedTensor2d(new Tensor2d("2*L(ia,jb) - L(ib,ja)", naoccA, navirA, naoccA, navirA));
    U->sort(1432, l2, 1.0, 0.0);
    U->scale(-1.0);
    U->axpy(l2, 2.0);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    EcorrL = U->vector_dot(K);
    U.reset();
    K.reset();
    EccdL = Escf + EcorrL;

    // print
    //l2->print();

}// end ccdl_l2_amps

}} // End Namespaces
