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

void DFOCC::ccsd_t2_amps_low()
{

    // defs
    SharedTensor2d K, I, T, Tnew, U, Tau, W, X, Y;

    // Read old amplitudes
    T = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    T->read_symm(psio_, PSIF_DFOCC_AMPS);

    // t_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = \sum_{e} t_ij^ae Ft_be = \sum_{e} T(ia,je) Ft_be
    X = SharedTensor2d(new Tensor2d("X (IA|JB)", naoccA, navirA, naoccA, navirA));
    X->contract(false, true, naoccA * navirA * naoccA, navirA, navirA, T, FtabA, 1.0, 0.0);

    // t_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = -\sum_{m} t_mj^ab Ft_mi = -\sum_{m} Ft(m,i) T(ma,jb)
    X->contract(true, false, naoccA, naoccA * navirA * navirA, naoccA, FtijA, T, -1.0, 1.0);
    T.reset();

    // t_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = \sum_{Q} T'_ia^Q b_jb^Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA));
    K->read(psio_, PSIF_DFOCC_INTS);
    T = SharedTensor2d(new Tensor2d("T1p (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(true, false, T, K, 1.0, 1.0);
    K.reset();
    T.reset();

    // t_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = -\sum_{Q} t_ai^Q t_jb^Q
    U = SharedTensor2d(new Tensor2d("T1 (Q|AI)", nQ, navirA, naoccA));
    U->read(psio_, PSIF_DFOCC_AMPS);
    K = SharedTensor2d(new Tensor2d("Temp (Q|IA)", nQ, naoccA, navirA));
    K->swap_3index_col(U);
    U.reset();
    T = SharedTensor2d(new Tensor2d("T1 (Q|IA)", nQ, naoccA, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(true, false, K, T, -1.0, 1.0);
    T.reset();
    X->symmetrize();

    // t_ij^ab <= <ij|ab>
    Tnew = SharedTensor2d(new Tensor2d("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(Tnew);

    // Contributions of X
    Tnew->axpy(X, 2.0);
    X.reset();

    // Write and close
    Tnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    // WmnijT2
    ccsd_WmnijT2_low();

    // WmbejT2
    ccsd_WmbejT2_low();

    // WijamT2
    //if (itr_occ > 1) ccsd_WijamT2_low();

    // WabefT2
    //ccsd_WabefT2_low();
    ccsd_Wabef2T2_low();

    // Denom
    Tnew = SharedTensor2d(new Tensor2d("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->apply_denom_chem(nfrzc, noccA, FockA);

    // Reset T1
    rms_t1 = t1newA->rms(t1A);
    SharedTensor2d Rt1A = SharedTensor2d(new Tensor2d("RT1 <I|A>", naoccA, navirA));
    Rt1A->copy(t1newA);
    Rt1A->subtract(t1A);
    t1A->copy(t1newA);

    // Reset T2
    T = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    T->read_symm(psio_, PSIF_DFOCC_AMPS);
    rms_t2 = Tnew->rms(T);
    // Error vector
    Tau = SharedTensor2d(new Tensor2d("RT2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tau->copy(Tnew);
    Tau->subtract(T);
    T->copy(Tnew);
    Tnew.reset();
    if (do_diis_ == 0) T->write_symm(psio_, PSIF_DFOCC_AMPS);

    // DIIS
    std::shared_ptr<Matrix> RT2(new Matrix("RT2", naoccA*navirA, naoccA*navirA));
    Tau->to_matrix(RT2);
    Tau.reset();
    std::shared_ptr<Matrix> T2(new Matrix("T2", naoccA*navirA, naoccA*navirA));
    T->to_matrix(T2);
    T.reset();
    std::shared_ptr<Matrix> RT1(new Matrix("RT1", naoccA, navirA));
    Rt1A->to_matrix(RT1);
    Rt1A.reset();
    std::shared_ptr<Matrix> T1(new Matrix("T1", naoccA, navirA));
    t1A->to_matrix(T1);

    // add entry
    if (do_diis_ == 1) ccsdDiisManager->add_entry(4, RT2.get(), RT1.get(), T2.get(), T1.get());
    RT2.reset();
    RT1.reset();

    // extrapolate
    if (do_diis_ == 1) {
        if (ccsdDiisManager->subspace_size() >= cc_mindiis_) ccsdDiisManager->extrapolate(2, T2.get(), T1.get());
        T = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
        T->set2(T2);
        t1A->set2(T1);
        T->write_symm(psio_, PSIF_DFOCC_AMPS);
    }
    T1.reset();
    T2.reset();

    if (do_diis_ == 0) {
        T = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
        T->read_symm(psio_, PSIF_DFOCC_AMPS);
    }

    // Form T'(ib,ja) = T(ia,jb)
    U = SharedTensor2d(new Tensor2d("T2p (IA|JB)", naoccA, navirA, naoccA, navirA));
    U->sort(1432, T, 1.0, 0.0);
    U->write_symm(psio_, PSIF_DFOCC_AMPS);
    U.reset();

    // Form U(ia,jb) = 2*T(ia,jb) - T (ib,ja)
    U = SharedTensor2d(new Tensor2d("U2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    U->sort(1432, T, 1.0, 0.0);
    U->scale(-1.0);
    U->axpy(T, 2.0);
    U->write_symm(psio_, PSIF_DFOCC_AMPS);
    U.reset();

    // Form Tau(ia,jb) = T(ia,jb) + t(ia) * t(jb)
    Tau = SharedTensor2d(new Tensor2d("Tau (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tau->dirprd224(t1A, t1A);
    Tau->add(T);
    T.reset();
    Tau->write_symm(psio_, PSIF_DFOCC_AMPS);

    // Energy
    U = SharedTensor2d(new Tensor2d("2*Tau(ia,jb) - Tau(ib,ja)", naoccA, navirA, naoccA, navirA));
    U->sort(1432, Tau, 1.0, 0.0);
    U->scale(-1.0);
    U->axpy(Tau, 2.0);
    Tau.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(K);
    Ecorr = U->vector_dot(K);
    U.reset();
    K.reset();
    Eccsd = Escf + Ecorr;

}// end ccsd_t2_amps_low
}} // End Namespaces
