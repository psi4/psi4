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
  
void DFOCC::ccd_t2_amps_low()
{

    // defs
    SharedTensor2d K, I, T, Tnew, U, Tau, W, X, Y;

    // Read old amplitudes
    T = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    T->read_symm(psio_, PSIF_DFOCC_AMPS);

    // t_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = \sum_{e} t_ij^ae F_be = \sum_{e} T(ia,je) F_be
    X = SharedTensor2d(new Tensor2d("X (IA|JB)", naoccA, navirA, naoccA, navirA));
    X->contract(false, true, naoccA * navirA * naoccA, navirA, navirA, T, FabA, 1.0, 0.0);

    // t_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = -\sum_{m} t_mj^ab F_mi = -\sum_{m} F(m,i) T(ma,jb)
    X->contract(true, false, naoccA, naoccA * navirA * navirA, naoccA, FijA, T, -1.0, 1.0);
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
    ccd_WmnijT2_low();

    // WmbejT2
    ccd_WmbejT2_low();

    // WabefT2
    ccd_WabefT2_low();

    // Denom
    Tnew = SharedTensor2d(new Tensor2d("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->apply_denom_chem(nfrzc, noccA, FockA);

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
    boost::shared_ptr<Matrix> RT2(new Matrix("RT2", naoccA*navirA, naoccA*navirA));
    Tau->to_matrix(RT2);
    Tau.reset();
    boost::shared_ptr<Matrix> T2(new Matrix("T2", naoccA*navirA, naoccA*navirA));
    T->to_matrix(T2);
    T.reset();

    // add entry
    if (do_diis_ == 1) ccsdDiisManager->add_entry(2, RT2.get(), T2.get());
    RT2.reset();

    // extrapolate
    if (do_diis_ == 1) {
        if (ccsdDiisManager->subspace_size() >= cc_mindiis_) ccsdDiisManager->extrapolate(1, T2.get());
        T = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
        T->set2(T2);
        T->write_symm(psio_, PSIF_DFOCC_AMPS);
    }
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
    T.reset();
    U->write_symm(psio_, PSIF_DFOCC_AMPS);

    // Energy
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(K);
    Ecorr = U->vector_dot(K);
    U.reset();
    K.reset();
    Eccd = Escf + Ecorr;

}// end ccd_t2_amps_low
}} // End Namespaces
