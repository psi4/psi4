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
  
void DFOCC::ccd_t2_amps()
{

    // defs
    SharedTensor2d K, I, T, Tnew, U, Tau, W, X, Y;

    // t_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = \sum_{e} t_ij^ae F_be = \sum_{e} T(ia,je) F_be
    X = SharedTensor2d(new Tensor2d("X (IA|JB)", naoccA, navirA, naoccA, navirA));
    X->contract(false, true, naoccA * navirA * naoccA, navirA, navirA, t2, FabA, 1.0, 0.0);

    // t_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = -\sum_{m} t_mj^ab F_mi = -\sum_{m} F(m,i) T(ma,jb)
    X->contract(true, false, naoccA, naoccA * navirA * navirA, naoccA, FijA, t2, -1.0, 1.0);
    X->symmetrize();

    // t_ij^ab <= <ij|ab>
    Tnew = SharedTensor2d(new Tensor2d("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tnew->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0); 

    // Contributions of X
    Tnew->axpy(X, 2.0);
    X.reset();

    // Write and close
    Tnew->write_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    // WmnijT2
    ccd_WmnijT2();

    // WmbejT2
    ccd_WmbejT2();

    // WabefT2
    ccd_WabefT2();

    // Denom
    Tnew = SharedTensor2d(new Tensor2d("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->apply_denom_chem(nfrzc, noccA, FockA);

    // Reset T2
    rms_t2 = Tnew->rms(t2); 

    // Error vector
    Tau = SharedTensor2d(new Tensor2d("RT2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tau->copy(Tnew);
    Tau->subtract(t2);
    t2->copy(Tnew);
    Tnew.reset();

    // DIIS
    boost::shared_ptr<Matrix> RT2(new Matrix("RT2", naoccA*navirA, naoccA*navirA));
    Tau->to_matrix(RT2);
    Tau.reset();
    boost::shared_ptr<Matrix> T2(new Matrix("T2", naoccA*navirA, naoccA*navirA));
    t2->to_matrix(T2);

    // add entry
    if (do_diis_ == 1) ccsdDiisManager->add_entry(2, RT2.get(), T2.get());
    RT2.reset();

    // extrapolate
    if (do_diis_ == 1) {
        if (ccsdDiisManager->subspace_size() >= cc_mindiis_) ccsdDiisManager->extrapolate(1, T2.get());
        t2->set2(T2);
    }
    T2.reset();

    // Form U(ia,jb) = 2*T(ia,jb) - T (ib,ja)
    U = SharedTensor2d(new Tensor2d("U2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_u2_amps(U,t2);

    // Energy
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    Ecorr = U->vector_dot(K);
    U.reset();
    K.reset();
    Eccd = Escf + Ecorr;

}// end ccd_t2_amps

}} // End Namespaces

