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
  
void DFOCC::ccsd_t2_amps()
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
    ccsd_WmnijT2();

    // WmbejT2
    ccsd_WmbejT2();

    // WmbjeT2
    ccsd_WmbjeT2();

    // WijamT2
    if (itr_occ > 1) ccsd_WijamT2();

    // WabefT2
    ccsd_WabefT2();

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
    T->copy(Tnew);
    Tnew.reset();
    T->write_symm(psio_, PSIF_DFOCC_AMPS);

    /*
    // DIIS
    Tau = SharedTensor2d(new Tensor2d("RT2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tau->copy(Tnew);
    Tau->subtract(T);
    double** RT2 = Tau->to_block_matrix();
    Tau.reset();
    double** T2 = T->to_block_matrix();
    T.reset();
    double** RT1 = Rt1A->to_block_matrix();
    Rt1A.reset();
    double** T1 = t1A->to_block_matrix();
    t2DiisManager->add_entry(4, RT2[0], RT1[0], T2[0], T1[0]);
    free_block(RT2);
    free_block(RT1);
    if (t2DiisManager->subspace_size() >= cc_mindiis_) t2DiisManager->extrapolate(2, T2[0], T1[0]);
    T = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    T->set(T2);
    free_block(T2);
    T->write_symm(psio_, PSIF_DFOCC_AMPS);
    t1A->set(T1);
    free_block(T1);
    */

    /*
    // print
    U = SharedTensor2d(new Tensor2d("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    U->sort(1324, T, 1.0, 0.0);
    U->print();
    U.reset();
    */

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

}// end ccsd_t2_amps
}} // End Namespaces

