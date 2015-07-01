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
  
void DFOCC::t2_2nd_sc()
{

    // defs
    SharedTensor2d K, L, M, I, T, Tnew, T1, T2, U, Tau, W, X, Y;

if (reference_ == "RESTRICTED") {
    // Read DF integrals
    bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA));
    bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA));
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQijA->read(psio_, PSIF_DFOCC_INTS);
    bQiaA->read(psio_, PSIF_DFOCC_INTS);
    bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);

    /*
    // t_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = \sum_{e} t_ij^ae F_be = \sum_{e} T(ia,je) F_be
    X = SharedTensor2d(new Tensor2d("X (IA|JB)", naoccA, navirA, naoccA, navirA));
    X->contract(false, true, naoccA * navirA * naoccA, navirA, navirA, t2, FabA, 1.0, 0.0);

    // t_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = -\sum_{m} t_mj^ab F_mi = -\sum_{m} F(m,i) T(ma,jb)
    X->contract(true, false, naoccA, naoccA * navirA * navirA, naoccA, FijA, t2, -1.0, 1.0);
    X->symmetrize();

    // Contributions of X
    Tnew->axpy(X, 2.0);
    X.reset();
    */

    // Read T2_1
    t2 = SharedTensor2d(new Tensor2d("T2_1 (ia|jb)", naoccA, navirA, naoccA, navirA));
    t2->read_symm(psio_, PSIF_DFOCC_AMPS);

    // WmnijT2
    mp3_WmnijT2();

    // WmbejT2
    mp3_WmbejT2();

    // WabefT2
    mp3_WabefT2();

    // Denom
    Tnew = SharedTensor2d(new Tensor2d("New T2_2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->apply_denom_chem(nfrzc, noccA, FockA);

    // Form T2 = T2(1) + T2(2)
    t2->axpy(Tnew, 1.0);
    Tnew.reset();

    // Form U(ia,jb) = 2*T(ia,jb) - T (ib,ja)
    U = SharedTensor2d(new Tensor2d("U2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_u2_amps(U,t2);
    t2.reset();

    // Energy
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    Ecorr = U->vector_dot(K);
    U.reset();
    K.reset();
    Emp3 = Escf + Ecorr;

    // Free ints
    bQijA.reset();
    bQiaA.reset();
    bQabA.reset();
}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // Read DF integrals
    bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA));
    bQijB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ij)", nQ, naoccB, naoccB));
    bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA));
    bQiaB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ia)", nQ, naoccB, navirB));
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQabB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ab)", nQ, navirB, navirB));
    bQijA->read(psio_, PSIF_DFOCC_INTS);
    bQijB->read(psio_, PSIF_DFOCC_INTS);
    bQiaA->read(psio_, PSIF_DFOCC_INTS);
    bQiaB->read(psio_, PSIF_DFOCC_INTS);
    bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);
    bQabB->read(psio_, PSIF_DFOCC_INTS, true, true);

    //=========================
    // T2AA
    //=========================

    // WmnijT2
    mp3_WmnijT2AA();

    // WmbejT2
    mp3_WmbejT2AA();

    // WabefT2
    mp3_WabefT2AA();

    // Denom
    Tnew = SharedTensor2d(new Tensor2d("New T2_2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    Tnew->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->apply_denom(nfrzc, noccA, FockA);

    // Form T2 = T2(1) + T2(2)
    U = SharedTensor2d(new Tensor2d("T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    U->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    U->axpy(Tnew, 1.0);
    Tnew.reset();

    // Energy
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    L->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ||AB>", naoccA, naoccA, navirA, navirA));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    Emp3AA = 0.25 * U->vector_dot(K);
    U.reset();
    K.reset();

    //=========================
    // T2BB
    //=========================

    // WmnijT2
    mp3_WmnijT2BB();

    // WmbejT2
    mp3_WmbejT2BB();

    // WabefT2
    mp3_WabefT2BB();

    // Denom
    Tnew = SharedTensor2d(new Tensor2d("New T2_2 <ij|ab>", naoccB, naoccB, navirB, navirB));
    Tnew->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->apply_denom(nfrzc, noccB, FockB);

    // Form T2 = T2(1) + T2(2)
    U = SharedTensor2d(new Tensor2d("T2_1 <ij|ab>", naoccB, naoccB, navirB, navirB));
    U->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    U->axpy(Tnew, 1.0);
    Tnew.reset();

    // Energy BB
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB));
    L->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|ab>", naoccB, naoccB, navirB, navirB));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij||ab>", naoccB, naoccB, navirB, navirB));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    Emp3BB = 0.25 * U->vector_dot(K);
    U.reset();
    K.reset();

    //=========================
    // T2AB
    //=========================

    // WmnijT2
    mp3_WmnijT2AB();

    // WmbejT2
    mp3_WmbejT2AB();

    // WabefT2
    mp3_WabefT2AB();

    // Denom
    Tnew = SharedTensor2d(new Tensor2d("New T2_2 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    Tnew->read(psio_, PSIF_DFOCC_AMPS);
    Tnew->apply_denom_os(nfrzc, noccA, noccB, FockA, FockB);

    // Form T2 = T2(1) + T2(2)
    U = SharedTensor2d(new Tensor2d("T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    U->read(psio_, PSIF_DFOCC_AMPS);
    U->axpy(Tnew, 1.0);
    Tnew.reset();

    // Energy
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|jb)", naoccA, navirA, naoccB, navirB));
    L->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    Emp3AB = U->vector_dot(K);
    U.reset();
    K.reset();

    // Overall energy
    Ecorr = Emp3AA + Emp3BB + Emp3AB;
    Emp3 = Escf + Ecorr;

    // Free ints
    bQijA.reset();
    bQijB.reset();
    bQiaA.reset();
    bQiaB.reset();
    bQabA.reset();
    bQabB.reset();

}// else if (reference_ == "UNRESTRICTED")

}// end t2_2nd_sc

}} // End Namespaces

