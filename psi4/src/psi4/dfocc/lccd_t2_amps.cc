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

void DFOCC::lccd_t2_amps()
{

    // defs
    SharedTensor2d K, L, M, I, T, Tnew, T1, T2, U, Tau, W, X, Y;
    SharedTensor2d R, RAA, RBB, RAB, TAA, TBB, TAB;

if (reference_ == "RESTRICTED") {
    // Read DF integrals
    bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA));
    bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA));
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQijA->read(psio_, PSIF_DFOCC_INTS);
    bQiaA->read(psio_, PSIF_DFOCC_INTS);
    bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);

    // Read T2 amps
    t2 = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    t2->read_symm(psio_, PSIF_DFOCC_AMPS);

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
    lccd_WmnijT2();

    // WmbejT2
    lccd_WmbejT2();

    // WabefT2
    lccd_WabefT2();

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
    if (orb_opt_ == "FALSE" || mo_optimized == 1) {
        std::shared_ptr<Matrix> RT2(new Matrix("RT2", naoccA*navirA, naoccA*navirA));
        Tau->to_matrix(RT2);
        Tau.reset();
        std::shared_ptr<Matrix> T2(new Matrix("T2", naoccA*navirA, naoccA*navirA));
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
    }

    // Write
    t2->write_symm(psio_, PSIF_DFOCC_AMPS);

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
    Elccd = Eref + Ecorr;

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

    // Read T2
    U = SharedTensor2d(new Tensor2d("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    U->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    T->sort(1324, U, 1.0, 0.0);
    U.reset();

    // t_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = \sum_{e} t_ij^ae F_be = \sum_{e} T(ia,je) F_be
    X = SharedTensor2d(new Tensor2d("X (IA|JB)", naoccA, navirA, naoccA, navirA));
    X->contract(false, true, naoccA * navirA * naoccA, navirA, navirA, T, FabA, 1.0, 0.0);

    // t_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = -\sum_{m} t_mj^ab F_mi = -\sum_{m} F(m,i) T(ma,jb)
    X->contract(true, false, naoccA, naoccA * navirA * navirA, naoccA, FijA, T, -1.0, 1.0);
    T.reset();
    X->symmetrize();
    X->scale(2.0);

    // t_ij^ab = <ij||ab>
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    L->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ||AB>", naoccA, naoccA, navirA, navirA));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    Tnew = SharedTensor2d(new Tensor2d("New T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    Tnew->copy(K);
    K.reset();

    // Contributions of X
    Tnew->sort(1324, X, 1.0, 1.0);
    X.reset();

    // Write and close
    Tnew->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    // WmnijT2
    lccd_WmnijT2AA();

    // WmbejT2
    lccd_WmbejT2AA();

    // WabefT2
    lccd_WabefT2AA();

    // Denom
    Tnew = SharedTensor2d(new Tensor2d("New T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    Tnew->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->apply_denom(nfrzc, noccA, FockA);
    Tnew->write_anti_symm(psio_, PSIF_DFOCC_AMPS);

    // Error vector
    T = SharedTensor2d(new Tensor2d("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    T->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    R = SharedTensor2d(new Tensor2d("RT2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    R->copy(Tnew);
    Tnew.reset();
    R->subtract(T);
    T.reset();
    R->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    R.reset();

    //=========================
    // T2BB
    //=========================

    // Read T2
    U = SharedTensor2d(new Tensor2d("T2 <ij|ab>", naoccB, naoccB, navirB, navirB));
    U->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T = SharedTensor2d(new Tensor2d("T2 (ia|jb)", naoccB, navirB, naoccB, navirB));
    T->sort(1324, U, 1.0, 0.0);
    U.reset();

    // t_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = \sum_{e} t_ij^ae F_be = \sum_{e} T(ia,je) F_be
    X = SharedTensor2d(new Tensor2d("X (ia|jb)", naoccB, navirB, naoccB, navirB));
    X->contract(false, true, naoccB * navirB * naoccB, navirB, navirB, T, FabB, 1.0, 0.0);

    // t_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
    // X(ia,jb) = -\sum_{m} t_mj^ab F_mi = -\sum_{m} F(m,i) T(ma,jb)
    X->contract(true, false, naoccB, naoccB * navirB * navirB, naoccB, FijB, T, -1.0, 1.0);
    T.reset();
    X->symmetrize();
    X->scale(2.0);

    // t_ij^ab = <ij||ab>
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB));
    L->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|ab>", naoccB, naoccB, navirB, navirB));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij||ab>", naoccB, naoccB, navirB, navirB));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    Tnew = SharedTensor2d(new Tensor2d("New T2 <ij|ab>", naoccB, naoccB, navirB, navirB));
    Tnew->copy(K);
    K.reset();

    // Contributions of X
    Tnew->sort(1324, X, 1.0, 1.0);
    X.reset();

    // write and close
    Tnew->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    // WmnijT2
    lccd_WmnijT2BB();

    // WmbejT2
    lccd_WmbejT2BB();

    // WabefT2
    lccd_WabefT2BB();

    // Denom
    Tnew = SharedTensor2d(new Tensor2d("New T2 <ij|ab>", naoccB, naoccB, navirB, navirB));
    Tnew->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tnew->apply_denom(nfrzc, noccB, FockB);
    Tnew->write_anti_symm(psio_, PSIF_DFOCC_AMPS);

    // Error vector
    T = SharedTensor2d(new Tensor2d("T2 <ij|ab>", naoccB, naoccB, navirB, navirB));
    T->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    R = SharedTensor2d(new Tensor2d("RT2 <ij|ab>", naoccB, naoccB, navirB, navirB));
    R->copy(Tnew);
    Tnew.reset();
    R->subtract(T);
    T.reset();
    R->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    R.reset();

    //=========================
    // T2AB
    //=========================

    // t_ij^ab = <ij|ab>
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|jb)", naoccA, navirA, naoccB, navirB));
    L->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    Tnew = SharedTensor2d(new Tensor2d("New T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    Tnew->copy(K);
    K.reset();

    // Read T2
    T = SharedTensor2d(new Tensor2d("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    T->read(psio_, PSIF_DFOCC_AMPS);

    // T(Ij,Ab) += \sum_{e} T_Ij^Ae F_be
    Tnew->contract424(4, 2, T, FabB, 1.0, 1.0);

    // T(Ij,Ab) += \sum_{E} T_Ij^Eb F_AE
    Tnew->contract424(3, 2, T, FabA, 1.0, 1.0);

    // T(Ij,Ab) -= \sum_{m} T_Im^Ab F_mj
    Tnew->contract424(2, 1, T, FijB, -1.0, 1.0);

    // T(Ij,Ab) -= \sum_{M} T_Mj^Ab F_MI
    Tnew->contract424(1, 1, T, FijA, -1.0, 1.0);
    T.reset();

    // write and close
    Tnew->write(psio_, PSIF_DFOCC_AMPS);
    Tnew.reset();

    // WmnijT2
    lccd_WmnijT2AB();

    // WmbejT2
    lccd_WmbejT2AB();

    // WabefT2
    lccd_WabefT2AB();

    // Denom
    Tnew = SharedTensor2d(new Tensor2d("New T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    Tnew->read(psio_, PSIF_DFOCC_AMPS);
    Tnew->apply_denom_os(nfrzc, noccA, noccB, FockA, FockB);
    Tnew->write(psio_, PSIF_DFOCC_AMPS);

    // Error vector
    T = SharedTensor2d(new Tensor2d("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    T->read(psio_, PSIF_DFOCC_AMPS);
    R = SharedTensor2d(new Tensor2d("RT2 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    R->copy(Tnew);
    Tnew.reset();
    R->subtract(T);
    T.reset();
    R->write(psio_, PSIF_DFOCC_AMPS);
    R.reset();

    //=========================
    // DIIS
    //=========================

 if (orb_opt_ == "FALSE" || mo_optimized == 1) {
    // RAA
    RAA = SharedTensor2d(new Tensor2d("RT2 <IJ|AB>", ntri_anti_ijAA, ntri_anti_abAA));
    RAA->read(psio_, PSIF_DFOCC_AMPS);
    std::shared_ptr<Matrix> RT2AA(new Matrix("RT2AA", ntri_anti_ijAA, ntri_anti_abAA));
    RAA->to_matrix(RT2AA);
    RAA.reset();
    // TAA
    TAA = SharedTensor2d(new Tensor2d("New T2 <IJ|AB>", ntri_anti_ijAA, ntri_anti_abAA));
    TAA->read(psio_, PSIF_DFOCC_AMPS);
    std::shared_ptr<Matrix> T2AA(new Matrix("T2AA", ntri_anti_ijAA, ntri_anti_abAA));
    TAA->to_matrix(T2AA);

    // RBB
    RBB = SharedTensor2d(new Tensor2d("RT2 <ij|ab>", ntri_anti_ijBB, ntri_anti_abBB));
    RBB->read(psio_, PSIF_DFOCC_AMPS);
    std::shared_ptr<Matrix> RT2BB(new Matrix("RT2BB", ntri_anti_ijBB, ntri_anti_abBB));
    RBB->to_matrix(RT2BB);
    RBB.reset();
    // TBB
    TBB = SharedTensor2d(new Tensor2d("New T2 <ij|ab>", ntri_anti_ijBB, ntri_anti_abBB));
    TBB->read(psio_, PSIF_DFOCC_AMPS);
    std::shared_ptr<Matrix> T2BB(new Matrix("T2BB", ntri_anti_ijBB, ntri_anti_abBB));
    TBB->to_matrix(T2BB);

    // RAB
    RAB = SharedTensor2d(new Tensor2d("RT2 <Ij|Ab>", naoccA*naoccB, navirA*navirB));
    RAB->read(psio_, PSIF_DFOCC_AMPS);
    std::shared_ptr<Matrix> RT2AB(new Matrix("RT2AB", naoccA*naoccB, navirA*navirB));
    RAB->to_matrix(RT2AB);
    RAB.reset();
    // TAB
    TAB = SharedTensor2d(new Tensor2d("New T2 <Ij|Ab>", naoccA*naoccB, navirA*navirB));
    TAB->read(psio_, PSIF_DFOCC_AMPS);
    std::shared_ptr<Matrix> T2AB(new Matrix("T2AB", naoccA*naoccB, navirA*navirB));
    TAB->to_matrix(T2AB);

    // add entry
    if (do_diis_ == 1) ccsdDiisManager->add_entry(6, RT2AA.get(), RT2BB.get(), RT2AB.get(), T2AA.get(), T2BB.get(), T2AB.get());
    RT2AA.reset();
    RT2BB.reset();
    RT2BB.reset();

    // extrapolate
    if (do_diis_ == 1) {
        if (ccsdDiisManager->subspace_size() >= cc_mindiis_) ccsdDiisManager->extrapolate(3, T2AA.get(), T2BB.get(), T2AB.get());
        TAA->set2(T2AA);
        TBB->set2(T2BB);
        TAB->set2(T2AB);
    }
    T2AA.reset();
    T2BB.reset();
    T2AB.reset();
    TAA->write(psio_, PSIF_DFOCC_AMPS);
    TBB->write(psio_, PSIF_DFOCC_AMPS);
    TAB->write(psio_, PSIF_DFOCC_AMPS);
 }// end if diis

    //=========================
    // Reset & Energy
    //=========================

    // Reset T2AA
    Tnew = SharedTensor2d(new Tensor2d("New T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    Tnew->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T = SharedTensor2d(new Tensor2d("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    T->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    rms_t2AA = Tnew->rms(T);
    T->copy(Tnew);
    Tnew.reset();
    T->write_anti_symm(psio_, PSIF_DFOCC_AMPS);

    // AA Energy
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    L->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ||AB>", naoccA, naoccA, navirA, navirA));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    ElccdAA = 0.25 * T->vector_dot(K);
    K.reset();

    // Form T2(IA|JB)
    U = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    U->sort(1324, T, 1.0, 0.0);
    T.reset();
    U->write_symm(psio_, PSIF_DFOCC_AMPS);
    U.reset();

    // Reset T2BB
    Tnew = SharedTensor2d(new Tensor2d("New T2 <ij|ab>", naoccB, naoccB, navirB, navirB));
    Tnew->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T = SharedTensor2d(new Tensor2d("T2 <ij|ab>", naoccB, naoccB, navirB, navirB));
    T->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    rms_t2BB = Tnew->rms(T);
    T->copy(Tnew);
    Tnew.reset();
    T->write_anti_symm(psio_, PSIF_DFOCC_AMPS);

    // BB Energy
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB));
    L->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|ab>", naoccB, naoccB, navirB, navirB));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij||ab>", naoccB, naoccB, navirB, navirB));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    ElccdBB = 0.25 * T->vector_dot(K);
    K.reset();

    // Form T2(ia|jb)
    U = SharedTensor2d(new Tensor2d("T2 (ia|jb)", naoccB, navirB, naoccB, navirB));
    U->sort(1324, T, 1.0, 0.0);
    T.reset();
    U->write_symm(psio_, PSIF_DFOCC_AMPS);
    U.reset();

    // Reset T2AB
    Tnew = SharedTensor2d(new Tensor2d("New T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    Tnew->read(psio_, PSIF_DFOCC_AMPS);
    T = SharedTensor2d(new Tensor2d("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    T->read(psio_, PSIF_DFOCC_AMPS);
    rms_t2AB = Tnew->rms(T);
    T->copy(Tnew);
    Tnew.reset();
    T->write(psio_, PSIF_DFOCC_AMPS);

    // AB Energy
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|jb)", naoccA, navirA, naoccB, navirB));
    L->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    ElccdAB = T->vector_dot(K);
    K.reset();

    // Overall energy
    Ecorr = ElccdAA + ElccdBB + ElccdAB;
    Elccd = Eref + Ecorr;

    // Form T2(IA|jb)
    U = SharedTensor2d(new Tensor2d("T2 (IA|jb)", naoccA, navirA, naoccB, navirB));
    U->sort(1324, T, 1.0, 0.0);
    T.reset();
    U->write(psio_, PSIF_DFOCC_AMPS);
    U.reset();

    // Free ints
    bQijA.reset();
    bQijB.reset();
    bQiaA.reset();
    bQiaB.reset();
    bQabA.reset();
    bQabB.reset();

    // combined rms
    double rms_ss = MAX0(rms_t2AA, rms_t2BB);
    rms_t2 = MAX0(rms_ss, rms_t2AB);

}// else if (reference_ == "UNRESTRICTED")

}// end lccd_t2_amps

}} // End Namespaces
