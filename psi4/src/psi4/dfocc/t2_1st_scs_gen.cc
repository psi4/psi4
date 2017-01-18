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

/** Standard library includes */
#include "psi4/libqt/qt.h"
#include "defines.h"
#include "dfocc.h"

using namespace psi;
using namespace std;


namespace psi{ namespace dfoccwave{

void DFOCC::t2_1st_scs_gen()
{

    SharedTensor2d K, L, M;
    timer_on("1st-order T2");
    Fint_zero();

if (reference_ == "RESTRICTED") {
    // Build amplitudes in Mulliken order
    t2p_1new = SharedTensor2d(new Tensor2d("New T2_1 (ia|jb)", naoccA, navirA, naoccA, navirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(K);
    t2p_1new->copy(K);
    K.reset();

    // Read old amps
    t2p_1 = SharedTensor2d(new Tensor2d("T2_1 (ia|jb)", naoccA, navirA, naoccA, navirA));
    t2p_1->read(psio_, PSIF_DFOCC_AMPS);

    // Fint contributions
    // T'(ia,jb) += \sum_{e} T_ij^ae F_be = \sum_{e} T'(ia,je) F_be
    t2p_1new->contract424(4, 2, t2p_1, FabA, 1.0, 1.0);

    // T'(ia,jb) += \sum_{e} T_ij^eb F_ae = \sum_{e} T'(ie,jb) F_ae
    t2p_1new->contract424(2, 2, t2p_1, FabA, 1.0, 1.0);

    // T'(ia,jb) -= \sum_{m} T_im^ab F_mj = \sum_{e} T'(ia,mb) F_mj
    t2p_1new->contract424(3, 2, t2p_1, FijA, -1.0, 1.0);

    // T'(ia,jb) -= \sum_{m} T_mj^ab F_mi = \sum_{e} T'(ma,jb) F_mi
    t2p_1new->contract424(1, 2, t2p_1, FijA, -1.0, 1.0);

    // Aplly denominators
    if (regularization == "FALSE") t2p_1new->apply_denom_chem(nfrzc, noccA, FockA);
    else if (regularization == "TRUE") t2p_1new->reg_denom_chem(nfrzc, noccA, FockA, reg_param);

    // rms
    rms_t2 = 0.0;
    rms_t2 = t2p_1new->rms(t2p_1);
    //t2p_1new->print();

    // reset
    t2p_1->copy(t2p_1new);
    t2p_1new.reset();
    t2p_1->write(psio_, PSIF_DFOCC_AMPS);

    // Sort amplitudes to Dirac order
    t2_1 = SharedTensor2d(new Tensor2d("T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    t2_1->sort(1324, t2p_1, 1.0, 0.0);
    t2_1->write(psio_, PSIF_DFOCC_AMPS);
    //t2_1->print();
    t2_1.reset();

    SharedTensor2d temp = SharedTensor2d(new Tensor2d("U2_1 (ia|jb)", naoccA, navirA, naoccA, navirA));
    temp->sort(1432, t2p_1, 1.0, 0.0);
    temp->scale(-1.0);
    temp->add(t2p_1);
    temp->write(psio_, PSIF_DFOCC_AMPS);

    u2p_1 = SharedTensor2d(new Tensor2d("U2_1 (ia|jb)", naoccA, navirA, naoccA, navirA));
    u2p_1->copy(temp);
    temp.reset();
    u2p_1->add(t2p_1);
    t2p_1.reset();
    u2p_1->write(psio_, PSIF_DFOCC_AMPS);
    u2p_1.reset();
}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // T2AA
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(L);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ||AB>", naoccA, naoccA, navirA, navirA));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    t2_1newAA = SharedTensor2d(new Tensor2d("New T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    t2_1newAA->copy(K);
    K.reset();

    // Fint contributions
    t2_1AA = SharedTensor2d(new Tensor2d("T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    t2_1AA->read(psio_, PSIF_DFOCC_AMPS);

    // T(IJ,AB) += \sum_{E} T_IJ^AE F_BE
    t2_1newAA->contract424(4, 2, t2_1AA, FabA, 1.0, 1.0);

    // T(IJ,AB) += \sum_{E} T_IJ^EB F_AE
    t2_1newAA->contract424(3, 2, t2_1AA, FabA, 1.0, 1.0);

    // T(IJ,AB) -= \sum_{M} T_IM^AB F_MJ
    t2_1newAA->contract424(2, 1, t2_1AA, FijA, -1.0, 1.0);

    // T(IJ,AB) -= \sum_{M} T_MJ^AB F_MI
    t2_1newAA->contract424(1, 1, t2_1AA, FijA, -1.0, 1.0);

    // apply denom
    if (regularization == "FALSE") t2_1newAA->apply_denom(nfrzc, noccA, FockA);
    else if (regularization == "TRUE") t2_1newAA->reg_denom(nfrzc, noccA, FockA, reg_param);
    if (print_ > 2) t2_1newAA->print();

    // rms
    rms_t2AA = 0.0;
    rms_t2AA = t2_1newAA->rms(t2_1AA);

    // reset
    t2_1AA->copy(t2_1newAA);
    t2_1newAA.reset();
    t2_1AA->write(psio_, PSIF_DFOCC_AMPS);
    t2p_1 = SharedTensor2d(new Tensor2d("T2_1 (IA|JB)", naoccA, navirA, naoccA, navirA));
    t2p_1->sort(1324, t2_1AA, 1.0, 0.0);
    t2_1AA.reset();
    t2p_1->write(psio_, PSIF_DFOCC_AMPS);
    t2p_1.reset();

    // T2BB
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB));
    tei_iajb_chem_directBB(L);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|ab>", naoccB, naoccB, navirB, navirB));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij||ab>", naoccB, naoccB, navirB, navirB));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    t2_1newBB = SharedTensor2d(new Tensor2d("New T2_1 <ij|ab>", naoccB, naoccB, navirB, navirB));
    t2_1newBB->copy(K);
    K.reset();

    // Fint contributions
    t2_1BB = SharedTensor2d(new Tensor2d("T2_1 <ij|ab>", naoccB, naoccB, navirB, navirB));
    t2_1BB->read(psio_, PSIF_DFOCC_AMPS);

    // T(ij,ab) += \sum_{e} T_ij^ae F_be
    t2_1newBB->contract424(4, 2, t2_1BB, FabB, 1.0, 1.0);

    // T(ij,ab) += \sum_{e} T_ij^eb F_ae
    t2_1newBB->contract424(3, 2, t2_1BB, FabB, 1.0, 1.0);

    // T(ij,ab) -= \sum_{m} T_im^ab F_mj
    t2_1newBB->contract424(2, 1, t2_1BB, FijB, -1.0, 1.0);

    // T(ij,ab) -= \sum_{m} T_mj^ab F_mi
    t2_1newBB->contract424(1, 1, t2_1BB, FijB, -1.0, 1.0);

    // apply denom
    if (regularization == "FALSE") t2_1newBB->apply_denom(nfrzc, noccB, FockB);
    else if (regularization == "TRUE") t2_1newBB->reg_denom(nfrzc, noccB, FockB, reg_param);
    if (print_ > 2) t2_1newBB->print();

    // rms
    rms_t2BB = 0.0;
    rms_t2BB = t2_1newBB->rms(t2_1BB);

    // reset
    t2_1BB->copy(t2_1newBB);
    t2_1newBB.reset();
    t2_1BB->write(psio_, PSIF_DFOCC_AMPS);
    t2p_1 = SharedTensor2d(new Tensor2d("T2_1 (ia|jb)", naoccB, navirB, naoccB, navirB));
    t2p_1->sort(1324, t2_1BB, 1.0, 0.0);
    t2_1BB.reset();
    t2p_1->write(psio_, PSIF_DFOCC_AMPS);
    t2p_1.reset();

    // T2AB
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|jb)", naoccA, navirA, naoccB, navirB));
    tei_iajb_chem_directAB(L);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    t2_1newAB = SharedTensor2d(new Tensor2d("New T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    t2_1newAB->copy(K);
    K.reset();

    // Fint contributions
    t2_1AB = SharedTensor2d(new Tensor2d("T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    t2_1AB->read(psio_, PSIF_DFOCC_AMPS);

    // T(Ij,Ab) += \sum_{e} T_Ij^Ae F_be
    t2_1newAB->contract424(4, 2, t2_1AB, FabB, 1.0, 1.0);

    // T(Ij,Ab) += \sum_{E} T_Ij^Eb F_AE
    t2_1newAB->contract424(3, 2, t2_1AB, FabA, 1.0, 1.0);

    // T(Ij,Ab) -= \sum_{m} T_Im^Ab F_mj
    t2_1newAB->contract424(2, 1, t2_1AB, FijB, -1.0, 1.0);

    // T(Ij,Ab) -= \sum_{M} T_Mj^Ab F_MI
    t2_1newAB->contract424(1, 1, t2_1AB, FijA, -1.0, 1.0);

    // apply denom
    if (regularization == "FALSE") t2_1newAB->apply_denom_os(nfrzc, noccA, noccB, FockA, FockB);
    else if (regularization == "TRUE") t2_1newAB->reg_denom_os(nfrzc, noccA, noccB, FockA, FockB, reg_param);
    if (print_ > 2) t2_1newAB->print();

    // rms
    rms_t2AB = 0.0;
    rms_t2AB = t2_1newAB->rms(t2_1AB);

    // reset
    t2_1AB->copy(t2_1newAB);
    t2_1newAB.reset();
    t2_1AB->write(psio_, PSIF_DFOCC_AMPS);
    t2p_1 = SharedTensor2d(new Tensor2d("T2_1 (IA|jb)", naoccA, navirA, naoccB, navirB));
    t2p_1->sort(1324, t2_1AB, 1.0, 0.0);
    t2_1AB.reset();
    t2p_1->write(psio_, PSIF_DFOCC_AMPS);
    t2p_1.reset();

    // combined rms
    double rms_ss = MAX0(rms_t2AA, rms_t2BB);
    rms_t2 = MAX0(rms_ss, rms_t2AB);

}// else if (reference_ == "UNRESTRICTED")
    timer_off("1st-order T2");
} // end t2_1st_sc

}} // End Namespaces
