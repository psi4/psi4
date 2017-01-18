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

void DFOCC::t2_1st_scs_sc()
{
    SharedTensor2d K, L, M;
    timer_on("1st-order T2");
if (reference_ == "RESTRICTED") {
    // Example: init from file
    //Array2d *temp = new Array2d(psio_, PSIF_DFOCC_INTS, "(ia|jb)", naoccA, navirA, naoccA, navirA);

    // Build amplitudes in Mulliken order
    t2p_1 = SharedTensor2d(new Tensor2d("T2_1 (ia|jb)", naoccA, navirA, naoccA, navirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(K);
    t2p_1->copy(K);
    K.reset();
    if (regularization == "FALSE") t2p_1->apply_denom_chem(nfrzc, noccA, FockA);
    else if (regularization == "TRUE") t2p_1->reg_denom_chem(nfrzc, noccA, FockA, reg_param);
    t2p_1->write(psio_, PSIF_DFOCC_AMPS);

    // Sort amplitudes to Dirac order
    t2_1 = SharedTensor2d(new Tensor2d("T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    t2_1->sort(1324, t2p_1, 1.0, 0.0);
    t2_1->write(psio_, PSIF_DFOCC_AMPS);
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
    t2_1AA = SharedTensor2d(new Tensor2d("T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    t2_1AA->copy(K);
    K.reset();
    if (regularization == "FALSE") t2_1AA->apply_denom(nfrzc, noccA, FockA);
    else if (regularization == "TRUE") t2_1AA->reg_denom(nfrzc, noccA, FockA, reg_param);
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
    t2_1BB = SharedTensor2d(new Tensor2d("T2_1 <ij|ab>", naoccB, naoccB, navirB, navirB));
    t2_1BB->copy(K);
    K.reset();
    if (regularization == "FALSE") t2_1BB->apply_denom(nfrzc, noccB, FockB);
    else if (regularization == "TRUE") t2_1BB->reg_denom(nfrzc, noccB, FockB, reg_param);
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
    t2_1AB = SharedTensor2d(new Tensor2d("T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    t2_1AB->copy(K);
    K.reset();
    if (regularization == "FALSE") t2_1AB->apply_denom_os(nfrzc, noccA, noccB, FockA, FockB);
    else if (regularization == "TRUE") t2_1AB->reg_denom_os(nfrzc, noccA, noccB, FockA, FockB, reg_param);
    t2_1AB->write(psio_, PSIF_DFOCC_AMPS);
    t2p_1 = SharedTensor2d(new Tensor2d("T2_1 (IA|jb)", naoccA, navirA, naoccB, navirB));
    t2p_1->sort(1324, t2_1AB, 1.0, 0.0);
    t2_1AB.reset();
    t2p_1->write(psio_, PSIF_DFOCC_AMPS);
    t2p_1.reset();
}// else if (reference_ == "UNRESTRICTED")
    timer_off("1st-order T2");
} // end t2_1st_sc

}} // End Namespaces
