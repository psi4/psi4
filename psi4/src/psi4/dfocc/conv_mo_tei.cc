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

void DFOCC::tei_ijkl_chem()
{
    timer_on("Build (oo|oo)");
    // AA spin case
    JijklAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|KL)", naoccA, naoccA, naoccA, naoccA));
    bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA * naoccA));
    bQijA->read(psio_, PSIF_DFOCC_INTS);
    JijklAA->gemm(true, false, bQijA, bQijA, 1.0, 0.0);
    if (reference_ == "RESTRICTED") bQijA.reset();
    JijklAA->write(psio_, PSIF_DFOCC_INTS);
    JijklAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    JijklBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ij|kl)", naoccB, naoccB, naoccB, naoccB));
    bQijB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ij)", nQ, naoccB * naoccB));
    bQijB->read(psio_, PSIF_DFOCC_INTS);
    JijklBB->gemm(true, false, bQijB, bQijB, 1.0, 0.0);
    JijklBB->write(psio_, PSIF_DFOCC_INTS);
    JijklBB.reset();

    // AB spin case
    JijklAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|kl)", naoccA, naoccA, naoccB, naoccB));
    JijklAB->gemm(true, false, bQijA, bQijB, 1.0, 0.0);
    bQijA.reset();
    bQijB.reset();
    JijklAB->write(psio_, PSIF_DFOCC_INTS);
    JijklAB.reset();
 }
    timer_off("Build (oo|oo)");
}// end tei_ijkl_chem

//=======================================================
//          (oo|oo)
//=======================================================
void DFOCC::tei_oooo_chem()
{
    timer_on("Build (oo|oo)");
    // AA spin case
    JooooAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|OO)", noccA, noccA, noccA, noccA));
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA * noccA));
    bQooA->read(psio_, PSIF_DFOCC_INTS);
    JooooAA->gemm(true, false, bQooA, bQooA, 1.0, 0.0);
    if (reference_ == "RESTRICTED") bQooA.reset();
    JooooAA->write(psio_, PSIF_DFOCC_INTS);
    JooooAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    JooooBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (oo|oo)", noccB, noccB, noccB, noccB));
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|oo)", nQ, noccB * noccB));
    bQooB->read(psio_, PSIF_DFOCC_INTS);
    JooooBB->gemm(true, false, bQooB, bQooB, 1.0, 0.0);
    JooooBB->write(psio_, PSIF_DFOCC_INTS);
    JooooBB.reset();

    // AB spin case
    JooooAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|oo)", noccA, noccA, noccB, noccB));
    JooooAB->gemm(true, false, bQooA, bQooB, 1.0, 0.0);
    bQooA.reset();
    bQooB.reset();
    JooooAB->write(psio_, PSIF_DFOCC_INTS);
    JooooAB.reset();
 }
    timer_off("Build (oo|oo)");
}// end tei_oooo_chem

//=======================================================
//          (IJ|KA)
//=======================================================
void DFOCC::tei_ijka_chem()
{
    timer_on("Build (oo|ov)");
    // AA spin case
    JijkaAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|KA)", naoccA, naoccA, naoccA, navirA));
    bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA * naoccA));
    bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA * navirA));
    bQijA->read(psio_, PSIF_DFOCC_INTS);
    bQiaA->read(psio_, PSIF_DFOCC_INTS);
    JijkaAA->gemm(true, false, bQijA, bQiaA, 1.0, 0.0);
    if (reference_ == "RESTRICTED") bQijA.reset();
    if (reference_ == "RESTRICTED") bQiaA.reset();
    JijkaAA->write(psio_, PSIF_DFOCC_INTS);
    JijkaAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    JijkaBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ij|ka)", naoccB, naoccB, naoccB, navirB));
    bQijB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ij)", nQ, naoccB * naoccB));
    bQiaB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ia)", nQ, naoccB * navirB));
    bQijB->read(psio_, PSIF_DFOCC_INTS);
    bQiaB->read(psio_, PSIF_DFOCC_INTS);
    JijkaBB->gemm(true, false, bQijB, bQiaB, 1.0, 0.0);
    JijkaBB->write(psio_, PSIF_DFOCC_INTS);
    JijkaBB.reset();

    // AB spin case
    JijkaAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|ka)", naoccA, naoccA, naoccB, navirB));
    JijkaAB->gemm(true, false, bQijA, bQiaB, 1.0, 0.0);
    bQijA.reset();
    bQiaB.reset();
    JijkaAB->write(psio_, PSIF_DFOCC_INTS);
    JijkaAB.reset();

    // BA spin case
    JiajkAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|jk)", naoccA, navirA, naoccB, naoccB));
    JiajkAB->gemm(true, false, bQiaA, bQijB, 1.0, 0.0);
    bQijB.reset();
    bQiaA.reset();
    JiajkAB->write(psio_, PSIF_DFOCC_INTS);
    JiajkAB.reset();
 }
    timer_off("Build (oo|ov)");
}// end tei_ijka_chem

//=======================================================
//          (OO|OV)
//=======================================================
void DFOCC::tei_ooov_chem()
{
    timer_on("Build (oo|ov)");
    // AA spin case
    JooovAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|OV)", noccA, noccA, noccA, nvirA));
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA * noccA));
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA * nvirA));
    bQooA->read(psio_, PSIF_DFOCC_INTS);
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    JooovAA->gemm(true, false, bQooA, bQovA, 1.0, 0.0);
    if (reference_ == "RESTRICTED") bQooA.reset();
    if (reference_ == "RESTRICTED") bQovA.reset();
    JooovAA->write(psio_, PSIF_DFOCC_INTS);
    JooovAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    JooovBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (oo|ov)", noccB, noccB, noccB, nvirB));
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|oo)", nQ, noccB * noccB));
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ov)", nQ, noccB * nvirB));
    bQooB->read(psio_, PSIF_DFOCC_INTS);
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    JooovBB->gemm(true, false, bQooB, bQovB, 1.0, 0.0);
    JooovBB->write(psio_, PSIF_DFOCC_INTS);
    JooovBB.reset();

    // AB spin case
    JooovAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|ov)", noccA, noccA, noccB, nvirB));
    JooovAB->gemm(true, false, bQooA, bQovB, 1.0, 0.0);
    bQooA.reset();
    bQovB.reset();
    JooovAB->write(psio_, PSIF_DFOCC_INTS);
    JooovAB.reset();

    // BA spin case
    JovooAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OV|oo)", noccA, nvirA, noccB, noccB));
    JovooAB->gemm(true, false, bQovA, bQooB, 1.0, 0.0);
    bQooB.reset();
    bQovA.reset();
    JovooAB->write(psio_, PSIF_DFOCC_INTS);
    JovooAB.reset();
 }
    timer_off("Build (oo|ov)");
}// end tei_ooov_chem

//=======================================================
//          (IJ|AB)
//=======================================================
void DFOCC::tei_ijab_chem()
{
    timer_on("Build (oo|vv)");
    // AA spin case
    JijabAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|AB)", naoccA, naoccA, navirA, navirA));
    bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA));
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQijA->read(psio_, PSIF_DFOCC_INTS);
    bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);
    JijabAA->gemm(true, false, bQijA, bQabA, 1.0, 0.0);
    if (reference_ == "RESTRICTED") bQijA.reset();
    if (reference_ == "RESTRICTED") bQabA.reset();
    JijabAA->write(psio_, PSIF_DFOCC_INTS);
    JijabAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    JijabBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ij|ab)", naoccB, naoccB, navirB, navirB));
    bQijB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ij)", nQ, naoccB, naoccB));
    bQabB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ab)", nQ, navirB, navirB));
    bQijB->read(psio_, PSIF_DFOCC_INTS);
    bQabB->read(psio_, PSIF_DFOCC_INTS, true, true);
    JijabBB->gemm(true, false, bQijB, bQabB, 1.0, 0.0);
    JijabBB->write(psio_, PSIF_DFOCC_INTS);
    JijabBB.reset();

    // AB spin case
    JijabAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|ab)", naoccA, naoccA, navirB, navirB));
    JijabAB->gemm(true, false, bQijA, bQabB, 1.0, 0.0);
    bQijA.reset();
    bQabB.reset();
    JijabAB->write(psio_, PSIF_DFOCC_INTS);
    JijabAB.reset();

    // BA spin case
    JabijAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (AB|ij)",navirA, navirA, naoccB, naoccB));
    JabijAB->gemm(true, false, bQabA, bQijB, 1.0, 0.0);
    bQijB.reset();
    bQabA.reset();
    JabijAB->write(psio_, PSIF_DFOCC_INTS);
    JabijAB.reset();
 }
    timer_off("Build (oo|vv)");
}// end tei_ijab_chem

//=======================================================
//          (OO|VV)
//=======================================================
void DFOCC::tei_oovv_chem()
{
    timer_on("Build (oo|vv)");
    // AA spin case
    JoovvAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|VV)", noccA, noccA, nvirA, nvirA));
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA * noccA));
    bQvvA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA));
    bQooA->read(psio_, PSIF_DFOCC_INTS);
    bQvvA->read(psio_, PSIF_DFOCC_INTS, true, true);
    JoovvAA->gemm(true, false, bQooA, bQvvA, 1.0, 0.0);
    if (reference_ == "RESTRICTED") bQooA.reset();
    if (reference_ == "RESTRICTED") bQvvA.reset();
    JoovvAA->write(psio_, PSIF_DFOCC_INTS);
    JoovvAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    JoovvBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (oo|vv)", noccB, noccB, nvirB, nvirB));
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|oo)", nQ, noccB * noccB));
    bQvvB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|vv)", nQ, nvirB, nvirB));
    bQooB->read(psio_, PSIF_DFOCC_INTS);
    bQvvB->read(psio_, PSIF_DFOCC_INTS, true, true);
    JoovvBB->gemm(true, false, bQooB, bQvvB, 1.0, 0.0);
    JoovvBB->write(psio_, PSIF_DFOCC_INTS);
    JoovvBB.reset();

    // AB spin case
    JoovvAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|vv)", noccA, noccA, nvirB, nvirB));
    JoovvAB->gemm(true, false, bQooA, bQvvB, 1.0, 0.0);
    bQooA.reset();
    bQvvB.reset();
    JoovvAB->write(psio_, PSIF_DFOCC_INTS);
    JoovvAB.reset();

    // BA spin case
    JvvooAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (VV|oo)",nvirA, nvirA, noccB, noccB));
    JvvooAB->gemm(true, false, bQvvA, bQooB, 1.0, 0.0);
    bQooB.reset();
    bQvvA.reset();
    JvvooAB->write(psio_, PSIF_DFOCC_INTS);
    JvvooAB.reset();
 }
    timer_off("Build (oo|vv)");
}// end tei_oovv_chem

//=======================================================
//          (ia|jb)
//=======================================================
void DFOCC::tei_iajb_chem()
{
    timer_on("Build (ia|jb)");

    // AA spin case
    JiajbAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA * navirA));
    bQiaA->read(psio_, PSIF_DFOCC_INTS);
    JiajbAA->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    if (reference_ == "RESTRICTED") bQiaA.reset();
    JiajbAA->write(psio_, PSIF_DFOCC_INTS);
    JiajbAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    JiajbBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB));
    bQiaB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ia)", nQ, naoccB * navirB));
    bQiaB->read(psio_, PSIF_DFOCC_INTS);
    JiajbBB->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    JiajbBB->write(psio_, PSIF_DFOCC_INTS);
    JiajbBB.reset();

    // AB spin case
    JiajbAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|jb)", naoccA, navirA, naoccB, navirB));
    JiajbAB->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    bQiaA.reset();
    bQiaB.reset();
    JiajbAB->write(psio_, PSIF_DFOCC_INTS);
    JiajbAB.reset();
 }
    timer_off("Build (ia|jb)");
}// end tei_iajb_chem

//=======================================================
//          (OV|OV)
//=======================================================
void DFOCC::tei_ovov_chem()
{
    timer_on("Build (ov|ov)");
    // AA spin case
    JovovAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OV|OV)", noccA, nvirA, noccA, nvirA));
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA * nvirA));
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    JovovAA->gemm(true, false, bQovA, bQovA, 1.0, 0.0);
    if (reference_ == "RESTRICTED") bQovA.reset();
    JovovAA->write(psio_, PSIF_DFOCC_INTS);
    JovovAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    JovovBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ov|ov)", noccB, nvirB, noccB, nvirB));
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ov)", nQ, noccB * nvirB));
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    JovovBB->gemm(true, false, bQovB, bQovB, 1.0, 0.0);
    JovovBB->write(psio_, PSIF_DFOCC_INTS);
    JovovBB.reset();

    // AB spin case
    JovovAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OV|ov)", noccA, nvirA, noccB, nvirB));
    JovovAB->gemm(true, false, bQovA, bQovB, 1.0, 0.0);
    bQovA.reset();
    bQovB.reset();
    JovovAB->write(psio_, PSIF_DFOCC_INTS);
    JovovAB.reset();
 }
    timer_off("Build (ov|ov)");
}// end tei_ovov_chem

//=======================================================
//          <IJ|KL>
//=======================================================
void DFOCC::tei_ijkl_phys()
{
    timer_on("Build <ij|kl>");
    // AA spin case
    IijklAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|KL>", naoccA, naoccA, naoccA, naoccA));
    JijklAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|KL)", naoccA, naoccA, naoccA, naoccA));
    JijklAA->read(psio_, PSIF_DFOCC_INTS);
    IijklAA->sort(1324, JijklAA, 1.0, 0.0);
    JijklAA.reset();
    IijklAA->write(psio_, PSIF_DFOCC_INTS);
    IijklAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    IijklBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|kl>", naoccB, naoccB, naoccB, naoccB));
    JijklBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ij|kl)", naoccB, naoccB, naoccB, naoccB));
    JijklBB->read(psio_, PSIF_DFOCC_INTS);
    IijklBB->sort(1324, JijklBB, 1.0, 0.0);
    JijklBB.reset();
    IijklBB->write(psio_, PSIF_DFOCC_INTS);
    IijklBB.reset();

    // AB spin case
    IijklAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ij|Kl>", naoccA, naoccB, naoccA, naoccB));
    JijklAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|kl)", naoccA, naoccA, naoccB, naoccB));
    JijklAB->read(psio_, PSIF_DFOCC_INTS);
    IijklAB->sort(1324, JijklAB, 1.0, 0.0);
    JijklAB.reset();
    IijklAB->write(psio_, PSIF_DFOCC_INTS);
    IijklAB.reset();
 }
    timer_off("Build <ij|kl>");
}// end tei_ijkl_phys

//=======================================================
//          <OO|OO>
//=======================================================
void DFOCC::tei_oooo_phys()
{
    timer_on("Build <ij|kl>");
    // AA spin case
    IooooAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <OO|OO>", noccA, noccA, noccA, noccA));
    JooooAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|OO)", noccA, noccA, noccA, noccA));
    JooooAA->read(psio_, PSIF_DFOCC_INTS);
    IooooAA->sort(1324, JooooAA, 1.0, 0.0);
    JooooAA.reset();
    IooooAA->write(psio_, PSIF_DFOCC_INTS);
    IooooAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    IooooBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <oo|oo>", noccB, noccB, noccB, noccB));
    JooooBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (oo|oo)", noccB, noccB, noccB, noccB));
    JooooBB->read(psio_, PSIF_DFOCC_INTS);
    IooooBB->sort(1324, JooooBB, 1.0, 0.0);
    JooooBB.reset();
    IooooBB->write(psio_, PSIF_DFOCC_INTS);
    IooooBB.reset();

    // AB spin case
    IooooAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Oo|Oo>", noccA, noccB, noccA, noccB));
    JooooAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|oo)", noccA, noccA, noccB, noccB));
    JooooAB->read(psio_, PSIF_DFOCC_INTS);
    IooooAB->sort(1324, JooooAB, 1.0, 0.0);
    JooooAB.reset();
    IooooAB->write(psio_, PSIF_DFOCC_INTS);
    IooooAB.reset();
 }
    timer_off("Build <ij|kl>");
}// end tei_oooo_phys

//=======================================================
//          <IJ|KA>
//=======================================================
void DFOCC::tei_ijka_phys()
{
    timer_on("Build <ij|ka>");
    // AA spin case
    IijkaAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|KA>", naoccA, naoccA, naoccA, navirA));
    JijkaAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|KA)", naoccA, naoccA, naoccA, navirA));
    JijkaAA->read(psio_, PSIF_DFOCC_INTS);
    IijkaAA->sort(1324, JijkaAA, 1.0, 0.0);
    JijkaAA.reset();
    IijkaAA->write(psio_, PSIF_DFOCC_INTS);
    IijkaAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    IijkaBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|ka>", naoccB, naoccB, naoccB, navirB));
    JijkaBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ij|ka)", naoccB, naoccB, naoccB, navirB));
    JijkaBB->read(psio_, PSIF_DFOCC_INTS);
    IijkaBB->sort(1324, JijkaBB, 1.0, 0.0);
    JijkaBB.reset();
    IijkaBB->write(psio_, PSIF_DFOCC_INTS);
    IijkaBB.reset();

    // AB spin case
    IijkaAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ij|Ka>", naoccA, naoccB, naoccA, navirB));
    JijkaAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|ka)", naoccA, naoccA, naoccB, navirB));
    JijkaAB->read(psio_, PSIF_DFOCC_INTS);
    IijkaAB->sort(1324, JijkaAB, 1.0, 0.0);
    JijkaAB.reset();
    IijkaAB->write(psio_, PSIF_DFOCC_INTS);
    IijkaAB.reset();

    // BA spin case
    IijakAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ij|Ak>", naoccA, naoccB, navirA, naoccB));
    JiajkAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|jk)", naoccA, navirA, naoccB, naoccB));
    JiajkAB->read(psio_, PSIF_DFOCC_INTS);
    IijakAB->sort(1324, JiajkAB, 1.0, 0.0);
    JiajkAB.reset();
    IijakAB->write(psio_, PSIF_DFOCC_INTS);
    IijakAB.reset();
 }
    timer_off("Build <ij|ka>");
}// end tei_ijka_phys

//=======================================================
//          <OO|OV>
//=======================================================
void DFOCC::tei_ooov_phys()
{
    timer_on("Build <ij|ka>");
    // AA spin case
    IooovAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <OO|OV>", noccA, noccA, noccA, nvirA));
    JooovAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|OV)", noccA, noccA, noccA, nvirA));
    JooovAA->read(psio_, PSIF_DFOCC_INTS);
    IooovAA->sort(1324, JooovAA, 1.0, 0.0);
    JooovAA.reset();
    IooovAA->write(psio_, PSIF_DFOCC_INTS);
    IooovAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    IooovBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <oo|ov>", noccB, noccB, noccB, nvirB));
    JooovBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (oo|ov)", noccB, noccB, noccB, nvirB));
    JooovBB->read(psio_, PSIF_DFOCC_INTS);
    IooovBB->sort(1324, JooovBB, 1.0, 0.0);
    JooovBB.reset();
    IooovBB->write(psio_, PSIF_DFOCC_INTS);
    IooovBB.reset();

    // AB spin case
    IooovAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Oo|Ov>", noccA, noccB, noccA, nvirB));
    JooovAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|ov)", noccA, noccA, noccB, nvirB));
    JooovAB->read(psio_, PSIF_DFOCC_INTS);
    IooovAB->sort(1324, JooovAB, 1.0, 0.0);
    JooovAB.reset();
    IooovAB->write(psio_, PSIF_DFOCC_INTS);
    IooovAB.reset();

    // BA spin case
    IoovoAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Oo|Vo>", noccA, noccB, nvirA, noccB));
    JovooAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OV|oo)", noccA, nvirA, noccB, noccB));
    JovooAB->read(psio_, PSIF_DFOCC_INTS);
    IoovoAB->sort(1324, JovooAB, 1.0, 0.0);
    JovooAB.reset();
    IoovoAB->write(psio_, PSIF_DFOCC_INTS);
    IoovoAB.reset();
 }
    timer_off("Build <ij|ka>");
}// end tei_ooov_phys

//=======================================================
//          <IJ|AB>
//=======================================================
void DFOCC::tei_ijab_phys()
{
    timer_on("Build <ij|ab>");

    // AA spin case
    IijabAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA));
    JiajbAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    JiajbAA->read(psio_, PSIF_DFOCC_INTS);
    IijabAA->sort(1324, JiajbAA, 1.0, 0.0);
    JiajbAA.reset();
    IijabAA->write(psio_, PSIF_DFOCC_INTS);
    IijabAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    IijabBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|ab>", naoccB, naoccB, navirB, navirB));
    JiajbBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB));
    JiajbBB->read(psio_, PSIF_DFOCC_INTS);
    IijabBB->sort(1324, JiajbBB, 1.0, 0.0);
    JiajbBB.reset();
    IijabBB->write(psio_, PSIF_DFOCC_INTS);
    IijabBB.reset();

    // AB spin case
    IijabAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    JiajbAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|jb)", naoccA, navirA, naoccB, navirB));
    JiajbAB->read(psio_, PSIF_DFOCC_INTS);
    IijabAB->sort(1324, JiajbAB, 1.0, 0.0);
    JiajbAB.reset();
    IijabAB->write(psio_, PSIF_DFOCC_INTS);
    IijabAB.reset();
 }
    timer_off("Build <ij|ab>");
}// end tei_ijab_phys

//=======================================================
//          <OO|VV>
//=======================================================
void DFOCC::tei_oovv_phys()
{
    timer_on("Build <ij|ab>");
    // AA spin case
    IoovvAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <OO|VV>", noccA, noccA, nvirA, nvirA));
    JovovAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OV|OV)", noccA, nvirA, noccA, nvirA));
    JovovAA->read(psio_, PSIF_DFOCC_INTS);
    IoovvAA->sort(1324, JovovAA, 1.0, 0.0);
    JovovAA.reset();
    IoovvAA->write(psio_, PSIF_DFOCC_INTS);
    IoovvAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    IoovvBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <oo|vv>", noccB, noccB, nvirB, nvirB));
    JovovBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ov|ov)", noccB, nvirB, noccB, nvirB));
    JovovBB->read(psio_, PSIF_DFOCC_INTS);
    IoovvBB->sort(1324, JovovBB, 1.0, 0.0);
    JovovBB.reset();
    IoovvBB->write(psio_, PSIF_DFOCC_INTS);
    IoovvBB.reset();

    // AB spin case
    IoovvAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Oo|Vv>", noccA, noccB, nvirA, nvirB));
    JovovAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OV|ov)", noccA, nvirA, noccB, nvirB));
    JovovAB->read(psio_, PSIF_DFOCC_INTS);
    IoovvAB->sort(1324, JovovAB, 1.0, 0.0);
    JovovAB.reset();
    IoovvAB->write(psio_, PSIF_DFOCC_INTS);
    IoovvAB.reset();
 }
    timer_off("Build <ij|ab>");
}// end tei_oovv_phys

//=======================================================
//          <IA|JB>
//=======================================================
void DFOCC::tei_iajb_phys()
{
    timer_on("Build <ia|jb>");
    // AA spin case
    IiajbAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IA|JB>", naoccA, navirA, naoccA, navirA));
    JijabAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|AB)", naoccA, naoccA, navirA, navirA));
    JijabAA->read(psio_, PSIF_DFOCC_INTS);
    IiajbAA->sort(1324, JijabAA, 1.0, 0.0);
    JijabAA.reset();
    IiajbAA->write(psio_, PSIF_DFOCC_INTS);
    IiajbAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    IiajbBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ia|jb>", naoccB, navirB, naoccB, navirB));
    JijabBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ij|ab)", naoccB, naoccB, navirB, navirB));
    JijabBB->read(psio_, PSIF_DFOCC_INTS);
    IiajbBB->sort(1324, JijabBB, 1.0, 0.0);
    JijabBB.reset();
    IiajbBB->write(psio_, PSIF_DFOCC_INTS);
    IiajbBB.reset();

    // AB spin case
    IiajbAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ia|Jb>", naoccA, navirB, naoccA, navirB));
    JijabAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|ab)", naoccA, naoccA, navirB, navirB));
    JijabAB->read(psio_, PSIF_DFOCC_INTS);
    IiajbAB->sort(1324, JijabAB, 1.0, 0.0);
    JijabAB.reset();
    IiajbAB->write(psio_, PSIF_DFOCC_INTS);
    IiajbAB.reset();

    // BA spin case
    IaibjAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ai|Bj>", navirA, naoccB, navirA, naoccB));
    JabijAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (AB|ij)", navirA, navirA, naoccB, naoccB));
    JabijAB->read(psio_, PSIF_DFOCC_INTS);
    IaibjAB->sort(1324, JabijAB, 1.0, 0.0);
    JabijAB.reset();
    IaibjAB->write(psio_, PSIF_DFOCC_INTS);
    IaibjAB.reset();
 }
    timer_off("Build <ia|jb>");
}// end tei_iajb_phys

//=======================================================
//          <OV|OV>
//=======================================================
void DFOCC::tei_ovov_phys()
{
    timer_on("Build <ia|jb>");
    // AA spin case
    IovovAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <OV|OV>", noccA, nvirA, noccA, nvirA));
    JoovvAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|VV)", noccA, noccA, nvirA, nvirA));
    JoovvAA->read(psio_, PSIF_DFOCC_INTS);
    IovovAA->sort(1324, JoovvAA, 1.0, 0.0);
    JoovvAA.reset();
    IovovAA->write(psio_, PSIF_DFOCC_INTS);
    IovovAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    IovovBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ov|ov>", noccB, nvirB, noccB, nvirB));
    JoovvBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (oo|vv)", noccB, noccB, nvirB, nvirB));
    JoovvBB->read(psio_, PSIF_DFOCC_INTS);
    IovovBB->sort(1324, JoovvBB, 1.0, 0.0);
    JoovvBB.reset();
    IovovBB->write(psio_, PSIF_DFOCC_INTS);
    IovovBB.reset();

    // AB spin case
    IovovAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ov|Ov>", noccA, nvirB, noccA, nvirB));
    JoovvAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|vv)", noccA, noccA, nvirB, nvirB));
    JoovvAB->read(psio_, PSIF_DFOCC_INTS);
    IovovAB->sort(1324, JoovvAB, 1.0, 0.0);
    JoovvAB.reset();
    IovovAB->write(psio_, PSIF_DFOCC_INTS);
    IovovAB.reset();

    // BA spin case
    IvovoAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Vo|Vo>", nvirA, noccB, nvirA, noccB));
    JvvooAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (VV|oo)", nvirA, nvirA, noccB, noccB));
    JvvooAB->read(psio_, PSIF_DFOCC_INTS);
    IvovoAB->sort(1324, JvvooAB, 1.0, 0.0);
    JvvooAB.reset();
    IvovoAB->write(psio_, PSIF_DFOCC_INTS);
    IvovoAB.reset();
 }
    timer_off("Build <ia|jb>");
}// end tei_ovov_phys

//=======================================================
//          <IJ||KL>
//=======================================================
void DFOCC::tei_ijkl_anti_symm()
{
    timer_on("Build <ij||kl>");
    // AA spin case
    AIijklAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ||KL>", naoccA, naoccA, naoccA, naoccA));
    IijklAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|KL>", naoccA, naoccA, naoccA, naoccA));
    IijklAA->read(psio_, PSIF_DFOCC_INTS);
    AIijklAA->sort(1243, IijklAA, 1.0, 0.0);
    AIijklAA->scale(-1.0);
    AIijklAA->add(IijklAA);
    IijklAA.reset();
    AIijklAA->write(psio_, PSIF_DFOCC_INTS);
    AIijklAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    AIijklBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij||kl>", naoccB, naoccB, naoccB, naoccB));
    IijklBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|kl>", naoccB, naoccB, naoccB, naoccB));
    IijklBB->read(psio_, PSIF_DFOCC_INTS);
    AIijklBB->sort(1243, IijklBB, 1.0, 0.0);
    AIijklBB->scale(-1.0);
    AIijklBB->add(IijklBB);
    IijklBB.reset();
    AIijklBB->write(psio_, PSIF_DFOCC_INTS);
    AIijklBB.reset();
 }
    timer_off("Build <ij||kl>");
}// end tei_ijkl_anti_symm

//=======================================================
//          <OO||OO>
//=======================================================
void DFOCC::tei_oooo_anti_symm()
{
    timer_on("Build <ij||kl>");
    // AA spin case
    AIooooAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <OO||OO>", noccA, noccA, noccA, noccA));
    IooooAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <OO|OO>", noccA, noccA, noccA, noccA));
    IooooAA->read(psio_, PSIF_DFOCC_INTS);
    AIooooAA->sort(1243, IooooAA, 1.0, 0.0);
    AIooooAA->scale(-1.0);
    AIooooAA->add(IooooAA);
    IooooAA.reset();
    AIooooAA->write(psio_, PSIF_DFOCC_INTS);
    AIooooAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    AIooooBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <oo||oo>", noccB, noccB, noccB, noccB));
    IooooBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <oo|oo>", noccB, noccB, noccB, noccB));
    IooooBB->read(psio_, PSIF_DFOCC_INTS);
    AIooooBB->sort(1243, IooooBB, 1.0, 0.0);
    AIooooBB->scale(-1.0);
    AIooooBB->add(IooooBB);
    IooooBB.reset();
    AIooooBB->write(psio_, PSIF_DFOCC_INTS);
    AIooooBB.reset();
 }
    timer_off("Build <ij||kl>");
}// end tei_oooo_anti_symm

//=======================================================
//          <IJ||KA>
//=======================================================
void DFOCC::tei_ijka_anti_symm()
{
    timer_on("Build <ij||ka>");
    // <ij||ka> = <ij|ka> - <ji|ka>
    // AA spin case
    AIijkaAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ||KA>", naoccA, naoccA, naoccA, navirA));
    IijkaAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|KA>", naoccA, naoccA, naoccA, navirA));
    IijkaAA->read(psio_, PSIF_DFOCC_INTS);
    AIijkaAA->sort(2134, IijkaAA, 1.0, 0.0);
    AIijkaAA->scale(-1.0);
    AIijkaAA->add(IijkaAA);
    IijkaAA.reset();
    AIijkaAA->write(psio_, PSIF_DFOCC_INTS);
    AIijkaAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    AIijkaBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij||ka>", naoccB, naoccB, naoccB, navirB));
    IijkaBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|ka>", naoccB, naoccB, naoccB, navirB));
    IijkaBB->read(psio_, PSIF_DFOCC_INTS);
    AIijkaBB->sort(2134, IijkaBB, 1.0, 0.0);
    AIijkaBB->scale(-1.0);
    AIijkaBB->add(IijkaBB);
    IijkaBB.reset();
    AIijkaBB->write(psio_, PSIF_DFOCC_INTS);
    AIijkaBB.reset();
 }
    timer_off("Build <ij||ka>");
}// end tei_ijka_anti_symm

//=======================================================
//          <OO||OV>
//=======================================================
void DFOCC::tei_ooov_anti_symm()
{
    timer_on("Build <ij||ka>");
    // <ij||ka> = <ij|ka> - <ji|ka>
    // AA spin case
    AIooovAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <OO||OV>", noccA, noccA, noccA, nvirA));
    IooovAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <OO|OV>", noccA, noccA, noccA, nvirA));
    IooovAA->read(psio_, PSIF_DFOCC_INTS);
    AIooovAA->sort(2134, IooovAA, 1.0, 0.0);
    AIooovAA->scale(-1.0);
    AIooovAA->add(IooovAA);
    IooovAA.reset();
    AIooovAA->write(psio_, PSIF_DFOCC_INTS);
    AIooovAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    AIooovBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <oo||ov>", noccB, noccB, noccB, nvirB));
    IooovBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <oo|ov>", noccB, noccB, noccB, nvirB));
    IooovBB->read(psio_, PSIF_DFOCC_INTS);
    AIooovBB->sort(2134, IooovBB, 1.0, 0.0);
    AIooovBB->scale(-1.0);
    AIooovBB->add(IooovBB);
    IooovBB.reset();
    AIooovBB->write(psio_, PSIF_DFOCC_INTS);
    AIooovBB.reset();
 }
    timer_off("Build <ij||ka>");
}// end tei_ooov_anti_symm

//=======================================================
//          <IJ||AB>
//=======================================================
void DFOCC::tei_ijab_anti_symm()
{
    timer_on("Build <ij||ab>");

    // AA spin case
    AIijabAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ||AB>", naoccA, naoccA, navirA, navirA));
    IijabAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA));
    IijabAA->read(psio_, PSIF_DFOCC_INTS);
    AIijabAA->sort(1243, IijabAA, 1.0, 0.0);
    AIijabAA->scale(-1.0);
    AIijabAA->add(IijabAA);
    IijabAA.reset();
    AIijabAA->write(psio_, PSIF_DFOCC_INTS);
    AIijabAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    AIijabBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij||ab>", naoccB, naoccB, navirB, navirB));
    IijabBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|ab>", naoccB, naoccB, navirB, navirB));
    IijabBB->read(psio_, PSIF_DFOCC_INTS);
    AIijabBB->sort(1243, IijabBB, 1.0, 0.0);
    AIijabBB->scale(-1.0);
    AIijabBB->add(IijabBB);
    IijabBB.reset();
    AIijabBB->write(psio_, PSIF_DFOCC_INTS);
    AIijabBB.reset();
 }
    timer_off("Build <ij||ab>");
}// end tei_ijab_anti_symm

//=======================================================
//          <OO||VV>
//=======================================================
void DFOCC::tei_oovv_anti_symm()
{
    timer_on("Build <ij||ab>");
    // AA spin case
    AIoovvAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <OO||VV>", noccA, noccA, nvirA, nvirA));
    IoovvAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <OO|VV>", noccA, noccA, nvirA, nvirA));
    IoovvAA->read(psio_, PSIF_DFOCC_INTS);
    AIoovvAA->sort(1243, IoovvAA, 1.0, 0.0);
    AIoovvAA->scale(-1.0);
    AIoovvAA->add(IoovvAA);
    IoovvAA.reset();
    AIoovvAA->write(psio_, PSIF_DFOCC_INTS);
    AIoovvAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    AIoovvBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <oo||vv>", noccB, noccB, nvirB, nvirB));
    IoovvBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <oo|vv>", noccB, noccB, nvirB, nvirB));
    IoovvBB->read(psio_, PSIF_DFOCC_INTS);
    AIoovvBB->sort(1243, IoovvBB, 1.0, 0.0);
    AIoovvBB->scale(-1.0);
    AIoovvBB->add(IoovvBB);
    IoovvBB.reset();
    AIoovvBB->write(psio_, PSIF_DFOCC_INTS);
    AIoovvBB.reset();
 }
    timer_off("Build <ij||ab>");
}// end tei_oovv_anti_symm

//=======================================================
//          <IA||JB>
//=======================================================
void DFOCC::tei_iajb_anti_symm()
{
    timer_on("Build <ia||jb>");
    // <ia||jb> = <ia|jb> - <ia|bj> = <ia|jb> - (ib|ja)
    // AA spin case
    AIiajbAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IA||JB>", naoccA, navirA, naoccA, navirA));
    JiajbAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    JiajbAA->read(psio_, PSIF_DFOCC_INTS);
    AIiajbAA->sort(1432, JiajbAA, 1.0, 0.0);
    JiajbAA.reset();
    AIiajbAA->scale(-1.0);
    IiajbAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IA|JB>", naoccA, navirA, naoccA, navirA));
    IiajbAA->read(psio_, PSIF_DFOCC_INTS);
    AIiajbAA->add(IiajbAA);
    IiajbAA.reset();
    AIiajbAA->write(psio_, PSIF_DFOCC_INTS);
    AIiajbAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    AIiajbBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ia||jb>", naoccB, navirB, naoccB, navirB));
    JiajbBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB));
    JiajbBB->read(psio_, PSIF_DFOCC_INTS);
    AIiajbBB->sort(1432, JiajbBB, 1.0, 0.0);
    JiajbBB.reset();
    AIiajbBB->scale(-1.0);
    IiajbBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ia|jb>", naoccB, navirB, naoccB, navirB));
    IiajbBB->read(psio_, PSIF_DFOCC_INTS);
    AIiajbBB->add(IiajbBB);
    IiajbBB.reset();
    AIiajbBB->write(psio_, PSIF_DFOCC_INTS);
    AIiajbBB.reset();
    //outfile->Printf("\tI am here.\n");
 }
    timer_off("Build <ia||jb>");
}// end tei_iajb_anti_symm

//=======================================================
//          <OV||OV>
//=======================================================
void DFOCC::tei_ovov_anti_symm()
{
    timer_on("Build <ia||jb>");
    // <ia||jb> = <ia|jb> - <ia|bj> = <ia|jb> - (ib|ja)
    // AA spin case
    AIovovAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <OV||OV>", noccA, nvirA, noccA, nvirA));
    JovovAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OV|OV)", noccA, nvirA, noccA, nvirA));
    JovovAA->read(psio_, PSIF_DFOCC_INTS);
    AIovovAA->sort(1432, JovovAA, 1.0, 0.0);
    JovovAA.reset();
    AIovovAA->scale(-1.0);
    IovovAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <OV|OV>", noccA, nvirA, noccA, nvirA));
    IovovAA->read(psio_, PSIF_DFOCC_INTS);
    AIovovAA->add(IovovAA);
    IovovAA.reset();
    AIovovAA->write(psio_, PSIF_DFOCC_INTS);
    AIovovAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    AIovovBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ov||ov>", noccB, nvirB, noccB, nvirB));
    JovovBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ov|ov)", noccB, nvirB, noccB, nvirB));
    JovovBB->read(psio_, PSIF_DFOCC_INTS);
    AIovovBB->sort(1432, JovovBB, 1.0, 0.0);
    JovovBB.reset();
    AIovovBB->scale(-1.0);
    IovovBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ov|ov>", noccB, nvirB, noccB, nvirB));
    IovovBB->read(psio_, PSIF_DFOCC_INTS);
    AIovovBB->add(IovovBB);
    IovovBB.reset();
    AIovovBB->write(psio_, PSIF_DFOCC_INTS);
    AIovovBB.reset();
    //outfile->Printf("\tI am here.\n");
 }
    timer_off("Build <ia||jb>");
}// end tei_ovov_anti_symm


}} // End Namespaces
