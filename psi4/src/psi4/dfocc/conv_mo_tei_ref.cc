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

void DFOCC::tei_oooo_chem_ref()
{
    timer_on("Build (oo|oo)");
    // AA spin case
    JooooAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|OO)", noccA, noccA, noccA, noccA));
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA * noccA));
    bQooA->read(psio_, PSIF_DFOCC_INTS);
    JooooAA->gemm(true, false, bQooA, bQooA, 1.0, 0.0);
    if (reference_ == "RESTRICTED") bQooA.reset();
    JooooAA->write(psio_, PSIF_DFOCC_INTS);
    JooooAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    JooooBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (oo|oo)", noccB, noccB, noccB, noccB));
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB * noccB));
    bQooB->read(psio_, PSIF_DFOCC_INTS);
    JooooBB->gemm(true, false, bQooB, bQooB, 1.0, 0.0);
    JooooBB->write(psio_, PSIF_DFOCC_INTS);
    JooooBB.reset();

    // AB spin case
    JooooAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|oo)", noccA, noccA, noccB, noccB));
    JooooAB->gemm(true, false, bQooA, bQooB, 1.0, 0.0);
    bQooA.reset();
    bQooB.reset();
    JooooAB->write(psio_, PSIF_DFOCC_INTS);
    JooooAB.reset();
 }
    timer_off("Build (oo|oo)");
}// end tei_oooo_chem_ref

//=======================================================
//          (OO|OV)
//=======================================================
void DFOCC::tei_ooov_chem_ref()
{
    timer_on("Build (oo|ov)");
    // AA spin case
    JooovAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|OV)", noccA, noccA, noccA, nvirA));
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA * noccA));
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA * nvirA));
    bQooA->read(psio_, PSIF_DFOCC_INTS);
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    JooovAA->gemm(true, false, bQooA, bQovA, 1.0, 0.0);
    if (reference_ == "RESTRICTED") bQooA.reset();
    if (reference_ == "RESTRICTED") bQovA.reset();
    JooovAA->write(psio_, PSIF_DFOCC_INTS);
    JooovAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    JooovBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (oo|ov)", noccB, noccB, noccB, nvirB));
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB * noccB));
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB * nvirB));
    bQooB->read(psio_, PSIF_DFOCC_INTS);
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    JooovBB->gemm(true, false, bQooB, bQovB, 1.0, 0.0);
    JooovBB->write(psio_, PSIF_DFOCC_INTS);
    JooovBB.reset();

    // AB spin case
    JooovAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|ov)", noccA, noccA, noccB, nvirB));
    JooovAB->gemm(true, false, bQooA, bQovB, 1.0, 0.0);
    bQooA.reset();
    bQovB.reset();
    JooovAB->write(psio_, PSIF_DFOCC_INTS);
    JooovAB.reset();

    // BA spin case
    JovooAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|oo)", noccA, nvirA, noccB, noccB));
    JovooAB->gemm(true, false, bQovA, bQooB, 1.0, 0.0);
    bQooB.reset();
    bQovA.reset();
    JovooAB->write(psio_, PSIF_DFOCC_INTS);
    JovooAB.reset();
 }
    timer_off("Build (oo|ov)");
}// end tei_ooov_chem_ref

//=======================================================
//          (OO|VV)
//=======================================================
void DFOCC::tei_oovv_chem_ref()
{
    timer_on("Build (oo|vv)");
    // AA spin case
    JoovvAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|VV)", noccA, noccA, nvirA, nvirA));
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA * noccA));
    bQvvA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    bQooA->read(psio_, PSIF_DFOCC_INTS);
    bQvvA->read(psio_, PSIF_DFOCC_INTS, true, true);
    JoovvAA->gemm(true, false, bQooA, bQvvA, 1.0, 0.0);
    if (reference_ == "RESTRICTED") bQooA.reset();
    if (reference_ == "RESTRICTED") bQvvA.reset();
    JoovvAA->write(psio_, PSIF_DFOCC_INTS);
    JoovvAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    JoovvBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (oo|vv)", noccB, noccB, nvirB, nvirB));
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB * noccB));
    bQvvB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vv)", nQ_ref, nvirB, nvirB));
    bQooB->read(psio_, PSIF_DFOCC_INTS);
    bQvvB->read(psio_, PSIF_DFOCC_INTS, true, true);
    JoovvBB->gemm(true, false, bQooB, bQvvB, 1.0, 0.0);
    JoovvBB->write(psio_, PSIF_DFOCC_INTS);
    JoovvBB.reset();

    // AB spin case
    JoovvAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|vv)", noccA, noccA, nvirB, nvirB));
    JoovvAB->gemm(true, false, bQooA, bQvvB, 1.0, 0.0);
    bQooA.reset();
    bQvvB.reset();
    JoovvAB->write(psio_, PSIF_DFOCC_INTS);
    JoovvAB.reset();

    // BA spin case
    JvvooAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (VV|oo)", nvirA, nvirA, noccB, noccB));
    JvvooAB->gemm(true, false, bQvvA, bQooB, 1.0, 0.0);
    bQooB.reset();
    bQvvA.reset();
    JvvooAB->write(psio_, PSIF_DFOCC_INTS);
    JvvooAB.reset();
 }
    timer_off("Build (oo|vv)");
}// end tei_oovv_chem_ref

//=======================================================
//          (OV|OV)
//=======================================================
void DFOCC::tei_ovov_chem_ref()
{
    timer_on("Build (ov|ov)");
    // AA spin case
    JovovAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|OV)", noccA, nvirA, noccA, nvirA));
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA * nvirA));
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    JovovAA->gemm(true, false, bQovA, bQovA, 1.0, 0.0);
    if (reference_ == "RESTRICTED") bQovA.reset();
    JovovAA->write(psio_, PSIF_DFOCC_INTS);
    JovovAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    JovovBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (ov|ov)", noccB, nvirB, noccB, nvirB));
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB * nvirB));
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    JovovBB->gemm(true, false, bQovB, bQovB, 1.0, 0.0);
    JovovBB->write(psio_, PSIF_DFOCC_INTS);
    JovovBB.reset();

    // AB spin case
    JovovAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|ov)", noccA, nvirA, noccB, nvirB));
    JovovAB->gemm(true, false, bQovA, bQovB, 1.0, 0.0);
    bQovA.reset();
    bQovB.reset();
    JovovAB->write(psio_, PSIF_DFOCC_INTS);
    JovovAB.reset();
 }
    timer_off("Build (ov|ov)");
}// end tei_ovov_chem_ref

//=======================================================
//          <OO|OO>
//=======================================================
void DFOCC::tei_oooo_phys_ref()
{
    timer_on("Build <ij|kl>");

    // AA spin case
    IooooAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OO|OO>", noccA, noccA, noccA, noccA));
    JooooAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|OO)", noccA, noccA, noccA, noccA));
    JooooAA->read(psio_, PSIF_DFOCC_INTS);
    IooooAA->sort(1324, JooooAA, 1.0, 0.0);
    JooooAA.reset();
    IooooAA->write(psio_, PSIF_DFOCC_INTS);
    IooooAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    IooooBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <oo|oo>", noccB, noccB, noccB, noccB));
    JooooBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (oo|oo)", noccB, noccB, noccB, noccB));
    JooooBB->read(psio_, PSIF_DFOCC_INTS);
    IooooBB->sort(1324, JooooBB, 1.0, 0.0);
    JooooBB.reset();
    IooooBB->write(psio_, PSIF_DFOCC_INTS);
    IooooBB.reset();

    // AB spin case
    IooooAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <Oo|Oo>", noccA, noccB, noccA, noccB));
    JooooAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|oo)", noccA, noccA, noccB, noccB));
    JooooAB->read(psio_, PSIF_DFOCC_INTS);
    IooooAB->sort(1324, JooooAB, 1.0, 0.0);
    JooooAB.reset();
    IooooAB->write(psio_, PSIF_DFOCC_INTS);
    IooooAB.reset();
 }
    timer_off("Build <ij|kl>");
}// end tei_oooo_phys_ref

//=======================================================
//          <OO|OV>
//=======================================================
void DFOCC::tei_ooov_phys_ref()
{
    timer_on("Build <ij|ka>");
    // AA spin case
    IooovAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OO|OV>", noccA, noccA, noccA, nvirA));
    JooovAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|OV)", noccA, noccA, noccA, nvirA));
    JooovAA->read(psio_, PSIF_DFOCC_INTS);
    IooovAA->sort(1324, JooovAA, 1.0, 0.0);
    JooovAA.reset();
    IooovAA->write(psio_, PSIF_DFOCC_INTS);
    IooovAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    IooovBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <oo|ov>", noccB, noccB, noccB, nvirB));
    JooovBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (oo|ov)", noccB, noccB, noccB, nvirB));
    JooovBB->read(psio_, PSIF_DFOCC_INTS);
    IooovBB->sort(1324, JooovBB, 1.0, 0.0);
    JooovBB.reset();
    IooovBB->write(psio_, PSIF_DFOCC_INTS);
    IooovBB.reset();

    // AB spin case
    IooovAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <Oo|Ov>", noccA, noccB, noccA, nvirB));
    JooovAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|ov)", noccA, noccA, noccB, nvirB));
    JooovAB->read(psio_, PSIF_DFOCC_INTS);
    IooovAB->sort(1324, JooovAB, 1.0, 0.0);
    JooovAB.reset();
    IooovAB->write(psio_, PSIF_DFOCC_INTS);
    IooovAB.reset();

    // BA spin case
    IoovoAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <Oo|Vo>", noccA, noccB, nvirA, noccB));
    JovooAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|oo)", noccA, nvirA, noccB, noccB));
    JovooAB->read(psio_, PSIF_DFOCC_INTS);
    IoovoAB->sort(1324, JovooAB, 1.0, 0.0);
    JovooAB.reset();
    IoovoAB->write(psio_, PSIF_DFOCC_INTS);
    IoovoAB.reset();
 }
    timer_off("Build <ij|ka>");
}// end tei_ooov_phys_ref

//=======================================================
//          <OO|VV>
//=======================================================
void DFOCC::tei_oovv_phys_ref()
{
    timer_on("Build <ij|ab>");

    // AA spin case
    IoovvAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OO|VV>", noccA, noccA, nvirA, nvirA));
    JovovAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|OV)", noccA, nvirA, noccA, nvirA));
    JovovAA->read(psio_, PSIF_DFOCC_INTS);
    IoovvAA->sort(1324, JovovAA, 1.0, 0.0);
    JovovAA.reset();
    IoovvAA->write(psio_, PSIF_DFOCC_INTS);
    IoovvAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    IoovvBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <oo|vv>", noccB, noccB, nvirB, nvirB));
    JovovBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (ov|ov)", noccB, nvirB, noccB, nvirB));
    JovovBB->read(psio_, PSIF_DFOCC_INTS);
    IoovvBB->sort(1324, JovovBB, 1.0, 0.0);
    JovovBB.reset();
    IoovvBB->write(psio_, PSIF_DFOCC_INTS);
    IoovvBB.reset();

    // AB spin case
    IoovvAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <Oo|Vv>", noccA, noccB, nvirA, nvirB));
    JovovAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|ov)", noccA, nvirA, noccB, nvirB));
    JovovAB->read(psio_, PSIF_DFOCC_INTS);
    IoovvAB->sort(1324, JovovAB, 1.0, 0.0);
    JovovAB.reset();
    IoovvAB->write(psio_, PSIF_DFOCC_INTS);
    IoovvAB.reset();
 }
    timer_off("Build <ij|ab>");
}// end tei_oovv_phys_ref

//=======================================================
//          <OV|OV>
//=======================================================
void DFOCC::tei_ovov_phys_ref()
{
    timer_on("Build <ia|jb>");
    // AA spin case
    IovovAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OV|OV>", noccA, nvirA, noccA, nvirA));
    JoovvAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|VV)", noccA, noccA, nvirA, nvirA));
    JoovvAA->read(psio_, PSIF_DFOCC_INTS);
    IovovAA->sort(1324, JoovvAA, 1.0, 0.0);
    JoovvAA.reset();
    IovovAA->write(psio_, PSIF_DFOCC_INTS);
    IovovAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    IovovBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <ov|ov>", noccB, nvirB, noccB, nvirB));
    JoovvBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (oo|vv)", noccB, noccB, nvirB, nvirB));
    JoovvBB->read(psio_, PSIF_DFOCC_INTS);
    IovovBB->sort(1324, JoovvBB, 1.0, 0.0);
    JoovvBB.reset();
    IovovBB->write(psio_, PSIF_DFOCC_INTS);
    IovovBB.reset();

    // AB spin case
    IovovAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <Ov|Ov>", noccA, nvirB, noccA, nvirB));
    JoovvAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|vv)", noccA, noccA, nvirB, nvirB));
    JoovvAB->read(psio_, PSIF_DFOCC_INTS);
    IovovAB->sort(1324, JoovvAB, 1.0, 0.0);
    JoovvAB.reset();
    IovovAB->write(psio_, PSIF_DFOCC_INTS);
    IovovAB.reset();

    // BA spin case
    IvovoAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <Vo|Vo>", nvirA, noccB, nvirA, noccB));
    JvvooAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (VV|oo)", nvirA, nvirA, noccB, noccB));
    JvvooAB->read(psio_, PSIF_DFOCC_INTS);
    IvovoAB->sort(1324, JvvooAB, 1.0, 0.0);
    JvvooAB.reset();
    IvovoAB->write(psio_, PSIF_DFOCC_INTS);
    IvovoAB.reset();
 }
    timer_off("Build <ia|jb>");
}// end tei_ovov_phys_ref

//=======================================================
//          <OO||OO>
//=======================================================
void DFOCC::tei_oooo_anti_symm_ref()
{
    timer_on("Build <ij||kl>");

    // AA spin case
    AIooooAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OO||OO>", noccA, noccA, noccA, noccA));
    IooooAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OO|OO>", noccA, noccA, noccA, noccA));
    IooooAA->read(psio_, PSIF_DFOCC_INTS);
    AIooooAA->sort(1243, IooooAA, 1.0, 0.0);
    AIooooAA->scale(-1.0);
    AIooooAA->add(IooooAA);
    IooooAA.reset();
    AIooooAA->write(psio_, PSIF_DFOCC_INTS);
    AIooooAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    AIooooBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <oo||oo>", noccB, noccB, noccB, noccB));
    IooooBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <oo|oo>", noccB, noccB, noccB, noccB));
    IooooBB->read(psio_, PSIF_DFOCC_INTS);
    AIooooBB->sort(1243, IooooBB, 1.0, 0.0);
    AIooooBB->scale(-1.0);
    AIooooBB->add(IooooBB);
    IooooBB.reset();
    AIooooBB->write(psio_, PSIF_DFOCC_INTS);
    AIooooBB.reset();
 }
    timer_off("Build <ij||kl>");
}// end tei_oooo_anti_symm_ref

//=======================================================
//          <OO||OV>
//=======================================================
void DFOCC::tei_ooov_anti_symm_ref()
{
    timer_on("Build <ij||ka>");
    // <ij||ka> = <ij|ka> - <ji|ka>
    // AA spin case
    AIooovAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OO||OV>", noccA, noccA, noccA, nvirA));
    IooovAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OO|OV>", noccA, noccA, noccA, nvirA));
    IooovAA->read(psio_, PSIF_DFOCC_INTS);
    AIooovAA->sort(2134, IooovAA, 1.0, 0.0);
    AIooovAA->scale(-1.0);
    AIooovAA->add(IooovAA);
    IooovAA.reset();
    AIooovAA->write(psio_, PSIF_DFOCC_INTS);
    AIooovAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    AIooovBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <oo||ov>", noccB, noccB, noccB, nvirB));
    IooovBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <oo|ov>", noccB, noccB, noccB, nvirB));
    IooovBB->read(psio_, PSIF_DFOCC_INTS);
    AIooovBB->sort(2134, IooovBB, 1.0, 0.0);
    AIooovBB->scale(-1.0);
    AIooovBB->add(IooovBB);
    IooovBB.reset();
    AIooovBB->write(psio_, PSIF_DFOCC_INTS);
    AIooovBB.reset();
 }
    timer_off("Build <ij||ka>");
}// end tei_ooov_anti_symm_ref

//=======================================================
//          <OO||VV>
//=======================================================
void DFOCC::tei_oovv_anti_symm_ref()
{
    timer_on("Build <ij||ab>");

    // AA spin case
    AIoovvAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OO||VV>", noccA, noccA, nvirA, nvirA));
    IoovvAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OO|VV>", noccA, noccA, nvirA, nvirA));
    IoovvAA->read(psio_, PSIF_DFOCC_INTS);
    AIoovvAA->sort(1243, IoovvAA, 1.0, 0.0);
    AIoovvAA->scale(-1.0);
    AIoovvAA->add(IoovvAA);
    IoovvAA.reset();
    AIoovvAA->write(psio_, PSIF_DFOCC_INTS);
    AIoovvAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    AIoovvBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <oo||vv>", noccB, noccB, nvirB, nvirB));
    IoovvBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <oo|vv>", noccB, noccB, nvirB, nvirB));
    IoovvBB->read(psio_, PSIF_DFOCC_INTS);
    AIoovvBB->sort(1243, IoovvBB, 1.0, 0.0);
    AIoovvBB->scale(-1.0);
    AIoovvBB->add(IoovvBB);
    IoovvBB.reset();
    AIoovvBB->write(psio_, PSIF_DFOCC_INTS);
    AIoovvBB.reset();
 }
    timer_off("Build <ij||ab>");
}// end tei_oovv_anti_symm_ref

//=======================================================
//          <OV||OV>
//=======================================================
void DFOCC::tei_ovov_anti_symm_ref()
{
    timer_on("Build <ia||jb>");
    // <ia||jb> = <ia|jb> - <ia|bj> = <ia|jb> - (ib|ja)
    // AA spin case
    AIovovAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OV||OV>", noccA, nvirA, noccA, nvirA));
    JovovAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|OV)", noccA, nvirA, noccA, nvirA));
    JovovAA->read(psio_, PSIF_DFOCC_INTS);
    AIovovAA->sort(1432, JovovAA, 1.0, 0.0);
    JovovAA.reset();
    AIovovAA->scale(-1.0);
    IovovAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OV|OV>", noccA, nvirA, noccA, nvirA));
    IovovAA->read(psio_, PSIF_DFOCC_INTS);
    AIovovAA->add(IovovAA);
    IovovAA.reset();
    AIovovAA->write(psio_, PSIF_DFOCC_INTS);
    AIovovAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    AIovovBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <ov||ov>", noccB, nvirB, noccB, nvirB));
    JovovBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (ov|ov)", noccB, nvirB, noccB, nvirB));
    JovovBB->read(psio_, PSIF_DFOCC_INTS);
    AIovovBB->sort(1432, JovovBB, 1.0, 0.0);
    JovovBB.reset();
    AIovovBB->scale(-1.0);
    IovovBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <ov|ov>", noccB, nvirB, noccB, nvirB));
    IovovBB->read(psio_, PSIF_DFOCC_INTS);
    AIovovBB->add(IovovBB);
    IovovBB.reset();
    AIovovBB->write(psio_, PSIF_DFOCC_INTS);
    AIovovBB.reset();
    //outfile->Printf("\tI am here.\n");
 }
    timer_off("Build <ia||jb>");
}// end tei_ovov_anti_symm_ref


}} // End Namespaces
