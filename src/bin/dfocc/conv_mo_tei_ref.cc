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

/** Standard library includes */
#include <libqt/qt.h>
#include "defines.h"
#include "dfocc.h"

using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace dfoccwave{
  
void DFOCC::tei_oooo_chem_ref()
{   
    timer_on("Build (oo|oo)");
    // AA spin case
    JooooAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|OO)", noccA, noccA, noccA, noccA));
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA * noccA));
    timer_on("I/O");
    bQooA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JooooAA->gemm(true, false, bQooA, bQooA, 1.0, 0.0);
    if (reference_ == "RESTRICTED") bQooA.reset();
    timer_on("I/O");
    JooooAA->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JooooAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    JooooBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (oo|oo)", noccB, noccB, noccB, noccB));
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB * noccB));
    timer_on("I/O");
    bQooB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JooooBB->gemm(true, false, bQooB, bQooB, 1.0, 0.0);
    timer_on("I/O");
    JooooBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JooooBB.reset();

    // AB spin case
    JooooAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|oo)", noccA, noccA, noccB, noccB));
    JooooAB->gemm(true, false, bQooA, bQooB, 1.0, 0.0);
    bQooA.reset();
    bQooB.reset();
    timer_on("I/O");
    JooooAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
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
    timer_on("I/O");
    bQooA->read(psio_, PSIF_DFOCC_INTS);
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JooovAA->gemm(true, false, bQooA, bQovA, 1.0, 0.0);
    if (reference_ == "RESTRICTED") bQooA.reset();
    if (reference_ == "RESTRICTED") bQovA.reset();
    timer_on("I/O");
    JooovAA->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JooovAA.reset();
    
 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    JooovBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (oo|ov)", noccB, noccB, noccB, nvirB));
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB * noccB));
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB * nvirB));
    timer_on("I/O");
    bQooB->read(psio_, PSIF_DFOCC_INTS);
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JooovBB->gemm(true, false, bQooB, bQovB, 1.0, 0.0);
    timer_on("I/O");
    JooovBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JooovBB.reset();

    // AB spin case
    JooovAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|ov)", noccA, noccA, noccB, nvirB));
    JooovAB->gemm(true, false, bQooA, bQovB, 1.0, 0.0);
    bQooA.reset();
    bQovB.reset();
    timer_on("I/O");
    JooovAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JooovAB.reset();

    // BA spin case
    JovooAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|oo)", noccA, nvirA, noccB, noccB));
    JovooAB->gemm(true, false, bQovA, bQooB, 1.0, 0.0);
    bQooB.reset();
    bQovA.reset();
    timer_on("I/O");
    JovooAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
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
    bQvvA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA * nvirA));
    timer_on("I/O");
    bQooA->read(psio_, PSIF_DFOCC_INTS);
    bQvvA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JoovvAA->gemm(true, false, bQooA, bQvvA, 1.0, 0.0);
    if (reference_ == "RESTRICTED") bQooA.reset();
    if (reference_ == "RESTRICTED") bQvvA.reset();
    timer_on("I/O");
    JoovvAA->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JoovvAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    JoovvBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (oo|vv)", noccB, noccB, nvirB, nvirB));
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB * noccB));
    bQvvB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vv)", nQ_ref, nvirB * nvirB));
    timer_on("I/O");
    bQooB->read(psio_, PSIF_DFOCC_INTS);
    bQvvB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JoovvBB->gemm(true, false, bQooB, bQvvB, 1.0, 0.0);
    timer_on("I/O");
    JoovvBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JoovvBB.reset();

    // AB spin case
    JoovvAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|vv)", noccA, noccA, nvirB, nvirB));
    JoovvAB->gemm(true, false, bQooA, bQvvB, 1.0, 0.0);
    bQooA.reset();
    bQvvB.reset();
    timer_on("I/O");
    JoovvAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JoovvAB.reset();

    // BA spin case
    JvvooAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (VV|oo)", nvirA, nvirA, noccB, noccB));
    JvvooAB->gemm(true, false, bQvvA, bQooB, 1.0, 0.0);
    bQooB.reset();
    bQvvA.reset();
    timer_on("I/O");
    JvvooAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
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
    timer_on("I/O");
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JovovAA->gemm(true, false, bQovA, bQovA, 1.0, 0.0);
    if (reference_ == "RESTRICTED") bQovA.reset();
    timer_on("I/O");
    JovovAA->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JovovAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    JovovBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (ov|ov)", noccB, nvirB, noccB, nvirB));
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB * nvirB));
    timer_on("I/O");
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JovovBB->gemm(true, false, bQovB, bQovB, 1.0, 0.0);
    timer_on("I/O");
    JovovBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JovovBB.reset();

    // AB spin case
    JovovAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|ov)", noccA, nvirA, noccB, nvirB));
    JovovAB->gemm(true, false, bQovA, bQovB, 1.0, 0.0);
    bQovA.reset();
    bQovB.reset();
    timer_on("I/O");
    JovovAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JovovAB.reset();
 }
    timer_off("Build (ov|ov)");
}// end tei_ovov_chem_ref 

//=======================================================
//          <OO|OO>
//=======================================================          
void DFOCC::tei_oooo_phys_ref()
{   
    timer_on("Build <ij|kl>)");

    // AA spin case
    IooooAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OO|OO>", noccA, noccA, noccA, noccA));
    JooooAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|OO)", noccA, noccA, noccA, noccA));
    timer_on("I/O");
    JooooAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IooooAA->sort(1324, JooooAA, 1.0, 0.0);
    JooooAA.reset();
    timer_on("I/O");
    IooooAA->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IooooAA.reset();
    
 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    IooooBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <oo|oo>", noccB, noccB, noccB, noccB));
    JooooBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (oo|oo)", noccB, noccB, noccB, noccB));
    timer_on("I/O");
    JooooBB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IooooBB->sort(1324, JooooBB, 1.0, 0.0);
    JooooBB.reset();
    timer_on("I/O");
    IooooBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IooooBB.reset();

    // AB spin case
    IooooAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <Oo|Oo>", noccA, noccB, noccA, noccB));
    JooooAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|oo)", noccA, noccA, noccB, noccB));
    timer_on("I/O");
    JooooAB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IooooAB->sort(1324, JooooAB, 1.0, 0.0);
    JooooAB.reset();
    timer_on("I/O");
    IooooAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IooooAB.reset();
 }
    timer_off("Build <ij|kl>)");
}// end tei_oooo_phys_ref

//=======================================================
//          <OO|OV>
//=======================================================          
void DFOCC::tei_ooov_phys_ref()
{   
    timer_on("Build <ij|ka>)");
    // AA spin case
    IooovAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OO|OV>", noccA, noccA, noccA, nvirA));
    JooovAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|OV)", noccA, noccA, noccA, nvirA));
    timer_on("I/O");
    JooovAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IooovAA->sort(1324, JooovAA, 1.0, 0.0);
    JooovAA.reset();
    timer_on("I/O");
    IooovAA->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IooovAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    IooovBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <oo|ov>", noccB, noccB, noccB, nvirB));
    JooovBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (oo|ov)", noccB, noccB, noccB, nvirB));
    timer_on("I/O");
    JooovBB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IooovBB->sort(1324, JooovBB, 1.0, 0.0);
    JooovBB.reset();
    timer_on("I/O");
    IooovBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IooovBB.reset();

    // AB spin case
    IooovAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <Oo|Ov>", noccA, noccB, noccA, nvirB));
    JooovAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|ov)", noccA, noccA, noccB, nvirB));
    timer_on("I/O");
    JooovAB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IooovAB->sort(1324, JooovAB, 1.0, 0.0);
    JooovAB.reset();
    timer_on("I/O");
    IooovAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IooovAB.reset();

    // BA spin case
    IoovoAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <Oo|Vo>", noccA, noccB, nvirA, noccB));
    JovooAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|oo)", noccA, nvirA, noccB, noccB));
    timer_on("I/O");
    JovooAB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IoovoAB->sort(1324, JovooAB, 1.0, 0.0);
    JovooAB.reset();
    timer_on("I/O");
    IoovoAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IoovoAB.reset();
 }
    timer_off("Build <ij|ka>)");
}// end tei_ooov_phys_ref

//=======================================================
//          <OO|VV>
//=======================================================          
void DFOCC::tei_oovv_phys_ref()
{   
    timer_on("Build <ij|ab>)");

    // AA spin case
    IoovvAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OO|VV>", noccA, noccA, nvirA, nvirA));
    JovovAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|OV)", noccA, nvirA, noccA, nvirA));
    timer_on("I/O");
    JovovAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IoovvAA->sort(1324, JovovAA, 1.0, 0.0);
    JovovAA.reset();
    timer_on("I/O");
    IoovvAA->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IoovvAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    IoovvBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <oo|vv>", noccB, noccB, nvirB, nvirB));
    JovovBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (ov|ov)", noccB, nvirB, noccB, nvirB));
    timer_on("I/O");
    JovovBB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IoovvBB->sort(1324, JovovBB, 1.0, 0.0);
    JovovBB.reset();
    timer_on("I/O");
    IoovvBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IoovvBB.reset();

    // AB spin case
    IoovvAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <Oo|Vv>", noccA, noccB, nvirA, nvirB));
    JovovAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|ov)", noccA, nvirA, noccB, nvirB));
    timer_on("I/O");
    JovovAB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IoovvAB->sort(1324, JovovAB, 1.0, 0.0);
    JovovAB.reset();
    timer_on("I/O");
    IoovvAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IoovvAB.reset();
 }
    timer_off("Build <ij|ab>)");
}// end tei_oovv_phys_ref

//=======================================================
//          <OV|OV>
//=======================================================          
void DFOCC::tei_ovov_phys_ref()
{   
    timer_on("Build <ia|jb>)");
    // AA spin case
    IovovAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OV|OV>", noccA, nvirA, noccA, nvirA));
    JoovvAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|VV)", noccA, noccA, nvirA, nvirA));
    timer_on("I/O");
    JoovvAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IovovAA->sort(1324, JoovvAA, 1.0, 0.0);
    JoovvAA.reset();
    timer_on("I/O");
    IovovAA->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IovovAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    IovovBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <ov|ov>", noccB, nvirB, noccB, nvirB));
    JoovvBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (oo|vv)", noccB, noccB, nvirB, nvirB));
    timer_on("I/O");
    JoovvBB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IovovBB->sort(1324, JoovvBB, 1.0, 0.0);
    JoovvBB.reset();
    timer_on("I/O");
    IovovBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IovovBB.reset();

    // AB spin case
    IovovAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <Ov|Ov>", noccA, nvirB, noccA, nvirB));
    JoovvAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|vv)", noccA, noccA, nvirB, nvirB));
    timer_on("I/O");
    JoovvAB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IovovAB->sort(1324, JoovvAB, 1.0, 0.0);
    JoovvAB.reset();
    timer_on("I/O");
    IovovAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IovovAB.reset();

    // BA spin case
    IvovoAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <Vo|Vo>", nvirA, noccB, nvirA, noccB));
    JvvooAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (VV|oo)", nvirA, nvirA, noccB, noccB));
    timer_on("I/O");
    JvvooAB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IvovoAB->sort(1324, JvvooAB, 1.0, 0.0);
    JvvooAB.reset();
    timer_on("I/O");
    IvovoAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IvovoAB.reset();
 }
    timer_off("Build <ia|jb>)");
}// end tei_ovov_phys_ref

//=======================================================
//          <OO||OO>
//=======================================================          
void DFOCC::tei_oooo_anti_symm_ref()
{   
    timer_on("Build <ij||kl>)");

    // AA spin case
    AIooooAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OO||OO>", noccA, noccA, noccA, noccA));
    IooooAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OO|OO>", noccA, noccA, noccA, noccA));
    timer_on("I/O");
    IooooAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIooooAA->sort(1243, IooooAA, 1.0, 0.0);
    AIooooAA->scale(-1.0);
    AIooooAA->add(IooooAA);
    IooooAA.reset();
    timer_on("I/O");
    AIooooAA->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIooooAA.reset();
    
 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    AIooooBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <oo||oo>", noccB, noccB, noccB, noccB));
    IooooBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <oo|oo>", noccB, noccB, noccB, noccB));
    timer_on("I/O");
    IooooBB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIooooBB->sort(1243, IooooBB, 1.0, 0.0);
    AIooooBB->scale(-1.0);
    AIooooBB->add(IooooBB);
    IooooBB.reset();
    timer_on("I/O");
    AIooooBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIooooBB.reset();
 }
    timer_off("Build <ij||kl>)");
}// end tei_oooo_anti_symm_ref

//=======================================================
//          <OO||OV>
//=======================================================          
void DFOCC::tei_ooov_anti_symm_ref()
{   
    timer_on("Build <ij||ka>)");
    // <ij||ka> = <ij|ka> - <ji|ka>
    // AA spin case
    AIooovAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OO||OV>", noccA, noccA, noccA, nvirA));
    IooovAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OO|OV>", noccA, noccA, noccA, nvirA));
    timer_on("I/O");
    IooovAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIooovAA->sort(2134, IooovAA, 1.0, 0.0);
    AIooovAA->scale(-1.0);
    AIooovAA->add(IooovAA);
    IooovAA.reset();
    timer_on("I/O");
    AIooovAA->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIooovAA.reset();
    
 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    AIooovBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <oo||ov>", noccB, noccB, noccB, nvirB));
    IooovBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <oo|ov>", noccB, noccB, noccB, nvirB));
    timer_on("I/O");
    IooovBB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIooovBB->sort(2134, IooovBB, 1.0, 0.0);
    AIooovBB->scale(-1.0);
    AIooovBB->add(IooovBB);
    IooovBB.reset();
    timer_on("I/O");
    AIooovBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIooovBB.reset();
 }
    timer_off("Build <ij||ka>)");
}// end tei_ooov_anti_symm_ref

//=======================================================
//          <OO||VV>
//=======================================================          
void DFOCC::tei_oovv_anti_symm_ref()
{   
    timer_on("Build <ij||ab>)");

    // AA spin case
    AIoovvAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OO||VV>", noccA, noccA, nvirA, nvirA));
    IoovvAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OO|VV>", noccA, noccA, nvirA, nvirA));
    timer_on("I/O");
    IoovvAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIoovvAA->sort(1243, IoovvAA, 1.0, 0.0);
    AIoovvAA->scale(-1.0);
    AIoovvAA->add(IoovvAA);
    IoovvAA.reset();
    timer_on("I/O");
    AIoovvAA->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIoovvAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    AIoovvBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <oo||vv>", noccB, noccB, nvirB, nvirB));
    IoovvBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <oo|vv>", noccB, noccB, nvirB, nvirB));
    timer_on("I/O");
    IoovvBB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIoovvBB->sort(1243, IoovvBB, 1.0, 0.0);
    AIoovvBB->scale(-1.0);
    AIoovvBB->add(IoovvBB);
    IoovvBB.reset();
    timer_on("I/O");
    AIoovvBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIoovvBB.reset();
 }
    timer_off("Build <ij||ab>)");
}// end tei_oovv_anti_symm_ref

//=======================================================
//          <OV||OV>
//=======================================================          
void DFOCC::tei_ovov_anti_symm_ref()
{   
    timer_on("Build <ia||jb>)");
    // <ia||jb> = <ia|jb> - <ia|bj> = <ia|jb> - (ib|ja)
    // AA spin case
    AIovovAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OV||OV>", noccA, nvirA, noccA, nvirA));
    JovovAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|OV)", noccA, nvirA, noccA, nvirA));
    timer_on("I/O");
    JovovAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIovovAA->sort(1432, JovovAA, 1.0, 0.0);
    JovovAA.reset();
    AIovovAA->scale(-1.0);
    IovovAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OV|OV>", noccA, nvirA, noccA, nvirA));
    timer_on("I/O");
    IovovAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIovovAA->add(IovovAA);
    IovovAA.reset();
    timer_on("I/O");
    AIovovAA->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIovovAA.reset();
    
 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    AIovovBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <ov||ov>", noccB, nvirB, noccB, nvirB));
    JovovBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (ov|ov)", noccB, nvirB, noccB, nvirB));
    timer_on("I/O");
    JovovBB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIovovBB->sort(1432, JovovBB, 1.0, 0.0);
    JovovBB.reset();
    AIovovBB->scale(-1.0);
    IovovBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <ov|ov>", noccB, nvirB, noccB, nvirB));
    timer_on("I/O");
    IovovBB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIovovBB->add(IovovBB);
    IovovBB.reset();
    timer_on("I/O");
    AIovovBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIovovBB.reset();
    //fprintf(outfile,"\tI am here.\n"); fflush(outfile);
 }
    timer_off("Build <ia||jb>)");
}// end tei_ovov_anti_symm_ref


}} // End Namespaces


