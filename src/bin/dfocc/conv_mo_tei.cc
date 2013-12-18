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

void DFOCC::tei_ijkl_chem()
{   
    timer_on("Build (oo|oo)");
    // AA spin case
    JijklAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|KL)", naoccA, naoccA, naoccA, naoccA));
    bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA * naoccA));
    timer_on("I/O");
    bQijA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JijklAA->gemm(true, false, bQijA, bQijA, 1.0, 0.0);
    if (reference_ == "RESTRICTED") bQijA.reset();
    timer_on("I/O");
    JijklAA->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JijklAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    JijklBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ij|kl)", naoccB, naoccB, naoccB, naoccB));
    bQijB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ij)", nQ, naoccB * naoccB));
    timer_on("I/O");
    bQijB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JijklBB->gemm(true, false, bQijB, bQijB, 1.0, 0.0);
    timer_on("I/O");
    JijklBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JijklBB.reset();

    // AB spin case
    JijklAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|kl)", naoccA, naoccA, naoccB, naoccB));
    JijklAB->gemm(true, false, bQijA, bQijB, 1.0, 0.0);
    bQijA.reset();
    bQijB.reset();
    timer_on("I/O");
    JijklAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
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
    JooooBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (oo|oo)", noccB, noccB, noccB, noccB));
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|oo)", nQ, noccB * noccB));
    timer_on("I/O");
    bQooB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JooooBB->gemm(true, false, bQooB, bQooB, 1.0, 0.0);
    timer_on("I/O");
    JooooBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JooooBB.reset();

    // AB spin case
    JooooAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|oo)", noccA, noccA, noccB, noccB));
    JooooAB->gemm(true, false, bQooA, bQooB, 1.0, 0.0);
    bQooA.reset();
    bQooB.reset();
    timer_on("I/O");
    JooooAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
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
    timer_on("I/O");
    bQijA->read(psio_, PSIF_DFOCC_INTS);
    bQiaA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JijkaAA->gemm(true, false, bQijA, bQiaA, 1.0, 0.0);
    if (reference_ == "RESTRICTED") bQijA.reset();
    if (reference_ == "RESTRICTED") bQiaA.reset();
    timer_on("I/O");
    JijkaAA->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JijkaAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    JijkaBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ij|ka)", naoccB, naoccB, naoccB, navirB));
    bQijB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ij)", nQ, naoccB * naoccB));
    bQiaB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ia)", nQ, naoccB * navirB));
    timer_on("I/O");
    bQijB->read(psio_, PSIF_DFOCC_INTS);
    bQiaB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JijkaBB->gemm(true, false, bQijB, bQiaB, 1.0, 0.0);
    timer_on("I/O");
    JijkaBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JijkaBB.reset();

    // AB spin case
    JijkaAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|ka)", naoccA, naoccA, naoccB, navirB));
    JijkaAB->gemm(true, false, bQijA, bQiaB, 1.0, 0.0);
    bQijA.reset();
    bQiaB.reset();
    timer_on("I/O");
    JijkaAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JijkaAB.reset();

    // BA spin case
    JiajkAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|jk)", naoccA, navirA, naoccB, naoccB));
    JiajkAB->gemm(true, false, bQiaA, bQijB, 1.0, 0.0);
    bQijB.reset();
    bQiaA.reset();
    timer_on("I/O");
    JiajkAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
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
    JooovBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (oo|ov)", noccB, noccB, noccB, nvirB));
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|oo)", nQ, noccB * noccB));
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ov)", nQ, noccB * nvirB));
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
    JooovAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|ov)", noccA, noccA, noccB, nvirB));
    JooovAB->gemm(true, false, bQooA, bQovB, 1.0, 0.0);
    bQooA.reset();
    bQovB.reset();
    timer_on("I/O");
    JooovAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JooovAB.reset();

    // BA spin case
    JovooAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OV|oo)", noccA, nvirA, noccB, noccB));
    JovooAB->gemm(true, false, bQovA, bQooB, 1.0, 0.0);
    bQooB.reset();
    bQovA.reset();
    timer_on("I/O");
    JovooAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
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
    bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA * naoccA));
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA * navirA));
    timer_on("I/O");
    bQijA->read(psio_, PSIF_DFOCC_INTS);
    bQabA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JijabAA->gemm(true, false, bQijA, bQabA, 1.0, 0.0);
    if (reference_ == "RESTRICTED") bQijA.reset();
    if (reference_ == "RESTRICTED") bQabA.reset();
    timer_on("I/O");
    JijabAA->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JijabAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    JijabBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ij|ab)", naoccB, naoccB, navirB, navirB));
    bQijB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ij)", nQ, naoccB * naoccB));
    bQabB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ab)", nQ, navirB * navirB));
    timer_on("I/O");
    bQijB->read(psio_, PSIF_DFOCC_INTS);
    bQabB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JijabBB->gemm(true, false, bQijB, bQabB, 1.0, 0.0);
    timer_on("I/O");
    JijabBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JijabBB.reset();

    // AB spin case
    JijabAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|ab)", naoccA, naoccA, navirB, navirB));
    JijabAB->gemm(true, false, bQijA, bQabB, 1.0, 0.0);
    bQijA.reset();
    bQabB.reset();
    timer_on("I/O");
    JijabAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JijabAB.reset();

    // BA spin case
    JabijAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (AB|ij)",navirA, navirA, naoccB, naoccB));
    JabijAB->gemm(true, false, bQabA, bQijB, 1.0, 0.0);
    bQijB.reset();
    bQabA.reset();
    timer_on("I/O");
    JabijAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
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
    bQvvA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA * nvirA));
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
    JoovvBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (oo|vv)", noccB, noccB, nvirB, nvirB));
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|oo)", nQ, noccB * noccB));
    bQvvB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|vv)", nQ, nvirB * nvirB));
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
    JoovvAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|vv)", noccA, noccA, nvirB, nvirB));
    JoovvAB->gemm(true, false, bQooA, bQvvB, 1.0, 0.0);
    bQooA.reset();
    bQvvB.reset();
    timer_on("I/O");
    JoovvAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JoovvAB.reset();

    // BA spin case
    JvvooAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (VV|oo)",nvirA, nvirA, noccB, noccB));
    JvvooAB->gemm(true, false, bQvvA, bQooB, 1.0, 0.0);
    bQooB.reset();
    bQvvA.reset();
    timer_on("I/O");
    JvvooAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
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
    timer_on("I/O");
    bQiaA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JiajbAA->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    if (reference_ == "RESTRICTED") bQiaA.reset();
    timer_on("I/O");
    JiajbAA->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JiajbAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    JiajbBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB));
    bQiaB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ia)", nQ, naoccB * navirB));
    timer_on("I/O");
    bQiaB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JiajbBB->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    timer_on("I/O");
    JiajbBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JiajbBB.reset();

    // AB spin case
    JiajbAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|jb)", naoccA, navirA, naoccB, navirB));
    JiajbAB->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    bQiaA.reset();
    bQiaB.reset();
    timer_on("I/O");
    JiajbAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
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
    JovovBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ov|ov)", noccB, nvirB, noccB, nvirB));
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ov)", nQ, noccB * nvirB));
    timer_on("I/O");
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JovovBB->gemm(true, false, bQovB, bQovB, 1.0, 0.0);
    timer_on("I/O");
    JovovBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JovovBB.reset();

    // AB spin case
    JovovAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OV|ov)", noccA, nvirA, noccB, nvirB));
    JovovAB->gemm(true, false, bQovA, bQovB, 1.0, 0.0);
    bQovA.reset();
    bQovB.reset();
    timer_on("I/O");
    JovovAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    JovovAB.reset();
 }
    timer_off("Build (ov|ov)");
}// end tei_ovov_chem

//=======================================================
//          <IJ|KL>
//=======================================================          
void DFOCC::tei_ijkl_phys()
{   
    timer_on("Build <ij|kl>)");
    // AA spin case
    IijklAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|KL>", naoccA, naoccA, naoccA, naoccA));
    JijklAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|KL)", naoccA, naoccA, naoccA, naoccA));
    timer_on("I/O");
    JijklAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IijklAA->sort(1324, JijklAA, 1.0, 0.0);
    JijklAA.reset();
    timer_on("I/O");
    IijklAA->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IijklAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    IijklBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|kl>", naoccB, naoccB, naoccB, naoccB));
    JijklBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ij|kl)", naoccB, naoccB, naoccB, naoccB));
    timer_on("I/O");
    JijklBB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IijklBB->sort(1324, JijklBB, 1.0, 0.0);
    JijklBB.reset();
    timer_on("I/O");
    IijklBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IijklBB.reset();

    // AB spin case
    IijklAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ij|Kl>", naoccA, naoccB, naoccA, naoccB));
    JijklAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|kl)", naoccA, naoccA, naoccB, naoccB));
    timer_on("I/O");
    JijklAB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IijklAB->sort(1324, JijklAB, 1.0, 0.0);
    JijklAB.reset();
    timer_on("I/O");
    IijklAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IijklAB.reset();
 }
    timer_off("Build <ij|kl>)");
}// end tei_ijkl_phys

//=======================================================
//          <OO|OO>
//=======================================================          
void DFOCC::tei_oooo_phys()
{   
    timer_on("Build <ij|kl>)");
    // AA spin case
    IooooAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <OO|OO>", noccA, noccA, noccA, noccA));
    JooooAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|OO)", noccA, noccA, noccA, noccA));
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
    IooooBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <oo|oo>", noccB, noccB, noccB, noccB));
    JooooBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (oo|oo)", noccB, noccB, noccB, noccB));
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
    IooooAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Oo|Oo>", noccA, noccB, noccA, noccB));
    JooooAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|oo)", noccA, noccA, noccB, noccB));
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
}// end tei_oooo_phys

//=======================================================
//          <IJ|KA>
//=======================================================          
void DFOCC::tei_ijka_phys()
{   
    timer_on("Build <ij|ka>)");
    // AA spin case
    IijkaAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|KA>", naoccA, naoccA, naoccA, navirA));
    JijkaAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|KA)", naoccA, naoccA, naoccA, navirA));
    timer_on("I/O");
    JijkaAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IijkaAA->sort(1324, JijkaAA, 1.0, 0.0);
    JijkaAA.reset();
    timer_on("I/O");
    IijkaAA->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IijkaAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    IijkaBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|ka>", naoccB, naoccB, naoccB, navirB));
    JijkaBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ij|ka)", naoccB, naoccB, naoccB, navirB));
    timer_on("I/O");
    JijkaBB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IijkaBB->sort(1324, JijkaBB, 1.0, 0.0);
    JijkaBB.reset();
    timer_on("I/O");
    IijkaBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IijkaBB.reset();

    // AB spin case
    IijkaAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ij|Ka>", naoccA, naoccB, naoccA, navirB));
    JijkaAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|ka)", naoccA, naoccA, naoccB, navirB));
    timer_on("I/O");
    JijkaAB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IijkaAB->sort(1324, JijkaAB, 1.0, 0.0);
    JijkaAB.reset();
    timer_on("I/O");
    IijkaAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IijkaAB.reset();

    // BA spin case
    IijakAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ij|Ak>", naoccA, naoccB, navirA, naoccB));
    JiajkAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|jk)", naoccA, navirA, naoccB, naoccB));
    timer_on("I/O");
    JiajkAB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IijakAB->sort(1324, JiajkAB, 1.0, 0.0);
    JiajkAB.reset();
    timer_on("I/O");
    IijakAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IijakAB.reset();
 }
    timer_off("Build <ij|ka>)");
}// end tei_ijka_phys

//=======================================================
//          <OO|OV>
//=======================================================          
void DFOCC::tei_ooov_phys()
{   
    timer_on("Build <ij|ka>)");
    // AA spin case
    IooovAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <OO|OV>", noccA, noccA, noccA, nvirA));
    JooovAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|OV)", noccA, noccA, noccA, nvirA));
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
    IooovBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <oo|ov>", noccB, noccB, noccB, nvirB));
    JooovBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (oo|ov)", noccB, noccB, noccB, nvirB));
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
    IooovAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Oo|Ov>", noccA, noccB, noccA, nvirB));
    JooovAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|ov)", noccA, noccA, noccB, nvirB));
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
    IoovoAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Oo|Vo>", noccA, noccB, nvirA, noccB));
    JovooAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OV|oo)", noccA, nvirA, noccB, noccB));
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
}// end tei_ooov_phys

//=======================================================
//          <IJ|AB>
//=======================================================          
void DFOCC::tei_ijab_phys()
{   
    timer_on("Build <ij|ab>)");

    // AA spin case
    IijabAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA));
    JiajbAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    timer_on("I/O");
    JiajbAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IijabAA->sort(1324, JiajbAA, 1.0, 0.0);
    JiajbAA.reset();
    timer_on("I/O");
    IijabAA->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IijabAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    IijabBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|ab>", naoccB, naoccB, navirB, navirB));
    JiajbBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB));
    timer_on("I/O");
    JiajbBB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IijabBB->sort(1324, JiajbBB, 1.0, 0.0);
    JiajbBB.reset();
    timer_on("I/O");
    IijabBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IijabBB.reset();

    // AB spin case
    IijabAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    JiajbAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|jb)", naoccA, navirA, naoccB, navirB));
    timer_on("I/O");
    JiajbAB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IijabAB->sort(1324, JiajbAB, 1.0, 0.0);
    JiajbAB.reset();
    timer_on("I/O");
    IijabAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IijabAB.reset();
 }
    timer_off("Build <ij|ab>)");
}// end tei_ijab_phys

//=======================================================
//          <OO|VV>
//=======================================================          
void DFOCC::tei_oovv_phys()
{   
    timer_on("Build <ij|ab>)");
    // AA spin case
    IoovvAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <OO|VV>", noccA, noccA, nvirA, nvirA));
    JovovAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OV|OV)", noccA, nvirA, noccA, nvirA));
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
    IoovvBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <oo|vv>", noccB, noccB, nvirB, nvirB));
    JovovBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ov|ov)", noccB, nvirB, noccB, nvirB));
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
    IoovvAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Oo|Vv>", noccA, noccB, nvirA, nvirB));
    JovovAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OV|ov)", noccA, nvirA, noccB, nvirB));
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
}// end tei_oovv_phys

//=======================================================
//          <IA|JB>
//=======================================================          
void DFOCC::tei_iajb_phys()
{   
    timer_on("Build <ia|jb>)");
    // AA spin case
    IiajbAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IA|JB>", naoccA, navirA, naoccA, navirA));
    JijabAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|AB)", naoccA, naoccA, navirA, navirA));
    timer_on("I/O");
    JijabAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IiajbAA->sort(1324, JijabAA, 1.0, 0.0);
    JijabAA.reset();
    timer_on("I/O");
    IiajbAA->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IiajbAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    IiajbBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ia|jb>", naoccB, navirB, naoccB, navirB));
    JijabBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ij|ab)", naoccB, naoccB, navirB, navirB));
    timer_on("I/O");
    JijabBB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IiajbBB->sort(1324, JijabBB, 1.0, 0.0);
    JijabBB.reset();
    timer_on("I/O");
    IiajbBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IiajbBB.reset();

    // AB spin case
    IiajbAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ia|Jb>", naoccA, navirB, naoccA, navirB));
    JijabAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IJ|ab)", naoccA, naoccA, navirB, navirB));
    timer_on("I/O");
    JijabAB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IiajbAB->sort(1324, JijabAB, 1.0, 0.0);
    JijabAB.reset();
    timer_on("I/O");
    IiajbAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IiajbAB.reset();

    // BA spin case
    IaibjAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ai|Bj>", navirA, naoccB, navirA, naoccB));
    JabijAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (AB|ij)", navirA, navirA, naoccB, naoccB));
    timer_on("I/O");
    JabijAB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IaibjAB->sort(1324, JabijAB, 1.0, 0.0);
    JabijAB.reset();
    timer_on("I/O");
    IaibjAB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    IaibjAB.reset();
 }
    timer_off("Build <ia|jb>)");
}// end tei_iajb_phys

//=======================================================
//          <OV|OV>
//=======================================================          
void DFOCC::tei_ovov_phys()
{   
    timer_on("Build <ia|jb>)");
    // AA spin case
    IovovAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <OV|OV>", noccA, nvirA, noccA, nvirA));
    JoovvAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|VV)", noccA, noccA, nvirA, nvirA));
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
    IovovBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ov|ov>", noccB, nvirB, noccB, nvirB));
    JoovvBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (oo|vv)", noccB, noccB, nvirB, nvirB));
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
    IovovAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ov|Ov>", noccA, nvirB, noccA, nvirB));
    JoovvAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OO|vv)", noccA, noccA, nvirB, nvirB));
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
    IvovoAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Vo|Vo>", nvirA, noccB, nvirA, noccB));
    JvvooAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (VV|oo)", nvirA, nvirA, noccB, noccB));
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
}// end tei_ovov_phys

//=======================================================
//          <IJ||KL>
//=======================================================          
void DFOCC::tei_ijkl_anti_symm()
{   
    timer_on("Build <ij||kl>)");
    // AA spin case
    AIijklAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ||KL>", naoccA, naoccA, naoccA, naoccA));
    IijklAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|KL>", naoccA, naoccA, naoccA, naoccA));
    timer_on("I/O");
    IijklAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIijklAA->sort(1243, IijklAA, 1.0, 0.0);
    AIijklAA->scale(-1.0);
    AIijklAA->add(IijklAA);
    IijklAA.reset();
    timer_on("I/O");
    AIijklAA->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIijklAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    AIijklBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij||kl>", naoccB, naoccB, naoccB, naoccB));
    IijklBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|kl>", naoccB, naoccB, naoccB, naoccB));
    timer_on("I/O");
    IijklBB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIijklBB->sort(1243, IijklBB, 1.0, 0.0);
    AIijklBB->scale(-1.0);
    AIijklBB->add(IijklBB);
    IijklBB.reset();
    timer_on("I/O");
    AIijklBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIijklBB.reset();
 }
    timer_off("Build <ij||kl>)");
}// end tei_ijkl_anti_symm

//=======================================================
//          <OO||OO>
//=======================================================          
void DFOCC::tei_oooo_anti_symm()
{   
    timer_on("Build <ij||kl>)");
    // AA spin case
    AIooooAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <OO||OO>", noccA, noccA, noccA, noccA));
    IooooAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <OO|OO>", noccA, noccA, noccA, noccA));
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
    AIooooBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <oo||oo>", noccB, noccB, noccB, noccB));
    IooooBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <oo|oo>", noccB, noccB, noccB, noccB));
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
}// end tei_oooo_anti_symm

//=======================================================
//          <IJ||KA>
//=======================================================          
void DFOCC::tei_ijka_anti_symm()
{   
    timer_on("Build <ij||ka>)");
    // <ij||ka> = <ij|ka> - <ji|ka>
    // AA spin case
    AIijkaAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ||KA>", naoccA, naoccA, naoccA, navirA));
    IijkaAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|KA>", naoccA, naoccA, naoccA, navirA));
    timer_on("I/O");
    IijkaAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIijkaAA->sort(2134, IijkaAA, 1.0, 0.0);
    AIijkaAA->scale(-1.0);
    AIijkaAA->add(IijkaAA);
    IijkaAA.reset();
    timer_on("I/O");
    AIijkaAA->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIijkaAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    AIijkaBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij||ka>", naoccB, naoccB, naoccB, navirB));
    IijkaBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|ka>", naoccB, naoccB, naoccB, navirB));
    timer_on("I/O");
    IijkaBB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIijkaBB->sort(2134, IijkaBB, 1.0, 0.0);
    AIijkaBB->scale(-1.0);
    AIijkaBB->add(IijkaBB);
    IijkaBB.reset();
    timer_on("I/O");
    AIijkaBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIijkaBB.reset();
 }
    timer_off("Build <ij||ka>)");
}// end tei_ijka_anti_symm

//=======================================================
//          <OO||OV>
//=======================================================          
void DFOCC::tei_ooov_anti_symm()
{   
    timer_on("Build <ij||ka>)");
    // <ij||ka> = <ij|ka> - <ji|ka>
    // AA spin case
    AIooovAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <OO||OV>", noccA, noccA, noccA, nvirA));
    IooovAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <OO|OV>", noccA, noccA, noccA, nvirA));
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
    AIooovBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <oo||ov>", noccB, noccB, noccB, nvirB));
    IooovBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <oo|ov>", noccB, noccB, noccB, nvirB));
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
}// end tei_ooov_anti_symm

//=======================================================
//          <IJ||AB>
//=======================================================          
void DFOCC::tei_ijab_anti_symm()
{   
    timer_on("Build <ij||ab>)");

    // AA spin case
    AIijabAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ||AB>", naoccA, naoccA, navirA, navirA));
    IijabAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA));
    timer_on("I/O");
    IijabAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIijabAA->sort(1243, IijabAA, 1.0, 0.0);
    AIijabAA->scale(-1.0);
    AIijabAA->add(IijabAA);
    IijabAA.reset();
    timer_on("I/O");
    AIijabAA->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIijabAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    AIijabBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij||ab>", naoccB, naoccB, navirB, navirB));
    IijabBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|ab>", naoccB, naoccB, navirB, navirB));
    timer_on("I/O");
    IijabBB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIijabBB->sort(1243, IijabBB, 1.0, 0.0);
    AIijabBB->scale(-1.0);
    AIijabBB->add(IijabBB);
    IijabBB.reset();
    timer_on("I/O");
    AIijabBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIijabBB.reset();
 }
    timer_off("Build <ij||ab>)");
}// end tei_ijab_anti_symm

//=======================================================
//          <OO||VV>
//=======================================================          
void DFOCC::tei_oovv_anti_symm()
{   
    timer_on("Build <ij||ab>)");
    // AA spin case
    AIoovvAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <OO||VV>", noccA, noccA, nvirA, nvirA));
    IoovvAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <OO|VV>", noccA, noccA, nvirA, nvirA));
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
    AIoovvBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <oo||vv>", noccB, noccB, nvirB, nvirB));
    IoovvBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <oo|vv>", noccB, noccB, nvirB, nvirB));
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
}// end tei_oovv_anti_symm

//=======================================================
//          <IA||JB>
//=======================================================          
void DFOCC::tei_iajb_anti_symm()
{   
    timer_on("Build <ia||jb>)");
    // <ia||jb> = <ia|jb> - <ia|bj> = <ia|jb> - (ib|ja)
    // AA spin case
    AIiajbAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IA||JB>", naoccA, navirA, naoccA, navirA));
    JiajbAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    timer_on("I/O");
    JiajbAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIiajbAA->sort(1432, JiajbAA, 1.0, 0.0);
    JiajbAA.reset();
    AIiajbAA->scale(-1.0);
    IiajbAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IA|JB>", naoccA, navirA, naoccA, navirA));
    timer_on("I/O");
    IiajbAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIiajbAA->add(IiajbAA);
    IiajbAA.reset();
    timer_on("I/O");
    AIiajbAA->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIiajbAA.reset();

 if (reference_ == "UNRESTRICTED") {
    // BB spin case
    AIiajbBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ia||jb>", naoccB, navirB, naoccB, navirB));
    JiajbBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB));
    timer_on("I/O");
    JiajbBB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIiajbBB->sort(1432, JiajbBB, 1.0, 0.0);
    JiajbBB.reset();
    AIiajbBB->scale(-1.0);
    IiajbBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ia|jb>", naoccB, navirB, naoccB, navirB));
    timer_on("I/O");
    IiajbBB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIiajbBB->add(IiajbBB);
    IiajbBB.reset();
    timer_on("I/O");
    AIiajbBB->write(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIiajbBB.reset();
    //fprintf(outfile,"\tI am here.\n"); fflush(outfile);
 }
    timer_off("Build <ia||jb>)");
}// end tei_iajb_anti_symm

//=======================================================
//          <OV||OV>
//=======================================================          
void DFOCC::tei_ovov_anti_symm()
{   
    timer_on("Build <ia||jb>)");
    // <ia||jb> = <ia|jb> - <ia|bj> = <ia|jb> - (ib|ja)
    // AA spin case
    AIovovAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <OV||OV>", noccA, nvirA, noccA, nvirA));
    JovovAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (OV|OV)", noccA, nvirA, noccA, nvirA));
    timer_on("I/O");
    JovovAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIovovAA->sort(1432, JovovAA, 1.0, 0.0);
    JovovAA.reset();
    AIovovAA->scale(-1.0);
    IovovAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <OV|OV>", noccA, nvirA, noccA, nvirA));
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
    AIovovBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ov||ov>", noccB, nvirB, noccB, nvirB));
    JovovBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ov|ov)", noccB, nvirB, noccB, nvirB));
    timer_on("I/O");
    JovovBB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    AIovovBB->sort(1432, JovovBB, 1.0, 0.0);
    JovovBB.reset();
    AIovovBB->scale(-1.0);
    IovovBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ov|ov>", noccB, nvirB, noccB, nvirB));
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
}// end tei_ovov_anti_symm


}} // End Namespaces


