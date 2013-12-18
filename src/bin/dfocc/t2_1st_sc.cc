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
  
void DFOCC::t2_1st_sc()
{   

    timer_on("Form 1st-order T2");
if (reference_ == "RESTRICTED") {
    // Example: init from file
    //Array2d *temp = new Array2d(psio_, PSIF_DFOCC_INTS, "(ia|jb)", naoccA, navirA, naoccA, navirA);

    // Build amplitudes in Mulliken order 
    t2p_1 = SharedTensor2d(new Tensor2d("T2_1(ia,jb)", naoccA, navirA, naoccA, navirA));
    JiajbAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    timer_on("I/O");
    JiajbAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    t2p_1->copy(JiajbAA);
    JiajbAA.reset();
    t2p_1->apply_denom_chem(nfrzc, noccA, FockA);
    timer_on("I/O");
    t2p_1->write(psio_, PSIF_DFOCC_AMPS);
    timer_off("I/O");
 
    // Sort amplitudes to Dirac order
    t2_1 = SharedTensor2d(new Tensor2d("T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    t2_1->sort(1324, t2p_1, 1.0, 0.0);
    timer_on("I/O");
    t2_1->write(psio_, PSIF_DFOCC_AMPS);
    timer_off("I/O");
    t2_1.reset();

    SharedTensor2d temp = SharedTensor2d(new Tensor2d("T2_1(ia,jb) - T2_1(ib,ja)", naoccA, navirA, naoccA, navirA));
    temp->sort(1432, t2p_1, 1.0, 0.0);
    temp->scale(-1.0);
    temp->add(t2p_1);
    timer_on("I/O");
    temp->write(psio_, PSIF_DFOCC_AMPS);
    timer_off("I/O");

    u2p_1 = SharedTensor2d(new Tensor2d("2*T2_1(ia,jb) - T2_1(ib,ja)", naoccA, navirA, naoccA, navirA));
    u2p_1->copy(temp);
    temp.reset();
    u2p_1->add(t2p_1);
    t2p_1.reset();
    timer_on("I/O");
    u2p_1->write(psio_, PSIF_DFOCC_AMPS);
    timer_off("I/O");
    u2p_1.reset();
}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // T2AA
    t2_1AA = SharedTensor2d(new Tensor2d("T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    AIijabAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ||AB>", naoccA, naoccA, navirA, navirA));
    timer_on("I/O");
    AIijabAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    t2_1AA->copy(AIijabAA);
    AIijabAA.reset();
    t2_1AA->apply_denom(nfrzc, noccA, FockA);
    timer_on("I/O");
    t2_1AA->write(psio_, PSIF_DFOCC_AMPS);
    timer_off("I/O");
    t2p_1 = SharedTensor2d(new Tensor2d("T2_1(IA,JB)", naoccA, navirA, naoccA, navirA));
    t2p_1->sort(1324, t2_1AA, 1.0, 0.0);
    t2_1AA.reset();
    timer_on("I/O");
    t2p_1->write(psio_, PSIF_DFOCC_AMPS);
    timer_off("I/O");
    t2p_1.reset();

    // T2BB
    t2_1BB = SharedTensor2d(new Tensor2d("T2_1 <ij|ab>", naoccB, naoccB, navirB, navirB));
    AIijabBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij||ab>", naoccB, naoccB, navirB, navirB));
    timer_on("I/O");
    AIijabBB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    t2_1BB->copy(AIijabBB);
    AIijabBB.reset();
    t2_1BB->apply_denom(nfrzc, noccB, FockB);
    timer_on("I/O");
    t2_1BB->write(psio_, PSIF_DFOCC_AMPS);
    timer_off("I/O");
    t2p_1 = SharedTensor2d(new Tensor2d("T2_1(ia,jb)", naoccB, navirB, naoccB, navirB));
    t2p_1->sort(1324, t2_1BB, 1.0, 0.0);
    t2_1BB.reset();
    timer_on("I/O");
    t2p_1->write(psio_, PSIF_DFOCC_AMPS);
    timer_off("I/O");
    t2p_1.reset();

    // T2AB
    t2_1AB = SharedTensor2d(new Tensor2d("T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    IijabAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    timer_on("I/O");
    IijabAB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    t2_1AB->copy(IijabAB);
    IijabAB.reset();
    t2_1AB->apply_denom_os(nfrzc, noccA, noccB, FockA, FockB);
    timer_on("I/O");
    t2_1AB->write(psio_, PSIF_DFOCC_AMPS);
    timer_off("I/O");
    t2p_1 = SharedTensor2d(new Tensor2d("T2_1(IA,jb)", naoccA, navirA, naoccB, navirB));
    t2p_1->sort(1324, t2_1AB, 1.0, 0.0);
    t2_1AB.reset();
    timer_on("I/O");
    t2p_1->write(psio_, PSIF_DFOCC_AMPS);
    timer_off("I/O");
    t2p_1.reset();
}// else if (reference_ == "UNRESTRICTED")
    timer_off("Form 1st-order T2");
} // end t2_1st_sc

}} // End Namespaces


