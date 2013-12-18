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
  
void DFOCC::omp2_tpdm()
{   

    timer_on("tpdm");
if (reference_ == "RESTRICTED") {
    // G_ia^Q = 2\sum_{m,e} b_me^Q (2t_im^ae - t_mi^ae)  
    G2c_ia = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|IA)", nQ, naoccA * navirA));
    u2p_1 = SharedTensor2d(new Tensor2d("2*T2_1(ia,jb) - T2_1(ib,ja)", naoccA, navirA, naoccA, navirA));
    bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA * navirA));
    timer_on("I/O");
    u2p_1->read(psio_, PSIF_DFOCC_AMPS);
    bQiaA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    G2c_ia->gemm(false, false, bQiaA, u2p_1, 2.0, 0.0);
    u2p_1.reset();
    bQiaA.reset();
    //G2c_ov = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA * nvirA));
    G2c_ov = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA, nvirA));
    G2c_ov->set3_act_ov(nfrzc, naoccA, navirA, nvirA, G2c_ia);
    G2c_ia.reset();
    timer_on("I/O");
    G2c_ov->write(psio_, PSIF_DFOCC_DENS);
    timer_off("I/O");
    if(print_ > 3) G2c_ov->print();

    // Form G_vo^Q
    G2c_vo = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VO)", nQ, nvirA, noccA));
    G2c_vo->swap_3index_col(G2c_ov);
    G2c_ov.reset();
    timer_on("I/O");
    G2c_vo->write(psio_, PSIF_DFOCC_DENS);
    timer_off("I/O");
    if(print_ > 3) G2c_vo->print();
    G2c_vo.reset();

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // G_IA^Q = \sum_{M,E} b_ME^Q t_IM^AE 
    G2c_iaA = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|IA)", nQ, naoccA * navirA));
    t2p_1 = SharedTensor2d(new Tensor2d("T2_1(IA,JB)", naoccA, navirA, naoccA, navirA));
    bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA * navirA));
    timer_on("I/O");
    t2p_1->read(psio_, PSIF_DFOCC_AMPS);
    bQiaA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    G2c_iaA->gemm(false, false, bQiaA, t2p_1, 1.0, 0.0);
    t2p_1.reset();
    bQiaA.reset();

    // G_IA^Q = \sum_{m,e} b_me^Q t_Im^Ae 
    t2p_1 = SharedTensor2d(new Tensor2d("T2_1(IA,jb)", naoccA, navirA, naoccB, navirB));
    bQiaB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ia)", nQ, naoccB * navirB));
    timer_on("I/O");
    t2p_1->read(psio_, PSIF_DFOCC_AMPS);
    bQiaB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    G2c_iaA->gemm(false, true, bQiaB, t2p_1, 1.0, 1.0);
    t2p_1.reset();
    bQiaB.reset();

    // G_IA^Q -> G_OV^Q
    //G2c_ovA = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA * nvirA));
    G2c_ovA = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA, nvirA));
    G2c_ovA->set3_act_ov(nfrzc, naoccA, navirA, nvirA, G2c_iaA);
    G2c_iaA.reset();
    timer_on("I/O");
    G2c_ovA->write(psio_, PSIF_DFOCC_DENS);
    timer_off("I/O");
    if(print_ > 3) G2c_ovA->print();
    //G2c_ovA.reset();

    // Form G_VO^Q
    G2c_voA = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VO)", nQ, nvirA, noccA));
    G2c_voA->swap_3index_col(G2c_ovA);
    G2c_ovA.reset();
    timer_on("I/O");
    G2c_voA->write(psio_, PSIF_DFOCC_DENS);
    timer_off("I/O");
    if(print_ > 3) G2c_voA->print();
    G2c_voA.reset();

    // G_ia^Q = \sum_{m,e} b_me^Q t_im^ae 
    G2c_iaB = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|ia)", nQ, naoccB * navirB));
    t2p_1 = SharedTensor2d(new Tensor2d("T2_1(ia,jb)", naoccB, navirB, naoccB, navirB));
    bQiaB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ia)", nQ, naoccB * navirB));
    timer_on("I/O");
    t2p_1->read(psio_, PSIF_DFOCC_AMPS);
    bQiaB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    G2c_iaB->gemm(false, false, bQiaB, t2p_1, 1.0, 0.0);
    t2p_1.reset();
    bQiaB.reset();
    //fprintf(outfile,"\tI am here.\n"); fflush(outfile);

    // G_ia^Q = \sum_{M,E} b_ME^Q t_Mi^Ea 
    t2p_1 = SharedTensor2d(new Tensor2d("T2_1(IA,jb)", naoccA, navirA, naoccB, navirB));
    bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA * navirA));
    timer_on("I/O");
    t2p_1->read(psio_, PSIF_DFOCC_AMPS);
    bQiaA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    G2c_iaB->gemm(false, false, bQiaA, t2p_1, 1.0, 1.0);
    t2p_1.reset();
    bQiaA.reset();

    // G_ia^Q -> G_ov^Q
    //G2c_ovB = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|ov)", nQ, noccB * nvirB));
    G2c_ovB = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|ov)", nQ, noccB, nvirB));
    G2c_ovB->set3_act_ov(nfrzc, naoccB, navirB, nvirB, G2c_iaB);
    G2c_iaB.reset();
    timer_on("I/O");
    G2c_ovB->write(psio_, PSIF_DFOCC_DENS);
    timer_off("I/O");
    if(print_ > 3) G2c_ovB->print();
    //G2c_ovB.reset();

    // Form G_vo^Q
    G2c_voB = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|vo)", nQ, nvirB, noccB));
    G2c_voB->swap_3index_col(G2c_ovB);
    G2c_ovB.reset();
    timer_on("I/O");
    G2c_voB->write(psio_, PSIF_DFOCC_DENS);
    timer_off("I/O");
    if(print_ > 3) G2c_voB->print();
    G2c_voB.reset();

}// else if (reference_ == "UNRESTRICTED")
    timer_off("tpdm");
} // end omp2_tpdm

}} // End Namespaces


