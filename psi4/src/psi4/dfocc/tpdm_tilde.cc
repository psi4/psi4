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

void DFOCC::tpdm_tilde()
{
    outfile->Printf("\tForming Gamma^tilde...\n");


    SharedTensor2d G, G2, Gref, Gsep, Gcorr;
    timer_on("tpdm_tilde");
if (reference_ == "RESTRICTED") {

    //=========================
    // Reference TPDM
    //=========================
    // Read J-1/2
    Jmhalf = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Jmhalf <P|Q>", nQ_ref, nQ_ref));
    Jmhalf->read(psio_, PSIF_DFOCC_INTS);

    G = SharedTensor2d(new Tensor2d("Reference 3-Index TPDM <Q|OO>", nQ_ref, noccA, noccA));
    Gref = SharedTensor2d(new Tensor2d("Reference 3-Index TPDM (Q|OO)", nQ_ref, noccA, noccA));
    Gref->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gref, 1.0, 0.0);
    Gref.reset();
    G->write(psio_, PSIF_DFOCC_DENS);
    G.reset();

    //=========================
    // Seprable TPDM
    //=========================
    // OO Block
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|OO>", nQ_ref, noccA, noccA));
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OO)", nQ_ref, noccA, noccA));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gsep, 1.0, 0.0);
    Gsep.reset();
    G->write(psio_, PSIF_DFOCC_DENS);
    G.reset();

    // OV Block
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|OV>", nQ_ref, noccA, nvirA));
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OV)", nQ_ref, noccA, nvirA));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gsep, 1.0, 0.0);
    Gsep.reset();
    G->write(psio_, PSIF_DFOCC_DENS);

    // VO Block
    G2 = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|VO>", nQ_ref, nvirA, noccA));
    G2->swap_3index_col(G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS);
    G2.reset();

    // VV Block
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|VV>", nQ_ref, nvirA, nvirA));
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VV)", nQ_ref, nvirA, nvirA));
    Gsep->read(psio_, PSIF_DFOCC_DENS, true, true);
    G->gemm(true, false, Jmhalf, Gsep, 1.0, 0.0);
    Gsep.reset();
    G->write(psio_, PSIF_DFOCC_DENS, true, true);
    G.reset();

    // Free J-1/2
    Jmhalf.reset();

    //=========================
    // Correlation TPDM
    //=========================
    // Read J-1/2
    Jmhalf = SharedTensor2d(new Tensor2d("DF_BASIS_CC Jmhalf <P|Q>", nQ, nQ));
    Jmhalf->read(psio_, PSIF_DFOCC_INTS);

    // OV Block
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|OV>", nQ, noccA, nvirA));
    Gcorr = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA, nvirA));
    Gcorr->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gcorr, 1.0, 0.0);
    Gcorr.reset();
    G->write(psio_, PSIF_DFOCC_DENS);

    // VO Block
    G2 = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|VO>", nQ, nvirA, noccA));
    G2->swap_3index_col(G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS);
    G2.reset();

    // Free J-1/2
    Jmhalf.reset();

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {

    //=========================
    // Reference TPDM
    //=========================
    // Read J-1/2
    Jmhalf = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Jmhalf <P|Q>", nQ_ref, nQ_ref));
    Jmhalf->read(psio_, PSIF_DFOCC_INTS);

    // OO Block
    G = SharedTensor2d(new Tensor2d("Reference 3-Index TPDM <Q|OO>", nQ_ref, noccA, noccA));
    Gref = SharedTensor2d(new Tensor2d("Reference 3-Index TPDM (Q|OO)", nQ_ref, noccA, noccA));
    Gref->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gref, 1.0, 0.0);
    Gref.reset();
    G->write(psio_, PSIF_DFOCC_DENS);
    G.reset();

    // oo Block
    G = SharedTensor2d(new Tensor2d("Reference 3-Index TPDM <Q|oo>", nQ_ref, noccB, noccB));
    Gref = SharedTensor2d(new Tensor2d("Reference 3-Index TPDM (Q|oo)", nQ_ref, noccB, noccB));
    Gref->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gref, 1.0, 0.0);
    Gref.reset();
    G->write(psio_, PSIF_DFOCC_DENS);
    G.reset();

    //=========================
    // Seprable TPDM
    //=========================
    // OO Block
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|OO>", nQ_ref, noccA, noccA));
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OO)", nQ_ref, noccA, noccA));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gsep, 1.0, 0.0);
    Gsep.reset();
    G->write(psio_, PSIF_DFOCC_DENS);
    G.reset();

    // oo Block
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|oo>", nQ_ref, noccB, noccB));
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|oo)", nQ_ref, noccB, noccB));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gsep, 1.0, 0.0);
    Gsep.reset();
    G->write(psio_, PSIF_DFOCC_DENS);
    G.reset();

    // OV Block
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|OV>", nQ_ref, noccA, nvirA));
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OV)", nQ_ref, noccA, nvirA));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gsep, 1.0, 0.0);
    Gsep.reset();
    G->write(psio_, PSIF_DFOCC_DENS);

    // VO Block
    G2 = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|VO>", nQ_ref, nvirA, noccA));
    G2->swap_3index_col(G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS);
    G2.reset();

    // ov Block
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|ov>", nQ_ref, noccB, nvirB));
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|ov)", nQ_ref, noccB, nvirB));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gsep, 1.0, 0.0);
    Gsep.reset();
    G->write(psio_, PSIF_DFOCC_DENS);

    // vo Block
    G2 = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|vo>", nQ_ref, nvirB, noccB));
    G2->swap_3index_col(G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS);
    G2.reset();

    // VV Block
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|VV>", nQ_ref, nvirA, nvirA));
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VV)", nQ_ref, nvirA, nvirA));
    Gsep->read(psio_, PSIF_DFOCC_DENS, true, true);
    G->gemm(true, false, Jmhalf, Gsep, 1.0, 0.0);
    Gsep.reset();
    G->write(psio_, PSIF_DFOCC_DENS, true, true);
    G.reset();

    // vv Block
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|vv>", nQ_ref, nvirB, nvirB));
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|vv)", nQ_ref, nvirB, nvirB));
    Gsep->read(psio_, PSIF_DFOCC_DENS, true, true);
    G->gemm(true, false, Jmhalf, Gsep, 1.0, 0.0);
    Gsep.reset();
    G->write(psio_, PSIF_DFOCC_DENS, true, true);
    G.reset();

    // Free J-1/2
    Jmhalf.reset();

    //=========================
    // Correlation TPDM
    //=========================
    // Read J-1/2
    Jmhalf = SharedTensor2d(new Tensor2d("DF_BASIS_CC Jmhalf <P|Q>", nQ, nQ));
    Jmhalf->read(psio_, PSIF_DFOCC_INTS);

    // OV Block
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|OV>", nQ, noccA, nvirA));
    Gcorr = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA, nvirA));
    Gcorr->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gcorr, 1.0, 0.0);
    Gcorr.reset();
    G->write(psio_, PSIF_DFOCC_DENS);

    // VO Block
    G2 = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|VO>", nQ, nvirA, noccA));
    G2->swap_3index_col(G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS);
    G2.reset();

    // ov Block
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|ov>", nQ, noccB, nvirB));
    Gcorr = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|ov)", nQ, noccB, nvirB));
    Gcorr->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gcorr, 1.0, 0.0);
    Gcorr.reset();
    G->write(psio_, PSIF_DFOCC_DENS);

    // vo Block
    G2 = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|vo>", nQ, nvirB, noccB));
    G2->swap_3index_col(G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS);
    G2.reset();

    // Free J-1/2
    Jmhalf.reset();

}// else if (reference_ == "UNRESTRICTED")
    timer_off("tpdm_tilde");
} // end tpdm_tilde

//======================================================================
//    CCSD: Gamma Tilde
//======================================================================
void DFOCC::tpdm_tilde_cc()
{
    outfile->Printf("\tForming Gamma^tilde...\n");


    SharedTensor2d G, G2, Gref, Gsep, Gcorr;
    timer_on("tpdm_tilde");
if (reference_ == "RESTRICTED") {

    //=========================
    // Reference TPDM
    //=========================
    // Read J-1/2
    Jmhalf = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Jmhalf <P|Q>", nQ_ref, nQ_ref));
    Jmhalf->read(psio_, PSIF_DFOCC_INTS);

    G = SharedTensor2d(new Tensor2d("Reference 3-Index TPDM <Q|OO>", nQ_ref, noccA, noccA));
    Gref = SharedTensor2d(new Tensor2d("Reference 3-Index TPDM (Q|OO)", nQ_ref, noccA, noccA));
    Gref->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gref, 1.0, 0.0);
    Gref.reset();
    G->write(psio_, PSIF_DFOCC_DENS);
    G.reset();

    //=========================
    // Seprable TPDM
    //=========================
    // OO Block
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|OO>", nQ_ref, noccA, noccA));
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OO)", nQ_ref, noccA, noccA));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gsep, 1.0, 0.0);
    Gsep.reset();
    G->write(psio_, PSIF_DFOCC_DENS);
    G.reset();

    // OV Block
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|OV>", nQ_ref, noccA, nvirA));
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OV)", nQ_ref, noccA, nvirA));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gsep, 1.0, 0.0);
    Gsep.reset();
    G->write(psio_, PSIF_DFOCC_DENS);

    // VO Block
    G2 = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|VO>", nQ_ref, nvirA, noccA));
    G2->swap_3index_col(G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS);
    G2.reset();

    // VV Block
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|VV>", nQ_ref, nvirA, nvirA));
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VV)", nQ_ref, nvirA, nvirA));
    Gsep->read(psio_, PSIF_DFOCC_DENS, true, true);
    G->gemm(true, false, Jmhalf, Gsep, 1.0, 0.0);
    Gsep.reset();
    G->write(psio_, PSIF_DFOCC_DENS, true, true);
    G.reset();

    // Free J-1/2
    Jmhalf.reset();

    //=========================
    // Correlation TPDM
    //=========================
    // Read J-1/2
    Jmhalf = SharedTensor2d(new Tensor2d("DF_BASIS_CC Jmhalf <P|Q>", nQ, nQ));
    Jmhalf->read(psio_, PSIF_DFOCC_INTS);

    // OO Block
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|OO>", nQ, noccA, noccA));
    Gcorr = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OO)", nQ, noccA, noccA));
    Gcorr->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gcorr, 1.0, 0.0);
    Gcorr.reset();
    G->write(psio_, PSIF_DFOCC_DENS);
    G.reset();

    // OV Block
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|OV>", nQ, noccA, nvirA));
    Gcorr = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA, nvirA));
    Gcorr->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gcorr, 1.0, 0.0);
    Gcorr.reset();
    G->write(psio_, PSIF_DFOCC_DENS);

    // VO Block
    G2 = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|VO>", nQ, nvirA, noccA));
    G2->swap_3index_col(G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS);
    G2.reset();

    // VV Block
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|VV>", nQ, nvirA, nvirA));
    Gcorr = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VV)", nQ, nvirA, nvirA));
    Gcorr->read(psio_, PSIF_DFOCC_DENS, true, true);
    G->gemm(true, false, Jmhalf, Gcorr, 1.0, 0.0);
    Gcorr.reset();
    G->write(psio_, PSIF_DFOCC_DENS, true, true);
    G.reset();

    // Free J-1/2
    Jmhalf.reset();

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {

    //=========================
    // Reference TPDM
    //=========================
    // Read J-1/2
    Jmhalf = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Jmhalf <P|Q>", nQ_ref, nQ_ref));
    Jmhalf->read(psio_, PSIF_DFOCC_INTS);

    // OO Block
    G = SharedTensor2d(new Tensor2d("Reference 3-Index TPDM <Q|OO>", nQ_ref, noccA, noccA));
    Gref = SharedTensor2d(new Tensor2d("Reference 3-Index TPDM (Q|OO)", nQ_ref, noccA, noccA));
    Gref->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gref, 1.0, 0.0);
    Gref.reset();
    G->write(psio_, PSIF_DFOCC_DENS);
    G.reset();

    // oo Block
    G = SharedTensor2d(new Tensor2d("Reference 3-Index TPDM <Q|oo>", nQ_ref, noccB, noccB));
    Gref = SharedTensor2d(new Tensor2d("Reference 3-Index TPDM (Q|oo)", nQ_ref, noccB, noccB));
    Gref->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gref, 1.0, 0.0);
    Gref.reset();
    G->write(psio_, PSIF_DFOCC_DENS);
    G.reset();

    //=========================
    // Seprable TPDM
    //=========================
    // OO Block
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|OO>", nQ_ref, noccA, noccA));
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OO)", nQ_ref, noccA, noccA));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gsep, 1.0, 0.0);
    Gsep.reset();
    G->write(psio_, PSIF_DFOCC_DENS);
    G.reset();

    // oo Block
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|oo>", nQ_ref, noccB, noccB));
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|oo)", nQ_ref, noccB, noccB));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gsep, 1.0, 0.0);
    Gsep.reset();
    G->write(psio_, PSIF_DFOCC_DENS);
    G.reset();

    // OV Block
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|OV>", nQ_ref, noccA, nvirA));
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OV)", nQ_ref, noccA, nvirA));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gsep, 1.0, 0.0);
    Gsep.reset();
    G->write(psio_, PSIF_DFOCC_DENS);

    // VO Block
    G2 = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|VO>", nQ_ref, nvirA, noccA));
    G2->swap_3index_col(G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS);
    G2.reset();

    // ov Block
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|ov>", nQ_ref, noccB, nvirB));
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|ov)", nQ_ref, noccB, nvirB));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gsep, 1.0, 0.0);
    Gsep.reset();
    G->write(psio_, PSIF_DFOCC_DENS);

    // vo Block
    G2 = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|vo>", nQ_ref, nvirB, noccB));
    G2->swap_3index_col(G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS);
    G2.reset();

    // VV Block
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|VV>", nQ_ref, nvirA, nvirA));
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VV)", nQ_ref, nvirA, nvirA));
    Gsep->read(psio_, PSIF_DFOCC_DENS, true, true);
    G->gemm(true, false, Jmhalf, Gsep, 1.0, 0.0);
    Gsep.reset();
    G->write(psio_, PSIF_DFOCC_DENS, true, true);
    G.reset();

    // vv Block
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM <Q|vv>", nQ_ref, nvirB, nvirB));
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|vv)", nQ_ref, nvirB, nvirB));
    Gsep->read(psio_, PSIF_DFOCC_DENS, true, true);
    G->gemm(true, false, Jmhalf, Gsep, 1.0, 0.0);
    Gsep.reset();
    G->write(psio_, PSIF_DFOCC_DENS, true, true);
    G.reset();

    // Free J-1/2
    Jmhalf.reset();

    //=========================
    // Correlation TPDM
    //=========================
    // Read J-1/2
    Jmhalf = SharedTensor2d(new Tensor2d("DF_BASIS_CC Jmhalf <P|Q>", nQ, nQ));
    Jmhalf->read(psio_, PSIF_DFOCC_INTS);

    // OO Block
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|OO>", nQ, noccA, noccA));
    Gcorr = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OO)", nQ, noccA, noccA));
    Gcorr->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gcorr, 1.0, 0.0);
    Gcorr.reset();
    G->write(psio_, PSIF_DFOCC_DENS);
    G.reset();

    // OV Block
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|OV>", nQ, noccA, nvirA));
    Gcorr = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA, nvirA));
    Gcorr->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gcorr, 1.0, 0.0);
    Gcorr.reset();
    G->write(psio_, PSIF_DFOCC_DENS);

    // VO Block
    G2 = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|VO>", nQ, nvirA, noccA));
    G2->swap_3index_col(G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS);
    G2.reset();

    // VV Block
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|VV>", nQ, nvirA, nvirA));
    Gcorr = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VV)", nQ, nvirA, nvirA));
    Gcorr->read(psio_, PSIF_DFOCC_DENS, true, true);
    G->gemm(true, false, Jmhalf, Gcorr, 1.0, 0.0);
    Gcorr.reset();
    G->write(psio_, PSIF_DFOCC_DENS, true, true);
    G.reset();

    // oo Block
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|oo>", nQ, noccB, noccB));
    Gcorr = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|oo)", nQ, noccB, noccB));
    Gcorr->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gcorr, 1.0, 0.0);
    Gcorr.reset();
    G->write(psio_, PSIF_DFOCC_DENS);
    G.reset();

    // ov Block
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|ov>", nQ, noccB, nvirB));
    Gcorr = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|ov)", nQ, noccB, nvirB));
    Gcorr->read(psio_, PSIF_DFOCC_DENS);
    G->gemm(true, false, Jmhalf, Gcorr, 1.0, 0.0);
    Gcorr.reset();
    G->write(psio_, PSIF_DFOCC_DENS);

    // vo Block
    G2 = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|vo>", nQ, nvirB, noccB));
    G2->swap_3index_col(G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS);
    G2.reset();

    // vv Block
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM <Q|vv>", nQ, nvirB, nvirB));
    Gcorr = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|vv)", nQ, nvirB, nvirB));
    Gcorr->read(psio_, PSIF_DFOCC_DENS, true, true);
    G->gemm(true, false, Jmhalf, Gcorr, 1.0, 0.0);
    Gcorr.reset();
    G->write(psio_, PSIF_DFOCC_DENS, true, true);
    G.reset();

    // Free J-1/2
    Jmhalf.reset();

}// else if (reference_ == "UNRESTRICTED")
    timer_off("tpdm_tilde");
} // end tpdm_tilde_cc


}} // End Namespaces
