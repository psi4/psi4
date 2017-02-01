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

void DFOCC::gfock_oo()
{

    timer_on("GFM OO");
    SharedTensor2d K;
    SharedTensor2d K2;
    SharedTensor2d G;

if (reference_ == "RESTRICTED") {
    //=========================
    // Reference Contribution
    //=========================
    GFoo->zero();
    GFoo->axpy(FooA, 2.0);

    //=========================
    // Correlation Contribution
    //=========================

    // Fij = \sum_{m} h_im G_mj
    GFoo->gemm(false, false, HooA, G1c_oo, 1.0, 1.0);

    // Fij += \sum_{Q} \sum_{e} G_je^Q b_ie^Q = \sum_{e} G_ej^Q b_ei^Q
    K2 = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA, nvirA));
    K2->read(psio_, PSIF_DFOCC_INTS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VO)", nQ, nvirA, noccA));
    K->swap_3index_col(K2);
    K2.reset();
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VO)", nQ, nvirA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    GFoo->contract(true, false, noccA, noccA, nQ * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    //=========================
    // Separable Part
    //=========================

    // Fij += \sum_{Q} \sum_{m} G_mj^Q b_mi^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OO)", nQ_ref, noccA * noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA * noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFoo->contract(true, false, noccA, noccA, nQ_ref * noccA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fia += \sum_{Q} \sum_{e} G_je^Q b_ie^Q = \sum_{e} G_ej^Q b_ei^Q
    K2 = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    K2->read(psio_, PSIF_DFOCC_INTS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VO)", nQ_ref, nvirA, noccA));
    K->swap_3index_col(K2);
    K2.reset();
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VO)", nQ_ref, nvirA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    GFoo->contract(true, false, noccA, noccA, nQ_ref * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Set global GF
    GF->set_oo(GFoo);

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    //=========================
    // Reference Contribution
    //=========================
    GFooA->copy(FooA);
    GFooB->copy(FooB);

    //=========================
    // Correlation Contribution
    //=========================

    // Fij = \sum_{m} h_im G_mj
    GFooA->gemm(false, false, HooA, G1c_ooA, 1.0, 1.0);
    GFooB->gemm(false, false, HooB, G1c_ooB, 1.0, 1.0);

 if (reference == "ROHF" && orb_opt_ == "FALSE") {
    // Fij = \sum_{e} h_ie G_ej
    GFooA->gemm(false, false, HovA, G1c_voA, 1.0, 1.0);
    GFooB->gemm(false, false, HovB, G1c_voB, 1.0, 1.0);
 }

    // F_IJ += \sum_{Q} \sum_{E} G_JE^Q b_IE^Q = \sum_{E} G_EJ^Q b_EI^Q
    K2 = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA, nvirA));
    K2->read(psio_, PSIF_DFOCC_INTS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VO)", nQ, nvirA, noccA));
    K->swap_3index_col(K2);
    K2.reset();
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VO)", nQ, nvirA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    GFooA->contract(true, false, noccA, noccA, nQ * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fij += \sum_{Q} \sum_{e} G_je^Q b_ie^Q = \sum_{e} G_ej^Q b_ei^Q
    K2 = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ov)", nQ, noccB, nvirB));
    K2->read(psio_, PSIF_DFOCC_INTS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|vo)", nQ, nvirB, noccB));
    K->swap_3index_col(K2);
    K2.reset();
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|vo)", nQ, nvirB, noccB));
    G->read(psio_, PSIF_DFOCC_DENS);
    GFooB->contract(true, false, noccB, noccB, nQ * nvirB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    //=========================
    // Separable Part
    //=========================

    // FIJ += \sum_{Q} \sum_{M} G_MJ^Q b_MI^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OO)", nQ_ref, noccA * noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA * noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFooA->contract(true, false, noccA, noccA, nQ_ref * noccA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fij += \sum_{Q} \sum_{m} G_mj^Q b_mi^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|oo)", nQ_ref, noccB * noccB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB * noccB));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFooB->contract(true, false, noccB, noccB, nQ_ref * noccB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // FIA += \sum_{Q} \sum_{E} G_JE^Q b_IE^Q = \sum_{E} G_EJ^Q b_EI^Q
    K2 = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    K2->read(psio_, PSIF_DFOCC_INTS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VO)", nQ_ref, nvirA, noccA));
    K->swap_3index_col(K2);
    K2.reset();
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VO)", nQ_ref, nvirA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    GFooA->contract(true, false, noccA, noccA, nQ_ref * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fia += \sum_{Q} \sum_{e} G_je^Q b_ie^Q = \sum_{e} G_ej^Q b_ei^Q
    K2 = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB, nvirB));
    K2->read(psio_, PSIF_DFOCC_INTS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vo)", nQ_ref, nvirB, noccB));
    K->swap_3index_col(K2);
    K2.reset();
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|vo)", nQ_ref, nvirB, noccB));
    G->read(psio_, PSIF_DFOCC_DENS);
    GFooB->contract(true, false, noccB, noccB, nQ_ref * nvirB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Set global GF
    GFA->set_oo(GFooA);
    GFB->set_oo(GFooB);
    //GFooA->print();
    //GFooB->print();

}// else if (reference_ == "UNRESTRICTED")
    timer_off("GFM OO");
} // end gfock_oo


//======================================================================
//    CCSD: GFOCK
//======================================================================
void DFOCC::gfock_cc_oo()
{

    timer_on("GFM OO");
    SharedTensor2d K, K2, X, G, G2;

if (reference_ == "RESTRICTED") {
    //=========================
    // Reference Contribution
    //=========================
    GFoo->zero();
    GFoo->axpy(FooA, 2.0);

    //=========================
    // Correlation Contribution
    //=========================

    // Read DF Integrals
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA, noccA));
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA, nvirA));
    bQooA->read(psio_, PSIF_DFOCC_INTS);
    bQovA->read(psio_, PSIF_DFOCC_INTS);

    // Fij = \sum_{m} h_im G_mj
    GFoo->gemm(false, false, HooA, G1c_oo, 1.0, 1.0);

    // Fij += \sum_{e} h_ie G_ej
    GFoo->gemm(false, true, HovA, G1c_ov, 1.0, 1.0);

    // Fij += \sum_{Q} \sum_{m} b_mi^Q G_mj^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OO)", nQ, noccA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    GFoo->contract(true, false, noccA, noccA, nQ * noccA, bQooA, G, 1.0, 1.0);
    G.reset();
    bQooA.reset();

    // Fij += \sum_{Q} \sum_{e} G_je^Q b_ie^Q = \sum_{e} G_ej^Q b_ei^Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VO)", nQ, nvirA, noccA));
    K->swap_3index_col(bQovA);
    bQovA.reset();
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VO)", nQ, nvirA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    GFoo->contract(true, false, noccA, noccA, nQ * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    //=========================
    // Separable Part
    //=========================

    // Read DF Integrals
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    bQooA->read(psio_, PSIF_DFOCC_INTS);
    bQovA->read(psio_, PSIF_DFOCC_INTS);

    // Fij += \sum_{Q} \sum_{m} G_mj^Q b_mi^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OO)", nQ_ref, noccA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    GFoo->contract(true, false, noccA, noccA, nQ_ref * noccA, bQooA, G, 1.0, 1.0);
    G.reset();
    bQooA.reset();

    // Fij += \sum_{Q} \sum_{e} G_je^Q b_ie^Q = \sum_{e} G_ej^Q b_ei^Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VO)", nQ_ref, nvirA, noccA));
    K->swap_3index_col(bQovA);
    bQovA.reset();
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VO)", nQ_ref, nvirA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    GFoo->contract(true, false, noccA, noccA, nQ_ref * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Set global GF
    GF->set_oo(GFoo);

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    //=========================
    // Reference Contribution
    //=========================
    GFooA->copy(FooA);
    GFooB->copy(FooB);

    //=========================
    // Correlation Contribution
    //=========================

    // Fij = \sum_{m} h_im G_mj
    GFooA->gemm(false, false, HooA, G1c_ooA, 1.0, 1.0);
    GFooB->gemm(false, false, HooB, G1c_ooB, 1.0, 1.0);

    // Fij = \sum_{e} h_ie G_ej : For CCSD
    GFooA->gemm(false, false, HovA, G1c_voA, 1.0, 1.0);
    GFooB->gemm(false, false, HovB, G1c_voB, 1.0, 1.0);

    // FIJ += \sum_{Q} \sum_{M} b_MI^Q G_MJ^Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OO)", nQ, noccA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    GFooA->contract(true, false, noccA, noccA, nQ * noccA, K, G, 1.0, 1.0);
    K.reset();
    G.reset();

    // Fij += \sum_{Q} \sum_{m} b_mi^Q G_mj^Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|oo)", nQ, noccB, noccB));
    K->read(psio_, PSIF_DFOCC_INTS);
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|oo)", nQ, noccB, noccB));
    G->read(psio_, PSIF_DFOCC_DENS);
    GFooB->contract(true, false, noccB, noccB, nQ * noccB, K, G, 1.0, 1.0);
    K.reset();
    G.reset();

    // F_IJ += \sum_{Q} \sum_{E} G_JE^Q b_IE^Q = \sum_{E} G_EJ^Q b_EI^Q
    K2 = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA, nvirA));
    K2->read(psio_, PSIF_DFOCC_INTS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VO)", nQ, nvirA, noccA));
    K->swap_3index_col(K2);
    K2.reset();
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VO)", nQ, nvirA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    GFooA->contract(true, false, noccA, noccA, nQ * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fij += \sum_{Q} \sum_{e} G_je^Q b_ie^Q = \sum_{e} G_ej^Q b_ei^Q
    K2 = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ov)", nQ, noccB, nvirB));
    K2->read(psio_, PSIF_DFOCC_INTS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|vo)", nQ, nvirB, noccB));
    K->swap_3index_col(K2);
    K2.reset();
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|vo)", nQ, nvirB, noccB));
    G->read(psio_, PSIF_DFOCC_DENS);
    GFooB->contract(true, false, noccB, noccB, nQ * nvirB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    //=========================
    // Separable Part
    //=========================

    // FIJ += \sum_{Q} \sum_{M} G_MJ^Q b_MI^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OO)", nQ_ref, noccA * noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA * noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFooA->contract(true, false, noccA, noccA, nQ_ref * noccA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fij += \sum_{Q} \sum_{m} G_mj^Q b_mi^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|oo)", nQ_ref, noccB * noccB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB * noccB));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFooB->contract(true, false, noccB, noccB, nQ_ref * noccB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // FIA += \sum_{Q} \sum_{E} G_JE^Q b_IE^Q = \sum_{E} G_EJ^Q b_EI^Q
    K2 = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    K2->read(psio_, PSIF_DFOCC_INTS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VO)", nQ_ref, nvirA, noccA));
    K->swap_3index_col(K2);
    K2.reset();
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VO)", nQ_ref, nvirA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    GFooA->contract(true, false, noccA, noccA, nQ_ref * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fia += \sum_{Q} \sum_{e} G_je^Q b_ie^Q = \sum_{e} G_ej^Q b_ei^Q
    K2 = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB, nvirB));
    K2->read(psio_, PSIF_DFOCC_INTS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vo)", nQ_ref, nvirB, noccB));
    K->swap_3index_col(K2);
    K2.reset();
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|vo)", nQ_ref, nvirB, noccB));
    G->read(psio_, PSIF_DFOCC_DENS);
    GFooB->contract(true, false, noccB, noccB, nQ_ref * nvirB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Set global GF
    GFA->set_oo(GFooA);
    GFB->set_oo(GFooB);
    //GFooA->print();
    //GFooB->print();

}// else if (reference_ == "UNRESTRICTED")
    timer_off("GFM OO");
} // end gfock_cc_oo


}} // End Namespaces
