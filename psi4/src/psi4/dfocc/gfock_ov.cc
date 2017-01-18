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

void DFOCC::gfock_ov()
{

    timer_on("GFM OV");
    SharedTensor2d K;
    SharedTensor2d K2;
    SharedTensor2d G;

if (reference_ == "RESTRICTED") {
    //=========================
    // Correlation Contribution
    //=========================

    // Fia = \sum_{e} h_ie G_ea
    GFov->gemm(true, false, HvoA, G1c_vv, 1.0, 0.0);

    // Fia += \sum_{Q} \sum_{m} G_ma^Q b_mi^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA, nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA * noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFov->contract(true, false, noccA, nvirA, nQ * noccA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    //=========================
    // Separable Part
    //=========================

    // Fia += \sum_{Q} \sum_{m} G_ma^Q b_mi^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OV)", nQ_ref, noccA * nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA * noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFov->contract(true, false, noccA, nvirA, nQ_ref * noccA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fia += \sum_{Q} \sum_{e} G_ea^Q b_ei^Q
    K2 = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    K2->read(psio_, PSIF_DFOCC_INTS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VO)", nQ_ref, nvirA, noccA));
    K->swap_3index_col(K2);
    K2.reset();
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VV)", nQ_ref, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    GFov->contract(true, false, noccA, nvirA, nQ_ref * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Set global GF
    GF->set_ov(GFov);

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    //=========================
    // Correlation Contribution
    //=========================

    // Fia = \sum_{e} h_ie G_ea
    GFovA->gemm(true, false, HvoA, G1c_vvA, 1.0, 0.0);
    GFovB->gemm(true, false, HvoB, G1c_vvB, 1.0, 0.0);

 if (reference == "ROHF" && orb_opt_ == "FALSE") {
    // Fia = \sum_{m} h_im G_ma
    GFovA->gemm(false, false, HooA, G1c_ovA, 1.0, 1.0);
    GFovB->gemm(false, false, HooB, G1c_ovB, 1.0, 1.0);
 }

    // F_IA += \sum_{Q} \sum_{M} G_MA^Q b_MI^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA, nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA * noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFovA->contract(true, false, noccA, nvirA, nQ * noccA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fia += \sum_{Q} \sum_{m} G_ma^Q b_mi^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|ov)", nQ, noccB, nvirB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|oo)", nQ, noccB * noccB));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFovB->contract(true, false, noccB, nvirB, nQ * noccB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    //=========================
    // Separable Part
    //=========================

    // FIA += \sum_{Q} \sum_{M} G_MA^Q b_MI^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OV)", nQ_ref, noccA * nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA * noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFovA->contract(true, false, noccA, nvirA, nQ_ref * noccA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fia += \sum_{Q} \sum_{m} G_ma^Q b_mi^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|ov)", nQ_ref, noccB * nvirB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB * noccB));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFovB->contract(true, false, noccB, nvirB, nQ_ref * noccB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // F_IA += \sum_{Q} \sum_{E} G_EA^Q b_EI^Q
    K2 = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    K2->read(psio_, PSIF_DFOCC_INTS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VO)", nQ_ref, nvirA, noccA));
    K->swap_3index_col(K2);
    K2.reset();
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VV)", nQ_ref, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    GFovA->contract(true, false, noccA, nvirA, nQ_ref * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fia += \sum_{Q} \sum_{e} G_ea^Q b_ei^Q
    K2 = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB, nvirB));
    K2->read(psio_, PSIF_DFOCC_INTS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vo)", nQ_ref, nvirB, noccB));
    K->swap_3index_col(K2);
    K2.reset();
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|vv)", nQ_ref, nvirB, nvirB));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    GFovB->contract(true, false, noccB, nvirB, nQ_ref * nvirB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Set global GF
    GFA->set_ov(GFovA);
    GFB->set_ov(GFovB);
    //GFovA->print();
    //GFovB->print();

}// else if (reference_ == "UNRESTRICTED")
    timer_off("GFM OV");
} // end gfock_ov

//======================================================================
//    CCSD: GFOCK
//======================================================================
void DFOCC::gfock_cc_ov()
{

    timer_on("GFM OV");
    SharedTensor2d K;
    SharedTensor2d K2;
    SharedTensor2d G;

if (reference_ == "RESTRICTED") {
    //=========================
    // Correlation Contribution
    //=========================

    // Fia = \sum_{e} h_ie G_ea
    GFov->gemm(false, false, HovA, G1c_vv, 1.0, 0.0);

    // Fia = \sum_{m} h_im G_ma
    GFov->gemm(false, false, HooA, G1c_ov, 1.0, 1.0);

    // Fia += \sum_{Q} \sum_{m} G_ma^Q b_mi^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA, nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFov->contract(true, false, noccA, nvirA, nQ * noccA, K, G, 1.0, 1.0);
    //GFov->cont332("IA", "MI", "MA", false, false, K, G, 1.0, 1.0); // it works
    G.reset();
    K.reset();

    // Fia += \sum_{Q} \sum_{e} G_ea^Q b_ei^Q
    K2 = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA, nvirA));
    K2->read(psio_, PSIF_DFOCC_INTS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VO)", nQ, nvirA, noccA));
    K->swap_3index_col(K2);
    K2.reset();
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VV)", nQ, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    GFov->contract(true, false, noccA, nvirA, nQ * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    //=========================
    // Separable Part
    //=========================

    // Fia += \sum_{Q} \sum_{m} G_ma^Q b_mi^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OV)", nQ_ref, noccA * nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFov->contract(true, false, noccA, nvirA, nQ_ref * noccA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fia += \sum_{Q} \sum_{e} G_ea^Q b_ei^Q
    K2 = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    K2->read(psio_, PSIF_DFOCC_INTS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VO)", nQ_ref, nvirA, noccA));
    K->swap_3index_col(K2);
    K2.reset();
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VV)", nQ_ref, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    GFov->contract(true, false, noccA, nvirA, nQ_ref * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Set global GF
    GF->set_ov(GFov);

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    //=========================
    // Correlation Contribution
    //=========================

    // Fia = \sum_{e} h_ie G_ea
    GFovA->gemm(true, false, HvoA, G1c_vvA, 1.0, 0.0);
    GFovB->gemm(true, false, HvoB, G1c_vvB, 1.0, 0.0);

    // Fia = \sum_{m} h_im G_ma
    GFovA->gemm(false, false, HooA, G1c_ovA, 1.0, 1.0);
    GFovB->gemm(false, false, HooB, G1c_ovB, 1.0, 1.0);

    // F_IA += \sum_{Q} \sum_{M} G_MA^Q b_MI^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA, nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA * noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFovA->contract(true, false, noccA, nvirA, nQ * noccA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fia += \sum_{Q} \sum_{m} G_ma^Q b_mi^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|ov)", nQ, noccB, nvirB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|oo)", nQ, noccB * noccB));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFovB->contract(true, false, noccB, nvirB, nQ * noccB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // F_IA += \sum_{Q} \sum_{E} G_EA^Q b_EI^Q
    K2 = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA, nvirA));
    K2->read(psio_, PSIF_DFOCC_INTS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VO)", nQ, nvirA, noccA));
    K->swap_3index_col(K2);
    K2.reset();
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VV)", nQ, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    GFovA->contract(true, false, noccA, nvirA, nQ * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fia += \sum_{Q} \sum_{e} G_ea^Q b_ei^Q
    K2 = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ov)", nQ, noccB, nvirB));
    K2->read(psio_, PSIF_DFOCC_INTS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|vo)", nQ, nvirB, noccB));
    K->swap_3index_col(K2);
    K2.reset();
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|vv)", nQ, nvirB, nvirB));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    GFovB->contract(true, false, noccB, nvirB, nQ * nvirB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    //=========================
    // Separable Part
    //=========================

    // FIA += \sum_{Q} \sum_{M} G_MA^Q b_MI^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OV)", nQ_ref, noccA * nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA * noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFovA->contract(true, false, noccA, nvirA, nQ_ref * noccA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fia += \sum_{Q} \sum_{m} G_ma^Q b_mi^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|ov)", nQ_ref, noccB * nvirB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB * noccB));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFovB->contract(true, false, noccB, nvirB, nQ_ref * noccB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // F_IA += \sum_{Q} \sum_{E} G_EA^Q b_EI^Q
    K2 = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    K2->read(psio_, PSIF_DFOCC_INTS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VO)", nQ_ref, nvirA, noccA));
    K->swap_3index_col(K2);
    K2.reset();
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VV)", nQ_ref, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    GFovA->contract(true, false, noccA, nvirA, nQ_ref * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fia += \sum_{Q} \sum_{e} G_ea^Q b_ei^Q
    K2 = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB, nvirB));
    K2->read(psio_, PSIF_DFOCC_INTS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vo)", nQ_ref, nvirB, noccB));
    K->swap_3index_col(K2);
    K2.reset();
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|vv)", nQ_ref, nvirB, nvirB));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    GFovB->contract(true, false, noccB, nvirB, nQ_ref * nvirB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Set global GF
    GFA->set_ov(GFovA);
    GFB->set_ov(GFovB);
    //GFovA->print();
    //GFovB->print();

}// else if (reference_ == "UNRESTRICTED")
    timer_off("GFM OV");
} // end gfock_cc_ov

}} // End Namespaces
