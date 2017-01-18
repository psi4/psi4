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

void DFOCC::gfock_vo()
{

    timer_on("GFM VO");
    SharedTensor2d K;
    SharedTensor2d G;

if (reference_ == "RESTRICTED") {
    //=========================
    // Reference Contribution
    //=========================
    GFvo->zero();
    GFvo->axpy(FvoA, 2.0);

    //=========================
    // Correlation Contribution
    //=========================

    // Fai = \sum_{m} h_am G_mi
    GFvo->gemm(false, false, HvoA, G1c_oo, 1.0, 1.0);

    // Fai += \sum_{Q} \sum_{e} G_ie^Q b_ae^Q = \sum_{e} G_ei^Q b_ea^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VO)", nQ, nvirA, noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    GFvo->contract(true, false, nvirA, noccA, nQ * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    //=========================
    // Separable Part
    //=========================

    // Fai += \sum_{Q} \sum_{m} G_mi^Q b_ma^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OO)", nQ_ref, noccA * noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA * nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFvo->contract(true, false, nvirA, noccA, nQ_ref * noccA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fai += \sum_{Q} \sum_{e} G_ie^Q b_ae^Q = \sum_{e} G_ei^Q b_ea^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VO)", nQ_ref, nvirA, noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    GFvo->contract(true, false, nvirA, noccA, nQ_ref * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Set global GF
    GF->set_vo(GFvo);

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    //=========================
    // Reference Contribution
    //=========================
    GFvoA->zero();
    GFvoB->zero();
    GFvoA->copy(FvoA);
    GFvoB->copy(FvoB);

    //=========================
    // Correlation Contribution
    //=========================

    // F_AI = \sum_{M} h_AM G_Mi
    GFvoA->gemm(false, false, HvoA, G1c_ooA, 1.0, 1.0);

 if (reference == "ROHF" && orb_opt_ == "FALSE") {
    // Fai = \sum_{e} h_ae G_ei
    GFvoA->gemm(false, false, HvvA, G1c_voA, 1.0, 1.0);
    GFvoB->gemm(false, false, HvvB, G1c_voB, 1.0, 1.0);
 }

    // F_AI += \sum_{Q} \sum_{E} G_IE^Q b_AE^Q = \sum_{E} G_EI^Q b_EA^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VO)", nQ, nvirA, noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    GFvoA->contract(true, false, nvirA, noccA, nQ * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fai = \sum_{m} h_am G_mi
    GFvoB->gemm(false, false, HvoB, G1c_ooB, 1.0, 1.0);

    // Fai += \sum_{Q} \sum_{e} G_ie^Q b_ae^Q = \sum_{e} G_ei^Q b_ea^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|vo)", nQ, nvirB, noccB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|vv)", nQ, nvirB, nvirB));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    GFvoB->contract(true, false, nvirB, noccB, nQ * nvirB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    //=========================
    // Separable Part
    //=========================

    // F_AI += \sum_{Q} \sum_{M} G_MI^Q b_MA^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OO)", nQ_ref, noccA * noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA * nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFvoA->contract(true, false, nvirA, noccA, nQ_ref * noccA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // F_AI += \sum_{Q} \sum_{E} G_IE^Q b_AE^Q = \sum_{E} G_EI^Q b_EA^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VO)", nQ_ref, nvirA, noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    GFvoA->contract(true, false, nvirA, noccA, nQ_ref * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fai += \sum_{Q} \sum_{m} G_mi^Q b_ma^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|oo)", nQ_ref, noccB * noccB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB * nvirB));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFvoB->contract(true, false, nvirB, noccB, nQ_ref * noccB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fai += \sum_{Q} \sum_{e} G_ie^Q b_ae^Q = \sum_{e} G_ei^Q b_ea^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|vo)", nQ_ref, nvirB, noccB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vv)", nQ_ref, nvirB, nvirB));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    GFvoB->contract(true, false, nvirB, noccB, nQ_ref * nvirB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Set global GF
    GFA->set_vo(GFvoA);
    GFB->set_vo(GFvoB);
    //GFvoA->print();
    //GFvoB->print();

}// else if (reference_ == "UNRESTRICTED")

    timer_off("GFM VO");
} // end gfock_vo

//======================================================================
//    CCSD: GFOCK
//======================================================================
void DFOCC::gfock_cc_vo()
{

    timer_on("GFM VO");
    SharedTensor2d K;
    SharedTensor2d G;

if (reference_ == "RESTRICTED") {
    //=========================
    // Reference Contribution
    //=========================
    GFvo->zero();
    GFvo->axpy(FvoA, 2.0);

    //=========================
    // Correlation Contribution
    //=========================

    // Fai = \sum_{m} h_am G_mi
    GFvo->gemm(false, false, HvoA, G1c_oo, 1.0, 1.0);

    // Fai = \sum_{m} h_ae G_ei
    GFvo->gemm(false, true, HvvA, G1c_ov, 1.0, 1.0);

    // Fai += \sum_{Q} \sum_{m} G_mi^Q b_ma^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OO)", nQ, noccA * noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA * nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFvo->contract(true, false, nvirA, noccA, nQ * noccA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fai += \sum_{Q} \sum_{e} G_ie^Q b_ae^Q = \sum_{e} G_ei^Q b_ea^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VO)", nQ, nvirA, noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    GFvo->contract(true, false, nvirA, noccA, nQ * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    //=========================
    // Separable Part
    //=========================

    // Fai += \sum_{Q} \sum_{m} G_mi^Q b_ma^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OO)", nQ_ref, noccA * noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA * nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFvo->contract(true, false, nvirA, noccA, nQ_ref * noccA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fai += \sum_{Q} \sum_{e} G_ie^Q b_ae^Q = \sum_{e} G_ei^Q b_ea^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VO)", nQ_ref, nvirA, noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    GFvo->contract(true, false, nvirA, noccA, nQ_ref * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Set global GF
    GF->set_vo(GFvo);

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    //=========================
    // Reference Contribution
    //=========================
    GFvoA->zero();
    GFvoB->zero();
    GFvoA->axpy(FvoA, 1.0);
    GFvoB->axpy(FvoB, 1.0);

    //=========================
    // Correlation Contribution
    //=========================

    // F_AI = \sum_{M} h_AM G_Mi
    GFvoA->gemm(false, false, HvoA, G1c_ooA, 1.0, 1.0);

    // F_AI = \sum_{E} h_AE G_EI
    GFvoA->gemm(false, false, HvvA, G1c_voA, 1.0, 1.0);

    // F_AI += \sum_{Q} \sum_{M} G_MI^Q b_MA^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OO)", nQ, noccA * noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA * nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFvoA->contract(true, false, nvirA, noccA, nQ * noccA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // F_AI += \sum_{Q} \sum_{E} G_IE^Q b_AE^Q = \sum_{E} G_EI^Q b_EA^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VO)", nQ, nvirA, noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    GFvoA->contract(true, false, nvirA, noccA, nQ * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fai = \sum_{m} h_am G_mi
    GFvoB->gemm(false, false, HvoB, G1c_ooB, 1.0, 1.0);

    // Fai = \sum_{e} h_ae G_ei
    GFvoB->gemm(false, false, HvvB, G1c_voB, 1.0, 1.0);

    // Fai += \sum_{Q} \sum_{m} G_mi^Q b_ma^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|oo)", nQ, noccB * noccB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ov)", nQ, noccB * nvirB));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFvoB->contract(true, false, nvirB, noccB, nQ * noccB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fai += \sum_{Q} \sum_{e} G_ie^Q b_ae^Q = \sum_{e} G_ei^Q b_ea^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|vo)", nQ, nvirB, noccB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|vv)", nQ, nvirB, nvirB));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    GFvoB->contract(true, false, nvirB, noccB, nQ * nvirB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    //=========================
    // Separable Part
    //=========================

    // F_AI += \sum_{Q} \sum_{M} G_MI^Q b_MA^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OO)", nQ_ref, noccA * noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA * nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFvoA->contract(true, false, nvirA, noccA, nQ_ref * noccA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // F_AI += \sum_{Q} \sum_{E} G_IE^Q b_AE^Q = \sum_{E} G_EI^Q b_EA^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VO)", nQ_ref, nvirA, noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    GFvoA->contract(true, false, nvirA, noccA, nQ_ref * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fai += \sum_{Q} \sum_{m} G_mi^Q b_ma^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|oo)", nQ_ref, noccB * noccB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB * nvirB));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFvoB->contract(true, false, nvirB, noccB, nQ_ref * noccB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fai += \sum_{Q} \sum_{e} G_ie^Q b_ae^Q = \sum_{e} G_ei^Q b_ea^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|vo)", nQ_ref, nvirB, noccB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vv)", nQ_ref, nvirB, nvirB));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    GFvoB->contract(true, false, nvirB, noccB, nQ_ref * nvirB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Set global GF
    GFA->set_vo(GFvoA);
    GFB->set_vo(GFvoB);
    //GFvoA->print();
    //GFvoB->print();

}// else if (reference_ == "UNRESTRICTED")

    timer_off("GFM VO");
} // end gfock_vo


}} // End Namespaces
