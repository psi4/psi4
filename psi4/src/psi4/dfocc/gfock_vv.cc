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

void DFOCC::gfock_vv()
{

    timer_on("GFM VV");
    SharedTensor2d K;
    SharedTensor2d G;

if (reference_ == "RESTRICTED") {
    //=========================
    // Correlation Contribution
    //=========================

    // Fab = \sum_{e} h_ae G_eb
    GFvv->gemm(true, false, HvvA, G1c_vv, 1.0, 0.0);

    // Fab += \sum_{Q} \sum_{m} G_mb^Q b_ma^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA, nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA * nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFvv->contract(true, false, nvirA, nvirA, nQ * noccA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    //=========================
    // Separable Part
    //=========================

    // Fab += \sum_{Q} \sum_{m} G_mb^Q b_ma^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OV)", nQ_ref, noccA * nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA * nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFvv->contract(true, false, nvirA, nvirA, nQ_ref * noccA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fab += \sum_{Q} \sum_{e} G_eb^Q b_ea^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VV)", nQ_ref, nvirA, nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    GFvv->contract(true, false, nvirA, nvirA, nQ_ref * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Set global GF
    GF->set_vv(noccA, GFvv);
    if (print_ > 2) GF->print();

    /*
    // Energy
    double Etemp = 0.0;
    G = SharedTensor2d(new Tensor2d("MO-basis correlation OPDM", nmo_, nmo_));
    G->copy(G1c);
    for (int i = 0; i < noccA; i++) G->add(i, i, 2.0);
    Etemp += 0.5 * G->vector_dot(HmoA);
    G.reset();
    for (int p = 0; p < nmo_; p++) {
         Etemp += 0.5*GF->get(p,p);
    }
    Etemp += Enuc;
    outfile->Printf("\tDF-MP2L Total Energy via GFM (a.u.): %20.14f\n", Etemp);
    */

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    //=========================
    // Correlation Contribution
    //=========================

    // Fab = \sum_{e} h_ae G_eb
    GFvvA->gemm(true, false, HvvA, G1c_vvA, 1.0, 0.0);
    GFvvB->gemm(true, false, HvvB, G1c_vvB, 1.0, 0.0);

 if (reference == "ROHF" && orb_opt_ == "FALSE") {
    // Fab = \sum_{m} h_am G_mb
    GFvvA->gemm(false, false, HvoA, G1c_ovA, 1.0, 1.0);
    GFvvB->gemm(false, false, HvoB, G1c_ovB, 1.0, 1.0);
 }

    // Fab += \sum_{Q} \sum_{m} G_mb^Q b_ma^Q
    // alpha spin
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA, nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA * nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFvvA->contract(true, false, nvirA, nvirA, nQ * noccA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // beta spin
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|ov)", nQ, noccB, nvirB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ov)", nQ, noccB * nvirB));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFvvB->contract(true, false, nvirB, nvirB, nQ * noccB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    //=========================
    // Separable Part
    //=========================

    // Fab += \sum_{Q} \sum_{m} G_mb^Q b_ma^Q
    // alpha spin
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OV)", nQ_ref, noccA * nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA * nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFvvA->contract(true, false, nvirA, nvirA, nQ_ref * noccA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // beta spin
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|ov)", nQ_ref, noccB * nvirB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB * nvirB));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFvvB->contract(true, false, nvirB, nvirB, nQ_ref * noccB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fab += \sum_{Q} \sum_{e} G_eb^Q b_ea^Q
    // alpha spin
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VV)", nQ_ref, nvirA, nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    GFvvA->contract(true, false, nvirA, nvirA, nQ_ref * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // beta spin
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|vv)", nQ_ref, nvirB, nvirB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vv)", nQ_ref, nvirB, nvirB));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    GFvvB->contract(true, false, nvirB, nvirB, nQ_ref * nvirB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Set global GF
    GFA->set_vv(noccA, GFvvA);
    GFB->set_vv(noccB, GFvvB);

    if (print_ > 2){
        GFA->print();
        GFB->print();
    }

    /*
    // Energy
    double Etemp = 0.0;
    G = SharedTensor2d(new Tensor2d("MO-basis correlation OPDM", nmo_, nmo_));
    G->copy(G1cA);
    for (int i = 0; i < noccA; i++) G->add(i, i, 1.0);
    Etemp += 0.5 * G->vector_dot(HmoA);
    G.reset();

    G = SharedTensor2d(new Tensor2d("MO-basis correlation OPDM", nmo_, nmo_));
    G->copy(G1cB);
    for (int i = 0; i < noccB; i++) G->add(i, i, 1.0);
    Etemp += 0.5 * G->vector_dot(HmoB);
    G.reset();

    for (int p = 0; p < nmo_; p++) {
         Etemp += 0.5 * GFA->get(p,p);
         Etemp += 0.5 * GFB->get(p,p);
    }
    Etemp += Enuc;
    outfile->Printf("\tDF-MP2L Total Energy via GFM (a.u.): %20.14f\n", Etemp);

    */

}// else if (reference_ == "UNRESTRICTED")
    timer_off("GFM VV");

    // Call Complementary GFM
    if (hess_type == "APPROX_DIAG_EKT") gftilde_vv();
} // end gfock_vv


//======================================================================
//    CCSD: GFOCK
//======================================================================
void DFOCC::gfock_cc_vv()
{

    timer_on("GFM VV");
    SharedTensor2d K;
    SharedTensor2d G;

if (reference_ == "RESTRICTED") {
    //=========================
    // Correlation Contribution
    //=========================

    // Fab = \sum_{e} h_ae G_eb
    GFvv->gemm(false, false, HvvA, G1c_vv, 1.0, 0.0);

    // Fab += \sum_{m} h_am G_mb
    GFvv->gemm(false, false, HvoA, G1c_ov, 1.0, 1.0);

    // Fab += \sum_{Q} \sum_{m} G_mb^Q b_ma^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA, nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFvv->contract(true, false, nvirA, nvirA, nQ * noccA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fab += \sum_{Q} \sum_{e} G_eb^Q b_ea^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VV)", nQ, nvirA, nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    GFvv->contract(true, false, nvirA, nvirA, nQ * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    //=========================
    // Separable Part
    //=========================

    // Fab += \sum_{Q} \sum_{m} G_mb^Q b_ma^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OV)", nQ_ref, noccA * nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA * nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFvv->contract(true, false, nvirA, nvirA, nQ_ref * noccA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fab += \sum_{Q} \sum_{e} G_eb^Q b_ea^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VV)", nQ_ref, nvirA, nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    GFvv->contract(true, false, nvirA, nvirA, nQ_ref * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Set global GF
    GF->set_vv(noccA, GFvv);
    if (print_ > 2) GF->print();

    /*
    // Energy
    double Etemp = 0.0;
    Etemp += 0.5 * G1->vector_dot(HmoA);
    for (int p = 0; p < nmo_; p++) {
         Etemp += 0.5*GF->get(p,p);
    }
    Etemp += Enuc;
    outfile->Printf("\tDF-CCL Total Energy via GFM (a.u.) : %20.14f\n", Etemp);
    */


}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    //=========================
    // Correlation Contribution
    //=========================

    // Fab = \sum_{e} h_ae G_eb
    GFvvA->gemm(true, false, HvvA, G1c_vvA, 1.0, 0.0);
    GFvvB->gemm(true, false, HvvB, G1c_vvB, 1.0, 0.0);

    // Fab = \sum_{m} h_am G_mb
    GFvvA->gemm(false, false, HvoA, G1c_ovA, 1.0, 1.0);
    GFvvB->gemm(false, false, HvoB, G1c_ovB, 1.0, 1.0);

    // Fab += \sum_{Q} \sum_{m} G_mb^Q b_ma^Q
    // alpha spin
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA, nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA * nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFvvA->contract(true, false, nvirA, nvirA, nQ * noccA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // beta spin
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|ov)", nQ, noccB, nvirB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ov)", nQ, noccB * nvirB));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFvvB->contract(true, false, nvirB, nvirB, nQ * noccB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fab += \sum_{Q} \sum_{e} G_eb^Q b_ea^Q
    // alpha spin
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VV)", nQ, nvirA, nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    GFvvA->contract(true, false, nvirA, nvirA, nQ * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // beta spin
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|vv)", nQ, nvirB, nvirB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|vv)", nQ, nvirB, nvirB));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    GFvvB->contract(true, false, nvirB, nvirB, nQ * nvirB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    //=========================
    // Separable Part
    //=========================

    // Fab += \sum_{Q} \sum_{m} G_mb^Q b_ma^Q
    // alpha spin
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OV)", nQ_ref, noccA * nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA * nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFvvA->contract(true, false, nvirA, nvirA, nQ_ref * noccA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // beta spin
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|ov)", nQ_ref, noccB * nvirB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB * nvirB));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    GFvvB->contract(true, false, nvirB, nvirB, nQ_ref * noccB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Fab += \sum_{Q} \sum_{e} G_eb^Q b_ea^Q
    // alpha spin
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VV)", nQ_ref, nvirA, nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    GFvvA->contract(true, false, nvirA, nvirA, nQ_ref * nvirA, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // beta spin
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|vv)", nQ_ref, nvirB, nvirB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vv)", nQ_ref, nvirB, nvirB));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    GFvvB->contract(true, false, nvirB, nvirB, nQ_ref * nvirB, K, G, 1.0, 1.0);
    G.reset();
    K.reset();

    // Set global GF
    GFA->set_vv(noccA, GFvvA);
    GFB->set_vv(noccB, GFvvB);

    if (print_ > 2){
        GFA->print();
        GFB->print();
    }

}// else if (reference_ == "UNRESTRICTED")
    timer_off("GFM VV");
} // end gfock_cc_vv



}} // End Namespaces
