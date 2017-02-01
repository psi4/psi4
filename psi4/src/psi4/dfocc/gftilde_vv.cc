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

void DFOCC::gftilde_vv()
{

    timer_on("GFtilde VV");
    SharedTensor2d K;
    SharedTensor2d G;

if (reference_ == "RESTRICTED") {
    // G_Q = \sum_{m,n} b_mn^Q G1c_mn
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA * noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    g1Q->gemv(false, nQ, noccA * noccA, K, G1c_oo, 1.0, 0.0);
    K.reset();

    // Gt_Q = \sum_{e,f} b_ef^Q G1c_ef
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    g1Qt2->gemv(false, nQ, nvirA * nvirA, K, G1c_vv, 1.0, 0.0);
    K.reset();

    // Ft_ab = 2 h_ab - F_ab
    GFtvv->copy(HvvA);
    GFtvv->scale(2.0);
    GFtvv->subtract(GFvv);

    // Ft_ab += 2.0*\sum_{Q} b_ab^Q J_Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    #pragma omp parallel for
    for (int a = 0; a < nvirA; a++) {
         for (int b = 0; b < nvirA; b++) {
              int ab = vv_idxAA->get(a,b);
              double sum = 0.0;
              for (int Q = 0; Q < nQ_ref; Q++) {
                   sum += K->get(Q,ab) * Jc->get(Q);
              }
              GFtvv->add(a, b, 2.0*sum);
         }
    }
    K.reset();

    // Ft_ab = 2.0 \sum_{Q} b_ab^Q (G_Q + Gt_Q)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    #pragma omp parallel for
    for (int a = 0; a < nvirA; a++) {
         for (int b = 0; b < nvirA; b++) {
              int ab = vv_idxAA->get(a,b);
              double sum = 0.0;
              for (int Q = 0; Q < nQ; Q++) {
                   sum += K->get(Q,ab) * ( g1Q->get(Q) + g1Qt2->get(Q) );
              }
              GFtvv->add(a, b, 2.0*sum);
         }
    }
    K.reset();

    // Ft_ab -= 2 \sum_{m} (ma|mb)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|OV)", noccA, nvirA, noccA, nvirA));
    tei_ovov_chem_ref_directAA(K);
    #pragma omp parallel for
    for (int a = 0; a < nvirA; a++) {
         for (int b = 0; b < nvirA; b++) {
              double sum = 0.0;
              for (int m = 0; m < noccA; m++) {
                   int ma = ov_idxAA->get(m,a);
                   int mb = ov_idxAA->get(m,b);
                   sum += K->get(ma,mb);
              }
              GFtvv->add(a, b, -2.0*sum);
         }
    }
    K.reset();

    // Ft_ab -= \sum_{Q} \sum_{m} b_ma^Q G_mb^Q
    G = SharedTensor2d(new Tensor2d("3-Index Corr OPDM (Q|OV)", nQ, noccA, nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA * nvirA));
    K->read(psio_, PSIF_DFOCC_INTS);
    // G_ia^Q = \sum_{m} G_im b_ma^Q
    G->contract233(false, false, noccA, nvirA, G1c_oo, K, 1.0, 0.0);
    GFtvv->contract(true, false, nvirA, nvirA, nQ * noccA, K, G, -1.0, 1.0);
    K.reset();
    G.reset();

    // Ft_ab -= \sum_{Q} \sum_{e} G_ea^Q b_eb^Q
    G = SharedTensor2d(new Tensor2d("3-Index Corr OPDM (Q|VV)", nQ, nvirA, nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    // G_ab^Q = \sum_{e} G_ae b_eb^Q
    G->contract233(false, false, nvirA, nvirA, G1c_vv, K, 1.0, 0.0);
    GFtvv->contract(true, false, nvirA, nvirA, nQ * nvirA, G, K, -1.0, 1.0);
    G.reset();
    K.reset();

    // Set global GF
    //GFt->set_vv(noccA, GFtvv);
    if (print_ > 2) GFtvv->print();

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // G_Q = \sum_{m,n} b_mn^Q G1c_mn
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA * noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    g1Q->gemv(false, nQ, noccA * noccA, K, G1c_ooA, 1.0, 0.0);
    K.reset();

    // beta
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|oo)", nQ, noccB * noccB));
    K->read(psio_, PSIF_DFOCC_INTS);
    g1Q->gemv(false, nQ, noccB * noccB, K, G1c_ooB, 1.0, 1.0);
    K.reset();

    // Gt_Q = \sum_{e,f} b_ef^Q G1c_ef
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    g1Qt2->gemv(false, nQ, nvirA * nvirA, K, G1c_vvA, 1.0, 0.0);
    K.reset();

    // beta
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|vv)", nQ, nvirB, nvirB));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    g1Qt2->gemv(false, nQ, nvirB * nvirB, K, G1c_vvB, 1.0, 1.0);
    K.reset();

    // Ft_ab =  h_ab - F_ab
    GFtvvA->copy(HvvA);
    GFtvvA->subtract(GFvvA);
    GFtvvB->copy(HvvB);
    GFtvvB->subtract(GFvvB);

    // Ft_AB += \sum_{Q} b_AB^Q J_Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    #pragma omp parallel for
    for (int a = 0; a < nvirA; a++) {
         for (int b = 0; b < nvirA; b++) {
              int ab = vv_idxAA->get(a,b);
              double sum = 0.0;
              for (int Q = 0; Q < nQ_ref; Q++) {
                   sum += K->get(Q,ab) * Jc->get(Q);
              }
              GFtvvA->add(a, b, sum);
         }
    }
    K.reset();

    // Ft_ab += \sum_{Q} b_ab^Q J_Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vv)", nQ_ref, nvirB, nvirB));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    #pragma omp parallel for
    for (int a = 0; a < nvirB; a++) {
         for (int b = 0; b < nvirB; b++) {
              int ab = vv_idxBB->get(a,b);
              double sum = 0.0;
              for (int Q = 0; Q < nQ_ref; Q++) {
                   sum += K->get(Q,ab) * Jc->get(Q);
              }
              GFtvvB->add(a, b, sum);
         }
    }
    K.reset();

    // Ft_AB =  \sum_{Q} b_AB^Q (G_Q + Gt_Q)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    #pragma omp parallel for
    for (int a = 0; a < nvirA; a++) {
         for (int b = 0; b < nvirA; b++) {
              int ab = vv_idxAA->get(a,b);
              double sum = 0.0;
              for (int Q = 0; Q < nQ; Q++) {
                   sum += K->get(Q,ab) * ( g1Q->get(Q) + g1Qt2->get(Q) );
              }
              GFtvvA->add(a, b, sum);
         }
    }
    K.reset();

    // Ft_ab =  \sum_{Q} b_ab^Q (G_Q + Gt_Q)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|vv)", nQ, nvirB, nvirB));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    #pragma omp parallel for
    for (int a = 0; a < nvirB; a++) {
         for (int b = 0; b < nvirB; b++) {
              int ab = vv_idxBB->get(a,b);
              double sum = 0.0;
              for (int Q = 0; Q < nQ; Q++) {
                   sum += K->get(Q,ab) * ( g1Q->get(Q) + g1Qt2->get(Q) );
              }
              GFtvvB->add(a, b, sum);
         }
    }
    K.reset();

    // Ft_AB -= \sum_{M} (MA|MB)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|OV)", noccA, nvirA, noccA, nvirA));
    tei_ovov_chem_ref_directAA(K);
    #pragma omp parallel for
    for (int a = 0; a < nvirA; a++) {
         for (int b = 0; b < nvirA; b++) {
              double sum = 0.0;
              for (int m = 0; m < noccA; m++) {
                   int ma = ov_idxAA->get(m,a);
                   int mb = ov_idxAA->get(m,b);
                   sum -= K->get(ma,mb);
              }
              GFtvvA->add(a, b, sum);
         }
    }
    K.reset();

    // Ft_ab -= \sum_{m} (ma|mb)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (ov|ov)", noccB, nvirB, noccB, nvirB));
    tei_ovov_chem_ref_directBB(K);
    #pragma omp parallel for
    for (int a = 0; a < nvirB; a++) {
         for (int b = 0; b < nvirB; b++) {
              double sum = 0.0;
              for (int m = 0; m < noccB; m++) {
                   int ma = ov_idxBB->get(m,a);
                   int mb = ov_idxBB->get(m,b);
                   sum -= K->get(ma,mb);
              }
              GFtvvB->add(a, b, sum);
         }
    }
    K.reset();

    // Ft_AB -= \sum_{Q} \sum_{M} b_MA^Q G_MB^Q
    G = SharedTensor2d(new Tensor2d("3-Index Corr OPDM (Q|OV)", nQ, noccA, nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA * nvirA));
    K->read(psio_, PSIF_DFOCC_INTS);
    // G_IA^Q = \sum_{M} G_IM b_MA^Q
    G->contract233(false, false, noccA, nvirA, G1c_ooA, K, 1.0, 0.0);
    GFtvvA->contract(true, false, nvirA, nvirA, nQ * noccA, K, G, -1.0, 1.0);
    K.reset();
    G.reset();

    // Ft_ab -= \sum_{Q} \sum_{m} b_ma^Q G_mb^Q
    G = SharedTensor2d(new Tensor2d("3-Index Corr OPDM (Q|ov)", nQ, noccB, nvirB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ov)", nQ, noccB * nvirB));
    K->read(psio_, PSIF_DFOCC_INTS);
    // G_ia^Q = \sum_{m} G_im b_ma^Q
    G->contract233(false, false, noccB, nvirB, G1c_ooB, K, 1.0, 0.0);
    GFtvvB->contract(true, false, nvirB, nvirB, nQ * noccB, K, G, -1.0, 1.0);
    K.reset();
    G.reset();

    // Ft_AB -= \sum_{Q} \sum_{E} G_EA^Q b_EB^Q
    G = SharedTensor2d(new Tensor2d("3-Index Corr OPDM (Q|VV)", nQ, nvirA, nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    // G_AB^Q = \sum_{E} G_AE b_EB^Q
    G->contract233(false, false, nvirA, nvirA, G1c_vvA, K, 1.0, 0.0);
    GFtvvA->contract(true, false, nvirA, nvirA, nQ * nvirA, G, K, -1.0, 1.0);
    G.reset();
    K.reset();

    // Ft_ab -= \sum_{Q} \sum_{e} G_ea^Q b_eb^Q
    G = SharedTensor2d(new Tensor2d("3-Index Corr OPDM (Q|vv)", nQ, nvirB, nvirB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|vv)", nQ, nvirB, nvirB));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    // G_ab^Q = \sum_{e} G_ae b_eb^Q
    G->contract233(false, false, nvirB, nvirB, G1c_vvB, K, 1.0, 0.0);
    GFtvvB->contract(true, false, nvirB, nvirB, nQ * nvirB, G, K, -1.0, 1.0);
    G.reset();
    K.reset();

    // Set global GF
    //GFtA->set_vv(noccA, GFtvvA);
    //GFtB->set_vv(noccB, GFtvvB);

    if (print_ > 2){
        GFtvvA->print();
        GFtvvB->print();
    }

}// else if (reference_ == "UNRESTRICTED")
    timer_off("GFtilde VV");
} // end gftilde_vv


}} // End Namespaces
