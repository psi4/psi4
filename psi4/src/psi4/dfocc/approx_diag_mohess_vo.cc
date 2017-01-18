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

void DFOCC::approx_diag_mohess_vo()
{

    timer_on("Diagonal MO Hessian VO");
    SharedTensor2d K;
    SharedTensor2d G;
    SharedTensor2d T;

if (reference_ == "RESTRICTED") {
    //=========================
    // Diagonal MO Hessian
    //=========================

    // A_ai = -2(Faa + Fii) + 2 h_aa G_ii + 2 h_ii G_aa
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              int ia = ov_idxAA->get(i,a);
              double value = -2.0 * (GF->get(a + noccA, a + noccA) + GF->get(i,i));
              value += 2.0 * HmoA->get(a + noccA, a + noccA) * G1->get(i,i);
              value += 2.0 * HmoA->get(i,i) * G1c_vv->get(a,a);
              AvoA->set(a, i, value);
         }
    }

    // A_ai += 2 G_ii \sum_{m} (ma|ma)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|OV)", noccA, nvirA, noccA, nvirA));
    tei_ovov_chem_ref_directAA(K);
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              int ia = ov_idxAA->get(i,a);
              double sum = 0.0;
              for (int m = 0; m < noccA; m++) {
                   int ma = ov_idxAA->get(m,a);
                   sum += K->get(ma,ma);
              }
              //double value = -2.0 * G1->get(i,i) * sum;
              double value = 2.0 * G1->get(i,i) * sum;
              AvoA->add(a, i, value);
         }
    }
    K.reset();

    //=========================
    // Reference Contribution
    //=========================

    // A_ai += 2\sum_{Q} b_aa^Q G_ii^Q
    G = SharedTensor2d(new Tensor2d("Reference 3-Index TPDM (Q|OO)", nQ_ref, noccA * noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    for (int a = 0; a < nvirA; a++) {
         int aa = vv_idxAA->get(a,a);
         for (int i = 0; i < noccA; i++) {
              int ii = oo_idxAA->get(i,i);
              double sum = 0.0;
              for (int Q = 0; Q < nQ_ref; Q++) {
                   sum += G->get(Q,ii) * K->get(Q,aa);
              }
              AvoA->add(a, i, 2.0*sum);
         }
    }
    G.reset();
    K.reset();

    //=========================
    // Separable Part
    //=========================

    // A_ai += 2\sum_{Q} b_aa^Q G_ii^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OO)", nQ_ref, noccA * noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    for (int a = 0; a < nvirA; a++) {
         int aa = vv_idxAA->get(a,a);
         for (int i = 0; i < noccA; i++) {
              int ii = oo_idxAA->get(i,i);
              double sum = 0.0;
              for (int Q = 0; Q < nQ_ref; Q++) {
                   sum += G->get(Q,ii) * K->get(Q,aa);
              }
              AvoA->add(a, i, 2.0*sum);
         }
    }
    G.reset();
    K.reset();

    // A_ai += 2\sum_{Q} b_ii^Q G_aa^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VV)", nQ_ref, nvirA, nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA * noccA));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    K->read(psio_, PSIF_DFOCC_INTS);
    for (int a = 0; a < nvirA; a++) {
         int aa = vv_idxAA->get(a,a);
         for (int i = 0; i < noccA; i++) {
              int ii = oo_idxAA->get(i,i);
              double sum = 0.0;
              for (int Q = 0; Q < nQ_ref; Q++) {
                   sum += G->get(Q,aa) * K->get(Q,ii);
              }
              AvoA->add(a, i, 2.0*sum);
         }
    }
    G.reset();
    K.reset();

    // The factor of 1/2 is necessary to convert spin-free Hessian to Alpha spin case
    AvoA->scale(0.5);
    if (print_ > 2) AvoA->print();

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    //=========================
    // Diagonal MO Hessian
    //=========================

    // A_AI = -2(F_AA + F_II) + 2 h_AA G_II + 2 h_II G_AA
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              int ia = ov_idxAA->get(i,a);
              double value = -2.0 * (GFA->get(a + noccA, a + noccA) + GFA->get(i,i));
              value += 2.0 * HmoA->get(a + noccA, a + noccA) * G1A->get(i,i);
              value += 2.0 * HmoA->get(i,i) * G1c_vvA->get(a,a);
              AvoA->set(a, i, value);
         }
    }

    // A_ai = -2(Faa + Fii) + 2 h_aa G_ii + 2 h_ii G_aa
    for (int a = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++) {
              int ia = ov_idxBB->get(i,a);
              double value = -2.0 * (GFB->get(a + noccB, a + noccB) + GFB->get(i,i));
              value += 2.0 * HmoB->get(a + noccB, a + noccB) * G1B->get(i,i);
              value += 2.0 * HmoB->get(i,i) * G1c_vvB->get(a,a);
              AvoB->set(a, i, value);
         }
    }

    // A_AI += 2 G_II \sum_{M} (MA|MA)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|OV)", noccA, nvirA, noccA, nvirA));
    tei_ovov_chem_ref_directAA(K);
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              int ia = ov_idxAA->get(i,a);
              double sum = 0.0;
              for (int m = 0; m < noccA; m++) {
                   int ma = ov_idxAA->get(m,a);
                   sum += K->get(ma,ma);
              }
              //double value = -2.0 * G1A->get(i,i) * sum;
              double value = 2.0 * G1A->get(i,i) * sum;
              AvoA->add(a, i, value);
         }
    }
    K.reset();

    // A_ai += -2 G_ii \sum_{m} (ma|ma)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (ov|ov)", noccB, nvirB, noccB, nvirB));
    tei_ovov_chem_ref_directBB(K);
    for (int a = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++) {
              int ia = ov_idxBB->get(i,a);
              double sum = 0.0;
              for (int m = 0; m < noccB; m++) {
                   int ma = ov_idxBB->get(m,a);
                   sum += K->get(ma,ma);
              }
              //double value = -2.0 * G1B->get(i,i) * sum;
              double value = 2.0 * G1B->get(i,i) * sum;
              AvoB->add(a, i, value);
         }
    }
    K.reset();

    //=========================
    // Reference Contribution
    //=========================

    // A_AI += 2\sum_{Q} b_AA^Q G_II^Q
    G = SharedTensor2d(new Tensor2d("Reference 3-Index TPDM (Q|OO)", nQ_ref, noccA * noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    for (int a = 0; a < nvirA; a++) {
         int aa = vv_idxAA->get(a,a);
         for (int i = 0; i < noccA; i++) {
              int ii = oo_idxAA->get(i,i);
              double sum = 0.0;
              for (int Q = 0; Q < nQ_ref; Q++) {
                   sum += G->get(Q,ii) * K->get(Q,aa);
              }
              AvoA->add(a, i, 2.0*sum);
         }
    }
    G.reset();
    K.reset();

    // A_ai += 2\sum_{Q} b_aa^Q G_ii^Q
    G = SharedTensor2d(new Tensor2d("Reference 3-Index TPDM (Q|oo)", nQ_ref, noccB * noccB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vv)", nQ_ref, nvirB, nvirB));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    for (int a = 0; a < nvirB; a++) {
         int aa = vv_idxBB->get(a,a);
         for (int i = 0; i < noccB; i++) {
              int ii = oo_idxBB->get(i,i);
              double sum = 0.0;
              for (int Q = 0; Q < nQ_ref; Q++) {
                   sum += G->get(Q,ii) * K->get(Q,aa);
              }
              AvoB->add(a, i, 2.0*sum);
         }
    }
    G.reset();
    K.reset();

    //=========================
    // Separable Part
    //=========================

    // A_AI += 2\sum_{Q} b_AA^Q G_II^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OO)", nQ_ref, noccA * noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    for (int a = 0; a < nvirA; a++) {
         int aa = vv_idxAA->get(a,a);
         for (int i = 0; i < noccA; i++) {
              int ii = oo_idxAA->get(i,i);
              double sum = 0.0;
              for (int Q = 0; Q < nQ_ref; Q++) {
                   sum += G->get(Q,ii) * K->get(Q,aa);
              }
              AvoA->add(a, i, 2.0*sum);
         }
    }
    G.reset();
    K.reset();

    // A_ai += 2\sum_{Q} b_aa^Q G_ii^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|oo)", nQ_ref, noccB * noccB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vv)", nQ_ref, nvirB, nvirB));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    for (int a = 0; a < nvirB; a++) {
         int aa = vv_idxBB->get(a,a);
         for (int i = 0; i < noccB; i++) {
              int ii = oo_idxBB->get(i,i);
              double sum = 0.0;
              for (int Q = 0; Q < nQ_ref; Q++) {
                   sum += G->get(Q,ii) * K->get(Q,aa);
              }
              AvoB->add(a, i, 2.0*sum);
         }
    }
    G.reset();
    K.reset();

    // A_AI += 2\sum_{Q} b_II^Q G_AA^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VV)", nQ_ref, nvirA, nvirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA * noccA));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    K->read(psio_, PSIF_DFOCC_INTS);
    for (int a = 0; a < nvirA; a++) {
         int aa = vv_idxAA->get(a,a);
         for (int i = 0; i < noccA; i++) {
              int ii = oo_idxAA->get(i,i);
              double sum = 0.0;
              for (int Q = 0; Q < nQ_ref; Q++) {
                   sum += G->get(Q,aa) * K->get(Q,ii);
              }
              AvoA->add(a, i, 2.0*sum);
         }
    }
    G.reset();
    K.reset();

    // A_ai += 2\sum_{Q} b_ii^Q G_aa^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|vv)", nQ_ref, nvirB, nvirB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB * noccB));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    K->read(psio_, PSIF_DFOCC_INTS);
    for (int a = 0; a < nvirB; a++) {
         int aa = vv_idxBB->get(a,a);
         for (int i = 0; i < noccB; i++) {
              int ii = oo_idxBB->get(i,i);
              double sum = 0.0;
              for (int Q = 0; Q < nQ_ref; Q++) {
                   sum += G->get(Q,aa) * K->get(Q,ii);
              }
              AvoB->add(a, i, 2.0*sum);
         }
    }
    G.reset();
    K.reset();

    if (print_ > 2) {
        AvoA->print();
        AvoB->print();
    }

}// else if (reference_ == "UNRESTRICTED")
    timer_off("Diagonal MO Hessian VO");
} // end diagonal_mohess_vo


}} // End Namespaces
