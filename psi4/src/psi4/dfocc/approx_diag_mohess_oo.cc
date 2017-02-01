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

void DFOCC::approx_diag_mohess_oo()
{

    timer_on("Diagonal MO Hessian OO");
    SharedTensor2d K;
    SharedTensor2d G;

if (reference_ == "RESTRICTED") {
    // A_ij = -2 (F_ii + F_jj) + 2 h_ii G_jj + 2 h_jj G_ii - 4 h_ij G_ij
    for (int i = 0; i < naoccA; i++) {
         for (int j = 0; j < nfrzc; j++) {
              double value = -2.0 * (GF->get(i + nfrzc, i + nfrzc) + GF->get(j,j));
              value += 2.0 * HmoA->get(i + nfrzc, i + nfrzc) * G1->get(j,j);
              value += 2.0 * HmoA->get(j, j) * G1->get(i + nfrzc, i + nfrzc);
              value -= 4.0 * HmoA->get(i + nfrzc, j) * G1->get(i + nfrzc, j);
              AooA->set(i, j, value);
         }
    }

    // A_ij += 4 \sum_{m} [(mi|mi) - (mj|mj)]
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|OO)", noccA, noccA, noccA, noccA));
    tei_oooo_chem_ref_directAA(K);
    for (int i = 0; i < naoccA; i++) {
         for (int j = 0; j < nfrzc; j++) {
              int ii = oo_idxAA->get(i + nfrzc, i + nfrzc);
              int jj = oo_idxAA->get(j,j);
              int ij = oo_idxAA->get(i + nfrzc,j);
              double value = 0.0;
              for (int m = 0; m < noccA; m++) {
                   int mi = oo_idxAA->get(m, i + nfrzc);
                   int mj = oo_idxAA->get(m, j);
                   value += K->get(mi,mi) - K->get(mj,mj);
              }
              AooA->add(i, j, 4.0*value);
         }
    }
    K.reset();

    //=========================
    // Reference Contribution
    //=========================

    // A_ij += \sum_{Q} 2 b_ii^Q G_jj^Q + 2 b_jj^Q G_ii^Q - 4 b_ij^Q G_ij^Q
    G = SharedTensor2d(new Tensor2d("Reference 3-Index TPDM (Q|OO)", nQ_ref, noccA * noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    for (int i = 0; i < naoccA; i++) {
         int ii = oo_idxAA->get(i + nfrzc,i + nfrzc);
         for (int j = 0; j < nfrzc; j++) {
              int jj = oo_idxAA->get(j,j);
              int ij = oo_idxAA->get(i + nfrzc,j);
              double sum = 0.0;
              for (int Q = 0; Q < nQ_ref; Q++) {
                   sum += 2.0 * G->get(Q,jj) * K->get(Q,ii);
                   sum += 2.0 * G->get(Q,ii) * K->get(Q,jj);
                   sum -= 4.0 * G->get(Q,ij) * K->get(Q,ij);
              }
              AooA->add(i, j, sum);
         }
    }
    G.reset();
    K.reset();

    //=========================
    // Separable Part
    //=========================

    // A_ij += \sum_{Q} 2 b_ii^Q G_jj^Q + 2 b_jj^Q G_ii^Q - 4 b_ij^Q G_ij^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OO)", nQ_ref, noccA * noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    for (int i = 0; i < naoccA; i++) {
         int ii = oo_idxAA->get(i + nfrzc,i + nfrzc);
         for (int j = 0; j < nfrzc; j++) {
              int jj = oo_idxAA->get(j,j);
              int ij = oo_idxAA->get(i + nfrzc,j);
              double sum = 0.0;
              for (int Q = 0; Q < nQ_ref; Q++) {
                   sum += 2.0 * G->get(Q,jj) * K->get(Q,ii);
                   sum += 2.0 * G->get(Q,ii) * K->get(Q,jj);
                   sum -= 4.0 * G->get(Q,ij) * K->get(Q,ij);
              }
              AooA->add(i, j, sum);
         }
    }
    G.reset();
    K.reset();

    AooA->scale(0.5);
    if (print_ > 2) AooA->print();

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // A_IJ = -2 (F_II + F_JJ) + 2 h_II G_JJ + 2 h_JJ G_II - 4 h_IJ G_IJ
    for (int i = 0; i < naoccA; i++) {
         for (int j = 0; j < nfrzc; j++) {
              double value = -2.0 * (GFA->get(i + nfrzc, i + nfrzc) + GFA->get(j,j));
              value += 2.0 * HmoA->get(i + nfrzc, i + nfrzc) * G1A->get(j,j);
              value += 2.0 * HmoA->get(j, j) * G1A->get(i + nfrzc, i + nfrzc);
              value -= 4.0 * HmoA->get(i + nfrzc, j) * G1A->get(i + nfrzc, j);
              AooA->set(i, j, value);
         }
    }

    // A_ij = -2 (F_ii + F_jj) + 2 h_ii G_jj + 2 h_jj G_ii - 4 h_ij G_ij
    for (int i = 0; i < naoccB; i++) {
         for (int j = 0; j < nfrzc; j++) {
              double value = -2.0 * (GFB->get(i + nfrzc, i + nfrzc) + GFB->get(j,j));
              value += 2.0 * HmoB->get(i + nfrzc, i + nfrzc) * G1B->get(j,j);
              value += 2.0 * HmoB->get(j, j) * G1B->get(i + nfrzc, i + nfrzc);
              value -= 4.0 * HmoB->get(i + nfrzc, j) * G1B->get(i + nfrzc, j);
              AooB->set(i, j, value);
         }
    }

    // A_IJ += 2 \sum_{M} [(MI|MI) - (MJ|MJ)]
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|OO)", noccA, noccA, noccA, noccA));
    tei_oooo_chem_ref_directAA(K);
    for (int i = 0; i < naoccA; i++) {
         for (int j = 0; j < nfrzc; j++) {
              int ii = oo_idxAA->get(i + nfrzc, i + nfrzc);
              int jj = oo_idxAA->get(j,j);
              int ij = oo_idxAA->get(i + nfrzc,j);
              double value = 0.0;
              for (int m = 0; m < noccA; m++) {
                   int mi = oo_idxAA->get(m, i + nfrzc);
                   int mj = oo_idxAA->get(m, j);
                   value += K->get(mi,mi) - K->get(mj,mj);
              }
              AooA->add(i, j, 2.0*value);
         }
    }
    K.reset();

    // A_ij += 2 \sum_{m} [(mi|mi) - (mj|mj)]
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (oo|oo)", noccB, noccB, noccB, noccB));
    tei_oooo_chem_ref_directBB(K);
    for (int i = 0; i < naoccB; i++) {
         for (int j = 0; j < nfrzc; j++) {
              int ii = oo_idxBB->get(i + nfrzc, i + nfrzc);
              int jj = oo_idxBB->get(j,j);
              int ij = oo_idxBB->get(i + nfrzc,j);
              double value = 0.0;
              for (int m = 0; m < noccB; m++) {
                   int mi = oo_idxBB->get(m, i + nfrzc);
                   int mj = oo_idxBB->get(m, j);
                   value += K->get(mi,mi) - K->get(mj,mj);
              }
              AooB->add(i, j, 2.0*value);
         }
    }
    K.reset();

    //=========================
    // Reference Contribution
    //=========================

    // A_IJ += \sum_{Q} 2 G_JJ^Q b_II^Q + 2 G_II^Q b_JJ^Q - 4 G_IJ^Q b_IJ^Q
    G = SharedTensor2d(new Tensor2d("Reference 3-Index TPDM (Q|OO)", nQ_ref, noccA * noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    for (int i = 0; i < naoccA; i++) {
         int ii = oo_idxAA->get(i + nfrzc,i + nfrzc);
         for (int j = 0; j < nfrzc; j++) {
              int jj = oo_idxAA->get(j,j);
              int ij = oo_idxAA->get(i + nfrzc,j);
              double sum = 0.0;
              for (int Q = 0; Q < nQ_ref; Q++) {
                   sum += 2.0 * G->get(Q,jj) * K->get(Q,ii);
                   sum += 2.0 * G->get(Q,ii) * K->get(Q,jj);
                   sum -= 4.0 * G->get(Q,ij) * K->get(Q,ij);
              }
              AooA->add(i, j, sum);
         }
    }
    G.reset();
    K.reset();

    // A_ij += \sum_{Q} 2 G_jj^Q b_ii^Q + 2 G_ii^Q b_jj^Q - 4 G_ij^Q b_ij^Q
    G = SharedTensor2d(new Tensor2d("Reference 3-Index TPDM (Q|oo)", nQ_ref, noccB * noccB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB, noccB));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    for (int i = 0; i < naoccB; i++) {
         int ii = oo_idxBB->get(i + nfrzc,i + nfrzc);
         for (int j = 0; j < nfrzc; j++) {
              int jj = oo_idxBB->get(j,j);
              int ij = oo_idxBB->get(i + nfrzc,j);
              double sum = 0.0;
              for (int Q = 0; Q < nQ_ref; Q++) {
                   sum += 2.0 * G->get(Q,jj) * K->get(Q,ii);
                   sum += 2.0 * G->get(Q,ii) * K->get(Q,jj);
                   sum -= 4.0 * G->get(Q,ij) * K->get(Q,ij);
              }
              AooB->add(i, j, sum);
         }
    }
    G.reset();
    K.reset();

    //=========================
    // Separable Part
    //=========================

    // A_IJ += \sum_{Q} 2 G_JJ^Q b_II^Q + 2 G_II^Q b_JJ^Q - 4 G_IJ^Q b_IJ^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OO)", nQ_ref, noccA * noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    for (int i = 0; i < naoccA; i++) {
         int ii = oo_idxAA->get(i + nfrzc,i + nfrzc);
         for (int j = 0; j < nfrzc; j++) {
              int jj = oo_idxAA->get(j,j);
              int ij = oo_idxAA->get(i + nfrzc,j);
              double sum = 0.0;
              for (int Q = 0; Q < nQ_ref; Q++) {
                   sum += 2.0 * G->get(Q,jj) * K->get(Q,ii);
                   sum += 2.0 * G->get(Q,ii) * K->get(Q,jj);
                   sum -= 4.0 * G->get(Q,ij) * K->get(Q,ij);
              }
              AooA->add(i, j, sum);
         }
    }
    G.reset();
    K.reset();

    // A_ij += \sum_{Q} 2 G_jj^Q b_ii^Q + 2 G_ii^Q b_jj^Q - 4 G_ij^Q b_ij^Q
    G = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|oo)", nQ_ref, noccB * noccB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB, noccB));
    G->read(psio_, PSIF_DFOCC_DENS);
    K->read(psio_, PSIF_DFOCC_INTS);
    for (int i = 0; i < naoccB; i++) {
         int ii = oo_idxBB->get(i + nfrzc,i + nfrzc);
         for (int j = 0; j < nfrzc; j++) {
              int jj = oo_idxBB->get(j,j);
              int ij = oo_idxBB->get(i + nfrzc,j);
              double sum = 0.0;
              for (int Q = 0; Q < nQ_ref; Q++) {
                   sum += 2.0 * G->get(Q,jj) * K->get(Q,ii);
                   sum += 2.0 * G->get(Q,ii) * K->get(Q,jj);
                   sum -= 4.0 * G->get(Q,ij) * K->get(Q,ij);
              }
              AooB->add(i, j, sum);
         }
    }
    G.reset();
    K.reset();

    if (print_ > 2) {
        AooA->print();
        AooB->print();
    }

}// else if (reference_ == "UNRESTRICTED")
    timer_off("Diagonal MO Hessian OO");
} // end diagonal_mohess_oo


}} // End Namespaces
