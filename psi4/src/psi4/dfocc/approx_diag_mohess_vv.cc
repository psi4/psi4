/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/** Standard library includes */
#include "psi4/libqt/qt.h"
#include "defines.h"
#include "dfocc.h"

using namespace psi;

namespace psi {
namespace dfoccwave {

void DFOCC::approx_diag_mohess_vv() {
    timer_on("Diagonal MO Hessian VV");
    SharedTensor2d K;
    SharedTensor2d G;
    SharedTensor2d T;

    //=========================
    // Diagonal MO Hessian
    //=========================

    if (reference_ == "RESTRICTED") {
        // A_ab = [\delta(ab) -1] 2(Faa + Fbb) + 2 h_aa G_bb - 4 h_ab G_ab + 2 h_bb G_aa
        AvvA->zero();
        for (int a = 0; a < nfrzv; a++) {
            for (int b = 0; b < navirA; b++) {

                 // A_ab = [\delta(ab) -1] 2(Faa + Fbb) 
                 double value = -2.0 * (GF->get(a + npop, a + npop) + GF->get(b + noccA, b + noccA));
                 if (a != b) AvvA->set(a,b,value);

                 // A_ab = 2 h_aa G_bb - 4 h_ab G_ab + 2 h_bb G_aa 
                 value = 2.0 * HmoA->get(a + npop, a + npop) * G1->get(b + noccA, b + noccA);
                 value -= 4.0 * HmoA->get(a + npop, b + noccA) * G1->get(a + npop, b + noccA);
                 value += 2.0 * HmoA->get(b + noccA, b + noccA) * G1->get(a + npop, a + npop);
                 //value -= 2*lshift_parameter;
                 AvvA->add(a, b, value);
            }
        }

        // Separable Part
        // A_ab += -4\sum_{Q} G_ab^Q b_ab^Q + 2\sum_{Q} G_aa^Q b_bb^Q + 2\sum_{Q} G_bb^Q b_aa^Q
        G = std::make_shared<Tensor2d>("3-Index Separable TPDM (Q|VV)", nQ_ref, nvirA,  nvirA);
        K = std::make_shared<Tensor2d>("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA);
        G->read(psio_, PSIF_DFOCC_DENS, true, true);
        K->read(psio_, PSIF_DFOCC_INTS, true, true);
        for (int a = 0; a < nfrzv; a++) {
             for (int b = 0; b < navirA; b++) {
                  int ab = ( (a+npop)*nvirA ) + b;
                  int aa = ( (a+npop)*nvirA ) + a + npop;
                  int bb = (b*nvirA) + b;
                  double sum = 0.0;
                  for (int Q = 0; Q < nQ_ref; Q++) {
                      sum = -2.0 * G->get(Q, ab) * K->get(Q, ab);
                      sum += G->get(Q, aa) * K->get(Q, bb);
                      sum += G->get(Q, bb) * K->get(Q, aa);
                  }
                  AvvA->add(a, b, 2.0 * sum);
             }
        }
        G.reset();
        K.reset();

        // Corr Part
        if (wfn_type_ != "DF-OMP2") {
            // A_ab += -4\sum_{Q} G_ab^Q b_ab^Q + 2\sum_{Q} G_aa^Q b_bb^Q + 2\sum_{Q} G_bb^Q b_aa^Q
            G = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|VV)", nQ, nvirA, nvirA);
            K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA);
            G->read(psio_, PSIF_DFOCC_DENS, true, true);
            K->read(psio_, PSIF_DFOCC_INTS, true, true);
            for (int a = 0; a < nfrzv; a++) {
                 for (int b = 0; b < navirA; b++) {
                      int ab = ( (a+npop)*nvirA ) + b;
                      int aa = ( (a+npop)*nvirA ) + a + npop;
                      int bb = (b*nvirA) + b;
                      double sum = 0.0;
                      for (int Q = 0; Q < nQ; Q++) {
                          sum = -2.0 * G->get(Q, ab) * K->get(Q, ab);
                          sum += G->get(Q, aa) * K->get(Q, bb);
                          sum += G->get(Q, bb) * K->get(Q, aa);
                      }
                      AvvA->add(a, b, 2.0 * sum);
                 }
            }
            G.reset();
            K.reset();
        }
        
        // The factor of 1/2 is necessary to convert spin-free Hessian to Alpha spin case
        AvvA->scale(0.5);
        // Level-shifting
        for (int a = 0; a < nfrzv; a++) {
            for (int b = 0; b < navirA; b++) {
                 //AvvA->subtract(a, b, lshift_parameter);
                 AvvA->add(a, b, lshift_parameter); // this the correct way
            }
        }
        // damping
        //AvvA->scale(1.0+lshift_parameter);

        if (print_ > 2) AvvA->print();
    }  // if (reference_ == "RESTRICTED")


    else if (reference_ == "UNRESTRICTED") {
        // Alpha part 
        // A_ab = [\delta(ab) -1] 2(Faa + Fbb) + 2 h_aa G_bb - 4 h_ab G_ab + 2 h_bb G_aa
        AvvA->zero();
        for (int a = 0; a < nfrzv; a++) {
            for (int b = 0; b < navirA; b++) {

                 // A_ab = [\delta(ab) -1] 2(Faa + Fbb) 
                 double value = -2.0 * (GFA->get(a + npop, a + npop) + GFA->get(b + noccA, b + noccA));
                 if (a != b) AvvA->set(a,b,value);

                 // A_ab = 2 h_aa G_bb - 4 h_ab G_ab + 2 h_bb G_aa 
                 value = 2.0 * HmoA->get(a + npop, a + npop) * G1A->get(b + noccA, b + noccA);
                 value -= 4.0 * HmoA->get(a + npop, b + noccA) * G1A->get(a + npop, b + noccA);
                 value += 2.0 * HmoA->get(b + noccA, b + noccA) * G1A->get(a + npop, a + npop);
                 //value -= 2*lshift_parameter;
                 AvvA->add(a, b, value);
            }
        }

        // Separable Part
        // A_ab += -4\sum_{Q} G_ab^Q b_ab^Q + 2\sum_{Q} G_aa^Q b_bb^Q + 2\sum_{Q} G_bb^Q b_aa^Q
        G = std::make_shared<Tensor2d>("3-Index Separable TPDM (Q|VV)", nQ_ref, nvirA,  nvirA);
        K = std::make_shared<Tensor2d>("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA);
        G->read(psio_, PSIF_DFOCC_DENS, true, true);
        K->read(psio_, PSIF_DFOCC_INTS, true, true);
        for (int a = 0; a < nfrzv; a++) {
             for (int b = 0; b < navirA; b++) {
                  int ab = ( (a+npop)*nvirA ) + b;
                  int aa = ( (a+npop)*nvirA ) + a + npop;
                  int bb = (b*nvirA) + b;
                  double sum = 0.0;
                  for (int Q = 0; Q < nQ_ref; Q++) {
                      sum = -2.0 * G->get(Q, ab) * K->get(Q, ab);
                      sum += G->get(Q, aa) * K->get(Q, bb);
                      sum += G->get(Q, bb) * K->get(Q, aa);
                  }
                  AvvA->add(a, b, 2.0 * sum);
             }
        }
        G.reset();
        K.reset();

        // Corr Part
        if (wfn_type_ != "DF-OMP2") {
            // A_ab += -4\sum_{Q} G_ab^Q b_ab^Q + 2\sum_{Q} G_aa^Q b_bb^Q + 2\sum_{Q} G_bb^Q b_aa^Q
            G = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|VV)", nQ, nvirA, nvirA);
            K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA);
            G->read(psio_, PSIF_DFOCC_DENS, true, true);
            K->read(psio_, PSIF_DFOCC_INTS, true, true);
            for (int a = 0; a < nfrzv; a++) {
                 for (int b = 0; b < navirA; b++) {
                      int ab = ( (a+npop)*nvirA ) + b;
                      int aa = ( (a+npop)*nvirA ) + a + npop;
                      int bb = (b*nvirA) + b;
                      double sum = 0.0;
                      for (int Q = 0; Q < nQ; Q++) {
                          sum = -2.0 * G->get(Q, ab) * K->get(Q, ab);
                          sum += G->get(Q, aa) * K->get(Q, bb);
                          sum += G->get(Q, bb) * K->get(Q, aa);
                      }
                      AvvA->add(a, b, 2.0 * sum);
                 }
            }
            G.reset();
            K.reset();
        }
        
        // Level-shifting
        for (int a = 0; a < nfrzv; a++) {
            for (int b = 0; b < navirA; b++) {
                 //AvvA->subtract(a, b, lshift_parameter);
                 AvvA->add(a, b, lshift_parameter); // this the correct way
            }
        }

        // apply damping
        //for (int a = 0; a < nfrzv; a++) {
        //    for (int b = 0; b < navirA; b++) {
        //         AvvA->set(a, b, AvvA->get(a,b) * (1 + lshift_parameter)); 
        //    }
        //}

        // just switch to steepest-descent
        //for (int a = 0; a < nfrzv; a++) {
        //    for (int b = 0; b < navirA; b++) {
        //         AvvA->set(a, b, 0.5); 
        //    }
        //}


        // Beta part 
        // A_ab = [\delta(ab) -1] 2(Faa + Fbb) + 2 h_aa G_bb - 4 h_ab G_ab + 2 h_bb G_aa
        AvvB->zero();
        for (int a = 0; a < nfrzv; a++) {
            for (int b = 0; b < navirB; b++) {

                 // A_ab = [\delta(ab) -1] 2(Faa + Fbb) 
                 double value = -2.0 * (GFB->get(a + npop, a + npop) + GFB->get(b + noccB, b + noccB));
                 if (a != b) AvvB->set(a,b,value);

                 // A_ab = 2 h_aa G_bb - 4 h_ab G_ab + 2 h_bb G_aa 
                 value = 2.0 * HmoB->get(a + npop, a + npop) * G1B->get(b + noccB, b + noccB);
                 value -= 4.0 * HmoB->get(a + npop, b + noccB) * G1B->get(a + npop, b + noccB);
                 value += 2.0 * HmoB->get(b + noccB, b + noccB) * G1B->get(a + npop, a + npop);
                 //value -= 2*lshift_parameter;
                 AvvB->add(a, b, value);
            }
        }

        // Separable Part
        // A_ab += -4\sum_{Q} G_ab^Q b_ab^Q + 2\sum_{Q} G_aa^Q b_bb^Q + 2\sum_{Q} G_bb^Q b_aa^Q
        G = std::make_shared<Tensor2d>("3-Index Separable TPDM (Q|vv)", nQ_ref, nvirB,  nvirB);
        K = std::make_shared<Tensor2d>("DF_BASIS_SCF B (Q|vv)", nQ_ref, nvirB, nvirB);
        G->read(psio_, PSIF_DFOCC_DENS, true, true);
        K->read(psio_, PSIF_DFOCC_INTS, true, true);
        for (int a = 0; a < nfrzv; a++) {
             for (int b = 0; b < navirB; b++) {
                  int ab = ( (a+npop)*nvirB ) + b;
                  int aa = ( (a+npop)*nvirB ) + a + npop;
                  int bb = (b*nvirB) + b;
                  double sum = 0.0;
                  for (int Q = 0; Q < nQ_ref; Q++) {
                      sum = -2.0 * G->get(Q, ab) * K->get(Q, ab);
                      sum += G->get(Q, aa) * K->get(Q, bb);
                      sum += G->get(Q, bb) * K->get(Q, aa);
                  }
                  AvvB->add(a, b, 2.0 * sum);
             }
        }
        G.reset();
        K.reset();

        // Corr Part
        if (wfn_type_ != "DF-OMP2") {
            // A_ab += -4\sum_{Q} G_ab^Q b_ab^Q + 2\sum_{Q} G_aa^Q b_bb^Q + 2\sum_{Q} G_bb^Q b_aa^Q
            G = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|vv)", nQ, nvirB, nvirB);
            K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|vv)", nQ, nvirB, nvirB);
            G->read(psio_, PSIF_DFOCC_DENS, true, true);
            K->read(psio_, PSIF_DFOCC_INTS, true, true);
            for (int a = 0; a < nfrzv; a++) {
                 for (int b = 0; b < navirB; b++) {
                      int ab = ( (a+npop)*nvirB ) + b;
                      int aa = ( (a+npop)*nvirB ) + a + npop;
                      int bb = (b*nvirB) + b;
                      double sum = 0.0;
                      for (int Q = 0; Q < nQ; Q++) {
                          sum = -2.0 * G->get(Q, ab) * K->get(Q, ab);
                          sum += G->get(Q, aa) * K->get(Q, bb);
                          sum += G->get(Q, bb) * K->get(Q, aa);
                      }
                      AvvB->add(a, b, 2.0 * sum);
                 }
            }
            G.reset();
            K.reset();
        }
        
        // Level-shifting
        for (int a = 0; a < nfrzv; a++) {
            for (int b = 0; b < navirB; b++) {
                 //AvvB->subtract(a, b, lshift_parameter);
                 AvvB->add(a, b, lshift_parameter); // this the correct way
            }
        }

        // apply damping
        //for (int a = 0; a < nfrzv; a++) {
        //    for (int b = 0; b < navirB; b++) {
        //         AvvB->set(a, b, AvvB->get(a,b) * (1 + lshift_parameter)); 
        //    }
        //}

        // just switch to steepest-descent
        //for (int a = 0; a < nfrzv; a++) {
        //    for (int b = 0; b < navirB; b++) {
        //         AvvB->set(a, b, 0.5); 
        //    }
        //}

        //AvvA->scale(4.0);
        //AvvB->scale(4.0);
        if (print_ > 2) {
            AvvA->print();
            AvvB->print();
        }

    }  // else if (reference_ == "UNRESTRICTED")
    timer_off("Diagonal MO Hessian VV");
}  // end diagonal_mohess_vo

}  // namespace dfoccwave
}  // namespace psi
