/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#include "psi4/libqt/qt.h"
#include "defines.h"
#include "dfocc.h"

using namespace psi;

namespace psi {
namespace dfoccwave {

void DFOCC::ccd_F_intr() {

    // RHF
    if (reference_ == "RESTRICTED") {
        // defs
        SharedTensor2d K, T, U, Tau;

        // OO block
        // F_mi = (1-\delta_{mi}) f_mi
        FijA->zero();
        for(int m=0; m<naoccA; ++m) {
             for(int i=0; i<naoccA; ++i) {
                 int mi = (m * naoccA) + i;
                 double value = 0.0;
                 if (m != i) value = FockA->get(m+nfrzc, i+nfrzc);
                 FijA->set(m, i, value);
             }
        }


        // read
        Tau = std::make_shared<Tensor2d>("T2 (Q|IA)", nQ, naoccA, navirA);
        Tau->read(psio_, PSIF_DFOCC_AMPS);

        // F_mi +=  \sum_{Q,e} Tau"_ie^Q b_me^Q
        FijA->contract332(false, true, navirA, bQiaA, Tau, 1.0, 1.0);

        // VV block
        // F_ae = (1-\delta_{ae}) f_ae
        FabA->zero();
        for(int a=0; a<navirA; ++a) {
             for(int e=0; e<navirA; ++e) {
                 int ae = (a * navirA) + e;
                 double value = 0.0;
                 if (a != e) value = FockA->get(a+noccA, e+noccA);
                 FabA->set(a, e, value);
             }
        }

        // F_ae -=  \sum_{Q,m} Tau'_ma^Q b_me^Q
        FabA->contract(true, false, navirA, navirA, nQ * naoccA, Tau, bQiaA, -1.0, 1.0);
    }// if (reference_ == "RESTRICTED")

    // UHF
    else if (reference_ == "UNRESTRICTED") {
        SharedTensor2d Tau, T;

        // F(M,I) OO Alpha Block
        // F(M,I) = (1-kronDelta(M,I)) f(M,I) + \sum_(E) 0.5*t(I,E)*f(M,E) + \sum_(Q) tQ*b(Q,MI) + \sum_(Q,E) Tau"(Q,IE)*b(Q,ME)
        // F(M,I) += (1-kronDelta(M,I)) f(M,I)
        for(int m=0; m<naoccA; ++m) {
             for(int i=0; i<naoccA; ++i) {
                 int mi = (m * naoccA) + i;
                 double value = 0.0;
                 if (m != i) value = FockA->get(m+nfrzc, i+nfrzc);
                 FijA->set(m, i, value);
             }
        }

        // F(M,I) +=  \sum_(Q,E) Tau"(IE,Q) * b (ME,Q)
        //Tau = std::make_shared<Tensor2d>("Tau2pp (Q|IA)", nQ, naoccA, navirA);
        Tau = std::make_shared<Tensor2d>("T2 (Q|IA)", nQ, naoccA, navirA);
        Tau->read(psio_, PSIF_DFOCC_AMPS);
        FijA->contract332(false, true, navirA, bQiaA, Tau, 1.0, 1.0);
        Tau.reset();

        // F(m,i) oo Beta Block
        // F(m,i) = (1-kronDelta(m,i)) f(m,i) + \sum_(e) 0.5*t(i,e)*f(m,e) + \sum_(Q) tQ*b(Q,mi) + \sum_(Q,e) Tau"(Q,ie)*b(Q,me)
        // F(m,i) += (1-kronDelta(m,i)) f(m,i)
        for(int m=0; m<naoccB; ++m) {
             for(int i=0; i<naoccB; ++i) {
                 int mi = (m * naoccB) + i;
                 double value = 0.0;
                 if (m != i) value = FockB->get(m+nfrzc, i+nfrzc);
                 FijB->set(m, i, value);
             }
        }

        // Fmi +=  \sum_(Q,e) Tau"(ie,Q) * b (me,Q)
        //Tau = std::make_shared<Tensor2d>("Tau2pp (Q|ia)", nQ, naoccB, navirB);
        Tau = std::make_shared<Tensor2d>("T2 (Q|ia)", nQ, naoccB, navirB);
        Tau->read(psio_, PSIF_DFOCC_AMPS);
        FijB->contract332(false, true, navirB, bQiaB, Tau, 1.0, 1.0);
        Tau.reset();

        // F(A,E) VV Alpha Block
        // F(A,E) = (1-kronDelta(A,E)) f(A,E) - \sum_(M) 0.5*t(M,A)*f(M,E) +  \sum_(Q) b(Q,AE)*tQ - \sum_(Q,M) Tau'(Q,MA)*b(Q,ME)
        // F(A,E) = (1-kronDelta(A,E)) f(A,E)
        for(int a=0; a<navirA; ++a) {
             for(int e=0; e<navirA; ++e) {
                 int ae = (a * navirA) + e;
                 double value = 0.0;
                 if (a != e) value = FockA->get(a+noccA, e+noccA);
                 FabA->set(a, e, value);
             }
        }

        // F(A,E) -= \sum_(Q,m) Tau'(Q,ma) * b(Q,me)
        //Tau = std::make_shared<Tensor2d>("Tau2p (Q|IA)", nQ, naoccA, navirA);
        Tau = std::make_shared<Tensor2d>("T2 (Q|IA)", nQ, naoccA, navirA);
        Tau->read(psio_, PSIF_DFOCC_AMPS);
        FabA->contract(true, false, navirA, navirA, nQ*naoccA, Tau, bQiaA, -1.0, 1.0);
        Tau.reset();

        // F(a,e) vv Beta Block
        // F(a,e) = (1-kronDelta(a,e)) f(a,e) - \sum_(m) 0.5*t(m,a)*f(m,e) +  \sum_(Q) b(Q,ae)*tQ - \sum_(Q,m) Tau'(Q,ma)*b(Q,me)
        // F(a,e) = (1-kronDelta(a,e)) f(a,e)
        for(int a=0; a<navirB; ++a) {
             for(int e=0; e<navirB; ++e) {
                 int ae = (a * navirB) + e;
                 double value = 0.0;
                 if (a != e) value = FockB->get(a+noccB, e+noccB);
                 FabB->set(a, e, value);
             }
        }

        // F(a,e) -= \sum_(Q,m) Tau'(Q,ma) * b(Q,me)
        //Tau = std::make_shared<Tensor2d>("Tau2p (Q|ia)", nQ, naoccB, navirB);
        Tau = std::make_shared<Tensor2d>("T2 (Q|ia)", nQ, naoccB, navirB);
        Tau->read(psio_, PSIF_DFOCC_AMPS);
        FabB->contract(true, false, navirB, navirB, nQ*naoccB, Tau, bQiaB, -1.0, 1.0);

    }// else if (reference_ == "UNRESTRICTED")

    // outfile->Printf("\tF int done.\n");

}  // end ccd_F_intr
}  // namespace dfoccwave
}  // namespace psi
