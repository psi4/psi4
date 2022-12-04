/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

void DFOCC::ccsd_F_intr() {

    // RHF
    if (reference_ == "RESTRICTED") {
        // defs
        SharedTensor2d K, T, U, Tau;

        // OO block
        // F_mi =  \sum_{Q} t_Q b_mi^Q
        FijA->gemv(true, bQijA, T1c, 1.0, 0.0);

        // F_mi +=  \sum_{Q,e} Tau"_ie^Q b_me^Q
        Tau = std::make_shared<Tensor2d>("Tau2pp (Q|IA)", nQ, naoccA, navirA);
        Tau->read(psio_, PSIF_DFOCC_AMPS);
        FijA->contract332(false, true, navirA, bQiaA, Tau, 1.0, 1.0);
        Tau.reset();

        // VV block
        // F_ae =  \sum_{Q} t_Q b_ae^Q
        FabA->gemv(true, bQabA, T1c, 1.0, 0.0);

        // F_ae -=  \sum_{Q,m} Tau'_ma^Q b_me^Q
        Tau = std::make_shared<Tensor2d>("Tau2p (Q|IA)", nQ, naoccA, navirA);
        Tau->read(psio_, PSIF_DFOCC_AMPS);
        FabA->contract(true, false, navirA, navirA, nQ * naoccA, Tau, bQiaA, -1.0, 1.0);
        Tau.reset();

        // OV block
        // F_me +=  \sum_{Q} t_Q b_me^Q
        FiaA->gemv(true, bQiaA, T1c, 1.0, 0.0);

        // F_me -=  \sum_{Q,n} t_nm^Q b_ne^Q
        T = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
        T->read(psio_, PSIF_DFOCC_AMPS);
        FiaA->contract(true, false, naoccA, navirA, nQ * naoccA, T, bQiaA, -1.0, 1.0);
        T.reset();
        FiaA->write(psio_, PSIF_DFOCC_AMPS);

        // Ft_mi = F_mi + 1/2 \sum_{e} t_i^e F_me
        FtijA->gemm(false, true, FiaA, t1A, 0.5, 0.0);
        FtijA->add(FijA);

        // Ft_ae = F_ae - 1/2 \sum_{m} t_m^a F_me
        FtabA->gemm(true, false, t1A, FiaA, -0.5, 0.0);
        FtabA->add(FabA);

        // outfile->Printf("\tF int done.\n");

    }  // if (reference_ == "RESTRICTED")

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

        // F(M,I) +=  \sum_(E) 0.5 * t(I,E) * f(M,E)
        //FijA->gemm(false, true, t1A, FaovA, 0.5, 1.0);
        for(int m=0; m<naoccA; ++m) {
             for(int i=0; i<naoccA; ++i) {
                 int mi = (m * naoccA) + i;
                 double value = 0.0;
                 for(int e=0; e<navirA; ++e) {
                     value += 0.5 * t1A->get(i,e) * FockA->get(m+nfrzc, e+noccA);
                 }
                 FijA->add(m, i, value);
             }
        }

        // F(M,I) +=  \sum_(Q) tQ * b(Q,MI)
        FijA->gemv(true, bQijA, T1c, 1.0, 1.0);

        // F(M,I) +=  \sum_(Q,E) Tau"(IE,Q) * b (ME,Q)
        Tau = std::make_shared<Tensor2d>("Tau2pp (Q|IA)", nQ, naoccA, navirA);
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

        // F(m,i) +=  \sum_(e) 0.5 * t(i,e) * f(m,e)
        //FijB->gemm(false, true, t1B, FaovB, 0.5, 1.0);
        for(int m=0; m<naoccB; ++m) {
             for(int i=0; i<naoccB; ++i) {
                 int mi = (m * naoccB) + i;
                 double value = 0.0;
                 for(int e=0; e<navirB; ++e) {
                     value += 0.5 * t1B->get(i,e) * FockB->get(m+nfrzc, e+noccB);
                 }
                 FijB->add(m, i, value);
             }
        }

        // Fmi +=  \sum_(Q) tQ * b(mi,Q)
        FijB->gemv(true, bQijB, T1c, 1.0, 1.0);

        // Fmi +=  \sum_(Q,e) Tau"(ie,Q) * b (me,Q)
        Tau = std::make_shared<Tensor2d>("Tau2pp (Q|ia)", nQ, naoccB, navirB);
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

        // F(A,E) -= \sum_(M) 0.5 * t(M,A) * f(M,E)
        //FabA->gemm(true, false, t1A, FaovA, -0.5, 1.0);
        for(int a=0; a<navirA; ++a) {
             for(int e=0; e<navirA; ++e) {
                 int ae = (a * navirA) + e;
                 double value = 0.0;
                 for(int m=0; m<naoccA; ++m) {
                     value -= 0.5 * t1A->get(m,a) * FockA->get(m+nfrzc, e+noccA);
                 }
                 FabA->add(a, e, value);
             }
        }

        // F(A,E) += \sum_(Q) b(Q,ae) * tQ
        FabA->gemv(true, bQabA, T1c, 1.0, 1.0);

        // F(A,E) -= \sum_(Q,m) Tau'(Q,ma) * b(Q,me)
        Tau = std::make_shared<Tensor2d>("Tau2p (Q|IA)", nQ, naoccA, navirA);
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

        // F(a,e) -= \sum_(m) 0.5 * t(m,a) * f(m,e)
        //FabB->gemm(true, false, t1B, FaovB, -0.5, 1.0);
        for(int a=0; a<navirB; ++a) {
             for(int e=0; e<navirB; ++e) {
                 int ae = (a * navirB) + e;
                 double value = 0.0;
                 for(int m=0; m<naoccB; ++m) {
                     value -= 0.5 * t1B->get(m,a) * FockB->get(m+nfrzc, e+noccB);
                 }
                 FabB->add(a, e, value);
             }
        }

        // F(a,e) += \sum_(Q) b(Q,ae) * tQ
        FabB->gemv(true, bQabB, T1c, 1.0, 1.0);

        // F(a,e) -= \sum_(Q,m) Tau'(Q,ma) * b(Q,me)
        Tau = std::make_shared<Tensor2d>("Tau2p (Q|ia)", nQ, naoccB, navirB);
        Tau->read(psio_, PSIF_DFOCC_AMPS);
        FabB->contract(true, false, navirB, navirB, nQ*naoccB, Tau, bQiaB, -1.0, 1.0);
        Tau.reset();

        // F(M,E) OV Alpha Block
        // F(M,E) = f(M,E) + \sum_(Q) b(Q,ME)*tQ - \sum_(Q,N) t(NM,Q) * b(NE,Q)
        // F(M,E) = FockA(m,e)
        for(int m=0; m<naoccA; ++m) {
             for(int e=0; e<navirA; ++e) {
                 int me = (m * navirA) + e;
                 //FiaA->set(m, e, FaovA->get(m,e));
                 FiaA->set(m, e, FockA->get(m+nfrzc,e+noccA));
             }
        }
        // F(m,e) +=  \sum_(Q) tQ * b(me,Q)
        FiaA->gemv(true, bQiaA, T1c, 1.0, 1.0);
        // F(m,e) -=  \sum_(Q,n) t(nm,Q) * b(ne,Q)
        T = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
        T->read(psio_, PSIF_DFOCC_AMPS);
        FiaA->contract(true, false, naoccA, navirA, nQ * naoccA, T, bQiaA, -1.0, 1.0);
        FiaA->write(psio_, PSIF_DFOCC_INTS);
        T.reset();
        FiaA->write(psio_, PSIF_DFOCC_AMPS);

        // F(m,e) ov Beta Block
        // F(m,e) = f(m,e) + \sum_(Q) b(Q,me)*tQ - \sum_(Q,n) t(nm,Q) * b(ne,Q)
        // F(m,e) = FockA(m,e)
        for(int m=0; m<naoccB; ++m) {
             for(int e=0; e<navirB; ++e) {
                 int me = (m * navirB) + e;
                 //FiaB->set(m, e, FaovB->get(m,e));
                 FiaB->set(m, e, FockB->get(m+nfrzc,e+noccB));
             }
        }
        // F(m,e) +=  \sum_(Q) tQ * b(me,Q)
        FiaB->gemv(true, bQiaB, T1c, 1.0, 1.0);
        // F(m,e) -=  \sum_(Q,n) t(nm,Q) * b(ne,Q)
        T = std::make_shared<Tensor2d>("T1 (Q|ij)", nQ, naoccB, naoccB);
        T->read(psio_, PSIF_DFOCC_AMPS);
        FiaB->contract(true, false, naoccB, navirB, nQ * naoccB, T, bQiaB, -1.0, 1.0);
        FiaB->write(psio_, PSIF_DFOCC_INTS);
        T.reset();
        FiaB->write(psio_, PSIF_DFOCC_AMPS);

        ////////////// F_ intermediates ////////////////

        // Ft_MI = F(M,I) + 0.5 * \sum_{E} t(I,E) * F(M,E)
        //FtijA = std::make_shared<Tensor2d>("Ftilde <I|J>", naoccA, naoccA);
        FtijA->gemm(false, true, FiaA, t1A, 0.5, 0.0);
        FtijA->axpy(FijA, 1.0);
        FtijA->write(psio_, PSIF_DFOCC_INTS);

        // Ft_mi = F(m,i) + 0.5 * \sum_{e} t(i,e) * F(m,e)
        //FtijB = std::make_shared<Tensor2d>("Ftilde <i|j>", naoccB, naoccB);
        FtijB->gemm(false, true, FiaB, t1B, 0.5, 0.0);
        FtijB->axpy(FijB, 1.0);
        FtijB->write(psio_, PSIF_DFOCC_INTS);

        // Ft_AE = F(A,E) - 0.5 * \sum_{M} t(M,A) * F(M,E)
        //FtabA = std::make_shared<Tensor2d>("Ftilde <A|B>", navirA, navirA);
        FtabA->gemm(true, false, t1A, FiaA, -0.5, 0.0);
        FtabA->axpy(FabA, 1.0);
        FtabA->write(psio_, PSIF_DFOCC_INTS);

        // Ft_ae = F(a,e) - 0.5 * \sum_{m} t(m,a) * F(m,e)
        //FtabB = std::make_shared<Tensor2d>("Ftilde <a|b>", navirB, navirB);
        FtabB->gemm(true, false, t1B, FiaB, -0.5, 0.0);
        FtabB->axpy(FabB, 1.0);
        FtabB->write(psio_, PSIF_DFOCC_INTS);

    }  // else if (reference_ == "UNRESTRICTED")

}  // end ccsd_F_intr
}  // namespace dfoccwave
}  // namespace psi
