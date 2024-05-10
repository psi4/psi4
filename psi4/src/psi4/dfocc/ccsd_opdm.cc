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

/** Standard library includes */
#include "psi4/libqt/qt.h"
#include "defines.h"
#include "dfocc.h"

using namespace psi;

namespace psi {
namespace dfoccwave {

void DFOCC::ccsd_opdm() {
    timer_on("opdm");

     if (reference_ == "RESTRICTED") {
         SharedTensor2d T, U, X;

         // G1_ij = -(G_ij + G_ji)
         T = std::make_shared<Tensor2d>("T Intermediate <I|J>", naoccA, naoccA);
         U = std::make_shared<Tensor2d>("U Intermediate <I|J>", naoccA, naoccA);
         U->copy(GtijA);
         if (wfn_type_ == "DF-CCSD(T)") U->axpy(G1c_ij, 1.0);
         T->symmetrize(U);
         U.reset();
         T->scale(-2.0);
         G1c_oo->set_act_oo(nfrzc, naoccA, T);
         T.reset();

         //  G1_ab = -(G_ab + G_ba)
         T = std::make_shared<Tensor2d>("T Intermediate <A|B>", navirA, navirA);
         U = std::make_shared<Tensor2d>("U Intermediate <A|B>", navirA, navirA);
         U->copy(GtabA);
         if (wfn_type_ == "DF-CCSD(T)") U->axpy(G1c_ab, 1.0);
         T->symmetrize(U);
         U.reset();
         T->scale(-2.0);
         G1c_vv->set_act_vv(T);
         T.reset();
         // G1c_vv->print();

         // G1_ia = t_i^a + l_i^a
         T = std::make_shared<Tensor2d>("Corr OPDM <I|A>", naoccA, navirA);
         T->axpy(t1A, 1.0);
         T->axpy(l1A, 1.0);

         // G1_ia += \sum(me) U(ia,me) l_m^e
         U = std::make_shared<Tensor2d>("U2 (IA|JB)", naoccA, navirA, naoccA, navirA);
         U->read_symm(psio_, PSIF_DFOCC_AMPS);
         T->gemv(false, U, l1A, 1.0, 1.0);
         U.reset();

         // G1_ia -= \sum(me) t_m^a t_i^e l_m^e
         X = std::make_shared<Tensor2d>("X <I|M>", naoccA, naoccA);
         X->gemm(false, true, t1A, l1A, 1.0, 0.0);
         T->gemm(false, false, X, t1A, -1.0, 1.0);
         X.reset();

         // G1_ia -= \sum(m) t_m^a G_im
         T->gemm(false, false, GijA, t1A, -1.0, 1.0);

         // G1_ia += \sum(e) t_i^e G_ea
         T->gemm(false, false, t1A, GabA, 1.0, 1.0);
         //T->print();

         // (T) Contribution
         if (wfn_type_ == "DF-CCSD(T)") {
             T->axpy(G1c_ia, 1.0);
             G1c_ij.reset();
             G1c_ia.reset();
             G1c_ab.reset();
         }

         // set OV block
         G1c_ov->set_act_ov(nfrzc, T);
         T.reset();

         // Build G1_ai
         G1c_vo->trans(G1c_ov);

         // Build G1c
         G1c->set_oo(G1c_oo);
         G1c->set_ov(G1c_ov);
         G1c->set_vo(G1c_vo);
         G1c->set_vv(noccA, G1c_vv);
         // G1c->print();

         // Build G1
         G1->copy(G1c);
         for (int i = 0; i < noccA; i++) G1->add(i, i, 2.0);

         if (print_ > 2) {
             G1->print();
             double trace = G1->trace();
             outfile->Printf("\t trace: %12.12f \n", trace);
         }

    }// end if (reference_ == "RESTRICTED")

    else if (reference_ == "UNRESTRICTED") {
         SharedTensor2d T, T2, L2, Tau, X, Y, Z, V, U, L;

        // G1_IJ = -1/2(G_IJ + G_JI)
        T = std::make_shared<Tensor2d>("G Intermediate <I|J>", naoccA, naoccA);
        T->symmetrize(GtijA);
        T->scale(-1.0);
        G1c_ooA->set_act_oo(nfrzc, naoccA, T);
        T.reset();

        // G1_ij = -1/2(G_ij + G_ji)
        T = std::make_shared<Tensor2d>("G Intermediate <i|j>", naoccB, naoccB);
        T->symmetrize(GtijB);
        T->scale(-1.0);
        G1c_ooB->set_act_oo(nfrzc, naoccB, T);
        T.reset();

        //  G1_AB = -1/2(G_AB + G_AB)
        T = std::make_shared<Tensor2d>("G Intermediate <A|B>", navirA, navirA);
        T->symmetrize(GtabA);
        T->scale(-1.0);
        G1c_vvA->set_act_vv(T);
        T.reset();
        //G1c_vvA->print();

        //  G1_ab = -1/2(G_ab + G_ab)
        T = std::make_shared<Tensor2d>("G Intermediate <a|b>", navirB, navirB);
        T->symmetrize(GtabB);
        T->scale(-1.0);
        G1c_vvB->set_act_vv(T);
        T.reset();

        if (wfn_type_ != "DF-CCSD(T)") {
            // G1_IA = 0.5 * t_I^A + 0.5 * l_I^A
            T = std::make_shared<Tensor2d>("Corr OPDM <I|A>", naoccA, navirA);
            T->axpy(t1A, 0.5);
            T->axpy(l1A, 0.5);

            // G1_IA += 0.5 * \sum(ME) t(IM,AE) l_M^E
            U = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
            U->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
            Z = std::make_shared<Tensor2d>("Z (IA|JB)", naoccA, navirA, naoccA, navirA);
            Z->sort(1324, U, 1.0, 0.0);
            U.reset();
            T->gemv(false, Z, l1A, 0.5, 1.0);
            Z.reset();

            // G1_IA -= 0.5 * \sum(ME) t_M^A t_I^E l_M^E
            X = std::make_shared<Tensor2d>("X <I|M>", naoccA, naoccA);
            X->gemm(false, true, t1A, l1A, 1.0, 0.0);
            T->gemm(false, false, X, t1A, -0.5, 1.0);
            X.reset();

            // G1_IA += 0.5 * \sum(me) t(Im,Ae) l_m^e
            U = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
            U->read(psio_, PSIF_DFOCC_AMPS);
            Z = std::make_shared<Tensor2d>("Z (IA|me)", naoccA, navirA, naoccB, navirB);
            Z->sort(1324, U, 1.0, 0.0);
            U.reset();
            T->gemv(false, Z, l1B, 0.5, 1.0);
            Z.reset();

            // G1_IA -= 0.5 * \sum(M) t_M^A G_IM
            T->gemm(false, false, GijA, t1A, -0.5, 1.0);

            // G1_IA += 0.5 * \sum(E) t_I^E G_EA
            T->gemm(false, false, t1A, GabA, 0.5, 1.0);
            //T->print();

            // set OV block
            G1c_ovA->set_act_ov(nfrzc, T);
            T.reset();

            // Build G1_ai
            G1c_voA->trans(G1c_ovA);

            // G1_ia = 0.5 * t_i^a + 0.5 * l_i^a
            T = std::make_shared<Tensor2d>("Corr OPDM <i|a>", naoccB, navirB);
            T->axpy(t1B, 0.5);
            T->axpy(l1B, 0.5);

            // G1_ia += 0.5 * \sum(me) t(im,ae) l_m^e
            U = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
            U->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
            Z = std::make_shared<Tensor2d>("Z (ia|jb)", naoccB, navirB, naoccB, navirB);
            Z->sort(1324, U, 1.0, 0.0);
            U.reset();
            T->gemv(false, Z, l1B, 0.5, 1.0);
            Z.reset();

            // G1_ia -= 0.5 * \sum(me) t_m^a t_i^e l_m^e
            X = std::make_shared<Tensor2d>("X <i|m>", naoccB, naoccB);
            X->gemm(false, true, t1B, l1B, 1.0, 0.0);
            T->gemm(false, false, X, t1B, -0.5, 1.0);
            X.reset();

            // G1_ia += 0.5 * \sum(ME) t(Mi,Ea) l_M^E
            U = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
            U->read(psio_, PSIF_DFOCC_AMPS);
            Z = std::make_shared<Tensor2d>("Z (ia|ME)", naoccB, navirB, naoccA, navirA);
            Z->sort(2413, U, 1.0, 0.0);
            U.reset();
            T->gemv(false, Z, l1A, 0.5, 1.0);
            Z.reset();

            // G1_ia -= 0.5 * \sum(m) t_m^a G_im
            T->gemm(false, false, GijB, t1B, -0.5, 1.0);

            // G1_ia += 0.5 * \sum(e) t_i^e G_ea
            T->gemm(false, false, t1B, GabB, 0.5, 1.0);

            // set OV block
            G1c_ovB->set_act_ov(nfrzc, T);
            T.reset();

            // Build G1_ai
            G1c_voB->trans(G1c_ovB);
        }

        if (wfn_type_ == "DF-CCSD(T)") {
            // Build G1c

            G1c_ooA->zero_off_diagonal();
            for (int i = 0; i < naoccA; i++) {
                G1c_ooA->add(nfrzc + i, nfrzc + i, G1c_iiA->get(i));
            }
            G1c_iiA.reset();
            G1cA->set_oo(G1c_ooA);


            G1c_vvA->zero_off_diagonal();
            for (int i = 0; i < navirA; i++) {
                G1c_vvA->add(i, i, G1c_aaA->get(i));
            }
            G1c_aaA.reset();
            G1cA->set_vv(noccA, G1c_vvA);

            // Build G1c
            G1c_ooB->zero_off_diagonal();
            for (int i = 0; i < naoccB; i++) {
                G1c_ooB->add(nfrzc + i, nfrzc + i, G1c_iiB->get(i));
            }
            G1c_iiB.reset();
            G1cB->set_oo(G1c_ooB);

            G1c_vvB->zero_off_diagonal();
            for (int i = 0; i < navirB; i++) {
                G1c_vvB->add(i, i, G1c_aaB->get(i));
            }
            G1c_aaB.reset();
            G1cB->set_vv(noccB, G1c_vvB);

            // Build G1
            G1A->copy(G1cA);
            G1B->copy(G1cB);
            for (int i = 0; i < noccA; i++) G1A->add(i, i, 1.0);
            for (int i = 0; i < noccB; i++) G1B->add(i, i, 1.0);
            G1A->zero_off_diagonal();
            G1B->zero_off_diagonal();
        }
        else {
            // Build G1c
            G1cA->set_oo(G1c_ooA);
            G1cA->set_ov(G1c_ovA);
            G1cA->set_vo(G1c_voA);
            G1cA->set_vv(noccA, G1c_vvA);

            // Build G1c
            G1cB->set_oo(G1c_ooB);
            G1cB->set_ov(G1c_ovB);
            G1cB->set_vo(G1c_voB);
            G1cB->set_vv(noccB, G1c_vvB);

            // Build G1
            G1A->copy(G1cA);
            G1B->copy(G1cB);
            for (int i = 0; i < noccA; i++) G1A->add(i, i, 1.0);
            for (int i = 0; i < noccB; i++) G1B->add(i, i, 1.0);
        }

        // print
        if (print_ > 2) {
            G1A->print();
            G1B->print();
            double trace = G1A->trace();
            outfile->Printf("\t Alpha trace: %12.12f \n", trace);
            trace = G1B->trace();
            outfile->Printf("\t Beta trace: %12.12f \n", trace);
        }


    }// else if (reference_ == "UNRESTRICTED")
    timer_off("opdm");
}  // end ccsd_opdm

//=======================================================
//       Diagonal OPDM
//=======================================================
void DFOCC::ccsd_diagonal_opdm() {
    SharedTensor2d T, U, X;
    timer_on("opdm");
    // if (reference_ == "RESTRICTED") {

    // G1_ij = -(G_ij + G_ji)
    T = std::make_shared<Tensor2d>("T Intermediate <I|J>", naoccA, naoccA);
    U = std::make_shared<Tensor2d>("U Intermediate <I|J>", naoccA, naoccA);
    U->copy(GtijA);
    T->symmetrize(U);
    U.reset();
    T->scale(-2.0);
    G1c_oo->set_act_oo(nfrzc, naoccA, T);
    G1c_oo->zero_off_diagonal();
    T.reset();

    //  G1_ab = -(G_ab + G_ba)
    T = std::make_shared<Tensor2d>("T Intermediate <A|B>", navirA, navirA);
    U = std::make_shared<Tensor2d>("U Intermediate <A|B>", navirA, navirA);
    U->copy(GtabA);
    T->symmetrize(U);
    U.reset();
    T->scale(-2.0);
    G1c_vv->set_act_vv(T);
    G1c_vv->zero_off_diagonal();
    T.reset();

    // (T) Contribution
    if (wfn_type_ == "DF-CCSD(T)") {
        for (int i = 0; i < naoccA; i++) G1c_oo->add(i + nfrzc, i + nfrzc, G1c_ii->get(i));
        for (int a = 0; a < navirA; a++) G1c_vv->add(a, a, G1c_aa->get(a));
        G1c_ii.reset();
        G1c_aa.reset();
    }

    // Build G1c
    G1c->set_oo(G1c_oo);
    G1c->set_vv(noccA, G1c_vv);
    // G1c->print();

    // Build G1
    G1->copy(G1c);
    for (int i = 0; i < noccA; i++) G1->add(i, i, 2.0);

    if (print_ > 2) {
        G1->print();
        double trace = G1->trace();
        outfile->Printf("\t trace: %12.12f \n", trace);
    }

    //}// end if (reference_ == "RESTRICTED")

    // else if (reference_ == "UNRESTRICTED") {
    //}// else if (reference_ == "UNRESTRICTED")
    timer_off("opdm");
}  // end ccsd_diagonal_opdm

}  // namespace dfoccwave
}  // namespace psi
