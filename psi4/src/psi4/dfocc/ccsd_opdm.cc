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
    SharedTensor2d T, U, X;
    timer_on("opdm");
    // if (reference_ == "RESTRICTED") {

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
    // Jc->print();

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

    //}// end if (reference_ == "RESTRICTED")

    // else if (reference_ == "UNRESTRICTED") {
    //}// else if (reference_ == "UNRESTRICTED")
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
