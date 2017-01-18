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

void DFOCC::fock()
{

    SharedTensor2d K, L, M;
    timer_on("Fock");
if (reference_ == "RESTRICTED") {
    // Build J_Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA * noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    for (int Q = 0; Q < nQ_ref; Q++) {
         double value = 0.0;
         for (int i = 0; i < noccA; i++) {
              int ii = oo_idxAA->get(i,i);
              value += K->get(Q,ii);
         }
         Jc->set(Q, 2.0*value);
    }

    // F_ij = h_ij + \sum_{Q} b_ij^Q J^Q - \sum_{Q} \sum_{m} b_mi^Q b_mj^Q
    FooA->copy(HooA);
    FooA->gemv(true, K, Jc, 1.0, 1.0);
    FooA->contract(true, false, noccA, noccA, nQ_ref * noccA, K, K, -1.0, 1.0);

    // F_ia = h_ia + \sum_{Q} b_ia^Q J^Q - \sum_{Q} \sum_{m} b_mi^Q b_ma^Q
    FovA->copy(HovA);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA * nvirA));
    L->read(psio_, PSIF_DFOCC_INTS);
    FovA->gemv(true, L, Jc, 1.0, 1.0);
    FovA->contract(true, false, noccA, nvirA, nQ_ref * noccA, K, L, -1.0, 1.0);
    K.reset();

    // F_ai = F_ia
    FvoA = FovA->transpose();

    // F_ab = h_ab + \sum_{Q} b_ab^Q J^Q - \sum_{Q} \sum_{m} b_ma^Q b_mb^Q
    FvvA->copy(HvvA);
    FvvA->contract(true, false, nvirA, nvirA, nQ_ref * noccA, L, L, -1.0, 1.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    FvvA->gemv(true, K, Jc, 1.0, 1.0);
    K.reset();

    // Set Fock Matrix
    FockA->set_oo(FooA);
    FockA->set_ov(FovA);
    FockA->set_vo(FvoA);
    FockA->set_vv(noccA, FvvA);

    if (print_ > 2) FockA->print();

    /*
    // Diagonalize
    SharedTensor1d eigA = std::shared_ptr<Tensor1d>(new Tensor1d("epsilon <P|Q>", nmo_));
    SharedTensor2d UmoA = std::shared_ptr<Tensor2d>(new Tensor2d("UmoA", nmo_, nmo_));
    FockA->diagonalize(UmoA, eigA, cutoff);
    eigA.reset();

    // Get new MOs
    SharedTensor2d Ca_new = std::shared_ptr<Tensor2d>(new Tensor2d("New alpha MO coefficients", nso_, nmo_));
    Ca_new->gemm(false, false, CmoA, UmoA, 1.0, 0.0);
    UmoA.reset();
    CmoA->copy(Ca_new);
    Ca_new.reset();
    */
}// end if (reference_ == "RESTRICTED")



else if (reference_ == "UNRESTRICTED") {
    // Build J_Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA * noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    for (int Q = 0; Q < nQ_ref; Q++) {
         double value = 0.0;
         for (int i = 0; i < noccA; i++) {
              int ii = oo_idxAA->get(i,i);
              value += K->get(Q,ii);
         }
         Jc->set(Q, value);
    }
    K.reset();

    // beta
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB * noccB));
    K->read(psio_, PSIF_DFOCC_INTS);
    for (int Q = 0; Q < nQ_ref; Q++) {
         double value = 0.0;
         for (int i = 0; i < noccB; i++) {
              int ii = oo_idxBB->get(i,i);
              value += K->get(Q,ii);
         }
         Jc->add(Q, value);
    }
    K.reset();

    // Alpha part
    // F_ij = h_ij + \sum_{Q} b_ij^Q J^Q - \sum_{Q} \sum_{m} b_mi^Q b_mj^Q
    FooA->copy(HooA);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA * noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    FooA->gemv(true, K, Jc, 1.0, 1.0);
    FooA->contract(true, false, noccA, noccA, nQ_ref * noccA, K, K, -1.0, 1.0);

    // F_ia = h_ia + \sum_{Q} b_ia^Q J^Q - \sum_{Q} \sum_{m} b_mi^Q b_ma^Q
    FovA->copy(HovA);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA * nvirA));
    L->read(psio_, PSIF_DFOCC_INTS);
    FovA->gemv(true, L, Jc, 1.0, 1.0);
    FovA->contract(true, false, noccA, nvirA, nQ_ref * noccA, K, L, -1.0, 1.0);
    K.reset();

    // F_ai = F_ia
    FvoA = FovA->transpose();

    // F_ab = h_ab + \sum_{Q} b_ab^Q J^Q - \sum_{Q} \sum_{m} b_ma^Q b_mb^Q
    FvvA->copy(HvvA);
    FvvA->contract(true, false, nvirA, nvirA, nQ_ref * noccA, L, L, -1.0, 1.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    FvvA->gemv(true, K, Jc, 1.0, 1.0);
    K.reset();

    // Beta part
    // F_ij = h_ij + \sum_{Q} b_ij^Q J^Q - \sum_{Q} \sum_{m} b_mi^Q b_mj^Q
    FooB->copy(HooB);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB * noccB));
    K->read(psio_, PSIF_DFOCC_INTS);
    FooB->gemv(true, K, Jc, 1.0, 1.0);
    FooB->contract(true, false, noccB, noccB, nQ_ref * noccB, K, K, -1.0, 1.0);

    // F_ia = h_ia + \sum_{Q} b_ia^Q J^Q - \sum_{Q} \sum_{m} b_mi^Q b_ma^Q
    FovB->copy(HovB);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB * nvirB));
    L->read(psio_, PSIF_DFOCC_INTS);
    FovB->gemv(true, L, Jc, 1.0, 1.0);
    FovB->contract(true, false, noccB, nvirB, nQ_ref * noccB, K, L, -1.0, 1.0);
    K.reset();

    // F_ai = F_ia
    FvoB = FovB->transpose();

    // F_ab = h_ab + \sum_{Q} b_ab^Q J^Q - \sum_{Q} \sum_{m} b_ma^Q b_mb^Q
    FvvB->copy(HvvB);
    FvvB->contract(true, false, nvirB, nvirB, nQ_ref * noccB, L, L, -1.0, 1.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vv)", nQ_ref, nvirB, nvirB));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    FvvB->gemv(true, K, Jc, 1.0, 1.0);
    K.reset();

    // Set Alpha Fock Matrix
    FockA->set_oo(FooA);
    FockA->set_ov(FovA);
    FockA->set_vo(FvoA);
    FockA->set_vv(noccA, FvvA);

    // Set Beta Fock Matrix
    FockB->set_oo(FooB);
    FockB->set_ov(FovB);
    FockB->set_vo(FvoB);
    FockB->set_vv(noccB, FvvB);

    if (print_ > 2) {
        FockA->print();
        FockB->print();
    }

}// end else if (reference_ == "UNRESTRICTED")
    timer_off("Fock");
} // end fock

}} // End Namespaces
