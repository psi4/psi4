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

#include "psi4/libqt/qt.h"
#include "defines.h"
#include "dfocc.h"

using namespace psi;

namespace psi {
namespace dfoccwave {

void DFOCC::ccsd_t1_amps() {
    // defs
    SharedTensor2d K, T1, T, U, Tau;

    // t_i^a <= \sum_{e} t_i^e Fae
    t1newA->gemm(false, true, t1A, FabA, 1.0, 0.0);

    // t_i^a <= -\sum_{m} t_m^a Fmi
    t1newA->gemm(true, false, FijA, t1A, -1.0, 1.0);

    // t_i^a <= \sum_{m,e} u_im^ae Fme
    U = std::make_shared<Tensor2d>("U2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    ccsd_u2_amps(U, t2);
    t1newA->gemv(false, U, FiaA, 1.0, 1.0);
    U.reset();

    // t_i^a <= \sum_{Q} t_Q b_ia^Q
    t1newA->gemv(true, bQiaA, T1c, 1.0, 1.0);

    // t_i^a <= -\sum_{Q,m} (T_ma^Q + t_ma^Q) b_mi^Q -= \sum_{Q,m} B(Qm,i) X(Qm,a)
    T = std::make_shared<Tensor2d>("T2 (Q|IA)", nQ, naoccA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    T1 = std::make_shared<Tensor2d>("T1 (Q|IA)", nQ, naoccA, navirA);
    T1->read(psio_, PSIF_DFOCC_AMPS);
    U = std::make_shared<Tensor2d>("T1+T2 (Q|IA)", nQ, naoccA, navirA);
    U->copy(T);
    T.reset();
    U->add(T1);
    T1.reset();
    t1newA->contract(true, false, naoccA, navirA, nQ * naoccA, bQijA, U, -1.0, 1.0);
    U.reset();

    // t_i^a <= \sum_{Q,e} T_ie^Q b_ea^Q = \sum_{Qe} T^Q(i,e) B^Q(e,a)
    T = std::make_shared<Tensor2d>("T2 (Q|IA)", nQ, naoccA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Temp (Q|AI)", nQ, navirA, naoccA);
    Tau->swap_3index_col(T);
    T.reset();
    t1newA->contract(true, false, naoccA, navirA, nQ * navirA, Tau, bQabA, 1.0, 1.0);
    Tau.reset();

    // Denom
    for (int i = 0; i < naoccA; ++i) {
        for (int a = 0; a < navirA; ++a) {
            double value = FockA->get(i + nfrzc, i + nfrzc) - FockA->get(a + noccA, a + noccA);
            t1newA->set(i, a, t1newA->get(i, a) / value);
        }
    }
    // t1newA->print();

}  // end ccsd_t1_amps
}  // namespace dfoccwave
}  // namespace psi
