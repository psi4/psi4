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
#include "psi4/libdiis/diismanager.h"
#include "psi4/libmints/matrix.h"

namespace psi {
namespace dfoccwave {

void DFOCC::ccd_t2_amps() {

    // RHF
    if (reference_ == "RESTRICTED") {
        // defs
        SharedTensor2d K, I, T, Tnew, U, Tau, W, X, Y;
        SharedMatrix T2, RT2;

        // t_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
        // X(ia,jb) = \sum_{e} t_ij^ae F_be = \sum_{e} T(ia,je) F_be
        X = std::make_shared<Tensor2d>("X (IA|JB)", naoccA, navirA, naoccA, navirA);
        X->contract(false, true, naoccA * navirA * naoccA, navirA, navirA, t2, FabA, 1.0, 0.0);
        // X->cont424("IAJB", "IAJE", "BE", false, t2, FabA, 1.0, 0.0); // it works

        // t_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
        // X(ia,jb) = -\sum_{m} t_mj^ab F_mi = -\sum_{m} F(m,i) T(ma,jb)
        X->contract(true, false, naoccA, naoccA * navirA * navirA, naoccA, FijA, t2, -1.0, 1.0);
        // X->cont244("IAJB", "MI", "MAJB", false, FijA, t2, -1.0, 1.0); // it works
        X->symmetrize();

        // t_ij^ab <= <ij|ab>
        Tnew = std::make_shared<Tensor2d>("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        Tnew->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);

        // Contributions of X
        Tnew->axpy(X, 2.0);
        X.reset();

        // Write and close
        Tnew->write_symm(psio_, PSIF_DFOCC_AMPS);
        Tnew.reset();

        // WmnijT2
        ccd_WmnijT2();

        // WmbejT2
        ccd_WmbejT2();

        // WabefT2
        if (Wabef_type_ == "AUTO") {
            if (!do_ppl_hm)
                ccd_WabefT2();
            else
                ccd_WabefT2_high_mem();
        } else if (Wabef_type_ == "LOW_MEM")
            ccd_WabefT2();
        else if (Wabef_type_ == "HIGH_MEM")
            ccd_WabefT2_high_mem();

        // Denom
        Tnew = std::make_shared<Tensor2d>("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
        Tnew->apply_denom_chem(nfrzc, noccA, FockA);

        // Reset T2
        rms_t2 = Tnew->rms(t2);

        // Error vector
        Tau = std::make_shared<Tensor2d>("RT2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        Tau->copy(Tnew);
        Tau->subtract(t2);
        t2->copy(Tnew);
        Tnew.reset();

        // DIIS
        RT2 = std::make_shared<Matrix>("RT2", naoccA * navirA, naoccA * navirA);
        Tau->to_matrix(RT2);
        Tau.reset();
        T2 = std::make_shared<Matrix>("T2", naoccA * navirA, naoccA * navirA);
        t2->to_matrix(T2);

        // add entry
        if (do_diis_ == 1 && orb_opt_ == "FALSE") ccsdDiisManager->add_entry(RT2.get(), T2.get());
        RT2.reset();

        // extrapolate
        if (do_diis_ == 1 && orb_opt_ == "FALSE") {
            if (ccsdDiisManager->subspace_size() >= cc_mindiis_) ccsdDiisManager->extrapolate(T2.get());
            t2->set2(T2);
        }
        T2.reset();

        // Form U(ia,jb) = 2*T(ia,jb) - T (ib,ja)
        U = std::make_shared<Tensor2d>("U2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        ccsd_u2_amps(U, t2);

        // Energy
        K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA);
        K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
        Ecorr = U->vector_dot(K);
        U.reset();
        K.reset();
        Eccd = Eref + Ecorr;
    }  // if (reference_ == "RESTRICTED")

}  // end ccd_t2_amps

}  // namespace dfoccwave
}  // namespace psi
