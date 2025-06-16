/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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

void DFOCC::ccd_3index_intr() {

    // RHF
    if (reference_ == "RESTRICTED") {
        // defs
        SharedTensor2d T, U;

        // T(Q,ia) = \sum_{jb} b_jb^Q u_ij^ab = \sum_{jb} b(Q,jb) U(jb,ia)
        U = std::make_shared<Tensor2d>("U2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        ccsd_u2_amps(U, t2);
        T = std::make_shared<Tensor2d>("T2 (Q|IA)", nQ, naoccA, navirA);
        T->gemm(false, false, bQiaA, U, 1.0, 0.0);
        U.reset();
        T->write(psio_, PSIF_DFOCC_AMPS);
    }// if (reference_ == "RESTRICTED")

    // UHF
    else if (reference_ == "UNRESTRICTED") {
        //std::cout << "ccd_3index_intr is starting... \n";
        SharedTensor2d TQiaA, TQiaB, TautAA, TautBB, TautAB, TautQA, TautQB;
        SharedTensor2d tQovA, tQovB, tQooA, tQooB, tQvoA, tQvoB, tQvvA, tQvvB, ttQA, ttQB;
        SharedTensor2d L, M, N, X, Y, T, K, U, Tau, T2;

        // T(Q,IA) = \sum(J,B) (T(IJ,AB) * b(Q,JB)) + \sum(j,b) (T2AB(Ij,Ab) * b(Q,jb)) (45)
        // T(Q,IA) += \sum(J,B) (T(IJ,AB) * b(Q,JB))
        L = std::make_shared<Tensor2d>("L (IA|JB)", naoccA, navirA, naoccA, navirA);
        T = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        T->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        L->sort(1324, T, 1.0, 0.0);
        T.reset();
        TQiaA = std::make_shared<Tensor2d>("T2 (Q|IA)", nQ, naoccA, navirA);
        TQiaA->gemm(false, true, bQiaA, L, 1.0, 0.0);
        L.reset();

        // T(Q,IA) += \sum(j,b) (T2AB(Ij,Ab) * b(Q,jb))
        M = std::make_shared<Tensor2d>("M (IA|jb)", naoccA, navirA, naoccB, navirB);
        T = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T->read(psio_, PSIF_DFOCC_AMPS);
        M->sort(1324, T, 1.0, 0.0);
        T.reset();
        TQiaA->gemm(false, true, bQiaB, M, 1.0, 1.0);
        M.reset();
        TQiaA->write(psio_, PSIF_DFOCC_AMPS);
        TQiaA.reset();

        // T(Q,ia) = \sum(j,b) (T2BB(ij,ab) * b(Q,jb)) + \sum(J,B) (T2AB(Ji,Ba) * b(Q,JB)) (46)
        // T(Q,ia) +=  \sum(j,b) (T2BB(ij,ab) * b(Q,jb))
        N = std::make_shared<Tensor2d>("N (ia|jb)", naoccB, navirB, naoccB, navirB);
        T = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        T->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        N->sort(1324, T, 1.0, 0.0);
        T.reset();
        TQiaB = std::make_shared<Tensor2d>("T2 (Q|ia)", nQ, naoccB, navirB);
        TQiaB->gemm(false, true, bQiaB, N, 1.0, 0.0);
        N.reset();
        // T(Q,ia) +=  \sum(J,B) (T2AB(Ji,Ba) * b(Q,JB))
        M = std::make_shared<Tensor2d>("M (IA|jb)", naoccA, navirA, naoccB, navirB);
        T = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T->read(psio_, PSIF_DFOCC_AMPS);
        M->sort(1324, T, 1.0, 0.0);
        T.reset();
        TQiaB->gemm(false, false, bQiaA, M, 1.0, 1.0);
        M.reset();
        TQiaB->write(psio_, PSIF_DFOCC_AMPS);
    }// else if (reference_ == "UNRESTRICTED")

    // outfile->Printf("\t3indices done.\n");

}  // end ccd_3index_intr
}  // namespace dfoccwave
}  // namespace psi
