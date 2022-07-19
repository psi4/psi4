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

#include "dfocc.h"

namespace psi{
namespace dfoccwave {

void DFOCC::uccsd_energy()
{
    SharedTensor2d X, Y, K, L, T, T2, T2new, G, J, W, U, Tau;

//T2AB   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    T2new = std::make_shared<Tensor2d>("New T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2new->read(psio_, PSIF_DFOCC_AMPS);

    // reset T2->copy(T2new)
    rms_t2AB = T2new->rms(T2);
    T2->copy(T2new);
    T2new.reset();
    T2->write(psio_, PSIF_DFOCC_AMPS);
    T2.reset();

    // EccsdAB Energy
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    uccsd_tau_amps_OS(naoccA, naoccB, navirA, navirB, Tau, T2, t1A, t1B);
    T2.reset();

    J = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|jb)", naoccA, navirA, naoccB, navirB);
    J->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    L = std::make_shared<Tensor2d>("L <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    L->sort(1324, J, 1.0, 0.0);
    J.reset();
    Eab = Tau->vector_dot(L);
    L.reset();
    Tau.reset();
    //outfile->Printf("\tAlpha-beta contribution to Ecorr   : %20.14f\n", Eab);

//T2AA   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T2new = std::make_shared<Tensor2d>("New T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2new->read_anti_symm(psio_, PSIF_DFOCC_AMPS);

    // reset T2->copy(T2new)
    rms_t1A = t1newA->rms(t1A);
    t1A->copy(t1newA);
    rms_t2AA = T2new->rms(T2);
    T2->copy(T2new);
    T2new.reset();
    T2->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T2.reset();

    // EccsdAA Energy
    // Tau_AA
    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA);
    uccsd_tau_amps(naoccA, naoccA, navirA, navirA, Tau, T2, t1A, t1A);
    T2.reset();
    J = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA);
    J->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    G = std::make_shared<Tensor2d>("G <IJ||AB>", naoccA, naoccA, navirA, navirA);
    G->sort(1324, J, 1.0, 0.0);
    G->sort(1342, J, -1.0, 1.0);
    J.reset();
    Eaa = 0.25 * Tau->vector_dot(G);
    G.reset();
    Tau.reset();
    //outfile->Printf("\tAlpha-alpha contribution to Ecorr  : %20.14f\n", Eaa);

//T2BB   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T2new = std::make_shared<Tensor2d>("New T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2new->read_anti_symm(psio_, PSIF_DFOCC_AMPS);

    // reset T2->copy(T2new)
    rms_t1B = t1newB->rms(t1B);
    t1B->copy(t1newB);
    rms_t2BB = T2new->rms(T2);
    T2->copy(T2new);
    T2new.reset();
    T2->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T2.reset();

    // EccsdBB Energy
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <ij|ab>", naoccB, naoccB, navirB, navirB);
    uccsd_tau_amps(naoccB, naoccB, navirB, navirB, Tau, T2, t1B, t1B);
    T2.reset();
    J = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB);
    J->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    G = std::make_shared<Tensor2d>("G <ij||ab>", naoccB, naoccB, navirB, navirB);
    G->sort(1324, J, 1.0, 0.0);
    G->sort(1342, J, -1.0, 1.0);
    J.reset();
    Ebb = 0.25 * Tau->vector_dot(G);
    G.reset();
    Tau.reset();
    //outfile->Printf("\tBeta-beta contribution to Ecorr    : %20.14f\n", Ebb);
    Ecorr = Eaa + Ebb + Eab;
    Eccsd = Escf + Ecorr;

    double rms_ss = MAX0(rms_t2AA, rms_t2BB);
    rms_t2 = MAX0(rms_ss, rms_t2AB);
    rms_t1 = MAX0(rms_t1A, rms_t1B);

}// end ccsd_energy

}  // namespace dfoccwave
}  // namespace psi

