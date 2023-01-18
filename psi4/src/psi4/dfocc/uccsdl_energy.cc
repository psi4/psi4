/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

namespace psi {
namespace dfoccwave {

void DFOCC::uccsdl_energy()
{
    SharedTensor2d X, Y, K, L, T, T2, L2, L2new, G, J, W, U, Tau;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //========================================================================================================//
    //                                               L2AB                                                     //
    //========================================================================================================//
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    L2 = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    L2->read(psio_, PSIF_DFOCC_AMPS);
    L2new = std::make_shared<Tensor2d>("New L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    L2new->read(psio_, PSIF_DFOCC_AMPS);

    // reset L2->copy(L2new)
    rms_t2AB = L2new->rms(L2);
    L2->copy(L2new);
    L2new.reset();
    L2->write(psio_, PSIF_DFOCC_AMPS);

    // EccsdLAB Energy
    Tau = std::make_shared<Tensor2d>("Tau <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    uccsd_tau_amps_OS(naoccA, naoccB, navirA, navirB, Tau, L2, l1A, l1B);
//L2->print();
    L2.reset();

    J = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|jb)", naoccA, navirA, naoccB, navirB);
    J->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    L = std::make_shared<Tensor2d>("MO Int <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    L->sort(1324, J, 1.0, 0.0);
    J.reset();
    EabL =  Tau->vector_dot(L);
    L.reset();
    Tau.reset();
    //outfile->Printf("\tAlpha-beta contribution to EcorrL   : %20.14f\n", EabL);


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //========================================================================================================//
    //                                               L2AA                                                     //
    //========================================================================================================//
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    L2 = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    L2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    L2new = std::make_shared<Tensor2d>("New L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    L2new->read_anti_symm(psio_, PSIF_DFOCC_AMPS);

    // reset L2->copy(L2new)
    rms_t1A = l1newA->rms(l1A);
    l1A->copy(l1newA);
    rms_t2AA = L2new->rms(L2);
    L2->copy(L2new);
    L2new.reset();
    L2->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    L2.reset();

    // EccsdAA Energy
    L2 = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    L2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA);
    uccsd_tau_amps(naoccA, naoccA, navirA, navirA, Tau, L2, l1A, l1A);
//L2->print();
    L2.reset();

    J = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA);
    J->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    G = std::make_shared<Tensor2d>("G <IJ||AB>", naoccA, naoccA, navirA, navirA);
    G->sort(1324, J, 1.0, 0.0);
    G->sort(1342, J, -1.0, 1.0);
    J.reset();
    EaaL = 0.25 * Tau->vector_dot(G);
    G.reset();
    Tau.reset();
    //outfile->Printf("\tAlpha-alpha contribution to EcorrL  : %20.14f\n", EaaL);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //========================================================================================================//
    //                                               L2BB                                                     //
    //========================================================================================================//
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    L2 = std::make_shared<Tensor2d>("L2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    L2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    L2new = std::make_shared<Tensor2d>("New L2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    L2new->read_anti_symm(psio_, PSIF_DFOCC_AMPS);

    // reset T2->copy(T2new)
    rms_t1B = l1newB->rms(l1B);
    l1B->copy(l1newB);
    rms_t2BB = L2new->rms(L2);
    L2->copy(L2new);
    L2new.reset();
    L2->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    L2.reset();

    // EccsdBB Energy
    L2 = std::make_shared<Tensor2d>("L2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    L2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <ij|ab>", naoccB, naoccB, navirB, navirB);
    uccsd_tau_amps(naoccB, naoccB, navirB, navirB, Tau, L2, l1B, l1B);
//L2->print();
    L2.reset();

    J = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB);
    J->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    G = std::make_shared<Tensor2d>("G <ij||ab>", naoccB, naoccB, navirB, navirB);
    G->sort(1324, J, 1.0, 0.0);
    G->sort(1342, J, -1.0, 1.0);
    J.reset();
    EbbL = 0.25 * Tau->vector_dot(G);
    G.reset();
    Tau.reset();
    //outfile->Printf("\tBeta-beta contribution to EcorrL    : %20.14f\n", EbbL);

    //double El1A =  l1A->vector_dot(FiaA);
    //double El1B =  l1B->vector_dot(FiaB);
    //outfile->Printf("\tL1A contribution to EcorrL          : %20.14f\n", El1A);
    //outfile->Printf("\tL1B contribution to EcorrL          : %20.14f\n", El1B);
    //EcorrL = EaaL + EbbL + EabL + El1A + El1B;
    EcorrL = EaaL + EbbL + EabL;
    EccsdL = Escf + EcorrL;

    // combined rms
    double rms_ss = MAX0(rms_t2AA, rms_t2BB);
    rms_t2 = MAX0(rms_ss, rms_t2AB);
    rms_t1 = MAX0(rms_t1A, rms_t1B);


}// end ccsd_energy

} // dfoccwave
} // End Namespaces

