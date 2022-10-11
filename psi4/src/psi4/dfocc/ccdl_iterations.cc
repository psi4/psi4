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
#include "psi4/libmints/matrix.h"
#include "psi4/libdiis/diismanager.h"

#include <cmath>

namespace psi {
namespace dfoccwave {

void DFOCC::ccdl_iterations() {
    //==========================================================================================
    //========================= Before iterations ==============================================
    //==========================================================================================
    // 3-index intermediates
    timer_on("CCD 3-index intr");
    ccd_3index_intr();
    timer_off("CCD 3-index intr");

    // F intermediates
    timer_on("CCD F intr");
    ccd_F_intr();
    timer_off("CCD F intr");

    // W intermediates
    ccdl_Wmnij();
    ccdl_Wmbej();
    ccdl_Wmbje();

    // Write
    // Form U_ij^ab
    SharedTensor2d U, T;
    U = std::make_shared<Tensor2d>("U2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    ccsd_u2_amps(U, t2);
    U->write_symm(psio_, PSIF_DFOCC_AMPS);
    U.reset();
    // Form T'(ib,ja) = t_ij^ab
    U = std::make_shared<Tensor2d>("T2p (IA|JB)", naoccA, navirA, naoccA, navirA);
    ccsd_t2_prime_amps(U, t2);
    U->write_symm(psio_, PSIF_DFOCC_AMPS);
    U.reset();
    // Form T(ij,ab)
    // T = SharedTensor2d(new Tensor2d("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    // T->sort(1324, t2, 1.0, 0.0);
    // T->write(psio_, PSIF_DFOCC_AMPS);
    // T.reset();

    // Malloc and Free
    gQ = std::make_shared<Tensor1d>("CCDL G_Q", nQ);
    l2 = std::make_shared<Tensor2d>("L2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    l2->copy(t2);
    t2.reset();

    //==========================================================================================
    //========================= Title ==========================================================
    //==========================================================================================
    outfile->Printf("\n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ================ Performing DF-CCDL iterations... ============================ \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf("\n");
    outfile->Printf("  Iter       E_corr                  DE                 L2 RMS      \n");
    outfile->Printf("  ----   ----------------      ----------------       ----------    \n");

    //==========================================================================================
    //========================= CCDL iterations ================================================
    //==========================================================================================
    itr_occ = 0;
    conver = 1;  // Assuming that the iterations will converge
    EccdL_old = Eccd;

    // DIIS
    if (do_diis_ == 1) {
        if (reference_ == "RESTRICTED") {
            Matrix L2("L2", naoccA * navirA, naoccA * navirA);
            ccsdlDiisManager = std::make_shared<DIISManager>(
                cc_maxdiis_, "CCDL DIIS L2 Amps", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::OnDisk);
            ccsdlDiisManager->set_error_vector_size(L2);
            ccsdlDiisManager->set_vector_size(L2);
        }
    }  // if diis true

    // head of loop
    do {
        // iterate
        itr_occ++;

        // 3-index intermediates
        timer_on("CCDL 3-index intr");
        ccdl_3index_intr();
        timer_off("CCDL 3-index intr");

        // L2 amplitudes
        timer_on("L2 AMPS");
        ccdl_l2_amps();
        timer_off("L2 AMPS");

        DE = EccdL - EccdL_old;
        EccdL_old = EccdL;

        // RMS
        if (reference_ == "UNRESTRICTED") {
            rms_t2 = MAX0(rms_t2AA, rms_t2BB);
            rms_t2 = MAX0(rms_t2, rms_t2AB);
        }

        // print
        outfile->Printf(" %3d      %13.10f         %13.10f     %12.2e  \n", itr_occ, EcorrL, DE, rms_t2);

        if (itr_occ >= cc_maxiter) {
            conver = 0;  // means iterations were NOT converged
            break;
        }

        if (rms_t2 >= DIVERGE) {
            throw PSIEXCEPTION("CCD iterations are diverging");
        }

    } while (std::fabs(DE) >= tol_Eod || rms_t2 >= tol_t2);

    if (conver == 1) {
        outfile->Printf("\n");
        outfile->Printf(" ============================================================================== \n");
        outfile->Printf(" ===================== DF-CCDL ITERATIONS ARE CONVERGED ======================= \n");
        outfile->Printf(" ============================================================================== \n");

    }

    else if (conver == 0) {
        outfile->Printf("\n ====================== DF-CCDL IS NOT CONVERGED IN %2d ITERATIONS ============ \n",
                        cc_maxiter);
        throw PSIEXCEPTION("DF-CCDL iterations did not converge");
    }

    // delete
    if (do_diis_ == 1) ccsdlDiisManager->delete_diis_file();

    // For GRAD
    if (dertype == "FIRST") {
        outfile->Printf("\n\tComputing 3-index intermediates...\n");
        timer_on("CCD PDM 3-index intr");
        ccd_pdm_3index_intr();
        timer_off("CCD PDM 3-index intr");
    }

    // free l2 amps
    U = std::make_shared<Tensor2d>("Ut2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    ccsd_u2_amps(U, l2);
    U->write_symm(psio_, PSIF_DFOCC_AMPS);
    U.reset();
    l2->write_symm(psio_, PSIF_DFOCC_AMPS);
    l2.reset();

    // For GRAD
    if (dertype == "FIRST") {
        timer_on("CCD PDM yQia");
        ccd_pdm_yQia();
        timer_off("CCD PDM yQia");
    }

    // Mem free for DF ints
    // if (df_ints_incore) {
    bQijA.reset();
    bQiaA.reset();
    bQabA.reset();
    //}

}  // end ccdl_iterations
}  // namespace dfoccwave
}  // namespace psi
