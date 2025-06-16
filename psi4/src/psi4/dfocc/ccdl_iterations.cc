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
#include "psi4/libmints/matrix.h"
#include "psi4/libdiis/diismanager.h"

#include <cmath>

namespace psi {
namespace dfoccwave {

void DFOCC::ccdl_iterations() {
    //==========================================================================================
    //========================= Before iterations ==============================================
    //==========================================================================================
    SharedTensor2d U, T;

    // 3-index intermediates
    //timer_on("CCD 3-index intr");
    //ccd_3index_intr();
    //timer_off("CCD 3-index intr");

    //// F intermediates
    //timer_on("CCD F intr");
    //ccd_F_intr();
    //timer_off("CCD F intr");

    if (reference_ == "RESTRICTED") {
        // W intermediates
        ccdl_Wmnij();
        ccdl_Wmbej();
        ccdl_Wmbje();

        // Write
        // Form U_ij^ab
        U = std::make_shared<Tensor2d>("U2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        ccsd_u2_amps(U, t2);
        U->write_symm(psio_, PSIF_DFOCC_AMPS);
        U.reset();
        // Form T'(ib,ja) = t_ij^ab
        U = std::make_shared<Tensor2d>("T2p (IA|JB)", naoccA, navirA, naoccA, navirA);
        ccsd_t2_prime_amps(U, t2);
        U->write_symm(psio_, PSIF_DFOCC_AMPS);
        U.reset();

        // Malloc and Free
        //gQ = std::make_shared<Tensor1d>("CCDL G_Q", nQ);
        l2 = std::make_shared<Tensor2d>("L2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        l2->copy(t2);
        //if (orbs_already_opt == 1)
        //    l2->read_symm(psio_, PSIF_DFOCC_AMPS);
        //else l2->copy(t2);
        t2.reset();
    }// if (reference_ == "RESTRICTED")

    else if (reference_ == "UNRESTRICTED") {
        // W intermediates
        timer_on("CCSDL W intr");
        ccdl_WmnijAA();
        ccdl_WmnijBB();
        ccdl_WmnijAB();

        ccdl_WMBEJ_AAAA();
        ccdl_Wmbej_BBBB();
        ccdl_WMbEj_ABAB();
        ccdl_WmBeJ_BABA();
        ccdl_WMbeJ_ABBA();
        ccdl_WmBEj_BAAB();
        timer_off("CCSDL W intr");

        SharedTensor2d T2, L2;
        T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        L2 = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        L2->copy(T2);
        T2.reset();
        L2->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
        L2.reset();

        T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        L2 = std::make_shared<Tensor2d>("L2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        L2->copy(T2);
        T2.reset();
        L2->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
        L2.reset();

        T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T2->read(psio_, PSIF_DFOCC_AMPS);
        L2 = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        L2->copy(T2);
        T2.reset();
        L2->write(psio_, PSIF_DFOCC_AMPS);
    }// else if (reference_ == "UNRESTRICTED")

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
        else if (reference_ == "UNRESTRICTED") {
            Matrix L2AA("L2AA", ntri_anti_ijAA, ntri_anti_abAA);
            Matrix L2BB("L2BB", ntri_anti_ijBB, ntri_anti_abBB);
            Matrix L2AB("L2AB", naoccA * naoccB, navirA * navirB);

            ccsdlDiisManager = std::make_shared<DIISManager>(
                    cc_maxdiis_, "CCSDL DIIS L Amps", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::OnDisk);
            ccsdlDiisManager->set_error_vector_size(L2AA, L2BB, L2AB);
            ccsdlDiisManager->set_vector_size(L2AA, L2BB, L2AB);
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

        // print
        outfile->Printf(" %3d      %13.10f         %13.10f     %12.2e  \n", itr_occ, EcorrL, DE, rms_l2);

        if (itr_occ >= cc_maxiter) {
            conver = 0;  // means iterations were NOT converged
            break;
        }

        if (rms_t2 >= DIVERGE) {
            throw PSIEXCEPTION("CCD iterations are diverging");
        }

    } while (std::fabs(DE) >= tol_Eod || rms_l2 >= tol_t2);  // rms_l2 was rms_t2

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
        //outfile->Printf("\t3-index intermediates are done.\n");
        timer_off("CCD PDM 3-index intr");
    }

    if (reference_ == "RESTRICTED") {
        // free l2 amps
        U = std::make_shared<Tensor2d>("Ut2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        ccsd_u2_amps(U, l2);
        U->write_symm(psio_, PSIF_DFOCC_AMPS);
        U.reset();
        l2->write_symm(psio_, PSIF_DFOCC_AMPS);
        l2.reset();
    }

    // For GRAD
    if (dertype == "FIRST") {
        timer_on("CCD PDM yQia");
        ccd_pdm_yQia();
        timer_off("CCD PDM yQia");
    }

    // Mem free for DF ints
    // if (df_ints_incore) {
    reset_mo_df_ints();
    //}

}  // end ccdl_iterations

}  // namespace dfoccwave
}  // namespace psi
