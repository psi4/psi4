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

#include <cmath>

namespace psi {
namespace dfoccwave {

void DFOCC::ccsdl_iterations() {
    //==========================================================================================
    //========================= Before iterations ==============================================
    //==========================================================================================

    SharedTensor2d U, T;
    if (reference_ == "RESTRICTED") {
        // 3-index intermediates
        //timer_on("CCSD 3-index intr");
        //ccsd_3index_intr();
        //timer_off("CCSD 3-index intr");

        // F intermediates
        //timer_on("CCSD F intr");
        //ccsd_F_intr();
        //timer_off("CCSD F intr");

        // W intermediates
        ccsdl_Wmnij();
        ccsdl_Wmbej();
        ccsdl_Wmbje();
        ccsdl_Wmbij();
        // ccsdl_Wmnie();

        // Write
        // Form U_ij^ab
        U = std::make_shared<Tensor2d>("U2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        ccsd_u2_amps(U, t2);
        U->write_symm(psio_, PSIF_DFOCC_AMPS);
        U.reset();
        // Form Tau_ij^ab
        U = std::make_shared<Tensor2d>("Tau (IA|JB)", naoccA, navirA, naoccA, navirA);
        ccsd_tau_amps(U, t2);
        U->write_symm(psio_, PSIF_DFOCC_AMPS);
        T = std::make_shared<Tensor2d>("2*Tau(ia,jb) - Tau(ib,ja)", naoccA, navirA, naoccA, navirA);
        ccsd_u2_amps(T, U);
        U.reset();
        T->write_symm(psio_, PSIF_DFOCC_AMPS);
        T.reset();
        // Form T'(ib,ja) = t_ij^ab
        U = std::make_shared<Tensor2d>("T2p (IA|JB)", naoccA, navirA, naoccA, navirA);
        ccsd_t2_prime_amps(U, t2);
        U->write_symm(psio_, PSIF_DFOCC_AMPS);
        U.reset();
        // Form T(ij,ab)
        // T = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        // T->sort(1324, t2, 1.0, 0.0);
        // T->write(psio_, PSIF_DFOCC_AMPS);
        // T.reset();
        // t1A->print();

        // Malloc and Free
        gQ = std::make_shared<Tensor1d>("CCSDL G_Q", nQ);
        gQp = std::make_shared<Tensor1d>("CCSDL G_Qp", nQ);
        l1A = std::make_shared<Tensor2d>("L1 <I|A>", naoccA, navirA);
        l1newA = std::make_shared<Tensor2d>("New L1 <I|A>", naoccA, navirA);
        l2 = std::make_shared<Tensor2d>("L2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        l2->copy(t2);
        t2.reset();
        l1A->copy(t1A);

    } // if restricted
    else if (reference_ == "UNRESTRICTED") {
        // W intermediates
        timer_on("CCSDL W intr");
        uccsdl_ZMBEJ_AAAA();
        uccsdl_Zmbej_BBBB();
        uccsdl_ZMbEj_ABAB();
        uccsdl_ZmBeJ_BABA();
        uccsdl_ZMbeJ_ABBA();
        uccsdl_ZmBEj_BAAB();

        uccsdl_WMBEJ_AAAA();
        uccsdl_Wmbej_BBBB();
        uccsdl_WMbEj_ABAB();
        uccsdl_WmBeJ_BABA();
        uccsdl_WMbeJ_ABBA();
        uccsdl_WmBEj_BAAB();

        uccsdl_WMNIE_AAAA();
        uccsdl_Wmnie_BBBB();
        uccsdl_WMnIe_ABAB();
        uccsdl_WmNiE_BABA();

        uccsdl_WMNIJ_AAAA();
        uccsdl_Wmnij_BBBB();
        uccsdl_WMnIj_ABAB();

        uccsdl_WMBIJ_AAAA();
        uccsdl_Wmbij_BBBB();
        uccsdl_WMbIj_ABAB();
        uccsdl_WmBiJ_BABA();
        timer_off("CCSDL W intr");

        // Global Tensors
        GijA = std::make_shared<Tensor2d>("G Intermediate <I|J>", naoccA, naoccA);
        GijB = std::make_shared<Tensor2d>("G Intermediate <i|j>", naoccB, naoccB);
        GabA = std::make_shared<Tensor2d>("G Intermediate <A|B>", navirA, navirA);
        GabB = std::make_shared<Tensor2d>("G Intermediate <a|b>", navirB, navirB);

        gQ = std::make_shared<Tensor1d>("CCSDL G_Q", nQ);
        gQp = std::make_shared<Tensor1d>("CCSDL G_Qp", nQ);
        l1A = std::make_shared<Tensor2d>("L1 <I|A>", naoccA, navirA);
        l1B = std::make_shared<Tensor2d>("L1 <i|a>", naoccB, navirB);
        l1newA = std::make_shared<Tensor2d>("New L1 <I|A>", naoccA, navirA);
        l1newB = std::make_shared<Tensor2d>("New L1 <i|a>", naoccB, navirB);
        l1A->copy(t1A);
        l1B->copy(t1B);

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
        L2.reset();
        //=== END DFUCCSD ===
    } // else if unrestricted

    //==========================================================================================
    //========================= Title ==========================================================
    //==========================================================================================
    outfile->Printf("\n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ================ Performing DF-CCSDL iterations... =========================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf("\n");
    outfile->Printf("  Iter       E_corr                  DE                 L2 RMS        L1 RMS     \n");
    outfile->Printf("  ----   ----------------      ----------------       ----------    ----------   \n");

    //==========================================================================================
    //========================= CCSDL iterations ===============================================
    //==========================================================================================
    itr_occ = 0;
    conver = 1;  // Assuming that the iterations will converge
    EccsdL_old = Eccsd;

    // DIIS
    if (do_diis_ == 1) {
        if (reference_ == "RESTRICTED") {
            std::shared_ptr<Matrix> L2(new Matrix("L2", naoccA * navirA, naoccA * navirA));
            std::shared_ptr<Matrix> L1(new Matrix("L1", naoccA, navirA));
            ccsdlDiisManager = std::shared_ptr<DIISManager>(
                new DIISManager(cc_maxdiis_, "CCSDL DIIS L Amps", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::OnDisk));
            ccsdlDiisManager->set_error_vector_size(L2.get(), L1.get());
            ccsdlDiisManager->set_vector_size(L2.get(), L1.get());
            L2.reset();
            L1.reset();
        }
        else if (reference_ == "UNRESTRICTED") {
            //=== BEGIN DFUCCSD ===
            std::shared_ptr<Matrix> L2AA(new Matrix("L2AA", ntri_anti_ijAA, ntri_anti_abAA));
            std::shared_ptr<Matrix> L2BB(new Matrix("L2BB", ntri_anti_ijBB, ntri_anti_abBB));
            std::shared_ptr<Matrix> L2AB(new Matrix("L2AB", naoccA * naoccB, navirA * navirB));
            std::shared_ptr<Matrix> L1A(new Matrix("L1A", naoccA, navirA));
            std::shared_ptr<Matrix> L1B(new Matrix("L1B", naoccB, navirB));

            ccsdlDiisManager = std::shared_ptr<DIISManager>(
                    new DIISManager(cc_maxdiis_, "CCSDL DIIS L Amps", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::OnDisk));
            ccsdlDiisManager->set_error_vector_size(L2AA.get(), L2BB.get(), L2AB.get(), L1A.get(), L1B.get());
            ccsdlDiisManager->set_vector_size(L2AA.get(), L2BB.get(), L2AB.get(), L1A.get(), L1B.get());
            L2AA.reset();
            L2BB.reset();
            L2AB.reset();
            L1A.reset();
            L1B.reset();
            //=== END DFUCCSD ===
        }
    }  // if diis true

    // head of loop
    do {
        // iterate
        itr_occ++;

        // 3-index intermediates
        timer_on("CCSDL 3-index intr");
        ccsdl_3index_intr();
        timer_off("CCSDL 3-index intr");

        // L1 amplitudes
        timer_on("L1 AMPS");
        ccsdl_l1_amps();
        timer_off("L1 AMPS");

        // L2 amplitudes
        timer_on("L2 AMPS");
        ccsdl_l2_amps();
        timer_off("L2 AMPS");

        DE = EccsdL - EccsdL_old;
        EccsdL_old = EccsdL;

        // print
        outfile->Printf(" %3d      %13.10f         %13.10f     %12.2e  %12.2e \n", itr_occ, EcorrL, DE, rms_t2, rms_t1);

        if (itr_occ >= cc_maxiter) {
            conver = 0;  // means iterations were NOT converged
            break;
        }

        if (rms_t2 >= DIVERGE || rms_t1 >= DIVERGE) {
            throw PSIEXCEPTION("CCSD iterations are diverging");
        }

    } while (std::fabs(DE) >= tol_Eod || rms_t2 >= tol_t2 || rms_t1 >= tol_t2);

    if (conver == 1) {
        outfile->Printf("\n");
        outfile->Printf(" ============================================================================== \n");
        outfile->Printf(" ===================== DF-CCSDL ITERATIONS ARE CONVERGED ====================== \n");
        outfile->Printf(" ============================================================================== \n");
    }

    else if (conver == 0) {
        outfile->Printf("\n ====================== DF-CCSDL IS NOT CONVERGED IN %2d ITERATIONS =========== \n",
                        cc_maxiter);
        throw PSIEXCEPTION("DF-CCSDL iterations did not converge");
    }

    // delete
    if (do_diis_ == 1) ccsdlDiisManager->delete_diis_file();

    // For GRAD
    if (reference_ == "RESTRICTED") {
        // free l2 amps
        U = std::make_shared<Tensor2d>("Ut2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        ccsd_u2_amps(U, l2);
        U->write_symm(psio_, PSIF_DFOCC_AMPS);
        U.reset();
        l2->write_symm(psio_, PSIF_DFOCC_AMPS);

        // For GRAD
        if (dertype == "FIRST" && wfn_type_ != "DF-CCSD(AT)") {
            outfile->Printf("\n\tComputing 3-index intermediates...\n");
            timer_on("CCSD PDM 3-index intr");
            ccsd_pdm_3index_intr();
            timer_off("CCSD PDM 3-index intr");

            timer_on("CCSD PDM yQia");
            ccsd_pdm_yQia();
            timer_off("CCSD PDM yQia");
        }
        l2.reset();

        // Mem free for DF ints
        // if (df_ints_incore) {
        bQijA.reset();
        bQiaA.reset();
        bQabA.reset();
    } // if restricted
    else if (reference_ == "UNRESTRICTED") {

        // Global Tensors
        GtijA = std::make_shared<Tensor2d>("Gtilde Intermediate <I|J>", naoccA, naoccA);
        GtijB = std::make_shared<Tensor2d>("Gtilde Intermediate <i|j>", naoccB, naoccB);
        GtabA = std::make_shared<Tensor2d>("Gtilde Intermediate <A|B>", navirA, navirA);
        GtabB = std::make_shared<Tensor2d>("Gtilde Intermediate <a|b>", navirB, navirB);

        L1c = std::make_shared<Tensor1d>("DF_BASIS_CC L1_Q", nQ);

        // For GRAD
        if (dertype == "FIRST" && wfn_type_ != "DF-CCSD(AT)") {
            uccsd_pdm_3index_intr();
            uccsd_pdm_yQia();

        }
    } // else if unrestricted

}  // end ccsdl_iterations
}  // namespace dfoccwave
}  // namespace psi
