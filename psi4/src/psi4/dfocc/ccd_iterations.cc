/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
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
#include "psi4/libdiis/diismanager.h"
#include "dfocc.h"
#include "psi4/libmints/matrix.h"

#include <cmath>

namespace psi {
namespace dfoccwave {

void DFOCC::ccd_iterations() {
    outfile->Printf("\n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ================ Performing DF-CCD iterations... ============================= \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf("\n");
    outfile->Printf("  Iter       E_corr                  DE                 T2 RMS      \n");
    outfile->Printf("  ----   ----------------      ----------------       ----------    \n");

    //==========================================================================================
    //========================= CCD iterations =================================================
    //==========================================================================================
    itr_occ = 0;
    conver = 1;  // Assuming that the iterations will converge
    Eccd_old = Eccd;
    //outfile->Printf("Initial Eccd: %13.10f \n", Eccd);

    // DIIS
    if (do_diis_ == 1) {
        if (reference_ == "RESTRICTED") {
            std::shared_ptr<Matrix> T2(new Matrix("T2", naoccA * navirA, naoccA * navirA));
            if (reference_ == "RESTRICTED") {
                ccsdDiisManager = std::shared_ptr<DIISManager>(
                    new DIISManager(cc_maxdiis_, "CCSD DIIS T Amps", DIISManager::LargestError, DIISManager::OnDisk));
                ccsdDiisManager->set_error_vector_size(1, DIISEntry::Matrix, T2.get());
                ccsdDiisManager->set_vector_size(1, DIISEntry::Matrix, T2.get());
            }
            T2.reset();
        }
        else if (reference_ == "UNRESTRICTED") {
            std::shared_ptr<Matrix> T2AA(new Matrix("T2AA", ntri_anti_ijAA, ntri_anti_abAA));
            std::shared_ptr<Matrix> T2BB(new Matrix("T2BB", ntri_anti_ijBB, ntri_anti_abBB));
            std::shared_ptr<Matrix> T2AB(new Matrix("T2AB", naoccA * naoccB, navirA * navirB));

            ccsdDiisManager = std::shared_ptr<DIISManager>(
                new DIISManager(cc_maxdiis_, "CCSD DIIS T Amps", DIISManager::LargestError, DIISManager::OnDisk));
            ccsdDiisManager->set_error_vector_size(3, DIISEntry::Matrix, T2AA.get(), DIISEntry::Matrix, T2BB.get(),
                                                   DIISEntry::Matrix, T2AB.get());
            ccsdDiisManager->set_vector_size(3, DIISEntry::Matrix, T2AA.get(), DIISEntry::Matrix, T2BB.get(),
                                             DIISEntry::Matrix, T2AB.get());
            T2AA.reset();
            T2BB.reset();
            T2AB.reset();
        }
    }  // if diis true

    // head of loop
    do {
        // iterate
        itr_occ++;

        // 3-index intermediates
        timer_on("CCSD 3-index intr");
        ccd_3index_intr();
        timer_off("CCSD 3-index intr");

        // F intermediates
        timer_on("CCSD F intr");
        ccd_F_intr();
        timer_off("CCSD F intr");

        // T2 amplitudes
        timer_on("T2 AMPS");
        ccd_t2_amps();
        timer_off("T2 AMPS");
        //std::cout << "I am here \n";

        DE = Eccd - Eccd_old;
        Eccd_old = Eccd;

        // print
        outfile->Printf(" %3d      %13.10f         %13.10f     %12.2e  \n", itr_occ, Ecorr, DE, rms_t2);

        if (itr_occ >= cc_maxiter) {
            conver = 0;  // means iterations were NOT converged
            break;
        }

        if (rms_t2 >= DIVERGE) {
            throw PSIEXCEPTION("CCD iterations are diverging");
        }

    } while (std::fabs(DE) >= tol_Eod || rms_t2 >= tol_t2);

    // delete
    if (do_diis_ == 1) ccsdDiisManager->delete_diis_file();

    // Mem alloc for DF ints
    //if (df_ints_incore) {
        if (cc_lambda_ == "FALSE") {
            //if (reference_ == "RESTRICTED") reset_mo_df_ints();
            reset_mo_df_ints();
        }
    //}

    // free t2 amps
    if (t2_incore && reference_ == "RESTRICTED") {
        if (cc_lambda_ == "TRUE") {
            t2->write_symm(psio_, PSIF_DFOCC_AMPS);
        } else
            t2.reset();
    }

    if (conver == 1) {
        outfile->Printf("\n");
        outfile->Printf(" ============================================================================== \n");
        outfile->Printf(" ===================== DF-CCD ITERATIONS ARE CONVERGED ======================== \n");
        outfile->Printf(" ============================================================================== \n");
    }

    else if (conver == 0) {
        outfile->Printf("\n ====================== DF-CCD IS NOT CONVERGED IN %2d ITERATIONS ============= \n",
                        cc_maxiter);
        throw PSIEXCEPTION("DF-CCD iterations did not converge");
    }

    // Debug
    //if (reference_ == "RESTRICTED") {
    //    auto T = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    //    T->read_symm(psio_, PSIF_DFOCC_AMPS);
    //    auto U = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    //    U->sort(1324, T, 1.0, 0.0);
    //    T.reset();
    //    U->print();
    //    U.reset();
    //}
    //else if (reference_ == "UNRESTRICTED") {
    //    auto U = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    //    U->read(psio_, PSIF_DFOCC_AMPS);
    //    U->print();
    //    U.reset();
    //}

}  // end ccd_iterations

//=========================
// CCD Step
//=========================
void DFOCC::ccd_step() {
    //==========================================================================================
    //========================= CCD iterations =================================================
    //==========================================================================================
    // DIIS
    if (do_diis_ == 1) {
        if (reference_ == "RESTRICTED") {
            std::shared_ptr<Matrix> T2(new Matrix("T2", naoccA * navirA, naoccA * navirA));
            if (reference_ == "RESTRICTED") {
                ccsdDiisManager = std::shared_ptr<DIISManager>(
                    new DIISManager(cc_maxdiis_, "CCSD DIIS T Amps", DIISManager::LargestError, DIISManager::OnDisk));
                ccsdDiisManager->set_error_vector_size(1, DIISEntry::Matrix, T2.get());
                ccsdDiisManager->set_vector_size(1, DIISEntry::Matrix, T2.get());
            }
            T2.reset();
        }
        else if (reference_ == "UNRESTRICTED") {
            std::shared_ptr<Matrix> T2AA(new Matrix("T2AA", ntri_anti_ijAA, ntri_anti_abAA));
            std::shared_ptr<Matrix> T2BB(new Matrix("T2BB", ntri_anti_ijBB, ntri_anti_abBB));
            std::shared_ptr<Matrix> T2AB(new Matrix("T2AB", naoccA * naoccB, navirA * navirB));

            ccsdDiisManager = std::shared_ptr<DIISManager>(
                new DIISManager(cc_maxdiis_, "CCSD DIIS T Amps", DIISManager::LargestError, DIISManager::OnDisk));
            ccsdDiisManager->set_error_vector_size(3, DIISEntry::Matrix, T2AA.get(), DIISEntry::Matrix, T2BB.get(),
                                                   DIISEntry::Matrix, T2AB.get());
            ccsdDiisManager->set_vector_size(3, DIISEntry::Matrix, T2AA.get(), DIISEntry::Matrix, T2BB.get(),
                                             DIISEntry::Matrix, T2AB.get());
            T2AA.reset();
            T2BB.reset();
            T2AB.reset();
        }
    }  // if diis true


    // head of loop
        // 3-index intermediates
        timer_on("CCSD 3-index intr");
        ccd_3index_intr();
        timer_off("CCSD 3-index intr");

        // F intermediates
        timer_on("CCSD F intr");
        ccd_F_intr();
        //std::cout << "ccd_F_intr is done. \n";
        timer_off("CCSD F intr");

        // T2 amplitudes
        timer_on("T2 AMPS");
        ccd_t2_amps();
        timer_off("T2 AMPS");

    // delete
    if (do_diis_ == 1) ccsdDiisManager->delete_diis_file();

    // Mem alloc for DF ints
    //if (df_ints_incore) {
        if (cc_lambda_ == "FALSE") {
            reset_mo_df_ints();
        }
    //}

    // free t2 amps
    if (t2_incore && reference_ == "RESTRICTED") {
        if (cc_lambda_ == "TRUE") {
            t2->write_symm(psio_, PSIF_DFOCC_AMPS);
        } else
            t2.reset();
    }

}  // end ccd_step


void DFOCC::malloc_t2_rhf() {
     t2 = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
     t2->read_symm(psio_, PSIF_DFOCC_AMPS);
}//

void DFOCC::reset_t2_rhf() {
     t2->write_symm(psio_, PSIF_DFOCC_AMPS);
     t2.reset();
}//

void DFOCC::malloc_l2_rhf() {
     l2 = std::make_shared<Tensor2d>("L2 (IA|JB)", naoccA, navirA, naoccA, navirA);
     l2->read_symm(psio_, PSIF_DFOCC_AMPS);
}//

void DFOCC::reset_l2_rhf() {
     l2->write_symm(psio_, PSIF_DFOCC_AMPS);
     l2.reset();
}//

void DFOCC::malloc_mo_df_ints() {
    // Read DF integrals
    bQijA = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA);
    bQiaA = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
    bQabA = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA);
    bQijA->read(psio_, PSIF_DFOCC_INTS);
    bQiaA->read(psio_, PSIF_DFOCC_INTS);
    bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);

    if (reference_ == "UNRESTRICTED") {
        bQijB = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|ij)", nQ, naoccB, naoccB);
        bQiaB = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|ia)", nQ, naoccB, navirB);
        bQabB = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|ab)", nQ, navirB, navirB);
        bQijB->read(psio_, PSIF_DFOCC_INTS);
        bQiaB->read(psio_, PSIF_DFOCC_INTS);
        bQabB->read(psio_, PSIF_DFOCC_INTS, true, true);
    }
}//

void DFOCC::reset_mo_df_ints() {
    bQijA.reset();
    bQiaA.reset();
    bQabA.reset();

    if (reference_ == "UNRESTRICTED") {
        bQijB.reset();
        bQiaB.reset();
        bQabB.reset();
    }
}//

}  // namespace dfoccwave
}  // namespace psi
