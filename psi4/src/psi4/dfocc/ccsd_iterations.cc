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

void DFOCC::ccsd_iterations() {
    // CD-PPL
    if (Wabef_type_ == "CD" && reference_ == "RESTRICTED") cd_abcd_cints();

    outfile->Printf("\n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ================ Performing DF-CCSD iterations... ============================ \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf("\n");
    outfile->Printf("  Iter       E_corr                  DE                 T2 RMS        T1 RMS     \n");
    outfile->Printf("  ----   ----------------      ----------------       ----------    ----------   \n");

    //==========================================================================================
    //========================= CCSD iterations ================================================
    //==========================================================================================
    itr_occ = 0;
    conver = 1;  // Assuming that the iterations will converge
    Eccsd_old = Eccsd;

    // DIIS
    if (do_diis_ == 1) {
        if (reference_ == "RESTRICTED") {
            std::shared_ptr<Matrix> T2(new Matrix("T2", naoccA * navirA, naoccA * navirA));
            std::shared_ptr<Matrix> T1(new Matrix("T1", naoccA, navirA));
            ccsdDiisManager = std::shared_ptr<DIISManager>(
                new DIISManager(cc_maxdiis_, "CCSD DIIS T Amps", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::OnDisk));
            ccsdDiisManager->set_error_vector_size(T2.get(), T1.get());
            ccsdDiisManager->set_vector_size(T2.get(), T1.get());
            T2.reset();
            T1.reset();
        }
        else if (reference_ == "UNRESTRICTED") {
            //=== BEGIN DFUCCSD ===
            std::shared_ptr<Matrix> T2AA(new Matrix("T2AA", ntri_anti_ijAA, ntri_anti_abAA));
            std::shared_ptr<Matrix> T2BB(new Matrix("T2BB", ntri_anti_ijBB, ntri_anti_abBB));
            std::shared_ptr<Matrix> T2AB(new Matrix("T2AB", naoccA * naoccB, navirA * navirB));
            std::shared_ptr<Matrix> T1A(new Matrix("T1A", naoccA, navirA));
            std::shared_ptr<Matrix> T1B(new Matrix("T1B", naoccB, navirB));

            ccsdDiisManager = std::shared_ptr<DIISManager>(
                new DIISManager(cc_maxdiis_, "CCSD DIIS T Amps", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::OnDisk));
            ccsdDiisManager->set_error_vector_size(T2AA.get(), T2BB.get(), T2AB.get(), T1A.get(), T1B.get());
            ccsdDiisManager->set_vector_size(T2AA.get(), T2BB.get(), T2AB.get(), T1A.get(), T1B.get());
            T2AA.reset();
            T2BB.reset();
            T2AB.reset();
            T1A.reset();
            T1B.reset();
            //=== END DFUCCSD ===
        }
    }  // if diis true

    // head of loop
    do {
        // iterate
        itr_occ++;

        // 3-index intermediates
        timer_on("CCSD 3-index intr");
        ccsd_3index_intr();
        timer_off("CCSD 3-index intr");

        // F intermediates
        timer_on("CCSD F intr");
        ccsd_F_intr();
        timer_off("CCSD F intr");

        // T1 amplitudes
        timer_on("T1 AMPS");
        ccsd_t1_amps();
        timer_off("T1 AMPS");

        if (reference_ == "UNRESTRICTED") {
            timer_on("CCSD W intr");
            uccsd_W_MBEJAAAA();
            uccsd_W_mbejBBBB();
            uccsd_W_mBeJBABA();
            uccsd_W_MbEjABAB();
            uccsd_W_mBEjBAAB();
            uccsd_W_MbeJABBA();
            timer_off("CCSD W intr");
        }

        // T2 amplitudes
        timer_on("T2 AMPS");
        ccsd_t2_amps();
        timer_off("T2 AMPS");

        DE = Eccsd - Eccsd_old;
        Eccsd_old = Eccsd;

        // RMS
        //if (reference_ == "UNRESTRICTED") {
        //    rms_t2 = MAX0(rms_t2AA, rms_t2BB);
        //    rms_t2 = MAX0(rms_t2, rms_t2AB);
        //    rms_t1 = MAX0(rms_t1A, rms_t1B);
        //}

        // print
        outfile->Printf(" %3d      %13.10f         %13.10f     %12.2e  %12.2e \n", itr_occ, Ecorr, DE, rms_t2, rms_t1);

        if (itr_occ >= cc_maxiter) {
            conver = 0;  // means iterations were NOT converged
            break;
        }

        if (rms_t2 >= DIVERGE || rms_t1 >= DIVERGE) {
            throw PSIEXCEPTION("CCSD iterations are diverging");
        }

    } while (std::fabs(DE) >= tol_Eod || rms_t2 >= tol_t2 || rms_t1 >= tol_t2);

    // delete
    if (do_diis_ == 1) ccsdDiisManager->delete_diis_file();

    // Mem alloc for DF ints
    if (df_ints_incore && cc_lambda_ == "FALSE") {
        reset_mo_df_ints();
    }

    // free t2 amps
    if (t2_incore && reference_ == "RESTRICTED") {
        /*
        if (cc_lambda_ == "TRUE") {
           t2->write_symm(psio_, PSIF_DFOCC_AMPS);
        }
        else t2.reset();
        */
        t2->write_symm(psio_, PSIF_DFOCC_AMPS);
        if (cc_lambda_ == "FALSE") t2.reset();
    }

    if (conver == 1) {
        outfile->Printf("\n");
        outfile->Printf(" ============================================================================== \n");
        outfile->Printf(" ===================== DF-CCSD ITERATIONS ARE CONVERGED ======================= \n");
        outfile->Printf(" ============================================================================== \n");

        // T1 Diagnostic
        double t1diag, t1norm, t1_ref;
        t1_ref = 0.02;
        t1norm = t1A->norm();
        t1diag = t1norm / std::sqrt(2.0 * naoccA);
        outfile->Printf("\n\tT1 diagnostic reference value: %20.14f\n", t1_ref);
        outfile->Printf("\tT1 diagnostic                : %20.14f\n", t1diag);

        // write to disk
        t1A->write(psio_, PSIF_DFOCC_AMPS);
        FijA->write(psio_, PSIF_DFOCC_AMPS);
        FabA->write(psio_, PSIF_DFOCC_AMPS);
        FtijA->write(psio_, PSIF_DFOCC_AMPS);
        FtabA->write(psio_, PSIF_DFOCC_AMPS);
        if (reference_ == "UNRESTRICTED") {
            t1B->write(psio_, PSIF_DFOCC_AMPS);
            FijB->write(psio_, PSIF_DFOCC_AMPS);
            FabB->write(psio_, PSIF_DFOCC_AMPS);
            FtijB->write(psio_, PSIF_DFOCC_AMPS);
            FtabB->write(psio_, PSIF_DFOCC_AMPS);
        }

    }

    else if (conver == 0) {
        outfile->Printf("\n ====================== DF-CCSD IS NOT CONVERGED IN %2d ITERATIONS ============ \n",
                        cc_maxiter);
        throw PSIEXCEPTION("DF-CCSD iterations did not converge");
    }

}  // end ccsd_iterations
}  // namespace dfoccwave
}  // namespace psi
