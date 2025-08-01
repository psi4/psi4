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
#include "psi4/libdiis/diismanager.h"
#include "psi4/libmints/matrix.h"

#include <cmath>

namespace psi {
namespace dfoccwave {

void DFOCC::ccd_iterations_low() {
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

    // DIIS
    if (do_diis_ == 1) {
        std::shared_ptr<Matrix> T2(new Matrix("T2", naoccA * navirA, naoccA * navirA));
        if (reference_ == "RESTRICTED") {
            ccsdDiisManager = std::shared_ptr<DIISManager>(
                new DIISManager(cc_maxdiis_, "CCSD DIIS T Amps", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::OnDisk));
            ccsdDiisManager->set_error_vector_size(T2.get());
            ccsdDiisManager->set_vector_size(T2.get());
        }
        T2.reset();
    }  // if diis true

    // head of loop
    do {
        // iterate
        itr_occ++;

        // 3-index intermediates
        timer_on("CCD 3-index intr");
        ccd_3index_intr_low();
        timer_off("CCD 3-index intr");

        // F intermediates
        timer_on("CCD F intr");
        ccd_F_intr_low();
        timer_off("CCD F intr");

        // T2 amplitudes
        timer_on("T2 AMPS");
        ccd_t2_amps_low();
        timer_off("T2 AMPS");

        DE = Eccd - Eccd_old;
        Eccd_old = Eccd;

        // RMS
        if (reference_ == "UNRESTRICTED") {
            rms_t2 = MAX0(rms_t2AA, rms_t2BB);
            rms_t2 = MAX0(rms_t2, rms_t2AB);
        }

        // print
        outfile->Printf(" %3d      %13.10f         %13.10f     %12.2e  \n", itr_occ, Ecorr, DE, rms_t2);

        if (itr_occ >= cc_maxiter) {
            conver = 0;  // means iterations were NOT converged
            break;
        }

        if (rms_t2 >= DIVERGE || rms_t1 >= DIVERGE) {
            throw PSIEXCEPTION("CCD iterations are diverging");
        }

    } while (std::fabs(DE) >= tol_Eod || rms_t2 >= tol_t2);

    // delete
    if (do_diis_ == 1) ccsdDiisManager->delete_diis_file();

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

}  // end ccd_iterations_low
}  // namespace dfoccwave
}  // namespace psi
