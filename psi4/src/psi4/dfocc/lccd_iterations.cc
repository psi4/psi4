/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/libqt/qt.h"
#include "defines.h"
#include "dfocc.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libdiis/diismanager.h"
using namespace std;


namespace psi{ namespace dfoccwave{

void DFOCC::lccd_iterations()
{

outfile->Printf("\n");
outfile->Printf(" ============================================================================== \n");
outfile->Printf(" ================ Performing DF-LCCD iterations... ============================= \n");
outfile->Printf(" ============================================================================== \n");
outfile->Printf("\n");
outfile->Printf("  Iter       E_corr                  DE                 T2 RMS      \n");
outfile->Printf("  ----   ----------------      ----------------       ----------    \n");

//==========================================================================================
//========================= CCD iterations =================================================
//==========================================================================================
      itr_occ = 0;
      conver = 1; // Assuming that the iterations will converge
      Elccd_old = Elccd;

      // DIIS
      if (do_diis_ == 1) {
	  // RHF
          if (reference_ == "RESTRICTED") {
              std::shared_ptr<Matrix> T2(new Matrix("T2", naoccA*navirA, naoccA*navirA));
              ccsdDiisManager = std::shared_ptr<DIISManager>(new DIISManager(cc_maxdiis_, "CCSD DIIS T Amps", DIISManager::LargestError, DIISManager::OnDisk));
              ccsdDiisManager->set_error_vector_size(1, DIISEntry::Matrix, T2.get());
              ccsdDiisManager->set_vector_size(1, DIISEntry::Matrix, T2.get());
              T2.reset();
          }

	  // UHF
	  else if (reference_ == "UNRESTRICTED") {
              std::shared_ptr<Matrix> T2AA(new Matrix("T2AA", ntri_anti_ijAA, ntri_anti_abAA));
              std::shared_ptr<Matrix> T2BB(new Matrix("T2BB", ntri_anti_ijBB, ntri_anti_abBB));
              std::shared_ptr<Matrix> T2AB(new Matrix("T2AB", naoccA*naoccB, navirA*navirB));
              ccsdDiisManager = std::shared_ptr<DIISManager>(new DIISManager(cc_maxdiis_, "CCSD DIIS T Amps", DIISManager::LargestError, DIISManager::OnDisk));
              ccsdDiisManager->set_error_vector_size(3, DIISEntry::Matrix, T2AA.get(), DIISEntry::Matrix, T2BB.get(), DIISEntry::Matrix, T2AB.get());
              ccsdDiisManager->set_vector_size(3, DIISEntry::Matrix, T2AA.get(), DIISEntry::Matrix, T2BB.get(), DIISEntry::Matrix, T2AB.get());
              T2AA.reset();
              T2BB.reset();
              T2AB.reset();
          }

      }// if diis true

// head of loop
do
{
        // iterate
        itr_occ++;

        // F intermediates
        timer_on("CCD F intr");
	Fint_zero();
        timer_off("CCD F intr");

        // T2 amplitudes
        timer_on("T2 AMPS");
	lccd_t2_amps();
        timer_off("T2 AMPS");

        DE = Elccd - Elccd_old;
        Elccd_old = Elccd;


   // print
   outfile->Printf(" %3d      %12.10f         %12.10f      %12.2e  \n", itr_occ, Ecorr, DE, rms_t2);

    if (itr_occ >= cc_maxiter) {
      conver = 0; // means iterations were NOT converged
      break;
    }

    if (rms_t2 >= DIVERGE) {
        throw PSIEXCEPTION("LCCD iterations are diverging");
    }

}
while(fabs(DE) >= tol_Eod || rms_t2 >= tol_t2);

 //delete
 if (do_diis_ == 1)  ccsdDiisManager->delete_diis_file();

 /*
 // Mem alloc for DF ints
 if (df_ints_incore && reference_ == "RESTRICTED") {
     bQijA.reset();
     bQiaA.reset();
     bQabA.reset();
 }

 // free t2 amps
 if (t2_incore && reference_ == "RESTRICTED") {
     if (dertype == "FIRST" || orb_opt_ == "TRUE") {
        t2->write_symm(psio_, PSIF_DFOCC_AMPS);
     }
     t2.reset();
 }
 */

if (conver == 1) {
outfile->Printf("\n");
outfile->Printf(" ============================================================================== \n");
outfile->Printf(" ===================== DF-LCCD ITERATIONS ARE CONVERGED ======================= \n");
outfile->Printf(" ============================================================================== \n");
}

else if (conver == 0) {
  outfile->Printf("\n ===================== DF-LCCD IS NOT CONVERGED IN %2d ITERATIONS ============= \n", cc_maxiter);
  throw PSIEXCEPTION("DF-LCCD iterations did not converge");
}

}// end lccd_iterations
}} // End Namespaces
