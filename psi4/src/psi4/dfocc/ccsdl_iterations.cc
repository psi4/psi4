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
#include "psi4/libdiis/diismanager.h"
#include "psi4/libmints/matrix.h"
using namespace std;


namespace psi{ namespace dfoccwave{

void DFOCC::ccsdl_iterations()
{

//==========================================================================================
//========================= Before iterations ==============================================
//==========================================================================================
        // 3-index intermediates
        timer_on("CCSD 3-index intr");
        ccsd_3index_intr();
        timer_off("CCSD 3-index intr");

        // F intermediates
        timer_on("CCSD F intr");
        ccsd_F_intr();
        timer_off("CCSD F intr");

        // W intermediates
        ccsdl_Wmnij();
        ccsdl_Wmbej();
        ccsdl_Wmbje();
        ccsdl_Wmbij();
        //ccsdl_Wmnie();

        // Write
        // Form U_ij^ab
        SharedTensor2d U, T;
        U = SharedTensor2d(new Tensor2d("U2 (IA|JB)", naoccA, navirA, naoccA, navirA));
        ccsd_u2_amps(U,t2);
        U->write_symm(psio_, PSIF_DFOCC_AMPS);
        U.reset();
        // Form Tau_ij^ab
        U = SharedTensor2d(new Tensor2d("Tau (IA|JB)", naoccA, navirA, naoccA, navirA));
        ccsd_tau_amps(U,t2);
        U->write_symm(psio_, PSIF_DFOCC_AMPS);
	T = SharedTensor2d(new Tensor2d("2*Tau(ia,jb) - Tau(ib,ja)", naoccA, navirA, naoccA, navirA));
        ccsd_u2_amps(T,U);
        U.reset();
        T->write_symm(psio_, PSIF_DFOCC_AMPS);
        T.reset();
        // Form T'(ib,ja) = t_ij^ab
        U = SharedTensor2d(new Tensor2d("T2p (IA|JB)", naoccA, navirA, naoccA, navirA));
        ccsd_t2_prime_amps(U,t2);
        U->write_symm(psio_, PSIF_DFOCC_AMPS);
        U.reset();
        // Form T(ij,ab)
        //T = SharedTensor2d(new Tensor2d("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
        //T->sort(1324, t2, 1.0, 0.0);
        //T->write(psio_, PSIF_DFOCC_AMPS);
        //T.reset();
	//t1A->print();

        // Malloc and Free
        gQ = SharedTensor1d(new Tensor1d("CCSDL G_Q", nQ));
        gQp = SharedTensor1d(new Tensor1d("CCSDL G_Qp", nQ));
	l1A = SharedTensor2d(new Tensor2d("L1 <I|A>", naoccA, navirA));
	l1newA = SharedTensor2d(new Tensor2d("New L1 <I|A>", naoccA, navirA));
        l2 = SharedTensor2d(new Tensor2d("L2 (IA|JB)", naoccA, navirA, naoccA, navirA));
        l2->copy(t2);
        t2.reset();
        l1A->copy(t1A);

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
      conver = 1; // Assuming that the iterations will converge
      EccsdL_old = Eccsd;

      // DIIS
      if (do_diis_ == 1) {
          std::shared_ptr<Matrix> L2(new Matrix("L2", naoccA*navirA, naoccA*navirA));
          std::shared_ptr<Matrix> L1(new Matrix("L1", naoccA, navirA));
          if (reference_ == "RESTRICTED") {
              ccsdlDiisManager = std::shared_ptr<DIISManager>(new DIISManager(cc_maxdiis_, "CCSDL DIIS L Amps", DIISManager::LargestError, DIISManager::OnDisk));
              ccsdlDiisManager->set_error_vector_size(2, DIISEntry::Matrix, L2.get(), DIISEntry::Matrix, L1.get());
              ccsdlDiisManager->set_vector_size(2, DIISEntry::Matrix, L2.get(), DIISEntry::Matrix, L1.get());
          }
          L2.reset();
          L1.reset();
      }// if diis true

// head of loop
do
{
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

    // RMS
    if (reference_ == "UNRESTRICTED") {
	rms_t2=MAX0(rms_t2AA,rms_t2BB);
	rms_t2=MAX0(rms_t2,rms_t2AB);
	rms_t1=MAX0(rms_t1A,rms_t1B);
    }

   // print
   outfile->Printf(" %3d      %12.10f         %12.10f      %12.2e  %12.2e \n", itr_occ, EcorrL, DE, rms_t2, rms_t1);

    if (itr_occ >= cc_maxiter) {
      conver = 0; // means iterations were NOT converged
      break;
    }

    if (rms_t2 >= DIVERGE || rms_t1 >= DIVERGE) {
        throw PSIEXCEPTION("CCSD iterations are diverging");
    }

}
while(fabs(DE) >= tol_Eod || rms_t2 >= tol_t2 || rms_t1 >= tol_t2);


if (conver == 1) {
outfile->Printf("\n");
outfile->Printf(" ============================================================================== \n");
outfile->Printf(" ===================== DF-CCSDL ITERATIONS ARE CONVERGED ====================== \n");
outfile->Printf(" ============================================================================== \n");
}

else if (conver == 0) {
  outfile->Printf("\n ====================== DF-CCSDL IS NOT CONVERGED IN %2d ITERATIONS =========== \n", cc_maxiter);
  throw PSIEXCEPTION("DF-CCSDL iterations did not converge");
}

 //delete
 if (do_diis_ == 1) ccsdlDiisManager->delete_diis_file();

     // For GRAD
     if (dertype == "FIRST") {
	 outfile->Printf("\n\tComputing 3-index intermediates...\n");
         timer_on("CCSD PDM 3-index intr");
         ccsd_pdm_3index_intr();
         timer_off("CCSD PDM 3-index intr");
     }

     // free l2 amps
     U = SharedTensor2d(new Tensor2d("Ut2 (IA|JB)", naoccA, navirA, naoccA, navirA));
     ccsd_u2_amps(U,l2);
     U->write_symm(psio_, PSIF_DFOCC_AMPS);
     U.reset();
     l2->write_symm(psio_, PSIF_DFOCC_AMPS);
     l2.reset();

     // For GRAD
     if (dertype == "FIRST") {
         timer_on("CCSD PDM yQia");
	 ccsd_pdm_yQia();
         timer_off("CCSD PDM yQia");
     }

     // Mem free for DF ints
     //if (df_ints_incore) {
     bQijA.reset();
     bQiaA.reset();
     bQabA.reset();
     //}

}// end ccsdl_iterations
}} // End Namespaces
