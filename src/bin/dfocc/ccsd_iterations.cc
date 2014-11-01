/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#include <libqt/qt.h>
#include "defines.h"
#include "dfocc.h"

using namespace psi;
using namespace std;


namespace psi{ namespace dfoccwave{
  
void DFOCC::ccsd_iterations()
{

outfile->Printf("\n");
outfile->Printf(" ============================================================================== \n");
outfile->Printf(" ================ Performing DF-CCSD iterations... ============================ \n");
outfile->Printf(" ============================================================================== \n");
outfile->Printf("\n");
outfile->Printf("  Iter    E_corr           E_total            DE           T2 RMS        T1 RMS     \n");
outfile->Printf("  ----   -------------    ---------------    ----------   ----------    ---------   \n");
  
//==========================================================================================
//========================= CCSD iterations ================================================
//==========================================================================================
      itr_occ = 0;
      conver = 1; // Assuming that the iterations will converge
      Eccsd_old = Eccsd;

      // DIIS
      SharedTensor2d T;
      T = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
      T->read_symm(psio_, PSIF_DFOCC_AMPS);
      if (reference_ == "RESTRICTED") {
          t2DiisManager = new DIISManager(cc_maxdiis_, "CCSD DIIS T2 Amps", DIISManager::LargestError, DIISManager::InCore);
          t2DiisManager->set_error_vector_size(2, DIISEntry::Pointer, &T, 
                                                  DIISEntry::Pointer, &t1A);
          t2DiisManager->set_vector_size(2, DIISEntry::Pointer, &T, 
                                            DIISEntry::Pointer, &t1A);
      }
      T.reset();

      /*
      else if (reference_ == "UNRESTRICTED") {
          t2DiisManager = new DIISManager(cc_maxdiis_, "CEPA DIIS T2 Amps", DIISManager::LargestError, DIISManager::InCore);
          t2DiisManager->set_error_vector_size(3, DIISEntry::DPDBuf4, &Taa,
                                                  DIISEntry::DPDBuf4, &Tbb,
                                                  DIISEntry::DPDBuf4, &Tab);
          t2DiisManager->set_vector_size(3, DIISEntry::DPDBuf4, &Taa,
                                            DIISEntry::DPDBuf4, &Tbb,
                                            DIISEntry::DPDBuf4, &Tab);
      }
      */

// head of loop      
do
{
        itr_occ++;

        // 3-index intermediates
        timer_on("CCSD 3-index intr");
        ccsd_3index_intr();
        timer_off("CCSD 3-index intr");

        // F intermediates
        timer_on("CCSD F intr");
        ccsd_F_intr();
        timer_off("CCSD F intr");

        // W intermediates
        timer_on("CCSD W intr");
        ccsd_W_intr();
        timer_off("CCSD W intr");

        // T1 amplitudes
        timer_on("T1 AMPS");
	ccsd_t1_amps();  
        timer_off("T1 AMPS");
   
        // T2 amplitudes
        timer_on("T2 AMPS");
	ccsd_t2_amps();  
        timer_off("T2 AMPS");

        // CCSD energy
        timer_on("CCSD Energy");
        ccsd_energy();
        timer_off("CCSD Energy");

        DE = Eccsd - Eccsd_old;
        Eccsd_old = Eccsd;

    // RMS
    if (reference_ == "UNRESTRICTED") {
	rms_t2=MAX0(rms_t2AA,rms_t2BB);
	rms_t2=MAX0(rms_t2,rms_t2AB);
	rms_t1=MAX0(rms_t1A,rms_t1B);
    }
	
   // print
   outfile->Printf(" %3d     %12.10f    %12.10f  %12.2e %12.2e %12.2e \n", itr_occ, Ecorr, Eccsd, DE, rms_t2, rms_t1);


    if (itr_occ >= cc_maxiter) {
      conver = 0; // means iterations were NOT converged
      break;  
    }

    if (rms_t2 >= DIVERGE || rms_t1 >= DIVERGE) {
        throw PSIEXCEPTION("CCSD iterations are diverging");
    }

}
while(fabs(DE) >= tol_Eod || rms_t2 >= tol_t2 || rms_t1 >= tol_t2); 

//delete
delete t2DiisManager;

if (conver == 1) {
outfile->Printf("\n");
outfile->Printf(" ============================================================================== \n");
outfile->Printf(" ===================== DF-CCSD ITERATIONS ARE CONVERGED ======================= \n");
outfile->Printf(" ============================================================================== \n");

}

else if (conver == 0) {
  outfile->Printf("\n ====================== DF-CCSD IS NOT CONVERGED IN %2d ITERATIONS ============ \n", cc_maxiter);
  throw PSIEXCEPTION("DF-CCSD iterations did not converge");
}

}// end ccsd_iterations
}} // End Namespaces

