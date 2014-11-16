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
  
void DFOCC::occ_iterations()
{
   
outfile->Printf("\n");
outfile->Printf(" ============================================================================== \n");
if (wfn_type_ == "DF-OMP2") outfile->Printf(" ================ Performing DF-OMP2 iterations... ============================ \n");  
else if (wfn_type_ == "DF-OMP3") outfile->Printf(" ================ Performing DF-OMP3 iterations... ============================ \n");  
else if (wfn_type_ == "DF-OCEPA") outfile->Printf(" ================ Performing DF-OCEPA iterations... =========================== \n");  
else if (wfn_type_ == "DF-OMP2.5") outfile->Printf(" ================ Performing DF-OMP2.5 iterations... ========================== \n");  
else if (wfn_type_ == "CD-OMP2") outfile->Printf(" ================ Performing CD-OMP2 iterations... ============================ \n");  
outfile->Printf(" ============================================================================== \n");
if (wfn_type_ == "DF-OMP2") outfile->Printf( "\t            Minimizing DF-MP2-L Functional \n");
else if (wfn_type_ == "DF-OMP3") outfile->Printf( "\t            Minimizing DF-MP3-L Functional \n");
else if (wfn_type_ == "DF-OCEPA") outfile->Printf( "\t            Minimizing DF-CEPA-L Functional \n");
else if (wfn_type_ == "DF-OMP2.5") outfile->Printf( "\t            Minimizing DF-MP2.5-L Functional \n");
else if (wfn_type_ == "CD-OMP2") outfile->Printf( "\t            Minimizing CD-MP2-L Functional \n");
outfile->Printf( "\t            ------------------------------ \n");
outfile->Printf( " Iter       E_total           DE           RMS MO Grad      MAX MO Grad      RMS T2    \n");
outfile->Printf( " ----    ---------------    ----------     -----------      -----------     ---------- \n");


  
//==========================================================================================
//========================= NR iterations ==================================================
//==========================================================================================
      itr_occ = 0;
      mu_ls = 0;
      conver = 1; // Assuming that the MOs will be optimized.
      mo_optimized = 0; 
      itr_diis = 0;

      // If diis?
      //if (noccA + noccB != 1) {
          if (do_diis_ == 1) {
              nvar = num_vecs +1;
              vecsA = SharedTensor2d(new Tensor2d("Alpha MO DIIS Vectors", num_vecs, nidpA));
              errvecsA = SharedTensor2d(new Tensor2d("Alpha MO DIIS Error Vectors", num_vecs, nidpA));

              if (reference_ == "UNRESTRICTED") {
                  vecsB = SharedTensor2d(new Tensor2d("Beta MO DIIS Vectors", num_vecs, nidpB));
                  errvecsB = SharedTensor2d(new Tensor2d("Beta MO DIIS Error Vectors", num_vecs, nidpB));
              }
          }
      //}
      
//==========================================================================================
//========================= Head of the Loop ===============================================
//==========================================================================================
do
{
       itr_occ++;
	
//==========================================================================================
//========================= New orbital step ===============================================
//==========================================================================================
        timer_on("kappa orb rot");
        if (hess_type == "HF") { 
           if (orb_resp_solver_ == "LINEQ") kappa_orb_resp();
           else if (orb_resp_solver_ == "PCG") kappa_orb_resp_pcg();
        }
        else kappa_diag_hess();
        timer_off("kappa orb rot");

//==========================================================================================
//========================= update mo coefficients =========================================
//==========================================================================================
        timer_on("update_mo");
        update_mo();
        timer_off("update_mo");

//==========================================================================================
//========================= Trans TEI ======================================================
//==========================================================================================
    // DF
    if (do_cd == "FALSE") { 
        timer_on("DF CC Integrals");
        trans_corr();
        timer_off("DF CC Integrals");

        timer_on("DF REF Integrals");
        trans_ref();
        timer_off("DF REF Integrals");
    }// end if (do_cd == "FALSE")  

    // CD
    else if (do_cd == "TRUE") { 
        timer_on("CD Trans");
        trans_cd();
        timer_off("CD Trans");
    }// end if (do_cd == "TRUE")  

        // Fock
        fock();

//==========================================================================================
//========================= New Amplitudes =================================================
//==========================================================================================
     if (wfn_type_ == "DF-OMP2" || wfn_type_ == "CD-OMP2") t2_1st_gen();

    /*
     else if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") {
        timer_on("T2(1)");
        omp3_t2_1st_general();  
        timer_off("T2(1)");
        timer_on("T2(2)");
	t2_2nd_general();  
        timer_off("T2(2)");
     }

     else if (wfn_type_ == "OCEPA") {
        timer_on("T2");
	t2_amps();  
        timer_off("T2");
     }
     */

//==========================================================================================
//========================= PDMs ===========================================================
//==========================================================================================
	if (wfn_type_ == "DF-OMP2"  || wfn_type_ == "CD-OMP2") {
	    omp2_opdm();
	    omp2_tpdm();
            separable_tpdm();
        }
        
        //else if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") omp3_response_pdms();
	//else if (wfn_type_ == "OCEPA") ocepa_response_pdms();

//==========================================================================================
//========================= GFM ============================================================
//==========================================================================================
        gfock_vo();
        gfock_ov();
        gfock_oo();
        gfock_vv();
        //if (nfrzc > 0) gfock_oo();
        //if (nfrzv > 0) gfock_vv();
	
//==========================================================================================
//========================= CCL ============================================================
//==========================================================================================
        // reference energy
        ref_energy();

     if (wfn_type_ == "DF-OMP2" || wfn_type_ == "CD-OMP2") mp2l_energy();

     /*
     else if (wfn_type_ == "DF-OMP3" || wfn_type_ == "DF-OMP2.5") { 
        if (compute_ccl == "TRUE") {
           timer_on("MP3L Energy");
	   ccl_energy();
           timer_off("MP3L Energy");
        }
       
        else {
           timer_on("REF Energy");
           ref_energy();
           timer_off("REF Energy");
           timer_on("MP3 Energy");
           mp3_energy();
           timer_off("MP3 Energy");
           Emp3L = Emp3;
           EcorrL = Emp3L-Escf;
           DE = Emp3L - Emp3L_old;
           Emp3L_old = Emp3L;
        }
     }

     else if (wfn_type_ == "DF-OCEPA") { 
        if (compute_ccl == "TRUE") {
           timer_on("CEPAL Energy");
	   ccl_energy();
           timer_off("CEPAL Energy");
        }
       
        else {
           timer_on("REF Energy");
           ref_energy();
           timer_off("REF Energy");
           timer_on("CEPA Energy");
           cepa_energy();
           timer_off("CEPA Energy");
           EcepaL = Ecepa;
           EcorrL = EcepaL-Escf;
           DE = EcepaL - EcepaL_old;
           EcepaL_old = EcepaL;
        }
     }
    */

//==========================================================================================
//========================= MO Grad ========================================================
//==========================================================================================
        timer_on("MO Grad");
	mograd();
        timer_off("MO Grad");
      
//==========================================================================================
//========================= Print ==========================================================
//==========================================================================================
    if (reference_ == "RESTRICTED") {
	nidp=nidpA;
	rms_wog=rms_wogA;
	biggest_mograd=biggest_mogradA;
	rms_kappa=rms_kappaA;
	biggest_kappa=biggest_kappaA;
    }

    else if (reference_ == "UNRESTRICTED") {
	nidp=MAX0(nidpA,nidpB);
	rms_wog=MAX0(rms_wogA,rms_wogB);
	biggest_mograd=MAX0(biggest_mogradA,biggest_mogradB);
	rms_kappa=MAX0(rms_kappaA,rms_kappaB);
	biggest_kappa=MAX0(biggest_kappaA,biggest_kappaB);
	rms_t2=MAX0(rms_t2AA,rms_t2BB);
	rms_t2=MAX0(rms_t2,rms_t2AB);
    }
	
if(wfn_type_ == "DF-OMP2") outfile->Printf(" %3d     %12.10f  %12.2e   %12.2e     %12.2e    %12.2e \n",itr_occ,Emp2L,DE,rms_wog,biggest_mograd,rms_t2);
else if(wfn_type_ == "DF-OMP3") outfile->Printf(" %3d     %12.10f  %12.2e   %12.2e     %12.2e    %12.2e \n",itr_occ,Emp3L,DE,rms_wog,biggest_mograd,rms_t2);
else if(wfn_type_ == "DF-OCEPA") outfile->Printf(" %3d     %12.10f  %12.2e   %12.2e     %12.2e    %12.2e \n",itr_occ,EcepaL,DE,rms_wog,biggest_mograd,rms_t2);
else if(wfn_type_ == "DF-OMP2.5") outfile->Printf(" %3d     %12.10f  %12.2e   %12.2e     %12.2e    %12.2e \n",itr_occ,Emp3L,DE,rms_wog,biggest_mograd,rms_t2);
else if(wfn_type_ == "CD-OMP2") outfile->Printf(" %3d     %12.10f  %12.2e   %12.2e     %12.2e    %12.2e \n",itr_occ,Emp2L,DE,rms_wog,biggest_mograd,rms_t2);


//==========================================================================================
//========================= Convergence? ===================================================
//==========================================================================================
    if (itr_occ >= mo_maxiter) {
      conver = 0; // means MOs are NOT optimized
      break;  
    }

    if (rms_wog < tol_grad && biggest_mograd < mograd_max) break;
    if (fabs(DE) <= tol_Eod) break;

    if (rms_wog >= DIVERGE) {
        throw PSIEXCEPTION("DF-OCC iterations are diverging");
    }

}
while(rms_wog >= tol_grad || biggest_mograd >= mograd_max); 

if (conver == 1) {
mo_optimized = 1; 
outfile->Printf("\n");
outfile->Printf(" ============================================================================== \n");
if (wfn_type_ == "DF-OMP2") outfile->Printf(" ======================== DF-OMP2 ITERATIONS ARE CONVERGED ==================== \n");
else if (wfn_type_ == "DF-OMP3") outfile->Printf(" ======================== DF-OMP3 ITERATIONS ARE CONVERGED ==================== \n");
else if (wfn_type_ == "DF-OCEPA") outfile->Printf(" ======================== DF-OCEPA ITERATIONS ARE CONVERGED =================== \n");
else if (wfn_type_ == "DF-OMP2.5") outfile->Printf(" ======================== DF-OMP2.5 ITERATIONS ARE CONVERGED ================== \n");
outfile->Printf(" ============================================================================== \n");

}

else if (conver == 0) {
  if (wfn_type_ == "DF-OMP2") outfile->Printf("\n ======================== DF-OMP2 IS NOT CONVERGED IN %2d ITERATIONS ========== \n", mo_maxiter);
  else if (wfn_type_ == "DF-OMP3") outfile->Printf("\n ======================== DF-OMP3 IS NOT CONVERGED IN %2d ITERATIONS ========== \n", mo_maxiter);
  else if (wfn_type_ == "DF-OCEPA") outfile->Printf("\n ======================== DF-OCEPA IS NOT CONVERGED IN %2d ITERATIONS ========= \n", mo_maxiter);
  else if (wfn_type_ == "DF-OMP2.5") outfile->Printf("\n ======================== DF-OMP2.5 IS NOT CONVERGED IN %2d ITERATIONS ======== \n", mo_maxiter);
  
  throw PSIEXCEPTION("DF-OCC iterations did not converge");
}

}
}} // End Namespaces

