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

#include "occwave.h"
#include "defines.h"


using namespace psi;


namespace psi{ namespace occwave{

void OCCWave::occ_iterations()
{

outfile->Printf("\n");
outfile->Printf(" ============================================================================== \n");
if (wfn_type_ == "OMP2") outfile->Printf(" ================ Performing OMP2 iterations... =============================== \n");
else if (wfn_type_ == "OMP3") outfile->Printf(" ================ Performing OMP3 iterations... =============================== \n");
else if (wfn_type_ == "OCEPA") outfile->Printf(" ================ Performing OCEPA iterations... ============================== \n");
else if (wfn_type_ == "OMP2.5") outfile->Printf(" ================ Performing OMP2.5 iterations... ============================= \n");
outfile->Printf(" ============================================================================== \n");
if (wfn_type_ == "OMP2") outfile->Printf( "\t            Minimizing MP2-L Functional \n");
else if (wfn_type_ == "OMP3") outfile->Printf( "\t            Minimizing MP3-L Functional \n");
else if (wfn_type_ == "OCEPA") outfile->Printf( "\t            Minimizing CEPA-L Functional \n");
else if (wfn_type_ == "OMP2.5") outfile->Printf( "\t            Minimizing MP2.5-L Functional \n");
outfile->Printf( "\t            --------------------------- \n");
outfile->Printf( " Iter       E_total           DE           RMS MO Grad      MAX MO Grad      RMS T2    \n");
outfile->Printf( " ----    ---------------    ----------     -----------      -----------     ---------- \n");



/********************************************************************************************/
/************************** NR iterations **************************************************/
/********************************************************************************************/
      itr_occ = 0;
      mu_ls = 0;
      conver = 1; // Assuming that the MOs will be optimized.
      mo_optimized = 0;
      itr_diis = 0;

      // If diis?
      //if (nooA + nooB != 1) {
          if (do_diis_ == 1) {
              nvar = num_vecs +1;
              vecsA = new Array2d("Alpha MO DIIS Vectors", num_vecs, nidpA);
              errvecsA = new Array2d("Alpha MO DIIS Error Vectors", num_vecs, nidpA);
              vecsA->zero();
              errvecsA->zero();

              if (reference_ == "UNRESTRICTED") {
                  vecsB = new Array2d("Beta MO DIIS Vectors", num_vecs, nidpB);
                  errvecsB = new Array2d("Beta MO DIIS Vectors", num_vecs, nidpB);
                  vecsB->zero();
                  errvecsB->zero();
              }
          }
      //}

      // Set up the orb-resp algorithm
      if (opt_method == "ORB_RESP") {

       // Lineq
       if (orb_resp_solver_ == "LINEQ") {
         if (reference_ == "RESTRICTED") Aorb = new Array2d("MO Hessian Matrix", nidpA, nidpA);
         else if (reference_ == "UNRESTRICTED") {
             nidp_tot = nidpA + nidpB;
             kappa = new Array1d("Total orb rot params vector of current step", nidp_tot);
         }
       }

       // PCG
       else if (orb_resp_solver_ == "PCG") {
          r_pcgA = new Array1d("Alpha PCG r vector", nidpA);
          z_pcgA = new Array1d("Alpha PCG z vector", nidpA);
          p_pcgA = new Array1d("Alpha PCG p vector", nidpA);
          r_pcg_newA = new Array1d("Alpha New PCG r vector", nidpA);
          z_pcg_newA = new Array1d("Alpha New PCG z vector", nidpA);
          p_pcg_newA = new Array1d("Alpha New PCG p vector", nidpA);
          sigma_pcgA = new Array1d("Alpha PCG sigma vector", nidpA);
          Minv_pcgA = new Array1d("Alpha PCG inverse of M matrix", nidpA);
          r_pcgA->zero();
          z_pcgA->zero();
          sigma_pcgA->zero();
          p_pcgA->zero();
          Minv_pcgA->zero();

        if (pcg_beta_type_ == "POLAK_RIBIERE") {
          dr_pcgA = new Array1d("Alpha PCG dr vector", nidpA);
          r_pcgA->zero();
        }

        if (reference_ == "UNRESTRICTED") {
            r_pcgB = new Array1d("Beta PCG r vector", nidpB);
            z_pcgB = new Array1d("Beta PCG z vector", nidpB);
            p_pcgB = new Array1d("Beta PCG p vector", nidpB);
            r_pcg_newB = new Array1d("Beta New PCG r vector", nidpB);
            z_pcg_newB = new Array1d("Beta New PCG z vector", nidpB);
            p_pcg_newB = new Array1d("Beta New PCG p vector", nidpB);
            sigma_pcgB = new Array1d("Beta PCG sigma vector", nidpB);
            Minv_pcgB = new Array1d("Beta PCG inverse of M matrix", nidpB);
            r_pcgB->zero();
            z_pcgB->zero();
            sigma_pcgB->zero();
            p_pcgB->zero();
            Minv_pcgB->zero();
            if (pcg_beta_type_ == "POLAK_RIBIERE") {
                dr_pcgB = new Array1d("Alpha PCG dr vector", nidpB);
                r_pcgB->zero();
            }
        }
       }// pcg if
      }// orb_resp if

/********************************************************************************************/
/************************** Head of the Loop ************************************************/
/********************************************************************************************/
do
{
       itr_occ++;

/********************************************************************************************/
/************************** New orbital step ************************************************/
/********************************************************************************************/
        timer_on("kappa orb rot");
        if (opt_method == "ORB_RESP") {
           if (orb_resp_solver_ == "LINEQ") kappa_orb_resp();
           else if (orb_resp_solver_ == "PCG") kappa_orb_resp_iter();
        }
        else if (opt_method == "MSD") kappa_msd();
        timer_off("kappa orb rot");

/********************************************************************************************/
/************************** update mo coefficients ******************************************/
/********************************************************************************************/
        timer_on("update_mo");
        update_mo();
        timer_off("update_mo");

/********************************************************************************************/
/************************** Transform TEI from SO to MO space *******************************/
/********************************************************************************************/
        timer_on("trans_ints");
	if (reference_ == "RESTRICTED") trans_ints_rhf();
	else if (reference_ == "UNRESTRICTED") trans_ints_uhf();
        timer_off("trans_ints");

/********************************************************************************************/
/************************** NEW amplitudes **************************************************/
/********************************************************************************************/
     if (wfn_type_ == "OMP2") {
        timer_on("T2(1)");
        omp2_t2_1st_general();
        timer_off("T2(1)");
     }

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


/********************************************************************************************/
/************************** One-particle and two-particle density matrices ******************/
/********************************************************************************************/
        timer_on("Response PDMs");
	if (wfn_type_ == "OMP2") omp2_response_pdms();
	else if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") omp3_response_pdms();
	else if (wfn_type_ == "OCEPA") ocepa_response_pdms();
        timer_off("Response PDMs");

/********************************************************************************************/
/************************** Asymmetric Generalizd-Fock matrix *******************************/
/********************************************************************************************/
        timer_on("Generalized-Fock");
	gfock();
        timer_off("Generalized-Fock");

/********************************************************************************************/
/************************** Compute Lagrangian Energy ***************************************/
/********************************************************************************************/
     if (wfn_type_ == "OMP2") {
        timer_on("MP2L Energy");
	ccl_energy();
        timer_off("MP2L Energy");
     }

     else if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") {
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

     else if (wfn_type_ == "OCEPA") {
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

/********************************************************************************************/
/************************** new orbital gradient ********************************************/
/********************************************************************************************/
        timer_on("MO Grad");
	mograd();
        timer_off("MO Grad");

/********************************************************************************************/
/************************** Print ***********************************************************/
/********************************************************************************************/
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

if(wfn_type_ == "OMP2") outfile->Printf(" %3d     %12.10f  %12.2e   %12.2e     %12.2e    %12.2e \n",itr_occ,Emp2L,DE,rms_wog,biggest_mograd,rms_t2);
else if(wfn_type_ == "OMP3") outfile->Printf(" %3d     %12.10f  %12.2e   %12.2e     %12.2e    %12.2e \n",itr_occ,Emp3L,DE,rms_wog,biggest_mograd,rms_t2);
else if(wfn_type_ == "OCEPA") outfile->Printf(" %3d     %12.10f  %12.2e   %12.2e     %12.2e    %12.2e \n",itr_occ,EcepaL,DE,rms_wog,biggest_mograd,rms_t2);
else if(wfn_type_ == "OMP2.5") outfile->Printf(" %3d     %12.10f  %12.2e   %12.2e     %12.2e    %12.2e \n",itr_occ,Emp3L,DE,rms_wog,biggest_mograd,rms_t2);
//outfile->Printf(" %3d     %12.10f  %12.2e   %12.2e     %12.2e    %12.2e  %12.2e  %12.2e \n",
//	           itr_occ,Emp2L,DE,rms_wog,biggest_mograd,rms_kappa,biggest_kappa,rms_t2);


/********************************************************************************************/
/********************************************************************************************/
    if (itr_occ >= mo_maxiter) {
      conver = 0; // means MOs are NOT optimized
      break;
    }

    if (rms_wog < tol_grad && biggest_mograd < mograd_max) break;

    if (rms_wog >= DIVERGE) {
        throw PSIEXCEPTION("OCC iterations are diverging");
    }

}
while(rms_wog >= tol_grad || biggest_mograd >= mograd_max);
//while(fabs(DE) >= tol_Eod || rms_wog >= tol_grad || rms_kappa >= tol_grad || biggest_mograd >= mograd_max ||
//      biggest_kappa >= mograd_max || rms_t2 >= tol_t2);

if (conver == 1) {
mo_optimized = 1;
outfile->Printf("\n");
outfile->Printf(" ============================================================================== \n");
if (wfn_type_ == "OMP2") outfile->Printf(" ======================== OMP2 ITERATIONS ARE CONVERGED ======================= \n");
else if (wfn_type_ == "OMP3") outfile->Printf(" ======================== OMP3 ITERATIONS ARE CONVERGED ======================= \n");
else if (wfn_type_ == "OCEPA") outfile->Printf(" ======================== OCEPA ITERATIONS ARE CONVERGED ====================== \n");
else if (wfn_type_ == "OMP2.5") outfile->Printf(" ======================== OMP2.5 ITERATIONS ARE CONVERGED ===================== \n");
outfile->Printf(" ============================================================================== \n");

}

else if (conver == 0) {
  if (wfn_type_ == "OMP2") outfile->Printf("\n ======================== OMP2 IS NOT CONVERGED IN %2d ITERATIONS ============= \n", mo_maxiter);
  else if (wfn_type_ == "OMP3") outfile->Printf("\n ======================== OMP3 IS NOT CONVERGED IN %2d ITERATIONS ============= \n", mo_maxiter);
  else if (wfn_type_ == "OCEPA") outfile->Printf("\n ======================== OCEPA IS NOT CONVERGED IN %2d ITERATIONS ============ \n", mo_maxiter);
  else if (wfn_type_ == "OMP2.5") outfile->Printf("\n ======================== OMP2.5 IS NOT CONVERGED IN %2d ITERATIONS =========== \n", mo_maxiter);

  throw PSIEXCEPTION("OCC iterations did not converge");
}
        // Clean up!
	delete [] idprowA;
	delete [] idpcolA;
	delete [] idpirrA;
        delete wogA;
        delete wog_intA;
	delete kappaA;
	delete kappa_newA;
	delete kappa_barA;

        if (reference_ == "UNRESTRICTED") {
	delete [] idprowB;
	delete [] idpcolB;
	delete [] idpirrB;
        delete wogB;
        delete wog_intB;
	delete kappaB;
	delete kappa_newB;
	delete kappa_barB;
        }

	if (do_diis_ == 1) {
          delete vecsA;
          delete errvecsA;
          if (reference_ == "UNRESTRICTED") delete vecsB;
          if (reference_ == "UNRESTRICTED") delete errvecsB;
	}

        // Clean up the mess of ORB-RESP
        if (opt_method == "ORB_RESP") {
            if (orb_resp_solver_ == "LINEQ") {
                if (reference_ == "UNRESTRICTED") {
	            delete kappa;
                 }
            }

         // PCG
         else if (orb_resp_solver_ == "PCG") {
          delete r_pcgA;
          delete z_pcgA;
          delete p_pcgA;
          delete sigma_pcgA;
          delete Minv_pcgA;
          delete r_pcg_newA;
          delete z_pcg_newA;
          delete p_pcg_newA;
          if (pcg_beta_type_ == "POLAK_RIBIERE") delete dr_pcgA;
          if(reference_ == "UNRESTRICTED") {
             delete r_pcgB;
             delete z_pcgB;
             delete p_pcgB;
             delete sigma_pcgB;
             delete Minv_pcgB;
             delete r_pcg_newB;
             delete z_pcg_newB;
             delete p_pcg_newB;
             if (pcg_beta_type_ == "POLAK_RIBIERE") delete dr_pcgB;
          }
	 }
        }// end orb resp if

}
}} // End Namespaces
