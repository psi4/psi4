/** Standard library includes */
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string> 
#include <iomanip>
#include <vector> 


/** Required PSI4 includes */
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <libtrans/mospace.h>
#include <libtrans/integraltransform.h>
#include <libdiis/diismanager.h>


/** Required libmints includes */
#include <libmints/mints.h>
#include <libmints/factory.h>
#include <libmints/wavefunction.h>

#include "occwave.h"
#include "defines.h"
#include "arrays.h"

using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace occwave{
  
void OCCWave::omp3_iterations()
{
  
fprintf(outfile,"\n  \n");      
fprintf(outfile," ============================================================================== \n");    
fprintf(outfile," ================ Performing OMP3 iterations... =============================== \n");  
fprintf(outfile," ============================================================================== \n");
fprintf(outfile, "\t            Minimizing MP3-L Functional \n");
fprintf(outfile, "\t            --------------------------- \n");
fprintf(outfile, " Iter       E_total           DE           RMS MO Grad      MAX MO Grad      RMS T2    \n");
fprintf(outfile, " ----    ---------------    ----------     -----------      -----------     ---------- \n");
fflush(outfile);

  
/********************************************************************************************/
/************************** NR iterations **************************************************/
/********************************************************************************************/
      itr_occ = 0;
      mu_ls = -lshift_parameter;
      conver = 1; // Assuming that the MOs will be optimized.
      mo_optimized = 0; 
      itr_diis = 0;

      if (do_diis_ == 1) {
	nvar = num_vecs +1;
        vecsA = new Array2d(num_vecs, nidpA, "Alpha MO DIIS Vectors");
        errvecsA = new Array2d(num_vecs, nidpA, "Alpha MO DIIS Error Vectors");
        vecsA->zero();
        errvecsA->zero();

        if (reference_ == "UNRESTRICTED") {
            vecsB = new Array2d(num_vecs, nidpB, "Beta MO DIIS Vectors");
            errvecsB = new Array2d(num_vecs, nidpB, "Beta MO DIIS Error Vectors");
            vecsB->zero();
            errvecsB->zero();
        }
      }
      
      if (opt_method == "NR") {
          r_pcgA = new Array1d(nidpA, "Alpha PCG r vector");
          z_pcgA = new Array1d(nidpA, "Alpha PCG z vector");
          p_pcgA = new Array1d(nidpA, "Alpha PCG p vector");
          r_pcg_newA = new Array1d(nidpA, "Alpha New PCG r vector");
          z_pcg_newA = new Array1d(nidpA, "Alpha New PCG z vector");
          p_pcg_newA = new Array1d(nidpA, "Alpha New PCG p vector");
          sigma_pcgA = new Array1d(nidpA, "Alpha PCG sigma vector");
          Minv_pcgA = new Array1d(nidpA, "Alpha PCG inverse of M matrix");
          r_pcgA->zero();
          z_pcgA->zero();
          sigma_pcgA->zero();
          p_pcgA->zero();
          Minv_pcgA->zero();

        if (pcg_beta_type_ == "POLAK_RIBIERE") {
          dr_pcgA = new Array1d(nidpA, "Alpha PCG dr vector");
          r_pcgA->zero();
        }

        if (reference_ == "UNRESTRICTED") {
            r_pcgB = new Array1d(nidpB, "Beta PCG r vector");
            z_pcgB = new Array1d(nidpB, "Beta PCG z vector");
            p_pcgB = new Array1d(nidpB, "Beta PCG p vector");
            r_pcg_newB = new Array1d(nidpB, "Beta New PCG r vector");
            z_pcg_newB = new Array1d(nidpB, "Beta New PCG z vector");
            p_pcg_newB = new Array1d(nidpB, "Beta New PCG p vector");
            sigma_pcgB = new Array1d(nidpB, "Beta PCG sigma vector");
            Minv_pcgB = new Array1d(nidpB, "Beta PCG inverse of M matrix");
            r_pcgB->zero();
            z_pcgB->zero();
            sigma_pcgB->zero();
            p_pcgB->zero();
            Minv_pcgB->zero();
            if (pcg_beta_type_ == "POLAK_RIBIERE") {
                dr_pcgB = new Array1d(nidpB, "Alpha PCG dr vector");
                r_pcgB->zero();
            }
        }
      }
 
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
        if (opt_method == "NR") kappa_orb_resp();
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
        timer_on("T2(1)");
        omp3_t2_1st_general();  
        timer_off("T2(1)");
        timer_on("T2(2)");
	t2_2nd_general();  
        timer_off("T2(2)");

/********************************************************************************************/
/************************** One-particle and two-particle density matrices ******************/
/********************************************************************************************/
        timer_on("Response PDMs");
	omp3_response_pdms();
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
	
fprintf(outfile," %3d     %12.10f  %12.2e   %12.2e     %12.2e    %12.2e \n",
	           itr_occ,Emp3L,DE,rms_wog,biggest_mograd,rms_t2);
//fprintf(outfile," %3d     %12.10f  %12.2e   %12.2e     %12.2e    %12.2e  %12.2e  %12.2e \n",
//	           itr_occ,Emp3L,DE,rms_wog,biggest_mograd,rms_kappa,biggest_kappa,rms_t2);
fflush(outfile);

/********************************************************************************************/
/********************************************************************************************/   
    if (itr_occ >= mo_maxiter) {
      conver = 0; // means MOs are NOT optimized
      break;  
    }

    if (rms_wog < tol_grad && biggest_mograd < mograd_max) break;

}
while(rms_wog >= tol_grad || biggest_mograd >= mograd_max); 
//while(fabs(DE) >= tol_Eod || rms_wog >= tol_grad || rms_kappa >= tol_grad || biggest_mograd >= mograd_max || 
//      biggest_kappa >= mograd_max || rms_t2 >= tol_t2); 

if (conver == 1) {
mo_optimized = 1; 
fprintf(outfile,"\n");
fprintf(outfile," ============================================================================== \n");
fprintf(outfile," ======================== OMP3 ITERATIONS ARE CONVERGED ======================= \n");
fprintf(outfile," ============================================================================== \n");
fflush(outfile);
}

else if (conver == 0) {
  fprintf(outfile,"\n ======================== OMP3 IS NOT CONVERGED IN %2d ITERATIONS ============= \n", mo_maxiter);
  fflush(outfile);
}

        // Clean up!
	delete [] idprowA;
	delete [] idpcolA;
	delete [] idpirrA;
        delete wogA;
	delete kappaA;
	delete kappa_newA;
	delete kappa_barA;

        if (reference_ == "UNRESTRICTED") {
	delete [] idprowB;
	delete [] idpcolB;
	delete [] idpirrB;
	delete wogB;
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

	if (opt_method == "NR") {
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

/********************************************************************************************/
/********************************************************************************************/

}
}} // End Namespaces

