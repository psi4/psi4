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

#include "omp3wave.h"
#include "defines.h"

using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace omp3wave{
  
void OMP3Wave::occ_iterations()
{
  
fprintf(outfile,"\n  \n");      
fprintf(outfile," ============================================================================== \n");    
fprintf(outfile," ================ Performing OMP3 iterations... =============================== \n");  
fprintf(outfile," ============================================================================== \n");
fprintf(outfile,"\n");
fprintf(outfile, "\t            Minimizing MP3-L Functional \n");
fprintf(outfile, "\t            --------------------------- \n");
fprintf(outfile, " Iter       E_total           DE           MO Grad RMS      MAX MO Grad      Korb RMS      MAX Korb      T2 RMS    \n");
fprintf(outfile, " ----    ---------------    ----------     -----------      -----------     ----------    -----------   ----------  \n");
fflush(outfile);

  
/********************************************************************************************/
/************************** NR iterations **************************************************/
/********************************************************************************************/
      itr_occ = 0;
      mu_ls = 0;
      conver = 1; // Assuming that the MOs will be optimized.
      mo_optimized = 0; 
      
      if (opt_method == "DIIS") {
	nvar = num_vecs +1;
	vecsA = block_matrix(num_vecs, nidpA);
	errvecsA = block_matrix(num_vecs, nidpA);
	memset(vecsA[0], 0, sizeof(double)*num_vecs*nidpA);  
	memset(errvecsA[0], 0, sizeof(double)*num_vecs*nidpA); 
      
	vecsB = block_matrix(num_vecs, nidpB);
	errvecsB = block_matrix(num_vecs, nidpB);
	memset(vecsB[0], 0, sizeof(double)*num_vecs*nidpB);  
	memset(errvecsB[0], 0, sizeof(double)*num_vecs*nidpB); 
      }
      
      
// head of loop      
do
{
       itr_occ++;
	
/********************************************************************************************/
/************************** New orbital step ************************************************/
/********************************************************************************************/
        timer_on("kappa orb rot");
        korbrot_sd();
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
        trans_ints(); 
        timer_off("trans_ints");
 
/********************************************************************************************/
/************************** NEW amplitudes **************************************************/
/********************************************************************************************/
        timer_on("T2(1)");
        t2_1st_general();  
        timer_off("T2(1)");
        timer_on("T2(2)");
	t2_2nd_general();  
        timer_off("T2(2)");

/********************************************************************************************/
/************************** One-particle and two-particle density matrices ******************/
/********************************************************************************************/
        timer_on("Response PDMs");
	response_pdms();
        timer_off("Response PDMs");

/********************************************************************************************/
/************************** Asymmetric Generalizd-Fock matrix *******************************/
/********************************************************************************************/
        timer_on("Generalized-Fock");
	GFockmo();
        timer_off("Generalized-Fock");
	
/********************************************************************************************/
/************************** Compute Lagrangian Energy ***************************************/
/********************************************************************************************/
        /*
        timer_on("MP3L Energy");
	mp3l_energy();
        timer_off("MP3L Energy");
        */

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

/********************************************************************************************/
/************************** new orbital gradient ********************************************/
/********************************************************************************************/	 
        timer_on("MO Grad");
	mograd();
        timer_off("MO Grad");
      
/********************************************************************************************/
/************************** Print ***********************************************************/
/********************************************************************************************/
	nidp=MAX0(nidpA,nidpB);
	rms_wog=MAX0(rms_wogA,rms_wogB);
	biggest_mograd=MAX0(biggest_mogradA,biggest_mogradB);
	rms_kappa=MAX0(rms_kappaA,rms_kappaB);
	biggest_kappa=MAX0(biggest_kappaA,biggest_kappaB);
	rms_t2=MAX0(rms_t2AA,rms_t2BB);
	rms_t2=MAX0(rms_t2,rms_t2AB);
	
fprintf(outfile," %3d     %12.10f  %12.2e   %12.2e     %12.2e    %12.2e  %12.2e  %12.2e \n",
	           itr_occ,Emp3L,DE,rms_wog,biggest_mograd,rms_kappa,biggest_kappa,rms_t2);
fflush(outfile);

/********************************************************************************************/
/********************************************************************************************/   
    if (itr_occ >= mo_maxiter) {
      break;  
      conver = 0; // means MOs are NOT optimized
    }

    if (rms_wog == 0.0) break;

}
while(fabs(DE) >= tol_Eod || rms_wog >= tol_grad || rms_kappa >= tol_grad || biggest_mograd >= mograd_max || 
      biggest_kappa >= mograd_max || rms_t2 >= tol_t2); 

if (conver == 1) {
timer_on("MP3L Energy");
mp3l_energy();
timer_off("MP3L Energy");

mo_optimized = 1; 
fprintf(outfile,"\n");
fprintf(outfile," ============================================================================== \n");
fprintf(outfile," ======================== OMP3 ITERATIONS ARE CONVERGED ======================= \n");
fprintf(outfile," ============================================================================== \n");
fprintf(outfile,"\n");
fflush(outfile);
}

else if (conver == 0) {
  fprintf(outfile,"\n ======================== OMP3 IS NOT CONVERGED IN %2d ITERATIONS ============= \n", mo_maxiter);
  fflush(outfile);
}

/********************************************************************************************/
/********************************************************************************************/

}
}} // End Namespaces

