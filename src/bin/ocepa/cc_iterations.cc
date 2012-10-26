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

#include "ocepawave.h"
#include "defines.h"
#include "arrays.h"

using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace ocepawave{
  
void OCEPAWave::cc_iterations()
{
  
fprintf(outfile,"\n  \n");      
fprintf(outfile," ============================================================================== \n");    
fprintf(outfile," ================ Performing CEPA iterations... =============================== \n");  
fprintf(outfile," ============================================================================== \n");
fprintf(outfile,"\n");
fprintf(outfile, "  Iter    E_corr           E_total            DE           T2 RMS        \n");
fprintf(outfile, "  ----   -------------    ---------------    ----------   ----------    \n");
fflush(outfile);

  
/********************************************************************************************/
/************************** NR iterations **************************************************/
/********************************************************************************************/
      itr_occ = 0;
      conver = 1; // Assuming that the iterations will converge
      

// head of loop      
do
{
        itr_occ++;
        timer_on("T2");
	t2_amps();  
        timer_off("T2");
        timer_on("CEPA Energy");
        cepa_energy();
        timer_off("CEPA Energy");
        EcepaL = Ecepa;
        EcorrL = EcepaL-Escf;
        DE = EcepaL - EcepaL_old;
        EcepaL_old = EcepaL;

    if (reference_ == "UNRESTRICTED") {
	rms_t2=MAX0(rms_t2AA,rms_t2BB);
	rms_t2=MAX0(rms_t2,rms_t2AB);
    }
	
fprintf(outfile," %3d     %12.10f    %12.10f  %12.2e %12.2e \n", itr_occ, Ecorr, Ecepa, DE, rms_t2);
fflush(outfile);

    if (itr_occ >= cc_maxiter) {
      conver = 0; // means iterations were NOT converged
      break;  
    }


}
while(fabs(DE) >= tol_Eod || rms_t2 >= tol_t2); 

if (conver == 1) {
fprintf(outfile,"\n");
fprintf(outfile," ============================================================================== \n");
fprintf(outfile," ======================== CEPA ITERATIONS ARE CONVERGED ======================= \n");
fprintf(outfile," ============================================================================== \n");
fflush(outfile);
}

else if (conver == 0) {
  fprintf(outfile,"\n ======================= CEPA IS NOT CONVERGED IN %2d ITERATIONS ============ \n", cc_maxiter);
  fflush(outfile);
}

}// end main
}} // End Namespaces

