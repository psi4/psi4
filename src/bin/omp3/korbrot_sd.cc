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


/** Required libmints includes */
#include <libmints/mints.h>
#include <libmints/factory.h>
#include <libmints/wavefunction.h>

#include "omp3wave.h"
#include "defines.h"
#include "arrays.h"


using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace omp3wave{ 

void OMP3Wave::korbrot_sd()
{ 
      //fprintf(outfile,"\n korbrot_sd is starting... \n"); fflush(outfile);

/********************************************************************************************/
/************************** initialize to -grad *********************************************/
/********************************************************************************************/ 
	// alpha
	for(int x = 0; x < nidpA; x++) {
	  int a = idprowA[x];
	  int i = idpcolA[x];
	  int h = idpirrA[x];
	  double value = FockA->get(h, a + occpiA[h], a + occpiA[h]) - FockA->get(h, i, i);  
	  kappaA->set(x, -wogA->get(x) / (2*value));
	}
	
	// beta
	for(int x = 0; x < nidpB; x++) {
	  int a = idprowB[x];
	  int i = idpcolB[x];
	  int h = idpirrB[x];
	  double value = FockB->get(h, a + occpiB[h], a + occpiB[h]) - FockB->get(h, i, i);  
	  kappaB->set(x, -wogB->get(x) / (2*value));
	}


/********************************************************************************************/
/************************** find biggest_kappa ***********************************************/
/********************************************************************************************/
	biggest_kappaA=0;            
	for (int i=0; i<nidpA;i++) { 
	    if (fabs(kappaA->get(i)) > biggest_kappaA) biggest_kappaA=fabs(kappaA->get(i));
	}
	
	biggest_kappaB=0;            
	for (int i=0; i<nidpB;i++){ 
	    if (fabs(kappaB->get(i)) > biggest_kappaB) biggest_kappaB=fabs(kappaB->get(i));
	}
	
/********************************************************************************************/
/************************** scale array *****************************************************/
/********************************************************************************************/	 	
	 if (biggest_kappaA > step_max) {   
	    for (int i=0; i<nidpA;i++) kappaA->set(i, kappaA->get(i) *(step_max/biggest_kappaA));
	 }
	 
	 if (biggest_kappaB > step_max) {   
	    for (int i=0; i<nidpB;i++) kappaB->set(i, kappaB->get(i) *(step_max/biggest_kappaB));
	 }
	 
/********************************************************************************************/
/************************** find biggest_kappa after scaling *********************************/
/********************************************************************************************/
	if (biggest_kappaA > step_max)
	{
	  biggest_kappaA=0;            
	  for (int i=0; i<nidpA;i++) 
	  { 
	      if (fabs(kappaA->get(i)) > biggest_kappaA)
	      {
		  biggest_kappaA = fabs(kappaA->get(i));
	      }
	  }
	}
	
	if (biggest_kappaB > step_max)
	{
	  biggest_kappaB=0;            
	  for (int i=0; i<nidpB;i++) 
	  { 
	      if (fabs(kappaB->get(i)) > biggest_kappaB)
	      {
		  biggest_kappaB=fabs(kappaB->get(i));
	      }
	  }
	}

/********************************************************************************************/
/************************** norm ************************************************************/
/********************************************************************************************/	 	
	rms_kappaA=0;
	rms_kappaB=0;
	rms_kappaA = kappaA->rms();
	rms_kappaB = kappaB->rms();
	
      if(print_ > 1){
       kappaA->print();
       kappaB->print();
      }
      
    //fprintf(outfile,"\n korbrot_sd done. \n"); fflush(outfile);

	
/********************************************************************************************/	
/********************************************************************************************/	
}
}} // End Namespaces

