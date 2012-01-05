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

#include "omp2wave.h"
#include "defines.h"


using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace omp2wave{ 

void OMP2Wave::korbrot_sd()
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
	  kappaA[x] = -wogA[x] / (2*value);
	}
	
	// beta
	for(int x = 0; x < nidpB; x++) {
	  int a = idprowB[x];
	  int i = idpcolB[x];
	  int h = idpirrB[x];
	  double value = FockB->get(h, a + occpiB[h], a + occpiB[h]) - FockB->get(h, i, i);  
	  kappaB[x] = -wogB[x] / (2*value);
	}

/********************************************************************************************/
/************************** find biggest_kappa ***********************************************/
/********************************************************************************************/
	biggest_kappaA=0;            
	for (int i=0; i<nidpA;i++) 
	{ 
	      if (fabs(kappaA[i]) > biggest_kappaA)
	      {
		  biggest_kappaA=fabs(kappaA[i]);
	      }
	}
	
	biggest_kappaB=0;            
	for (int i=0; i<nidpB;i++) 
	{ 
	      if (fabs(kappaB[i]) > biggest_kappaB)
	      {
		  biggest_kappaB=fabs(kappaB[i]);
	      }
	}
	
/********************************************************************************************/
/************************** scale array *****************************************************/
/********************************************************************************************/	 	
	 if (biggest_kappaA > step_max)
	 {   
	      for (int i=0; i<nidpA;i++) 
	      { 
		  kappaA[i]=kappaA[i]*(step_max/biggest_kappaA);
	      }
	 }
	 
	 if (biggest_kappaB > step_max)
	 {   
	      for (int i=0; i<nidpB;i++) 
	      { 
		  kappaB[i]=kappaB[i]*(step_max/biggest_kappaB);
	      }
	 }
	 
/********************************************************************************************/
/************************** find biggest_kappa after scaling *********************************/
/********************************************************************************************/
	if (biggest_kappaA > step_max)
	{
	  biggest_kappaA=0;            
	  for (int i=0; i<nidpA;i++) 
	  { 
	      if (fabs(kappaA[i]) > biggest_kappaA)
	      {
		  biggest_kappaA=fabs(kappaA[i]);
	      }
	  }
	}
	
	if (biggest_kappaB > step_max)
	{
	  biggest_kappaB=0;            
	  for (int i=0; i<nidpB;i++) 
	  { 
	      if (fabs(kappaB[i]) > biggest_kappaB)
	      {
		  biggest_kappaB=fabs(kappaB[i]);
	      }
	  }
	}

/********************************************************************************************/
/************************** norm ************************************************************/
/********************************************************************************************/	 	
	rms_kappaA=0;
	rms_kappaB=0;
	
	for (int i=0; i<nidpA;i++) rms_kappaA+=kappaA[i]*kappaA[i];
	for (int i=0; i<nidpB;i++) rms_kappaB+=kappaB[i]*kappaB[i];
	rms_kappaA=sqrt(rms_kappaA)/nidpA;
	rms_kappaB=sqrt(rms_kappaB)/nidpB;
	
	
      if(print_ > 1){
	for(int i = 0; i<nidpA; i++){
	  fprintf(outfile,"\n idpA, idprowA, idpcolA, kappaA: %3d %3d %3d %20.14f\n", i, idprowA[i],idpcolA[i],kappaA[i]);
	  fflush(outfile);
	}

	for(int i = 0; i<nidpB; i++){
	  fprintf(outfile,"\n idpB, idprowB, idpcolB, kappaB: %3d %3d %3d %20.14f\n", i, idprowB[i],idpcolB[i],kappaB[i]); 
	  fflush(outfile);
	}
      }
      
    //fprintf(outfile,"\n korbrot_sd done. \n"); fflush(outfile);

	
/********************************************************************************************/	
/********************************************************************************************/	
}
}} // End Namespaces

