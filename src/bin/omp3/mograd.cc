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
  

void OMP3Wave::mograd()
{
 
      double norm_wogA, norm_wogB;    
 
/********************************************************************************************/
/************************** memalloc ********************************************************/
/********************************************************************************************/
      WorbA->zero();
      WorbB->zero();

/********************************************************************************************/
/************************** wog0 ************************************************************/
/********************************************************************************************/
      SharedMatrix tempA(GFockA->transpose());
      WorbA->copy(GFockA);
      WorbA->subtract(tempA);       
      WorbA->scale(2.0);
      
      SharedMatrix tempB(GFockB->transpose());
      WorbB->copy(GFockB);
      WorbB->subtract(tempB);       
      WorbB->scale(2.0);

/********************************************************************************************/
/****************** set wog *****************************************************************/
/********************************************************************************************/
      // Alpha
      for(int x = 0; x < nidpA; x++) {
	int a = idprowA[x];
	int i = idpcolA[x];
	int h = idpirrA[x];
	wogA->set(x, WorbA->get(h, a + occpiA[h], i));
      }
      
      // Beta
      for(int x = 0; x < nidpB; x++) {
	int a = idprowB[x];
	int i = idpcolB[x];
	int h = idpirrB[x];
	wogB->set(x, WorbB->get(h, a + occpiB[h], i));
      }
      
    
/********************************************************************************************/
/************************** find biggest_mograd *********************************************/
/********************************************************************************************/          
    biggest_mogradA=0;
    for (int i=0; i<nidpA;i++){ 
      if (fabs(wogA->get(i)) > biggest_mogradA)  biggest_mogradA=fabs(wogA->get(i));
    }
    
    biggest_mogradB=0;
    for (int i=0; i<nidpB;i++){ 
      if (fabs(wogB->get(i)) > biggest_mogradB)  biggest_mogradB=fabs(wogB->get(i));
    }
      
/********************************************************************************************/
/************************** RMS *************************************************************/
/********************************************************************************************/	  
    rms_wogA=0;
    for (int i=0; i<nidpA;i++) rms_wogA += wogA->get(i) * wogA->get(i);
    norm_wogA=sqrt(rms_wogA);
    rms_wogA=sqrt(rms_wogA)/nidpA;  
    
    rms_wogB=0;
    for (int i=0; i<nidpB;i++) rms_wogB += wogB->get(i) * wogB->get(i);
    norm_wogB=sqrt(rms_wogB);
    rms_wogB=sqrt(rms_wogB)/nidpB;  
    
/********************************************************************************************/
/************************** print ***********************************************************/
/********************************************************************************************/
    if(print_ > 1){
      for(int i = 0; i < nidpA; i++){
        fprintf(outfile,"\n i, idpirrA, idprowA, idpcolA, wogA: %3d %3d %3d %3d %20.14f\n", i, idpirrA[i], idprowA[i],idpcolA[i],wogA->get(i)); 
	fflush(outfile);
      }
      
      for(int i = 0; i < nidpB; i++){
        fprintf(outfile,"\n i, idpirrB, idprowB, idpcolB, wogB: %3d %3d %3d %3d %20.14f\n", i, idpirrB[i], idprowB[i],idpcolB[i],wogB->get(i));
	fflush(outfile);
      }
    }

/********************************************************************************************/
/********************************************************************************************/

}// end of main
}} // End Namespaces

