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

void OMP3Wave::idp()
{
     int dim;
     
/********************************************************************************************/
/************************** Form IDPs *******************************************************/
/********************************************************************************************/
    nidpA=0;
    nidpB=0;

    // V-O: I exclude symmetry broken rotations from the list of IDPs since they already have zero gradient.
    for(int h = 0; h < nirreps; h++){
      nidpA += virtpiA[h] * occpiA[h]; 
      nidpB += virtpiB[h] * occpiB[h]; 
    }

    fprintf(outfile,"\n\tNumber of alpha independent-pairs:%3d\n", nidpA);
    fprintf(outfile,"\n\tNumber of beta independent-pairs :%3d\n", nidpB);
    fflush(outfile);  
    
    if (nidpA != 0) {
      idp_returnA = 1;
      wogA = new Array1d(nidpA, "Alpha MO grad vector");
      kappaA = new Array1d(nidpA, "Alpha orb rot params vector of current step");
      kappa_barA = new Array1d(nidpA, "Alpha orb rot params vector with respect to scf MOs");
      wogA->zero();
      kappaA->zero();
      kappa_barA->zero();
    }
    
    if (nidpB != 0) {
      idp_returnB = 1;
      wogB = new Array1d(nidpB, "Beta MO grad vector");
      kappaB = new Array1d(nidpB, "Beta orb rot params vector of current step");
      kappa_barB = new Array1d(nidpB, "Beta orb rot params vector with respect to scf MOs");
      wogB->zero();
      kappaB->zero();
      kappa_barB->zero();
    }
 
/********************************************************************************************/
/************************** form idprow & idpcol vectors ************************************/
/********************************************************************************************/
    // allocate memory 
    idprowA = new int[nidpA]; 
    idpcolA = new int[nidpA];
    idpirrA = new int[nidpA]; 
    idprowB = new int[nidpB]; 
    idpcolB = new int[nidpB]; 
    idpirrB = new int[nidpB]; 
    
    // initialize 
    memset(idprowA,0, sizeof(int)*nidpA);
    memset(idpcolA,0, sizeof(int)*nidpA);
    memset(idpirrA,0, sizeof(int)*nidpA);   
    memset(idprowB,0, sizeof(int)*nidpB);
    memset(idpcolB,0, sizeof(int)*nidpB);    
    memset(idpirrB,0, sizeof(int)*nidpB);

    // set idpA 
    dim=0;
    for(int h = 0; h < nirreps; h++){      
      for(int a = 0; a < virtpiA[h]; a++){
	for(int i = 0; i < occpiA[h]; i++){
	  idprowA[dim]=a;
	  idpcolA[dim]=i;
	  idpirrA[dim]=h;
	  dim++;  
	}
      }
    }
    
    // set idpB 
    dim=0;
    for(int h = 0; h < nirreps; h++){
      for(int a = 0; a < virtpiB[h]; a++){
	for(int i = 0; i < occpiB[h]; i++){
	  idprowB[dim]=a;
	  idpcolB[dim]=i;
	  idpirrB[dim]=h;
	  dim++;  
	}
      }
    }
    
    if(print_ > 1){
     for(int i = 0; i < nidpA; i++){
        fprintf(outfile,"\n i, idpirrA, idprowA, idpcolA: %3d %3d %3d %3d\n", i, idpirrA[i], idprowA[i],idpcolA[i]);
	fflush(outfile);
      }
      
      for(int i = 0; i < nidpB; i++){
        fprintf(outfile,"\n i, idpirrB, idprowB, idpcolB: %3d %3d %3d %3d\n", i, idpirrB[i], idprowB[i],idpcolB[i]); 
	fflush(outfile);
      }
    }
      

/********************************************************************************************/
/********************************************************************************************/  

}// end of main
}} // End Namespaces


