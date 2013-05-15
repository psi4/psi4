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

#include "occwave.h"

using namespace boost;
using namespace psi;
using namespace std;

namespace psi{ namespace occwave{

void OCCWave::idp()
{
     int dim;

if (reference_ == "RESTRICTED") {
    // Form IDPs
    nidpA=0;

    // V-O: I exclude symmetry broken rotations from the list of IDPs since they already have zero gradient.
    for(int h = 0; h < nirrep_; h++){
      nidpA += virtpiA[h] * occpiA[h]; 
    }

    fprintf(outfile,"\n\tNumber of independent-pairs: %3d\n", nidpA);
    fflush(outfile);  
    
    if (nidpA != 0) {
      idp_returnA = 1;
      wogA = new Array1d("Alpha MO grad vector", nidpA);
      kappaA = new Array1d("Alpha orb rot params vector of current step", nidpA);
      kappa_newA = new Array1d("Alpha New orb rot params vector of current step", nidpA);
      kappa_barA = new Array1d("Alpha orb rot params vector with respect to scf MOs", nidpA);
      wog_intA = new Array1d("Alpha Interpolated MO grad vector", nidpA);
      wogA->zero();
      kappaA->zero();
      kappa_barA->zero();
    }
    
    // allocate memory 
    idprowA = new int[nidpA]; 
    idpcolA = new int[nidpA];
    idpirrA = new int[nidpA]; 
    
    // initialize 
    memset(idprowA,0, sizeof(int)*nidpA);
    memset(idpcolA,0, sizeof(int)*nidpA);
    memset(idpirrA,0, sizeof(int)*nidpA);   

    // set idpA 
    dim=0;
    for(int h = 0; h < nirrep_; h++){      
      for(int a = 0; a < virtpiA[h]; a++){
	for(int i = 0; i < occpiA[h]; i++){
	  idprowA[dim]=a;
	  idpcolA[dim]=i;
	  idpirrA[dim]=h;
	  dim++;  
	}
      }
    }

    if(print_ > 2){
     for(int i = 0; i < nidpA; i++){
        fprintf(outfile,"\n i, idpirrA, idprowA, idpcolA: %3d %3d %3d %3d\n", i, idpirrA[i], idprowA[i],idpcolA[i]);
	fflush(outfile);
      }
    }
     
}// end if (reference_ == "RESTRICTED") 

else if (reference_ == "UNRESTRICTED") {
    // Form IDPs
    nidpA=0;
    nidpB=0;

    // V-O: I exclude symmetry broken rotations from the list of IDPs since they already have zero gradient.
    for(int h = 0; h < nirrep_; h++){
      nidpA += virtpiA[h] * occpiA[h]; 
      nidpB += virtpiB[h] * occpiB[h]; 
    }

    fprintf(outfile,"\n\tNumber of alpha independent-pairs:%3d\n", nidpA);
    fprintf(outfile,"\tNumber of beta independent-pairs :%3d\n", nidpB);
    fflush(outfile);  
    
    if (nidpA != 0) {
      idp_returnA = 1;
      wogA = new Array1d("Alpha MO grad vector", nidpA);
      kappaA = new Array1d("Alpha orb rot params vector of current step", nidpA);
      kappa_newA = new Array1d("Alpha New orb rot params vector of current step", nidpA);
      kappa_barA = new Array1d("Alpha orb rot params vector with respect to scf MOs", nidpA);
      wog_intA = new Array1d("Alpha Interpolated MO grad vector", nidpA);
      wogA->zero();
      kappaA->zero();
      kappa_barA->zero();
    }
    
    if (nidpB != 0) {
      idp_returnB = 1;
      wogB = new Array1d("Beta MO grad vector", nidpB);
      kappaB = new Array1d("Beta orb rot params vector of current step", nidpB);
      kappa_newB = new Array1d("Beta New orb rot params vector of current step", nidpB);
      kappa_barB = new Array1d("Beta orb rot params vector with respect to scf MOs", nidpB);
      wog_intB = new Array1d("Beta Interpolated MO grad vector", nidpB);
      wogB->zero();
      kappaB->zero();
      kappa_barB->zero();
    }
 
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
    for(int h = 0; h < nirrep_; h++){      
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
    for(int h = 0; h < nirrep_; h++){
      for(int a = 0; a < virtpiB[h]; a++){
	for(int i = 0; i < occpiB[h]; i++){
	  idprowB[dim]=a;
	  idpcolB[dim]=i;
	  idpirrB[dim]=h;
	  dim++;  
	}
      }
    }
    
    if(print_ > 2){
     for(int i = 0; i < nidpA; i++){
        fprintf(outfile,"\n i, idpirrA, idprowA, idpcolA: %3d %3d %3d %3d\n", i, idpirrA[i], idprowA[i],idpcolA[i]);
	fflush(outfile);
      }
      
      for(int i = 0; i < nidpB; i++){
        fprintf(outfile,"\n i, idpirrB, idprowB, idpcolB: %3d %3d %3d %3d\n", i, idpirrB[i], idprowB[i],idpcolB[i]); 
	fflush(outfile);
      }
    }
      
}// end if (reference_ == "UNRESTRICTED") 

}// end of main
}} // End Namespaces


