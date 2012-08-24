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
#include <boost/shared_ptr.hpp>
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h> 
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <libdiis/diismanager.h>


/** Required libmints includes */
#include <libmints/mints.h>
#include <libmints/factory.h>
#include <libmints/wavefunction.h>

#include "omp2wave.h"
#include "defines.h"
#include "arrays.h"

using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace omp2wave{

void OMP2Wave::update_mo()
{
      //fprintf(outfile,"\n update_mo is starting... \n"); fflush(outfile);     

/********************************************************************************************/
/************************** initialize array ************************************************/
/********************************************************************************************/	
	UorbA->zero();
	UorbB->zero();
	KorbA->zero();
	KorbB->zero();
	
/********************************************************************************************/
/************************** Build kappa_bar *************************************************/
/********************************************************************************************/ 
        kappa_barA->add(kappaA);        
        kappa_barB->add(kappaB);        

/********************************************************************************************/
/************************ DO DIIS ***********************************************************/
/********************************************************************************************/
if (opt_method == "DIIS") {
  
        // Form Diis Error Vector & Extrapolant Alpha Spin Case
	if (itr_occ <= num_vecs) {  
	  for(int i = 0; i < nidpA; i++){
	    errvecsA->set(itr_occ-1, i, wogA->get(i));
	    vecsA->set(itr_occ-1, i, kappa_barA->get(i));
	  }  
	}
	
	
	if (itr_occ > num_vecs) {  
	  for(int j = 0; j < (num_vecs-1); j++){
	    for(int i = 0; i < nidpA; i++){
	      errvecsA->set(j, i, errvecsA->get(j+1, i));
	      vecsA->set(j, i, vecsA->get(j+1, i));
	    }  
	  }
	  
	  for(int i = 0; i < nidpA; i++){
	    errvecsA->set(num_vecs-1, i, wogA->get(i));
	    vecsA->set(num_vecs-1, i, kappa_barA->get(i));
	  }    
	}
	
        // Form Diis Error Vector & Extrapolant Beta Spin Case
	if (itr_occ <= num_vecs) {  
	  for(int i = 0; i < nidpB; i++){
	    errvecsB->set(itr_occ-1, i, wogB->get(i));
	    vecsB->set(itr_occ-1, i, kappa_barB->get(i));
	  }  
	}
	
	
	if (itr_occ > num_vecs) {  
	  for(int j = 0; j < (num_vecs-1); j++){
	    for(int i = 0; i < nidpB; i++){
	      errvecsB->set(j, i, errvecsB->get(j+1, i));
	      vecsB->set(j, i, vecsB->get(j+1, i));
	    }  
	  }
	  
	  for(int i = 0; i < nidpB; i++){
	    errvecsB->set(num_vecs-1, i, wogB->get(i));
	    vecsB->set(num_vecs-1, i, kappa_barB->get(i));
	  }    
	}
	
	
        // Extrapolate 
        if (itr_occ >= num_vecs) {
	  diis(nidpA, vecsA, errvecsA, kappa_barA);
	  diis(nidpB, vecsB, errvecsB, kappa_barB);
	}
	
}// end if (opt_method == "DIIS") {

/********************************************************************************************/
/************************** Construct Korb **************************************************/
/********************************************************************************************/
	// alpha
	for(int x = 0; x < nidpA; x++) {
	  int a = idprowA[x];
	  int i = idpcolA[x];
	  int h = idpirrA[x];
	  KorbA->set(h, a + occpiA[h], i, kappa_barA->get(x));
	  KorbA->set(h, i, a + occpiA[h], -kappa_barA->get(x));
	}
	
	// beta
	for(int x = 0; x < nidpB; x++) {
	  int a = idprowB[x];
	  int i = idpcolB[x];
	  int h = idpirrB[x];
	  KorbB->set(h, a + occpiB[h], i, kappa_barB->get(x));
	  KorbB->set(h, i, a + occpiB[h], -kappa_barB->get(x));
	}
	
/********************************************************************************************/
/************************** Construct Uorb **************************************************/
/********************************************************************************************/	
	//set to identity
	UorbA->identity();
	UorbB->identity();
	
	// K contribution
	UorbA->add(KorbA);
	UorbB->add(KorbB);
	
	//form K^2
	KsqrA->gemm(false, false, 1.0, KorbA, KorbA, 0.0); 
	KsqrB->gemm(false, false, 1.0, KorbB, KorbB, 0.0); 
	KsqrA->scale(0.5);
	KsqrB->scale(0.5);
	
	// 0.5*K^2 contribution
	UorbA->add(KsqrA);
	UorbB->add(KsqrB);

/********************************************************************************************/
/************************** Orthogonalize U matrix ******************************************/
/********************************************************************************************/
if (orth_type == "MGS") {;
    double rmgs1a,rmgs2a,rmgs1b,rmgs2b;
    
    // loop-over nirreps
    for (int h=0; h<nirreps; h++) {
      
      // loop-1
      for (int k = 0; k < mopi[h]; k++) {
	rmgs1a=0.0;
	rmgs1b=0.0;
	
	// loop-1a
	for (int i=0; i < mopi[h]; i++) {  
	  rmgs1a += UorbA->get(h, i, k) * UorbA->get(h, i, k);
	  rmgs1b += UorbB->get(h, i, k) * UorbB->get(h, i, k);
	}// end 1a
	
	rmgs1a=sqrt(rmgs1a);
	rmgs1b=sqrt(rmgs1b);
	  
	// loop-1b
	for (int i=0; i < mopi[h]; i++) {  
	  UorbA->set(h, i, k, UorbA->get(h, i, k) / rmgs1a);
	  UorbB->set(h, i, k, UorbB->get(h, i, k) / rmgs1b);
	}// end 1b
	
	// loop-2
	for (int j=(k+1); j < mopi[h]; j++) {
	  rmgs2a=0; 
	  rmgs2b=0; 
	  
	  // loop-2a
	  for (int i=0; i < mopi[h]; i++) {  
	    rmgs2a += UorbA->get(h, i, k) * UorbA->get(h, i, j);
	    rmgs2b += UorbB->get(h, i, k) * UorbB->get(h, i, j);
	  }// end 2a
	  
	  // loop-2b
	  for (int i=0; i < mopi[h]; i++) {  
	    UorbA->set(h, i, j, UorbA->get(h, i, j) - (rmgs2a * UorbA->get(h, i, k)));
	    UorbB->set(h, i, j, UorbB->get(h, i, j) - (rmgs2b * UorbB->get(h, i, k)));
	  }// end 2b
	  
	}// end 2
      }// end 1
    }// end loop-over nirreps
}// end main if


else if (orth_type == "GS") {
    int rowA = UorbA->nrow();
    int colA = UorbA->ncol();
    
    double **AdumA = block_matrix(rowA, colA);
    memset(AdumA[0], 0, sizeof(double)*rowA*colA);
    AdumA = UorbA->to_block_matrix();    
    schmidt(AdumA, rowA, colA, outfile);  
    UorbA->set(AdumA);    
    free_block(AdumA);
    
    int rowB = UorbB->nrow();
    int colB = UorbB->ncol();
    
    double **AdumB = block_matrix(rowB, colB);
    memset(AdumB[0], 0, sizeof(double)*rowB*colB);
    AdumB = UorbB->to_block_matrix();    
    schmidt(AdumB, rowB, colB, outfile);  
    UorbB->set(AdumB);    
    free_block(AdumB);
}
   
/********************************************************************************************/
/************************** Build new MO coeff. *********************************************/
/********************************************************************************************/
	Ca_->gemm(false, false, 1.0, Ca_ref, UorbA, 0.0); 
	Cb_->gemm(false, false, 1.0, Cb_ref, UorbB, 0.0); 

       	if (print_ > 1) {
	  UorbA->print();
          UorbB->print();
	  Ca_->print();
	  Cb_->print();
	}

/********************************************************************************************/
/************************** Save MO coefficients to Chkpt file ******************************/
/********************************************************************************************/	 
	C_pitzerA = Ca_->to_block_matrix();    
	C_pitzerB = Cb_->to_block_matrix();    
	chkpt_->wt_alpha_scf(C_pitzerA);
	chkpt_->wt_beta_scf(C_pitzerB);
	free_block(C_pitzerA);
	free_block(C_pitzerB);

/********************************************************************************************/
/********************************************************************************************/	

}
}} // End Namespaces

