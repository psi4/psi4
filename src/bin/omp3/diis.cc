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


using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace omp3wave{ 

void OMP3Wave::diis(int dimvec, double **vecs, double **errvecs, double *vec_new)
{ 
  
       double **Bmat, *Cvec, *errvec_new, det, scale_factor;

/********************************************************************************************/
/************************** memalloc ********************************************************/
/********************************************************************************************/ 
        Bmat = block_matrix(nvar, nvar);
        Cvec = init_array(nvar);
	errvec_new = init_array(dimvec);
	memset(Bmat[0], 0, sizeof(double)*nvar*nvar); 
	memset(Cvec, 0, sizeof(double)*nvar); 
	memset(errvec_new, 0, sizeof(double)*dimvec); 

/********************************************************************************************/
/************************** Form B matrix ***************************************************/
/********************************************************************************************/
	// num_vecs/num_vecs part
	for(int i = 0; i < num_vecs; i++){
	  for(int j = 0; j < num_vecs; j++){
	    Bmat[i][j] = dot_product(dimvec, &(errvecs[i][0]), &(errvecs[j][0]));
	  }
	}
	
	for(int i = 0; i < num_vecs; i++){
	  Bmat[nvar-1][i] = -1.0;
	  Bmat[i][nvar-1] = -1.0;
	}
	
	Bmat[nvar-1][nvar-1] = 0.0;
	
	
	// scale Bmat
	scale_factor = 1 / Bmat[0][0];
	for(int i = 0; i < nvar; i++){
	  for(int j = 0; j < nvar; j++){
	    Bmat[i][j] = scale_factor * Bmat[i][j];
	  }
	}

	// level shift
	if (level_shift == "TRUE") {
	  //for(int i = 0; i < num_vecs; i++) Bmat[i][i] -= lshift_parameter; // this also works
	  for(int i = 0; i < num_vecs; i++) Bmat[i][i] *= (1 + lshift_parameter); 
	}
	
/********************************************************************************************/
/************************** Form C vector ***************************************************/
/********************************************************************************************/
	for(int i = 0; i < num_vecs; i++) Cvec[i] = 0.0;
	Cvec[nvar-1] = -1.0;
	
/********************************************************************************************/
/************************** Solve LINEQ *****************************************************/
/********************************************************************************************/
       flin(Bmat, Cvec, nvar, 1, &det);         
       if (fabs(det) < DIIS_MIN_DET) { 
	fprintf(outfile, "Warning: diis matrix near-singular\n");
	fprintf(outfile, "Determinant is %6.3E\n", det);
       }     
       
/********************************************************************************************/
/************************** Extrapolate *****************************************************/
/********************************************************************************************/
	for(int i = 0; i < num_vecs; i++){
	  for(int j = 0; j < dimvec; j++){    
	    vec_new[j] += Cvec[i] * vecs[i][j];
	    errvec_new[j] += Cvec[i] * errvecs[i][j];
	  }  

	}

	free(Cvec);
	free(errvec_new);
	free_block(Bmat);

/********************************************************************************************/
/********************************************************************************************/

}
}} // End Namespaces

