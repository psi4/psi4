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

void OMP3Wave::diis(int dimvec, Array2d *vecs, Array2d *errvecs, Array1d *vec_new)
{ 

/********************************************************************************************/
/************************** memalloc ********************************************************/
/********************************************************************************************/ 
        Array2d *Bmat = new Array2d(nvar, nvar, "DIIS B Matrix"); 
        Array1d *Cvec = new Array1d(nvar, "DIIS C Vector"); 
        Array1d *errvec_new = new Array1d(dimvec, "New Error Vector");
        int d2 = errvecs->dim2();
        Array1d *vrow = new Array1d(d2); 
        Array1d *vcol = new Array1d(d2); 
        Bmat->zero(); 
        Cvec->zero(); 
        errvec_new->zero(); 
        vrow->zero(); 
        vcol->zero(); 

/********************************************************************************************/
/************************** Form B matrix ***************************************************/
/********************************************************************************************/
        // num_vecs/num_vecs part
	for(int i = 0; i < num_vecs; i++){
          vrow->row_vector(errvecs, i); 
	  for(int j = 0; j < num_vecs; j++){
            vcol->row_vector(errvecs, j); 
            double value = vrow->dot(vcol); 
            Bmat->set(i, j, value);
	  }
	}

	for(int i = 0; i < num_vecs; i++){
          Bmat->set(nvar - 1, i, -1.0);
          Bmat->set(i, nvar - 1, -1.0);
	}
	
        Bmat->set(nvar - 1, nvar - 1, -1.0);
	
	// scale Bmat
	double scale_factor = 1 / Bmat->get(0, 0);
        Bmat->scale(scale_factor);
         
	// level shift
	if (level_shift == "TRUE") {
	  //for(int i = 0; i < num_vecs; i++) Bmat->set(i, i, -lshift_parameter);// this also works
	  for(int i = 0; i < num_vecs; i++) Bmat->set(i, i, Bmat->get(i, i) * (1 + lshift_parameter) ); 
	}
	
/********************************************************************************************/
/************************** Form C vector ***************************************************/
/********************************************************************************************/
	Cvec->set(nvar - 1,  -1.0);
	
/********************************************************************************************/
/************************** Solve LINEQ *****************************************************/
/********************************************************************************************/
         if (lineq == "CDGESV") Bmat->cdgesv(Cvec);
 
         else if (lineq == "FLIN") {
            double det = 0.0;      
            Bmat->lineq_flin(Cvec, &det);
            if (fabs(det) < DIIS_MIN_DET) { 
               fprintf(outfile, "Warning: diis matrix near-singular\n");
               fprintf(outfile, "Determinant is %6.3E\n", det);
            }
          }
 
          else if (lineq == "POPLE") Bmat->lineq_pople(Cvec, num_vecs, cutoff);

             
       
/********************************************************************************************/
/************************** Extrapolate *****************************************************/
/********************************************************************************************/
	for(int i = 0; i < num_vecs; i++){
	  for(int j = 0; j < dimvec; j++){    
	    vec_new->set(j, vec_new->get(j) + (Cvec->get(i) * vecs->get(i,j) ) );
	    errvec_new->set(j, errvec_new->get(j) + ( Cvec->get(i) * errvecs->get(i,j)));
	  }  

	}

	delete Bmat;
	delete Cvec;
	delete errvec_new;
	delete vrow;
	delete vcol;

/********************************************************************************************/
/********************************************************************************************/

}
}} // End Namespaces

