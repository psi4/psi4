#include "occwave.h"
#include "defines.h"


using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace occwave{ 

void OCCWave::diis(int dimvec, Array2d *vecs, Array2d *errvecs, Array1d *vec_new, Array1d *errvec_new)
{ 

/********************************************************************************************/
/************************** memalloc ********************************************************/
/********************************************************************************************/ 
        Array2d *Bmat = new Array2d("DIIS B Matrix", nvar, nvar); 
        Array1d *Cvec = new Array1d("DIIS C Vector", nvar); 
        Array1d *vrow = new Array1d(dimvec); 
        Array1d *vcol = new Array1d(dimvec); 
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
	
        Bmat->set(nvar - 1, nvar - 1, 0.0);
	
	// scale Bmat
	//double scale_factor = 1 / Bmat->get(0, 0);
        //Bmat->scale(scale_factor);
         
	// level shift
	if (level_shift == "TRUE") {
	  //for(int i = 0; i < num_vecs; i++) Bmat->set(i, i, -lshift_parameter);// this also works
	  for(int i = 0; i < num_vecs; i++) Bmat->set(i, i, Bmat->get(i, i) * (1 + lshift_parameter) ); 
	}
	
        // Form the c vector
	Cvec->set(nvar - 1,  -1.0);
	
/********************************************************************************************/
/************************** Solve LINEQ *****************************************************/
/********************************************************************************************/
         if (lineq == "CDGESV") Bmat->cdgesv(Cvec);
         else if (lineq == "FLIN") {
            double det = 0.0;      
            Bmat->lineq_flin(Cvec, &det);
            if (fabs(det) < DIIS_MIN_DET) { 
               fprintf(outfile, "Warning!!! Diis matrix is near-singular\n");
               fprintf(outfile, "Determinant is %6.3E\n", det);
               fflush(outfile);
            }
          }
          else if (lineq == "POPLE") Bmat->lineq_pople(Cvec, num_vecs, cutoff);
       
/********************************************************************************************/
/************************** Extrapolate *****************************************************/
/********************************************************************************************/
	for(int i = 0; i < dimvec; i++){
          double sum1 = 0.0;
          double sum2 = 0.0;
	  for(int j = 0; j < num_vecs; j++){    
	    sum1 += Cvec->get(j) * vecs->get(j,i);
            sum2 += Cvec->get(j) * errvecs->get(j,i);
	  }  
          vec_new->set(i, sum1);
          errvec_new->set(i, sum2);
	}
 

	delete Bmat;
	delete Cvec;
	delete vrow;
	delete vcol;

/********************************************************************************************/
/********************************************************************************************/

}
}} // End Namespaces

