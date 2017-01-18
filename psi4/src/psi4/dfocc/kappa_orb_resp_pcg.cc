/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

#include "defines.h"
#include "dfocc.h"
#include "psi4/psi4-dec.h"

using namespace psi;
using namespace std;

namespace psi{ namespace dfoccwave{

void DFOCC::kappa_orb_resp_pcg()
{ 
//outfile->Printf("\n kappa_orb_resp_pcg is starting... \n"); 

    SharedTensor2d K, L;

if (reference_ == "RESTRICTED") {
    // Memalloc
    zvectorA = SharedTensor1d(new Tensor1d("Alpha Z-Vector", noccA * nvirA));
    zvec_newA = SharedTensor1d(new Tensor1d("Alpha New Z-Vector", noccA * nvirA));
    Minv_pcgA = SharedTensor1d(new Tensor1d("Alpha PCG M inverse", noccA * nvirA));
    sigma_pcgA = SharedTensor1d(new Tensor1d("Alpha PCG sigma", noccA * nvirA));
    r_pcgA = SharedTensor1d(new Tensor1d("Alpha PCG r", noccA * nvirA));
    r_pcg_newA = SharedTensor1d(new Tensor1d("Alpha PCG new r", noccA * nvirA));
    z_pcgA = SharedTensor1d(new Tensor1d("Alpha PCG z", noccA * nvirA));
    z_pcg_newA = SharedTensor1d(new Tensor1d("Alpha PCG new z", noccA * nvirA));
    p_pcgA = SharedTensor1d(new Tensor1d("Alpha PCG p", noccA * nvirA));
    p_pcg_newA = SharedTensor1d(new Tensor1d("Alpha PCG new p", noccA * nvirA));
    dr_pcgA = SharedTensor1d(new Tensor1d("Alpha PCG dr", noccA * nvirA));
    residualA = SharedTensor1d(new Tensor1d("Alpha Residual Vector", noccA * nvirA));

    // Build kappa0 and M 
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
              double value = FockA->get(a + noccA, a + noccA) - FockA->get(i,i);
              zvectorA->set(ai, -WorbA->get(a + noccA, i) / (2.0*value));       
              Minv_pcgA->set(ai, 0.5/value);
         }
    }

    // Build S = A kappa_0
    sigma_rhf(sigma_pcgA, zvectorA);

    // Level Shift
    if (level_shift == "TRUE") sigma_pcgA->axpy(zvectorA, lshift_parameter);

    // Build r0
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
              residualA->set(ai, -WorbA->get(a + noccA, i) - sigma_pcgA->get(ai));       
         }
    }
    r_pcgA->copy(residualA);

    // Build z0
    z_pcgA->dirprd(Minv_pcgA, r_pcgA);

    // Build p0
    p_pcgA->copy(z_pcgA);

    // Call Orbital Response Solver
    orb_resp_pcg_rhf();

    // Memfree
    zvec_newA.reset(); 
    Minv_pcgA.reset(); 
    sigma_pcgA.reset();
    r_pcgA.reset();
    r_pcg_newA.reset();
    z_pcgA.reset();
    z_pcg_newA.reset();
    p_pcgA.reset();
    p_pcg_newA.reset();
    dr_pcgA.reset();
    residualA.reset();

    // Build kappa for VO block
    #pragma omp parallel for
    for (int x = 0; x < nidpA; x++) {
         int p = idprowA->get(x);
	 int q = idpcolA->get(x);
         if (p >= noccA && q < noccA) {
             int ai = vo_idxAA->get(p-noccA,q);
	     kappaA->set(x, zvectorA->get(ai));
         }
    }
    zvectorA.reset();

    // OO Block
    if (nfrzc > 0) {
      // Compute OO-Block orb rot params
      approx_diag_mohess_oo();
      #pragma omp parallel for
      for (int x = 0; x < nidpA; x++) {
	   int p = idprowA->get(x);
	   int q = idpcolA->get(x);
           if (p < noccA && q < noccA) {
               double value = AooA->get(p-nfrzc,q); 
	       kappaA->set(x, -wogA->get(x)/value);
           }
      }
    } // if (nfrzc > 0) 

    // If LINEQ FAILED!
    if (pcg_conver == 0) {
        //if (print_ > 1 ) {
        outfile->Printf("\tWarning!!! PCG did NOT converged in %2d iterations, switching to an approximately diagonal MO Hessian. \n", itr_pcg);
        
        //}

        // Compute VO-Block Hess
        approx_diag_mohess_vo();

        // Build kappa again
        #pragma omp parallel for
        for (int x = 0; x < nidpA; x++) {
	    int p = idprowA->get(x);
	    int q = idpcolA->get(x);
            double value = 0.0;
            if (p >= noccA && q < noccA) value = AvoA->get(p-noccA,q); 
            else if (p < noccA && q < noccA) value = AooA->get(p-nfrzc,q); 
	    kappaA->set(x, -wogA->get(x)/value);
        }
    } // end if pcg_conver = 0

        // find biggest_kappa 
	biggest_kappaA=0;            
	for (int i=0; i<nidpA;i++) { 
	    if (fabs(kappaA->get(i)) > biggest_kappaA) biggest_kappaA=fabs(kappaA->get(i));
	}

        // Scale
	if (biggest_kappaA > step_max) {   
	    for (int i=0; i<nidpA;i++) kappaA->set(i, kappaA->get(i) *(step_max/biggest_kappaA));
	}
	 
        // find biggest_kappa again 
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
	
        // norm
	rms_kappaA=0;
	rms_kappaA = kappaA->rms();
	
        // print
        if(print_ > 2) kappaA->print();

}// end if (reference_ == "RESTRICTED") 

else if (reference_ == "UNRESTRICTED") {
    nidp_tot = nidpA + nidpB;

    // Memalloc
    zvector = SharedTensor1d(new Tensor1d("UHF Z-Vector", nidp_tot));
    zvec_new = SharedTensor1d(new Tensor1d("New UHF Z-Vector", nidp_tot));
    Minv_pcg = SharedTensor1d(new Tensor1d("PCG M inverse", nidp_tot));
    sigma_pcg = SharedTensor1d(new Tensor1d("PCG sigma", nidp_tot));
    r_pcg = SharedTensor1d(new Tensor1d("PCG r", nidp_tot));
    r_pcg_new = SharedTensor1d(new Tensor1d("PCG new r", nidp_tot));
    z_pcg = SharedTensor1d(new Tensor1d("PCG z", nidp_tot));
    z_pcg_new = SharedTensor1d(new Tensor1d("PCG new z", nidp_tot));
    p_pcg = SharedTensor1d(new Tensor1d("PCG p", nidp_tot));
    p_pcg_new = SharedTensor1d(new Tensor1d("PCG new p", nidp_tot));
    dr_pcg = SharedTensor1d(new Tensor1d("PCG dr", nidp_tot));
    residual = SharedTensor1d(new Tensor1d("Residual Vector", nidp_tot));
    zvectorA = SharedTensor1d(new Tensor1d("Alpha Z-Vector", noccA * nvirA));
    zvectorB = SharedTensor1d(new Tensor1d("Beta Z-Vector", noccB * nvirB));
    sigma_pcgA = SharedTensor1d(new Tensor1d("Alpha PCG sigma", noccA * nvirA));
    sigma_pcgB = SharedTensor1d(new Tensor1d("Beta PCG sigma", noccB * nvirB));
    p_pcgA = SharedTensor1d(new Tensor1d("Alpha PCG p", noccA * nvirA));
    p_pcgB = SharedTensor1d(new Tensor1d("Beta PCG p", noccB * nvirB));

    // Build kappa0 and M
    // alpha
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
              double value = FockA->get(a + noccA, a + noccA) - FockA->get(i,i);
              zvectorA->set(ai, -WorbA->get(a + noccA, i) / (2.0*value));       
              Minv_pcg->set(ai, 0.5/value);
         }
    }

    // beta
    for (int a = 0, ai = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++, ai++) {
              double value = FockB->get(a + noccB, a + noccB) - FockB->get(i,i);
              zvectorB->set(ai, -WorbB->get(a + noccB, i) / (2.0*value));       
              Minv_pcg->set(ai + nidpA, 0.5/value);
         }
    }

    // Form initial zvector vector
    for (int ai = 0; ai < nidpA; ai++) zvector->set(ai, zvectorA->get(ai));
    for (int ai = 0; ai < nidpB; ai++) zvector->set(ai + nidpA, zvectorB->get(ai));

    // Build S = A kappa_0
    sigma_uhf(sigma_pcgA, sigma_pcgB, zvectorA, zvectorB);

    // Level Shift
    if (level_shift == "TRUE") {
        sigma_pcgA->axpy(zvectorA, lshift_parameter);
        sigma_pcgB->axpy(zvectorB, lshift_parameter);
    }

    // Form sigma vector
    for (int ai = 0; ai < nidpA; ai++) sigma_pcg->set(ai, sigma_pcgA->get(ai));
    for (int ai = 0; ai < nidpB; ai++) sigma_pcg->set(ai + nidpA, sigma_pcgB->get(ai));

    // Build r0
    // alpha
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
              residual->set(ai, -WorbA->get(a + noccA, i));       
         }
    }
    // beta
    for (int a = 0, ai = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++, ai++) {
              residual->set(ai + nidpA, -WorbB->get(a + noccB, i));       
         }
    }
    residual->subtract(sigma_pcg);
    r_pcg->copy(residual);

    // Build z0
    z_pcg->dirprd(Minv_pcg, r_pcg);

    // Build p0
    p_pcg->copy(z_pcg);

    // Form initial pA and pB vectors
    for (int ai = 0; ai < nidpA; ai++) p_pcgA->set(ai, p_pcg->get(ai));
    for (int ai = 0; ai < nidpB; ai++) p_pcgB->set(ai, p_pcg->get(ai + nidpA));

    // Call Orbital Response Solver
    orb_resp_pcg_uhf();

    // Memfree alpha
    zvec_new.reset(); 
    Minv_pcg.reset(); 
    sigma_pcg.reset();
    sigma_pcgA.reset();
    sigma_pcgB.reset();
    r_pcg.reset();
    r_pcg_new.reset();
    z_pcg.reset();
    z_pcg_new.reset();
    p_pcg.reset();
    p_pcgA.reset();
    p_pcgB.reset();
    p_pcg_new.reset();
    dr_pcg.reset();
    residual.reset();

    // Build zvector for VO block
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
	      zvectorA->set(ai, zvector->get(ai));
         }
    }

    // Build zvector for vo block
    for (int a = 0, ai = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++, ai++) {
	      zvectorB->set(ai, zvector->get(ai + nidpA));
         }
    }
    zvector.reset();

    // Build kappa for VO block
    // alpha
    #pragma omp parallel for
    for (int x = 0; x < nidpA; x++) {
         int p = idprowA->get(x);
	 int q = idpcolA->get(x);
         if (p >= noccA && q < noccA) {
             int ai = vo_idxAA->get(p-noccA,q);
	     kappaA->set(x, zvectorA->get(ai));
         }
    }
    zvectorA.reset();
    // beta
    #pragma omp parallel for
    for (int x = 0; x < nidpB; x++) {
         int p = idprowB->get(x);
	 int q = idpcolB->get(x);
         if (p >= noccB && q < noccB) {
             int ai = vo_idxBB->get(p-noccB,q);
	     kappaB->set(x, zvectorB->get(ai));
         }
    }
    zvectorB.reset();

    if (nfrzc > 0) {
      // Compute OO-Block orb rot params
      approx_diag_mohess_oo();
      #pragma omp parallel for
      for (int x = 0; x < nidpA; x++) {
	   int p = idprowA->get(x);
	   int q = idpcolA->get(x);
           if (p < noccA && q < noccA) {
               double value = AooA->get(p-nfrzc,q); 
	       kappaA->set(x, -wogA->get(x)/value);
           }
      }
      // Compute oo-Block orb rot params
      #pragma omp parallel for
      for (int x = 0; x < nidpB; x++) {
	   int p = idprowB->get(x);
	   int q = idpcolB->get(x);
           if (p < noccB && q < noccB) {
               double value = AooB->get(p-nfrzc,q); 
	       kappaB->set(x, -wogB->get(x)/value);
           }
      }
    } // if (nfrzc > 0) 

    // If LINEQ FAILED!
    if (pcg_conver == 0) {
        //if (print_ > 1 ) {
        outfile->Printf("\tWarning!!! PCG did NOT converged in %2d iterations, switching to an approximately diagonal MO Hessian. \n", itr_pcg);
        
        //}

        // Compute VO-Block Hess
        approx_diag_mohess_vo();

	// alpha
        #pragma omp parallel for
        for(int x = 0; x < nidpA; x++) {
	    int p = idprowA->get(x);
	    int q = idpcolA->get(x);
            double value = 0.0;
            if (p >= noccA && q < noccA) value = AvoA->get(p-noccA,q); 
            else if (p < noccA && q < noccA) value = AooA->get(p-nfrzc,q); 
	    kappaA->set(x, -wogA->get(x)/value);
        }
	
	// beta
        #pragma omp parallel for
        for(int x = 0; x < nidpB; x++) {
	    int p = idprowB->get(x);
	    int q = idpcolB->get(x);
            double value = 0.0;
            if (p >= noccB && q < noccB) value = AvoB->get(p-noccB,q); 
            else if (p < noccB && q < noccB) value = AooB->get(p-nfrzc,q); 
	    kappaB->set(x, -wogB->get(x)/value);
        }
    } // end if pcg_conver = 0

        // find biggest_kappa 
	biggest_kappaA=0;            
	for (int i=0; i<nidpA;i++) { 
	    if (fabs(kappaA->get(i)) > biggest_kappaA) biggest_kappaA=fabs(kappaA->get(i));
	}
	
	biggest_kappaB=0;            
	for (int i=0; i<nidpB;i++){ 
	    if (fabs(kappaB->get(i)) > biggest_kappaB) biggest_kappaB=fabs(kappaB->get(i));
	}
	
        // Scale
	if (biggest_kappaA > step_max) {   
	    for (int i=0; i<nidpA;i++) kappaA->set(i, kappaA->get(i) *(step_max/biggest_kappaA));
	}
	 
	if (biggest_kappaB > step_max) {   
	    for (int i=0; i<nidpB;i++) kappaB->set(i, kappaB->get(i) *(step_max/biggest_kappaB));
	}
	 
        // find biggest_kappa again 
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

        // norm
	rms_kappaA=0;
	rms_kappaB=0;
	rms_kappaA = kappaA->rms();
	rms_kappaB = kappaB->rms();
	
        // print
        if(print_ > 2){
          kappaA->print();
          kappaB->print();
        }
}// end if (reference_ == "UNRESTRICTED") 
 //outfile->Printf("\n kappa_orb_resp_pcg done. \n"); 
}// end kappa_orb_resp_pcg

//=======================================================
//          PCG (RHF)
//=======================================================          
void DFOCC::orb_resp_pcg_rhf()
{ 

    SharedTensor2d K, L;
    itr_pcg = 0;
    double rms_r_pcgA = 0.0;
    double rms_residual = 0.0;
    double rms_pcg = 0.0;
    double a_pcgA = 0.0;
    double b_pcgA = 0.0;
    pcg_conver = 1; // assuming pcg will converge

 // Head of the loop
 do
 {

    //outfile->Printf( "pcg iter: %3d \n", itr_pcg); 
    // Set 
    //SvoA->set(sigma_pcgA);
    //PvoA->set(p_pcgA);


    // Build sigma
    sigma_rhf(sigma_pcgA, p_pcgA);

    // Level Shift
    if (level_shift == "TRUE") sigma_pcgA->axpy(p_pcgA, lshift_parameter);

    // Compute line search parameter alpha
    a_pcgA = r_pcgA->dot(z_pcgA) / p_pcgA->dot(sigma_pcgA);

   // Build kappa-new
   zvec_newA->zero();
   zvec_newA->copy(p_pcgA);
   zvec_newA->scale(a_pcgA);
   zvec_newA->add(zvectorA);

   // Build r-new
   r_pcg_newA->zero();
   r_pcg_newA->copy(sigma_pcgA);
   r_pcg_newA->scale(-a_pcgA);
   r_pcg_newA->add(r_pcgA);

   // Build z-new
   z_pcg_newA->dirprd(Minv_pcgA, r_pcg_newA);

   // Build line search parameter beta
   if (pcg_beta_type_ == "FLETCHER_REEVES") {
       b_pcgA = r_pcg_newA->dot(z_pcg_newA) / r_pcgA->dot(z_pcgA);
   }

   else if (pcg_beta_type_ == "POLAK_RIBIERE") {
       dr_pcgA->copy(r_pcg_newA);
       dr_pcgA->subtract(r_pcgA);
       b_pcgA = z_pcg_newA->dot(dr_pcgA) / z_pcgA->dot(r_pcgA);
   }

   // Build p-new
   p_pcg_newA->copy(p_pcgA);
   p_pcg_newA->scale(b_pcgA);
   p_pcg_newA->add(z_pcg_newA);

   // RMS 
   rms_pcg = 0.0;
   rms_pcg = zvec_newA->rms(zvectorA);
   rms_r_pcgA = r_pcg_newA->rms();

   // Reset
   zvectorA->copy(zvec_newA);
   r_pcgA->copy(r_pcg_newA);
   z_pcgA->copy(z_pcg_newA);
   p_pcgA->copy(p_pcg_newA);

   // increment iteration index 
   itr_pcg++;

   // Print
   //fprintf(outfile,"\t%3d     %12.2e     %12.2e\n",itr_pcg,rms_pcg,rms_r_pcgA);
   //fflush(outfile);

   // If we exceed maximum number of iteration, break the loop
   if (itr_pcg >= pcg_maxiter) {
       pcg_conver = 0;
       break;
   }  

   if (rms_pcg < tol_pcg) break;  

 }
 while(fabs(rms_pcg) >= tol_pcg || fabs(rms_r_pcgA) >= tol_pcg);  

}// end orb_resp_pcg_rhf


//=======================================================
//          PCG (UHF)
//=======================================================          
void DFOCC::orb_resp_pcg_uhf()
{ 
    SharedTensor2d K, L;
    itr_pcg = 0;
    double rms_r_pcg = 0.0;
    double rms_pcg = 0.0;
    double a_pcg = 0.0;
    double b_pcg = 0.0;
    double rms_residual = 0.0;
    pcg_conver = 1; // assuming pcg will converge

 // Head of the loop
 do
 {

    //outfile->Printf( "pcg iter: %3d \n", itr_pcg); 
    // Build sigma
    sigma_uhf(sigma_pcgA, sigma_pcgB, p_pcgA, p_pcgB);


    // Level Shift
    if (level_shift == "TRUE") {
        sigma_pcgA->axpy(p_pcgA, lshift_parameter);
        sigma_pcgB->axpy(p_pcgB, lshift_parameter);
    }

    // Form sigma vector
    for (int ai = 0; ai < nidpA; ai++) sigma_pcg->set(ai, sigma_pcgA->get(ai));
    for (int ai = 0; ai < nidpB; ai++) sigma_pcg->set(ai + nidpA, sigma_pcgB->get(ai));

   // Build line search parameter alpha
   a_pcg = r_pcg->dot(z_pcg) / p_pcg->dot(sigma_pcg);

   // Build kappa-new
   zvec_new->copy(p_pcg);
   zvec_new->scale(a_pcg);
   zvec_new->add(zvector);

   // Build r-new
   r_pcg_new->copy(sigma_pcg);
   r_pcg_new->scale(-a_pcg);
   r_pcg_new->add(r_pcg);
   rms_r_pcg = r_pcg_new->rms();

   // Build z-new
   z_pcg_new->dirprd(Minv_pcg, r_pcg_new);

   // Build line search parameter beta
   if (pcg_beta_type_ == "FLETCHER_REEVES") {
       b_pcg = r_pcg_new->dot(z_pcg_new) / r_pcg->dot(z_pcg);
   }

   else if (pcg_beta_type_ == "POLAK_RIBIERE") {
       dr_pcg->copy(r_pcg_new);
       dr_pcg->subtract(r_pcg);
       b_pcg = z_pcg_new->dot(dr_pcg) / z_pcg->dot(r_pcg);
   }

   // Build p-new
   p_pcg_new->copy(p_pcg);
   p_pcg_new->scale(b_pcg);
   p_pcg_new->add(z_pcg_new);

   // RMS 
   rms_pcg = 0.0;
   rms_pcg = zvec_new->rms(zvector);
   rms_r_pcg = r_pcg_new->rms();

   // Reset
   zvector->copy(zvec_new);
   r_pcg->copy(r_pcg_new);
   z_pcg->copy(z_pcg_new);
   p_pcg->copy(p_pcg_new);

    // Form pA and pB vectors
    for (int ai = 0; ai < nidpA; ai++) p_pcgA->set(ai, p_pcg->get(ai));
    for (int ai = 0; ai < nidpB; ai++) p_pcgB->set(ai, p_pcg->get(ai + nidpA));

   // increment iteration index 
   itr_pcg++;

   // Print
   //fprintf(outfile,"\t%3d     %12.2e     %12.2e\n",itr_pcg,rms_pcg,rms_r_pcg);
   //fflush(outfile);

   // If we exceed maximum number of iteration, break the loop
   if (itr_pcg >= pcg_maxiter) {
       pcg_conver = 0;
       break;
   }  

   if (rms_pcg < tol_pcg) break;  

 }
 while(fabs(rms_pcg) >= tol_pcg || fabs(rms_r_pcg) >= tol_pcg);  
 
}// end orb_resp_pcg_uhf

}} // End Namespaces
