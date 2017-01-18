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

using namespace psi;
using namespace std;

namespace psi{ namespace dfoccwave{

void DFOCC::z_vector_solver()
{ 
//outfile->Printf("\n z_vector_solver is starting... \n"); 
    outfile->Printf("\tSolving orbital Z-vector equations...\n");
    

    SharedTensor2d K, L;

if (reference_ == "RESTRICTED") {
    // Memalloc
    zvectorA = SharedTensor1d(new Tensor1d("Alpha Z-Vector", noccA * nvirA));
    zvec_newA = SharedTensor1d(new Tensor1d("Alpha New Z-Vector", noccA * nvirA));
    Minv_pcgA = SharedTensor1d(new Tensor1d("Alpha PCG M inverse", noccA * nvirA));
    sigma_pcgA = SharedTensor1d(new Tensor1d("Alpha PCG sigma", noccA * nvirA));
    residualA = SharedTensor1d(new Tensor1d("Alpha PCG Residual", noccA * nvirA));
    Wvo_vecA = SharedTensor1d(new Tensor1d("Effective MO gradient <V|O>", noccA * nvirA));

    // DIIS
    nvar = num_vecs +1;
    itr_diis = 0;
    vecsA = SharedTensor2d(new Tensor2d("Alpha MO DIIS Vectors", num_vecs, nidpA));
    errvecsA = SharedTensor2d(new Tensor2d("Alpha MO DIIS Error Vectors", num_vecs, nidpA));
    wog_intA = SharedTensor1d(new Tensor1d("Alpha Interpolated MO grad vector", nidpA));

    // Build initial guess, Wvo, and M-1 
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
              double value = FockA->get(a + noccA, a + noccA) - FockA->get(i,i);
              zvectorA->set(ai, -WorbA->get(a + noccA, i) / (2.0*value));       
              Minv_pcgA->set(ai, 0.5/value);
              Wvo_vecA->set(ai, WorbA->get(a + noccA, i));       
         }
    }

    // Call Orbital Response Solver
    zvec_solver_rhf();

    // Memfree
    zvec_newA.reset(); 
    Minv_pcgA.reset(); 
    sigma_pcgA.reset();
    residualA.reset();
    Wvo_vecA.reset();
    vecsA.reset();
    errvecsA.reset();
    wog_intA.reset();

    // Build Zvo
    ZvoA = SharedTensor2d(new Tensor2d("Zvector <V|O>", nvirA, noccA));
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
              ZvoA->set(a, i, zvectorA->get(ai));
	 }
    }
    //zvectorA->print();

    // Build Z_ia = Z_ai
    ZovA = SharedTensor2d(new Tensor2d("Zvector <O|V>", noccA, nvirA));
    ZovA = ZvoA->transpose();
    //ZvoA->print();
    //ZovA->print();

    // If LINEQ FAILED!
    if (pcg_conver == 0) {
        outfile->Printf("\tWarning!!! Iterative solver did NOT converged in %2d iterations. \n", itr_pcg);
        outfile->Printf("\tI will solve the z-vector equation with a direct method.\n");
        
        z_vector();
    } // end if pcg_conver = 0


}// end if (reference_ == "RESTRICTED") 

else if (reference_ == "UNRESTRICTED") {
    // Alpha Memalloc
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

    // Beta Memalloc
    zvectorB = SharedTensor1d(new Tensor1d("Beta Z-Vector", noccB * nvirB));
    zvec_newB = SharedTensor1d(new Tensor1d("Beta New Z-Vector", noccB * nvirB));
    Minv_pcgB = SharedTensor1d(new Tensor1d("Beta PCG M inverse", noccB * nvirB));
    sigma_pcgB = SharedTensor1d(new Tensor1d("Beta PCG sigma", noccB * nvirB));
    r_pcgB = SharedTensor1d(new Tensor1d("Beta PCG r", noccB * nvirB));
    r_pcg_newB = SharedTensor1d(new Tensor1d("Beta PCG new r", noccB * nvirB));
    z_pcgB = SharedTensor1d(new Tensor1d("Beta PCG z", noccB * nvirB));
    z_pcg_newB = SharedTensor1d(new Tensor1d("Beta PCG new z", noccB * nvirB));
    p_pcgB = SharedTensor1d(new Tensor1d("Beta PCG p", noccB * nvirB));
    p_pcg_newB = SharedTensor1d(new Tensor1d("Beta PCG new p", noccB * nvirB));
    dr_pcgB = SharedTensor1d(new Tensor1d("Beta PCG dr", noccB * nvirB));

    // Build kappa0 and M0
    // alpha
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
              double value = FockA->get(a + noccA, a + noccA) - FockA->get(i,i);
              zvectorA->set(ai, -WorbA->get(a + noccA, i) / (2.0*value));       
              Minv_pcgA->set(ai, 0.5/value);
              p_pcgA->set(ai, zvectorA->get(ai));
         }
    }

    // beta
    for (int a = 0, ai = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++, ai++) {
              double value = FockB->get(a + noccB, a + noccB) - FockB->get(i,i);
              zvectorB->set(ai, -WorbB->get(a + noccB, i) / (2.0*value));       
              Minv_pcgB->set(ai, 0.5/value);
              p_pcgB->set(ai, zvectorB->get(ai));
         }
    }

    // Build sigma0
    // Memalloc
    SharedTensor2d SvoA = SharedTensor2d(new Tensor2d("PCG Sigma <V|O>", nvirA, noccA));
    SharedTensor2d SvoB = SharedTensor2d(new Tensor2d("PCG Sigma <v|o>", nvirB, noccB));
    SharedTensor2d PvoA = SharedTensor2d(new Tensor2d("PCG P <V|O>", nvirA, noccA));
    SharedTensor2d PvoB = SharedTensor2d(new Tensor2d("PCG P <v|o>", nvirB, noccB));
    SharedTensor1d pQ = SharedTensor1d(new Tensor1d("DF_BASIS_SCF p_Q", nQ_ref));

    // Set 
    SvoA->set(sigma_pcgA);
    SvoB->set(sigma_pcgB);
    PvoA->set(p_pcgA);
    PvoB->set(p_pcgB);

    // p_Q = \sum_{BJ} b_BJ^Q p_BJ + \sum_{bj} b_bj^Q p_bj
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB, nvirB));
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    SharedTensor2d bQvoB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vo)", nQ_ref, nvirB, noccB));
    bQvoB->swap_3index_col(bQovB);
    bQovB.reset();
    pQ->gemv(false, bQvoB, p_pcgB, 1.0, 0.0);
    bQvoB.reset();
    // alpha contribution
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    SharedTensor2d bQvoA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VO)", nQ_ref, nvirA, noccA));
    bQvoA->swap_3index_col(bQovA);
    bQovA.reset();
    pQ->gemv(false, bQvoA, p_pcgA, 1.0, 1.0);

    // Build alpha sigma
    // s_AI = 2 \sum_{B} f_AB p_BI - 2 \sum_{J} f_IJ p_AJ  
    SvoA->gemm(false, false, FvvA, PvoA, 2.0, 0.0);
    SvoA->gemm(false, false, PvoA, FooA, -2.0, 1.0);

    // s_AI += 4 \sum_{Q} b_AI^Q p^Q
    SvoA->gemv(true, bQvoA, pQ, 4.0, 1.0);

    // p_IJ^Q = \sum_{B} b_BI^Q p_BJ
    SharedTensor2d pQooA = SharedTensor2d(new Tensor2d("PCG P (Q|OO)", nQ_ref, noccA, noccA));
    pQooA->contract323(true, false, noccA, noccA, bQvoA, PvoA, 1.0, 0.0);
    // s_AI += -2 \sum_{Q} \sum_{J} b_AJ^Q p_IJ^Q
    SvoA->contract332(false, true, noccA, bQvoA, pQooA, -2.0, 1.0);
    pQooA.reset();
    bQvoA.reset();

    // p_AJ^Q = \sum_{B} b_BA^Q p_BJ
    bQvvA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    bQvvA->read(psio_, PSIF_DFOCC_INTS, true, true);
    SharedTensor2d pQvoA = SharedTensor2d(new Tensor2d("PCG P (Q|VO)", nQ_ref, nvirA, noccA));
    pQvoA->contract323(false, false, nvirA, noccA, bQvvA, PvoA, 1.0, 0.0);
    bQvvA.reset();
    // s_AI += -2 \sum_{Q} \sum_{J} b_IJ^Q p_AJ^Q
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA * noccA));
    bQooA->read(psio_, PSIF_DFOCC_INTS);
    SvoA->contract332(false, false, noccA, pQvoA, bQooA, -2.0, 1.0);
    bQooA.reset();
    pQvoA.reset();

    // Form sigma vector
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
              sigma_pcgA->set(ai, SvoA->get(a,i));
         }
    }

    // Build beta sigma
    // s_ai = 2 \sum_{b} f_ab p_bi - 2 \sum_{j} f_ij p_aj  
    SvoB->gemm(false, false, FvvB, PvoB, 2.0, 0.0);
    SvoB->gemm(false, false, PvoB, FooB, -2.0, 1.0);

    // s_ai += 4 \sum_{Q} bai^Q p^Q
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB, nvirB));
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    bQvoB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vo)", nQ_ref, nvirB, noccB));
    bQvoB->swap_3index_col(bQovB);
    bQovB.reset();
    SvoB->gemv(true, bQvoB, pQ, 4.0, 1.0);

    // p_ij^Q = \sum_{b} b_bi^Q p_bj
    SharedTensor2d pQooB = SharedTensor2d(new Tensor2d("PCG P (Q|oo)", nQ_ref, noccB, noccB));
    pQooB->contract323(true, false, noccB, noccB, bQvoB, PvoB, 1.0, 0.0);
    // s_ai += -2 \sum_{Q} \sum_{j} b_aj^Q p_ij^Q
    SvoB->contract332(false, true, noccB, bQvoB, pQooB, -2.0, 1.0);
    pQooB.reset();
    bQvoB.reset();

    // p_aj^Q = \sum_{b} b_ba^Q p_bj
    bQvvB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vv)", nQ_ref, nvirB, nvirB));
    bQvvB->read(psio_, PSIF_DFOCC_INTS, true, true);
    SharedTensor2d pQvoB = SharedTensor2d(new Tensor2d("PCG P (Q|vo)", nQ_ref, nvirB, noccB));
    pQvoB->contract323(false, false, nvirB, noccB, bQvvB, PvoB, 1.0, 0.0);
    bQvvB.reset();
    // s_ai += -2 \sum_{Q} \sum_{j} b_ij^Q p_aj^Q
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB * noccB));
    bQooB->read(psio_, PSIF_DFOCC_INTS);
    SvoB->contract332(false, false, noccB, pQvoB, bQooB, -2.0, 1.0);
    bQooB.reset();
    pQvoB.reset();

    // Form sigma vector
    for (int a = 0, ai = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++, ai++) {
              sigma_pcgB->set(ai, SvoB->get(a,i));
         }
    }

    // memfree
    PvoA.reset();
    PvoB.reset();
    SvoA.reset();
    SvoB.reset();
    pQ.reset();

    // Build r0
    // alpha
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
              r_pcgA->set(ai, -WorbA->get(a + noccA, i) - sigma_pcgA->get(ai));       
         }
    }
    // beta
    for (int a = 0, ai = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++, ai++) {
              r_pcgB->set(ai, -WorbB->get(a + noccB, i) - sigma_pcgB->get(ai));
         }
    }

    // Build z0
    z_pcgA->dirprd(Minv_pcgA, r_pcgA);
    z_pcgB->dirprd(Minv_pcgB, r_pcgB);

    // Build p0
    p_pcgA->copy(z_pcgA);
    p_pcgB->copy(z_pcgB);

    // Call Orbital Response Solver
    pcg_solver_uhf();

    // Memfree alpha
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

    // Memfree beta
    zvec_newB.reset(); 
    Minv_pcgB.reset(); 
    sigma_pcgB.reset();
    r_pcgB.reset();
    r_pcg_newB.reset();
    z_pcgB.reset();
    z_pcg_newB.reset();
    p_pcgB.reset();
    p_pcg_newB.reset();
    dr_pcgB.reset();

    // Build Zvo
    // Alpha
    ZvoA = SharedTensor2d(new Tensor2d("Zvector <V|O>", nvirA, noccA));
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
              ZvoA->set(a, i, zvectorA->get(ai));
	 }
    }

    // Build Z_ia = Z_ai
    ZovA = SharedTensor2d(new Tensor2d("Zvector <O|V>", noccA, nvirA));
    ZovA = ZvoA->transpose();

    // Beta
    ZvoB = SharedTensor2d(new Tensor2d("Zvector <v|o>", nvirB, noccB));
    for (int a = 0, ai = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++, ai++) {
              ZvoB->set(a, i, zvectorB->get(ai));
	 }
    }

    // Build Z_ia = Z_ai
    ZovB = SharedTensor2d(new Tensor2d("Zvector <o|v>", noccB, nvirB));
    ZovB = ZvoB->transpose();

    // If LINEQ FAILED!
    if (pcg_conver == 0) {
        outfile->Printf("\tWarning!!! PCG did NOT converged in %2d iterations. \n", itr_pcg);
        
    } 

}// end if (reference_ == "UNRESTRICTED") 
 //outfile->Printf("\n z_vector_pcg done. \n"); 
}// end z_vector_pcg

//=======================================================
//          Z-Vector (RHF)
//=======================================================          
void DFOCC::zvec_solver_rhf()
{ 

    SharedTensor2d K, L;
    itr_pcg = 0;
    double rms_residual = 0.0;
    double rms_pcg = 0.0;
    pcg_conver = 1; // assuming pcg will converge

outfile->Printf( "\n\t          Z-Vector Solver \n");
outfile->Printf( "\t   ------------------------------ \n");
outfile->Printf( "\tIter     RMS Z-Vector        RMS Residual  \n");
outfile->Printf( "\t----    ---------------    --------------\n");

   //Minv_pcgA->print();
   //zvectorA->print();

 // Head of the loop
 do
 {

   //outfile->Printf( "pcg iter: %3d \n", itr_pcg); 
   // Build sigma
   sigma_pcgA->zero();
   sigma_orb_resp_rhf(sigma_pcgA, zvectorA);

   // Build W + sigma
   sigma_pcgA->add(Wvo_vecA);
   //sigma_pcgA->print();

   // Build kappa-new
   zvec_newA->zero();
   zvec_newA->dirprd(Minv_pcgA, sigma_pcgA);
   zvec_newA->scale(-1.0);
   //zvec_newA->add(zvectorA);

    // Build sigma for kappa
    //sigma_rhf(sigma_pcgA, zvectorA);
    sigma_rhf(sigma_pcgA, zvec_newA);
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
              residualA->set(ai, -WorbA->get(a + noccA, i) - sigma_pcgA->get(ai));       
         }
    }
    rms_residual = residualA->rms();

   // RMS 
   rms_pcg = 0.0;
   rms_pcg = zvec_newA->rms(zvectorA);

/*
if (do_diis_ == 1) {
  
        // starting with itr = 1
        itr_diis++;
   
        // Form Diis Error Vector & Extrapolant Alpha Spin Case
	if (itr_diis <= num_vecs) {  
	  for(int i = 0; i < nidpA; i++){
	    //errvecsA->set(itr_diis-1, i, zvec_newA->get(i) - zvectorA->get(i));
	    errvecsA->set(itr_diis-1, i, residualA->get(i));
	    vecsA->set(itr_diis-1, i, zvec_newA->get(i));
	  }  
	}
	
	if (itr_diis > num_vecs) {  
	  for(int j = 0; j < (num_vecs-1); j++){
	    for(int i = 0; i < nidpA; i++){
	      errvecsA->set(j, i, errvecsA->get(j+1, i));
	      vecsA->set(j, i, vecsA->get(j+1, i));
	    }  
	  }
	  
	  for(int i = 0; i < nidpA; i++){
	    //errvecsA->set(num_vecs-1, i, zvec_newA->get(i) - zvectorA->get(i));
	    errvecsA->set(num_vecs-1, i, residualA->get(i));
	    vecsA->set(num_vecs-1, i, zvec_newA->get(i));
	  }    
	}
	
        // Extrapolate 
        if (itr_diis >= num_vecs) {
	  diis(nidpA, vecsA, errvecsA, zvec_newA, wog_intA);
	}
	
}// end if (do_diis_ == 1) 
*/

   // Reset
   //zvectorA->print();
   //zvec_newA->print();
   zvectorA->copy(zvec_newA);

   // increment iteration index 
   itr_pcg++;

   // Print
   outfile->Printf("\t%3d     %12.2e     %12.2e\n",itr_pcg,rms_pcg,rms_residual);
   

   // If we exceed maximum number of iteration, break the loop
   if (itr_pcg >= pcg_maxiter) {
       pcg_conver = 0;
       break;
   }  

   //if (rms_r_pcgA < tol_pcg || rms_pcg < tol_pcg) break;  
   //if (rms_residual < tol_pcg  || rms_pcg < tol_pcg) break;  
   if (rms_pcg < tol_pcg) break;  

 }
 while(fabs(rms_pcg) >= tol_pcg);  

    // Converged?
    //r_pcgA->print();
    //residualA->print();
    zvectorA->print();

}// end zvec_sover_rhf


//=======================================================
//          PCG (UHF)
//=======================================================          
void DFOCC::zvec_solver_uhf()
{ 
    SharedTensor2d K, L;
    itr_pcg = 0;
    double rms_r_pcgA = 0.0;
    double rms_r_pcgB = 0.0;
    double rms_r_pcg = 0.0;
    double rms_pcg = 0.0;
    double a_pcgA = 0.0;
    double a_pcgB = 0.0;
    double b_pcgA = 0.0;
    double b_pcgB = 0.0;
    pcg_conver = 1; // assuming pcg will converge

    // Memalloc
    SharedTensor2d SvoA = SharedTensor2d(new Tensor2d("PCG Sigma <V|O>", nvirA, noccA));
    SharedTensor2d SvoB = SharedTensor2d(new Tensor2d("PCG Sigma <v|o>", nvirB, noccB));
    SharedTensor2d PvoA = SharedTensor2d(new Tensor2d("PCG P <V|O>", nvirA, noccA));
    SharedTensor2d PvoB = SharedTensor2d(new Tensor2d("PCG P <v|o>", nvirB, noccB));
    SharedTensor1d pQ = SharedTensor1d(new Tensor1d("DF_BASIS_SCF p_Q", nQ_ref));

 // Head of the loop
 do
 {
    //outfile->Printf( "pcg iter: %3d \n", itr_pcg); 

    // Set 
    SvoA->set(sigma_pcgA);
    SvoB->set(sigma_pcgB);
    PvoA->set(p_pcgA);
    PvoB->set(p_pcgB);

    // p_Q = \sum_{BJ} b_BJ^Q p_BJ + \sum_{bj} b_bj^Q p_bj
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB, nvirB));
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    SharedTensor2d bQvoB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vo)", nQ_ref, nvirB, noccB));
    bQvoB->swap_3index_col(bQovB);
    bQovB.reset();
    pQ->gemv(false, bQvoB, p_pcgB, 1.0, 0.0);
    bQvoB.reset();
    // alpha contribution
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    SharedTensor2d bQvoA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VO)", nQ_ref, nvirA, noccA));
    bQvoA->swap_3index_col(bQovA);
    bQovA.reset();
    pQ->gemv(false, bQvoA, p_pcgA, 1.0, 1.0);

    // Build alpha sigma
    // s_AI = 2 \sum_{B} f_AB p_BI - 2 \sum_{J} f_IJ p_AJ  
    SvoA->gemm(false, false, FvvA, PvoA, 2.0, 0.0);
    SvoA->gemm(false, false, PvoA, FooA, -2.0, 1.0);

    // s_AI += 4 \sum_{Q} b_AI^Q p^Q
    SvoA->gemv(true, bQvoA, pQ, 4.0, 1.0);

    // p_IJ^Q = \sum_{B} b_BI^Q p_BJ
    SharedTensor2d pQooA = SharedTensor2d(new Tensor2d("PCG P (Q|OO)", nQ_ref, noccA, noccA));
    pQooA->contract323(true, false, noccA, noccA, bQvoA, PvoA, 1.0, 0.0);
    // s_AI += -2 \sum_{Q} \sum_{J} b_AJ^Q p_IJ^Q
    SvoA->contract332(false, true, noccA, bQvoA, pQooA, -2.0, 1.0);
    pQooA.reset();
    bQvoA.reset();

    // p_AJ^Q = \sum_{B} b_BA^Q p_BJ
    bQvvA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    bQvvA->read(psio_, PSIF_DFOCC_INTS, true, true);
    SharedTensor2d pQvoA = SharedTensor2d(new Tensor2d("PCG P (Q|VO)", nQ_ref, nvirA, noccA));
    pQvoA->contract323(false, false, nvirA, noccA, bQvvA, PvoA, 1.0, 0.0);
    bQvvA.reset();
    // s_AI += -2 \sum_{Q} \sum_{J} b_IJ^Q p_AJ^Q
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA * noccA));
    bQooA->read(psio_, PSIF_DFOCC_INTS);
    SvoA->contract332(false, false, noccA, pQvoA, bQooA, -2.0, 1.0);
    bQooA.reset();
    pQvoA.reset();

    // Form sigma vector
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
              sigma_pcgA->set(ai, SvoA->get(a,i));
         }
    }

    // Build beta sigma
    // s_ai = 2 \sum_{b} f_ab p_bi - 2 \sum_{j} f_ij p_aj  
    SvoB->gemm(false, false, FvvB, PvoB, 2.0, 0.0);
    SvoB->gemm(false, false, PvoB, FooB, -2.0, 1.0);

    // s_ai += 4 \sum_{Q} bai^Q p^Q
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB, nvirB));
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    bQvoB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vo)", nQ_ref, nvirB, noccB));
    bQvoB->swap_3index_col(bQovB);
    bQovB.reset();
    SvoB->gemv(true, bQvoB, pQ, 4.0, 1.0);

    // p_ij^Q = \sum_{b} b_bi^Q p_bj
    SharedTensor2d pQooB = SharedTensor2d(new Tensor2d("PCG P (Q|oo)", nQ_ref, noccB, noccB));
    pQooB->contract323(true, false, noccB, noccB, bQvoB, PvoB, 1.0, 0.0);
    // s_ai += -2 \sum_{Q} \sum_{j} b_aj^Q p_ij^Q
    SvoB->contract332(false, true, noccB, bQvoB, pQooB, -2.0, 1.0);
    pQooB.reset();
    bQvoB.reset();

    // p_aj^Q = \sum_{b} b_ba^Q p_bj
    bQvvB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vv)", nQ_ref, nvirB, nvirB));
    bQvvB->read(psio_, PSIF_DFOCC_INTS, true, true);
    SharedTensor2d pQvoB = SharedTensor2d(new Tensor2d("PCG P (Q|vo)", nQ_ref, nvirB, noccB));
    pQvoB->contract323(false, false, nvirB, noccB, bQvvB, PvoB, 1.0, 0.0);
    bQvvB.reset();
    // s_ai += -2 \sum_{Q} \sum_{j} b_ij^Q p_aj^Q
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB * noccB));
    bQooB->read(psio_, PSIF_DFOCC_INTS);
    SvoB->contract332(false, false, noccB, pQvoB, bQooB, -2.0, 1.0);
    bQooB.reset();
    pQvoB.reset();

    // Form sigma vector
    for (int a = 0, ai = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++, ai++) {
              sigma_pcgB->set(ai, SvoB->get(a,i));
         }
    }

   // Build line search parameter alpha
   a_pcgA = r_pcgA->dot(z_pcgA) / p_pcgA->dot(sigma_pcgA);
   a_pcgB = r_pcgB->dot(z_pcgB) / p_pcgB->dot(sigma_pcgB);

   // Build kappa-new
   zvec_newA->copy(p_pcgA);
   zvec_newA->scale(a_pcgA);
   zvec_newA->add(zvectorA);
   zvec_newB->copy(p_pcgB);
   zvec_newB->scale(a_pcgB);
   zvec_newB->add(zvectorB);

   // Build r-new
   r_pcg_newA->copy(sigma_pcgA);
   r_pcg_newA->scale(-a_pcgA);
   r_pcg_newA->add(r_pcgA);
   rms_r_pcgA = r_pcg_newA->rms();
   r_pcg_newB->copy(sigma_pcgB);
   r_pcg_newB->scale(-a_pcgB);
   r_pcg_newB->add(r_pcgB);
   rms_r_pcgB = r_pcg_newB->rms();
   rms_r_pcg = MAX0(rms_r_pcgA,rms_r_pcgB);

   // Build z-new
   z_pcg_newA->dirprd(Minv_pcgA, r_pcg_newA);
   z_pcg_newB->dirprd(Minv_pcgB, r_pcg_newB);

   // Build line search parameter beta
   if (pcg_beta_type_ == "FLETCHER_REEVES") {
       b_pcgA = r_pcg_newA->dot(z_pcg_newA) / r_pcgA->dot(z_pcgA);
       b_pcgB = r_pcg_newB->dot(z_pcg_newB) / r_pcgB->dot(z_pcgB);
   }

   else if (pcg_beta_type_ == "POLAK_RIBIERE") {
       dr_pcgA->copy(r_pcg_newA);
       dr_pcgA->subtract(r_pcgA);
       dr_pcgB->copy(r_pcg_newB);
       dr_pcgB->subtract(r_pcgB);
       b_pcgA = z_pcg_newA->dot(dr_pcgA) / z_pcgA->dot(r_pcgA);
       b_pcgB = z_pcg_newB->dot(dr_pcgB) / z_pcgB->dot(r_pcgB);
   }

   // Build p-new
   p_pcg_newA->copy(p_pcgA);
   p_pcg_newA->scale(b_pcgA);
   p_pcg_newA->add(z_pcg_newA);
   p_pcg_newB->copy(p_pcgB);
   p_pcg_newB->scale(b_pcgB);
   p_pcg_newB->add(z_pcg_newB);

   // Reset
   zvectorA->copy(zvec_newA);
   r_pcgA->copy(r_pcg_newA);
   z_pcgA->copy(z_pcg_newA);
   p_pcgA->copy(p_pcg_newA);
   zvectorB->copy(zvec_newB);
   r_pcgB->copy(r_pcg_newB);
   z_pcgB->copy(z_pcg_newB);
   p_pcgB->copy(p_pcg_newB);

   // RMS 
   rms_kappaA = zvectorA->rms();
   rms_kappaB = zvectorB->rms();
   rms_kappa=MAX0(rms_kappaA,rms_kappaB);
   rms_pcg = rms_kappa;

   // increment iteration index 
   itr_pcg++;

   // If we exceed maximum number of iteration, break the loop
   if (itr_pcg >= pcg_maxiter) {
       pcg_conver = 0;
       break;
   }  

   if (rms_r_pcg < tol_pcg || rms_pcg < tol_pcg) break;  

 }
 while(fabs(rms_pcg) >= tol_pcg);  

    // memfree
    PvoA.reset();
    PvoB.reset();
    SvoA.reset();
    SvoB.reset();
    pQ.reset();
}// end zvec_solver_uhf

//=======================================================
//          Sigma (RHF)
//=======================================================          
void DFOCC::sigma_orb_resp_rhf(SharedTensor1d& sigma, SharedTensor1d& p_vec)
{ 
    // Build sigma0
    // Memalloc
    SharedTensor2d SvoA = SharedTensor2d(new Tensor2d("PCG Sigma <V|O>", nvirA, noccA));
    SharedTensor1d pQ = SharedTensor1d(new Tensor1d("DF_BASIS_SCF p_Q", nQ_ref));
    SharedTensor2d PvoA = SharedTensor2d(new Tensor2d("PCG P <V|O>", nvirA, noccA));
    PvoA->set(p_vec);

    // Build sigma
    // s_ai += 4 \sum_{Q} bai^Q p^Q
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    SharedTensor2d bQvoA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VO)", nQ_ref, nvirA, noccA));
    bQvoA->swap_3index_col(bQovA);
    bQovA.reset();
    // p_Q = 2\sum_{bj} b_bj^Q p_bj
    pQ->gemv(false, bQvoA, p_vec, 2.0, 0.0);
    SvoA->gemv(true, bQvoA, pQ, 4.0, 0.0);

    // p_ij^Q = \sum_{b} b_bi^Q p_bj
    SharedTensor2d pQooA = SharedTensor2d(new Tensor2d("PCG P (Q|OO)", nQ_ref, noccA, noccA));
    pQooA->contract323(true, false, noccA, noccA, bQvoA, PvoA, 1.0, 0.0);
    // s_ai += -2 \sum_{Q} \sum_{j} b_aj^Q p_ij^Q
    SvoA->contract332(false, true, noccA, bQvoA, pQooA, -2.0, 1.0);
    pQooA.reset();
    bQvoA.reset();

    // p_aj^Q = \sum_{b} b_ba^Q p_bj
    bQvvA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    bQvvA->read(psio_, PSIF_DFOCC_INTS, true, true);
    SharedTensor2d pQvoA = SharedTensor2d(new Tensor2d("PCG P (Q|VO)", nQ_ref, nvirA, noccA));
    pQvoA->contract323(false, false, nvirA, noccA, bQvvA, PvoA, 1.0, 0.0);
    bQvvA.reset();
    // s_ai += -2 \sum_{Q} \sum_{j} b_ij^Q p_aj^Q
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    bQooA->read(psio_, PSIF_DFOCC_INTS);
    SvoA->contract332(false, false, noccA, pQvoA, bQooA, -2.0, 1.0);
    bQooA.reset();
    pQvoA.reset();

    // Form sigma vector
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
              sigma->set(ai, SvoA->get(a,i));
         }
    }

    // memfree
    SvoA.reset();
    PvoA.reset();
    pQ.reset();

}// end sigma_rhf


}} // End Namespaces
