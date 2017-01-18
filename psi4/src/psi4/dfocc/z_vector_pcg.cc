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
#include "psi4/libciomr/libciomr.h"

using namespace psi;
using namespace std;

namespace psi{ namespace dfoccwave{

void DFOCC::z_vector_pcg()
{ 
//outfile->Printf("\n z_vector_pcg is starting... \n"); 
    outfile->Printf("\tSolving orbital Z-vector equations...\n");
    

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
              //zvectorA->set(ai, -WorbA->get(a + noccA, i) / (2.0*value));       
              zvectorA->set(ai, -WvoA->get(a, i) / (2.0*value));       
              Minv_pcgA->set(ai, 0.5/value);
         }
    }

    // Build S = A kappa_0
    sigma_rhf(sigma_pcgA, zvectorA);

    // Build r0
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
              residualA->set(ai, -WvoA->get(a, i) - sigma_pcgA->get(ai));       
         }
    }
    r_pcgA->copy(residualA);

    // Build z0
    z_pcgA->dirprd(Minv_pcgA, r_pcgA);

    // Build p0
    p_pcgA->copy(z_pcgA);

    // Call Orbital Response Solver
    pcg_solver_rhf();

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

    // Build Zvo
    //zvectorA->print();
    ZvoA = SharedTensor2d(new Tensor2d("Zvector <V|O>", nvirA, noccA));
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
              ZvoA->set(a, i, zvectorA->get(ai));
	 }
    }

    // Build Z_ia = Z_ai
    ZovA = SharedTensor2d(new Tensor2d("Zvector <O|V>", noccA, nvirA));
    ZovA = ZvoA->transpose();

    // If LINEQ FAILED!
    if (pcg_conver == 0) {
        outfile->Printf("\tWarning!!! PCG did NOT converged in %2d iterations. \n", itr_pcg);
        outfile->Printf("\tI will solve the z-vector equation with a direct method.\n");
        
        z_vector();
        
    } // end if pcg_conver = 0


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
              //zvectorA->set(ai, -WorbA->get(a + noccA, i) / (2.0*value));       
              zvectorA->set(ai, -WvoA->get(a, i) / (2.0*value));       
              Minv_pcg->set(ai, 0.5/value);
         }
    }

    // beta
    for (int a = 0, ai = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++, ai++) {
              double value = FockB->get(a + noccB, a + noccB) - FockB->get(i,i);
              //zvectorB->set(ai, -WorbB->get(a + noccB, i) / (2.0*value));       
              zvectorB->set(ai, -WvoB->get(a, i) / (2.0*value));       
              Minv_pcg->set(ai + nidpA, 0.5/value);
         }
    }

    // Form initial zvector vector
    for (int ai = 0; ai < nidpA; ai++) zvector->set(ai, zvectorA->get(ai));
    for (int ai = 0; ai < nidpB; ai++) zvector->set(ai + nidpA, zvectorB->get(ai));

    // Build S = A kappa_0
    sigma_uhf(sigma_pcgA, sigma_pcgB, zvectorA, zvectorB);

    // Form sigma vector
    for (int ai = 0; ai < nidpA; ai++) sigma_pcg->set(ai, sigma_pcgA->get(ai));
    for (int ai = 0; ai < nidpB; ai++) sigma_pcg->set(ai + nidpA, sigma_pcgB->get(ai));

    // Build r0
    // alpha
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
              //residual->set(ai, -WorbA->get(a + noccA, i));       
              residual->set(ai, -WvoA->get(a, i));       
         }
    }
    // beta
    for (int a = 0, ai = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++, ai++) {
              //residual->set(ai + nidpA, -WorbB->get(a + noccB, i));       
              residual->set(ai + nidpA, -WvoB->get(a, i));       
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
    pcg_solver_uhf();

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
    //zvectorA->print();

    // Build zvector for vo block
    for (int a = 0, ai = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++, ai++) {
	      zvectorB->set(ai, zvector->get(ai + nidpA));
         }
    }
    //zvectorB->print();
    zvector.reset();

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
        outfile->Printf("\tI will solve the z-vector equation with a direct method.\n");
        
        z_vector();
        
    } 

}// end if (reference_ == "UNRESTRICTED") 
 //outfile->Printf("\n z_vector_pcg done. \n"); 
}// end z_vector_pcg

//=======================================================
//          PCG (RHF)
//=======================================================          
void DFOCC::pcg_solver_rhf()
{ 

    SharedTensor2d K, L;
    itr_pcg = 0;
    double rms_r_pcgA = 0.0;
    double rms_residual = 0.0;
    double rms_pcg = 0.0;
    double a_pcgA = 0.0;
    double b_pcgA = 0.0;
    pcg_conver = 1; // assuming pcg will converge

outfile->Printf( "\n\t            PCG Solver \n");
outfile->Printf( "\t   ------------------------------ \n");
outfile->Printf( "\tIter     RMS Z Vector        RMS Residual  \n");
outfile->Printf( "\t----    ---------------    ---------------\n");

 // Head of the loop
 do
 {

    //outfile->Printf( "pcg iter: %3d \n", itr_pcg); 
    // Build sigma
    sigma_rhf(sigma_pcgA, p_pcgA);

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

    /*
    // Build sigma for kappa
    sigma_rhf(sigma_pcgA, zvectorA);

    // Build r0
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
              residualA->set(ai, -WorbA->get(a + noccA, i) - sigma_pcgA->get(ai));       
         }
    }
    rms_residual = residualA->rms();
    */

   // increment iteration index 
   itr_pcg++;

   // Print
   outfile->Printf("\t%3d     %12.2e     %12.2e\n",itr_pcg,rms_pcg,rms_r_pcgA);
   

   // If we exceed maximum number of iteration, break the loop
   if (itr_pcg >= pcg_maxiter) {
       pcg_conver = 0;
       break;
   }  

   if (rms_pcg < tol_pcg) break;  

 }
 while(fabs(rms_pcg) >= tol_pcg || fabs(rms_r_pcgA) >= tol_pcg);  

    // Converged?
    outfile->Printf("\n");
    
    //residualA->print();

}// end pcg_sover_rhf


//=======================================================
//          PCG (UHF)
//=======================================================          
void DFOCC::pcg_solver_uhf()
{ 
    SharedTensor2d K, L;
    itr_pcg = 0;
    double rms_r_pcg = 0.0;
    double rms_pcg = 0.0;
    double a_pcg = 0.0;
    double b_pcg = 0.0;
    double rms_residual = 0.0;
    pcg_conver = 1; // assuming pcg will converge

outfile->Printf( "\n\t            PCG Solver \n");
outfile->Printf( "\t   ------------------------------ \n");
outfile->Printf( "\tIter     RMS Z Vector        RMS Residual  \n");
outfile->Printf( "\t----    ---------------    ---------------\n");


 // Head of the loop
 do
 {
    //outfile->Printf( "pcg iter: %3d \n", itr_pcg); 
    // Build sigma
    sigma_uhf(sigma_pcgA, sigma_pcgB, p_pcgA, p_pcgB);

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

    /*
    // Form zvectorA and zvectorB vectors
    for (int ai = 0; ai < nidpA; ai++) zvectorA->set(ai, zvector->get(ai));
    for (int ai = 0; ai < nidpB; ai++) zvectorB->set(ai, zvector->get(ai + nidpA));

    // Build sigma for kappa
    sigma_uhf(sigma_pcgA, sigma_pcgB, zvectorA, zvectorB);

    // Form sigma vector
    for (int ai = 0; ai < nidpA; ai++) sigma_pcg->set(ai, sigma_pcgA->get(ai));
    for (int ai = 0; ai < nidpB; ai++) sigma_pcg->set(ai + nidpA, sigma_pcgB->get(ai));

    // Build r0
    //#pragma omp parallel for
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
              residual->set(ai, -WorbA->get(a + noccA, i));       
         }
    }

    //#pragma omp parallel for
    for (int a = 0, ai = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++, ai++) {
              residual->set(ai + nidpA, -WorbB->get(a + noccB, i));       
         }
    }
    residual->subtract(sigma_pcg);
    rms_residual = residual->rms();
    */

   // increment iteration index 
   itr_pcg++;

   // Print
   outfile->Printf("\t%3d     %12.2e     %12.2e\n",itr_pcg,rms_pcg,rms_r_pcg);
   

   // If we exceed maximum number of iteration, break the loop
   if (itr_pcg >= pcg_maxiter) {
       pcg_conver = 0;
       break;
   }  

   if (rms_pcg < tol_pcg) break;  

 }
 while(fabs(rms_pcg) >= tol_pcg || fabs(rms_r_pcg) >= tol_pcg);  
    // Converged?
    outfile->Printf("\n");
    

}// end pcg_solver_uhf

//=======================================================
//          Sigma (RHF)
//=======================================================          
void DFOCC::sigma_rhf(SharedTensor1d& sigma, SharedTensor1d& p_vec)
{ 
    // Build sigma0
    // Memalloc
    SharedTensor2d SvoA = SharedTensor2d(new Tensor2d("PCG Sigma <V|O>", nvirA, noccA));
    SharedTensor1d pQ = SharedTensor1d(new Tensor1d("DF_BASIS_SCF p_Q", nQ_ref));
    SharedTensor2d PvoA = SharedTensor2d(new Tensor2d("PCG P <V|O>", nvirA, noccA));
    PvoA->set(p_vec);

    // Build sigma
    // s_ai = 2 \sum_{b} f_ab p_bi - 2 \sum_{j} f_ij p_aj  = 2 (f_aa - f_ii) p_ai 
    //SvoA->gemm(false, false, FvvA, PvoA, 2.0, 0.0);
    //SvoA->gemm(false, false, PvoA, FooA, -2.0, 1.0);
    #pragma omp parallel for
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              double value = FockA->get(a + noccA, a + noccA) - FockA->get(i,i);
              SvoA->set(a, i, 2.0 * value * PvoA->get(a, i));
         }
    }
    
    // p_Q = 2\sum_{bj} b_bj^Q p_bj
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    SharedTensor2d bQvoA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VO)", nQ_ref, nvirA, noccA));
    bQvoA->swap_3index_col(bQovA);
    bQovA.reset();
    pQ->gemv(false, bQvoA, p_vec, 2.0, 0.0);

    // s_ai += 4 \sum_{Q} bai^Q p^Q
    SvoA->gemv(true, bQvoA, pQ, 4.0, 1.0);
    
    // p_ij^Q = \sum_{b} b_bi^Q p_bj = \su_{b} b_ib^Q p_bj
    SharedTensor2d pQooA = SharedTensor2d(new Tensor2d("PCG P (Q|OO)", nQ_ref, noccA, noccA));
    pQooA->contract323(true, false, noccA, noccA, bQvoA, PvoA, 1.0, 0.0);

    // s_ai += -2 \sum_{Q} \sum_{j} b_aj^Q p_ij^Q = b[Q](a,j) p'[Q](j,i) 
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

//=======================================================
//          Sigma (UHF)
//=======================================================          
void DFOCC::sigma_uhf(SharedTensor1d& sigma_A, SharedTensor1d& sigma_B, SharedTensor1d& p_vecA, SharedTensor1d& p_vecB)
{ 
    // Memalloc
    SharedTensor2d SvoA = SharedTensor2d(new Tensor2d("PCG Sigma <V|O>", nvirA, noccA));
    SharedTensor2d SvoB = SharedTensor2d(new Tensor2d("PCG Sigma <v|o>", nvirB, noccB));
    SharedTensor2d PvoA = SharedTensor2d(new Tensor2d("PCG P <V|O>", nvirA, noccA));
    SharedTensor2d PvoB = SharedTensor2d(new Tensor2d("PCG P <v|o>", nvirB, noccB));
    SharedTensor1d pQ = SharedTensor1d(new Tensor1d("DF_BASIS_SCF p_Q", nQ_ref));

    // Build sigma
    PvoA->set(p_vecA);
    PvoB->set(p_vecB);

    // Build alpha sigma
    // s_AI = 2 \sum_{B} f_AB p_BI - 2 \sum_{J} f_IJ p_AJ  
    //SvoA->gemm(false, false, FvvA, PvoA, 2.0, 0.0);
    //SvoA->gemm(false, false, PvoA, FooA, -2.0, 1.0);
    #pragma omp parallel for
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              double value = FockA->get(a + noccA, a + noccA) - FockA->get(i,i);
              SvoA->set(a, i, 2.0 * value * PvoA->get(a, i));
         }
    }

    // p_Q = \sum_{BJ} b_BJ^Q p_BJ + \sum_{bj} b_bj^Q p_bj
    // beta contribution
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB, nvirB));
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    SharedTensor2d bQvoB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vo)", nQ_ref, nvirB, noccB));
    bQvoB->swap_3index_col(bQovB);
    bQovB.reset();
    pQ->gemv(false, bQvoB, p_vecB, 1.0, 0.0);
    bQvoB.reset();

    // alpha contribution
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    SharedTensor2d bQvoA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VO)", nQ_ref, nvirA, noccA));
    bQvoA->swap_3index_col(bQovA);
    bQovA.reset();
    pQ->gemv(false, bQvoA, p_vecA, 1.0, 1.0);

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
              sigma_A->set(ai, SvoA->get(a,i));
         }
    }
    PvoA.reset();
    SvoA.reset();

    // Build beta sigma
    // s_ai = 2 \sum_{b} f_ab p_bi - 2 \sum_{j} f_ij p_aj  
    //SvoB->gemm(false, false, FvvB, PvoB, 2.0, 0.0);
    //SvoB->gemm(false, false, PvoB, FooB, -2.0, 1.0);
    #pragma omp parallel for
    for (int a = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++) {
              double value = FockB->get(a + noccB, a + noccB) - FockB->get(i,i);
              SvoB->set(a, i, 2.0 * value * PvoB->get(a, i));
         }
    }

    // s_ai += 4 \sum_{Q} bai^Q p^Q
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB, nvirB));
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    bQvoB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vo)", nQ_ref, nvirB, noccB));
    bQvoB->swap_3index_col(bQovB);
    bQovB.reset();
    SvoB->gemv(true, bQvoB, pQ, 4.0, 1.0);
    pQ.reset();

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
              sigma_B->set(ai, SvoB->get(a,i));
         }
    }

    // memfree
    PvoB.reset();
    SvoB.reset();

}// end sigma_uhf

}} // End Namespaces
