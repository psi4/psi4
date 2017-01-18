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

void DFOCC::z_vector_cg()
{ 
//outfile->Printf("\n z_vector_cg is starting... \n"); 
    outfile->Printf("\tSolving orbital Z-vector equations...\n");
    

    SharedTensor2d K, L;

if (reference_ == "RESTRICTED") {
    // Memalloc
    zvectorA = SharedTensor1d(new Tensor1d("Alpha Z-Vector", noccA * nvirA));
    zvector = SharedTensor1d(new Tensor1d("Alpha Z-Vector", noccA * nvirA));
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
    Aorb = SharedTensor2d(new Tensor2d("MO Hessian Matrix", nvirA, noccA, nvirA, noccA));
    build_rhf_mohess(Aorb);
    sigma_pcgA->gemv(false, Aorb, zvectorA, 1.0, 0.0);
    Aorb.reset();

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
    zvector->copy(zvectorA);
    cg_solver();
    zvectorA->copy(zvector);

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
    zvector.reset();

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
    zvec_newA = SharedTensor1d(new Tensor1d("New UHF Z-Vector", nidp_tot));
    Minv_pcgA = SharedTensor1d(new Tensor1d("PCG M inverse", nidp_tot));
    sigma_pcgA = SharedTensor1d(new Tensor1d("PCG sigma", nidp_tot));
    r_pcgA = SharedTensor1d(new Tensor1d("PCG r", nidp_tot));
    r_pcg_newA = SharedTensor1d(new Tensor1d("PCG new r", nidp_tot));
    z_pcgA = SharedTensor1d(new Tensor1d("PCG z", nidp_tot));
    z_pcg_newA = SharedTensor1d(new Tensor1d("PCG new z", nidp_tot));
    p_pcgA = SharedTensor1d(new Tensor1d("PCG p", nidp_tot));
    p_pcg_newA = SharedTensor1d(new Tensor1d("PCG new p", nidp_tot));
    dr_pcgA = SharedTensor1d(new Tensor1d("PCG dr", nidp_tot));
    residualA = SharedTensor1d(new Tensor1d("Residual Vector", nidp_tot));
    zvectorA = SharedTensor1d(new Tensor1d("Alpha Z-Vector", noccA * nvirA));
    zvectorB = SharedTensor1d(new Tensor1d("Beta Z-Vector", noccB * nvirB));

    // Build kappa0 and M
    // alpha
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
              double value = FockA->get(a + noccA, a + noccA) - FockA->get(i,i);
              zvector->set(ai, -WorbA->get(a + noccA, i) / (2.0*value));       
              Minv_pcgA->set(ai, 0.5/value);
         }
    }

    // beta
    for (int a = 0, ai = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++, ai++) {
              double value = FockB->get(a + noccB, a + noccB) - FockB->get(i,i);
              zvector->set(ai + nidpA, -WorbB->get(a + noccB, i) / (2.0*value));       
              Minv_pcgA->set(ai + nidpA, 0.5/value);
         }
    }

    // Build sigma
    Aorb = SharedTensor2d(new Tensor2d("UHF MO Hessian Matrix", nidp_tot, nidp_tot));
    build_uhf_mohess(Aorb);
    sigma_pcgA->gemv(false, Aorb, zvector, 1.0, 0.0);
    Aorb.reset();

    // Build r0
    // alpha
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
              residualA->set(ai, -WorbA->get(a + noccA, i));       
         }
    }
    // beta
    for (int a = 0, ai = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++, ai++) {
              residualA->set(ai + nidpA, -WorbB->get(a + noccB, i));       
         }
    }
    residualA->subtract(sigma_pcgA);
    r_pcgA->copy(residualA);

    // Build z0
    z_pcgA->dirprd(Minv_pcgA, r_pcgA);

    // Build p0
    p_pcgA->copy(z_pcgA);

    // Call Orbital Response Solver
    cg_solver();

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
    residualA.reset();

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
        
    } 

}// end if (reference_ == "UNRESTRICTED") 
 //outfile->Printf("\n z_vector_cg done. \n"); 
}// end z_vector_cg

//=======================================================
//          CG (RHF)
//=======================================================          
void DFOCC::cg_solver()
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
    if (reference_ == "RESTRICTED") {
        Aorb = SharedTensor2d(new Tensor2d("MO Hessian Matrix", nvirA, noccA, nvirA, noccA));
        build_rhf_mohess(Aorb);
        sigma_pcgA->gemv(false, Aorb, p_pcgA, 1.0, 0.0);
    }

    else if (reference_ == "UNRESTRICTED") {
        Aorb = SharedTensor2d(new Tensor2d("UHF MO Hessian Matrix", nidp_tot, nidp_tot));
        build_uhf_mohess(Aorb);
        sigma_pcgA->gemv(false, Aorb, p_pcgA, 1.0, 0.0);
    }

    // Compute line search parameter alpha
    a_pcgA = r_pcgA->dot(z_pcgA) / p_pcgA->dot(sigma_pcgA);

   // Build kappa-new
   zvec_newA->zero();
   zvec_newA->copy(p_pcgA);
   zvec_newA->scale(a_pcgA);
   zvec_newA->add(zvector);

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
   rms_pcg = zvec_newA->rms(zvector);
   rms_r_pcgA = r_pcg_newA->rms();

   // Reset
   zvector->copy(zvec_newA);
   r_pcgA->copy(r_pcg_newA);
   z_pcgA->copy(z_pcg_newA);
   p_pcgA->copy(p_pcg_newA);

    // Build sigma for kappa
    sigma_pcgA->gemv(false, Aorb, zvector, 1.0, 0.0);
    Aorb.reset();

    // Build r0
    // alpha
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
              residualA->set(ai, -WorbA->get(a + noccA, i));       
         }
    }

    if (reference_ == "UNRESTRICTED") {
    // beta
    for (int a = 0, ai = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++, ai++) {
              residualA->set(ai + nidpA, -WorbB->get(a + noccB, i));       
         }
    }}
    residualA->subtract(sigma_pcgA);
    rms_residual = residualA->rms();

   // increment iteration index 
   itr_pcg++;

   // Print
   outfile->Printf("\t%3d     %12.2e     %12.2e\n",itr_pcg,rms_pcg,rms_residual);
   

   // If we exceed maximum number of iteration, break the loop
   if (itr_pcg >= pcg_maxiter) {
       pcg_conver = 0;
       break;
   }  

   //if (rms_residual < tol_pcg  || rms_pcg < tol_pcg) break;  

 }
 while(fabs(rms_pcg) >= tol_pcg || fabs(rms_residual) >= tol_pcg);  

    // Converged?
    outfile->Printf("\n");
    
}// cg_solver

}} // End Namespaces
