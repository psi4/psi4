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

void DFOCC::approx_diag_hf_mohess_vo()
{ 
     double value;
     if (reference_ == "RESTRICTED") {
         // VO Block
         for (int a = 0; a < nvirA; a++) {
              for (int i = 0; i < noccA; i++) {
                   value = 2.0 * (FockA->get(a + noccA, a + noccA) - FockA->get(i,i));
                   if (regularization == "TRUE") value += reg_param;
                   AvoA->set(a, i, value);
              }
         }
     } // end if (reference_ == "RESTRICTED") 

     else if (reference_ == "UNRESTRICTED") {
         // VO Block
         for (int a = 0; a < nvirA; a++) {
              for (int i = 0; i < noccA; i++) {
                   value = 2.0 * (FockA->get(a + noccA, a + noccA) - FockA->get(i,i));
                   if (regularization == "TRUE") value += reg_param;
                   AvoA->set(a, i, value);
              }
         }

         // vo Block
         for (int a = 0; a < nvirB; a++) {
              for (int i = 0; i < noccB; i++) {
                   value = 2.0 * (FockB->get(a + noccB, a + noccB) - FockB->get(i,i));
                   if (regularization == "TRUE") value += reg_param;
                   AvoB->set(a, i, value);
              }
         }
     }// end else if (reference_ == "UNRESTRICTED")
}//

//=========================
// APPROX_DIAG_HF OO
//=========================
void DFOCC::approx_diag_hf_mohess_oo()
{ 
     if (reference_ == "RESTRICTED") {
         // OO Block
         for (int i = 0; i < naoccA; i++) {
              for (int j = 0; j < nfrzc; j++) {
                   double value = 2.0 * msd_oo_scale * (FockA->get(i + nfrzc, i + nfrzc) - FockA->get(j,j));
                   AooA->set(i, j, value);
              }
         }
     } // end if (reference_ == "RESTRICTED") 

     else if (reference_ == "UNRESTRICTED") {
         // OO Block
         for (int i = 0; i < naoccA; i++) {
              for (int j = 0; j < nfrzc; j++) {
                   double value = 2.0 * msd_oo_scale * (FockA->get(i + nfrzc, i + nfrzc) - FockA->get(j,j));
                   AooA->set(i, j, value);
              }
         }

         // oo Block
         for (int i = 0; i < naoccB; i++) {
              for (int j = 0; j < nfrzc; j++) {
                   double value = 2.0 * msd_oo_scale * (FockB->get(i + nfrzc, i + nfrzc) - FockB->get(j,j));
                   AooB->set(i, j, value);
              }
         }
     }// end else if (reference_ == "UNRESTRICTED")
}//

//=========================
// APPROX_DIAG_EKT VO
//=========================
void DFOCC::approx_diag_ekt_mohess_vo()
{ 
     double value;
     if (reference_ == "RESTRICTED") {
         // VO Block
         for (int a = 0; a < nvirA; a++) {
              for (int i = 0; i < noccA; i++) {
                   value = GFtvv->get(a, a) - GF->get(i,i);
                   if (regularization == "TRUE") value += reg_param;
                   AvoA->set(a, i, value);
              }
         }
     } // end if (reference_ == "RESTRICTED") 

     else if (reference_ == "UNRESTRICTED") {
         // VO Block
         for (int a = 0; a < nvirA; a++) {
              for (int i = 0; i < noccA; i++) {
                   value = 2.0 * (GFtvvA->get(a, a) - GFA->get(i,i));
                   if (regularization == "TRUE") value += reg_param;
                   AvoA->set(a, i, value);
              }
         }

         // vo Block
         for (int a = 0; a < nvirB; a++) {
              for (int i = 0; i < noccB; i++) {
                   value = 2.0 * (GFtvvB->get(a, a) - GFB->get(i,i));
                   if (regularization == "TRUE") value += reg_param;
                   AvoB->set(a, i, value);
              }
         }
     }// end else if (reference_ == "UNRESTRICTED")
}//

//=========================
// APPROX_DIAG_EKT OO
//=========================
void DFOCC::approx_diag_ekt_mohess_oo()
{ 
     if (reference_ == "RESTRICTED") {
         // OO Block
         for (int i = 0; i < naoccA; i++) {
              for (int j = 0; j < nfrzc; j++) {
                   double value = msd_oo_scale * (GF->get(i + nfrzc, i + nfrzc) - GF->get(j,j));
                   AooA->set(i, j, value);
              }
         }
     } // end if (reference_ == "RESTRICTED") 

     else if (reference_ == "UNRESTRICTED") {
         // OO Block
         for (int i = 0; i < naoccA; i++) {
              for (int j = 0; j < nfrzc; j++) {
                   double value = 2.0 * msd_oo_scale * (GFA->get(i + nfrzc, i + nfrzc) - GFA->get(j,j));
                   AooA->set(i, j, value);
              }
         }

         // oo Block
         for (int i = 0; i < naoccB; i++) {
              for (int j = 0; j < nfrzc; j++) {
                   double value = 2.0 * msd_oo_scale * (GFB->get(i + nfrzc, i + nfrzc) - GFB->get(j,j));
                   AooB->set(i, j, value);
              }
         }
     }// end else if (reference_ == "UNRESTRICTED")
}//

//=========================
// kappa_diag_hess
//=========================
void DFOCC::kappa_diag_hess()
{ 
//outfile->Printf("\n kappa_diag_hess is starting... \n"); 
        double value;

 if (hess_type == "APPROX_DIAG") {
      approx_diag_mohess_vo();
      if (nfrzc > 0) approx_diag_mohess_oo();
 }

 else if (hess_type == "APPROX_DIAG_EKT") {
      approx_diag_ekt_mohess_vo();
      if (nfrzc > 0) approx_diag_ekt_mohess_oo();
 }

 else if (hess_type == "APPROX_DIAG_HF") {
      approx_diag_hf_mohess_vo();
      if (nfrzc > 0) approx_diag_hf_mohess_oo();
 }

 else if (hess_type == "DIAG") {
      diagonal_mohess_vo();
      if (nfrzc > 0) diagonal_mohess_oo();
 }

// Kappa
if (reference_ == "RESTRICTED") {
        // Get kappa
        for(int x = 0; x < nidpA; x++) {
	    int p = idprowA->get(x);
	    int q = idpcolA->get(x);
            if (p >= noccA && q < noccA) value = AvoA->get(p-noccA,q); 
            else if (p < noccA && q < noccA) value = AooA->get(p-nfrzc,q); 
	    kappaA->set(x, -wogA->get(x)/value);
        }

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
        // Get kappa
	// alpha
        for(int x = 0; x < nidpA; x++) {
	    int p = idprowA->get(x);
	    int q = idpcolA->get(x);
            if (p >= noccA && q < noccA) value = AvoA->get(p-noccA,q); 
            else if (p < noccA && q < noccA) value = AooA->get(p-nfrzc,q); 
	    kappaA->set(x, -wogA->get(x)/value);
        }
	
	// beta
        for(int x = 0; x < nidpB; x++) {
	    int p = idprowB->get(x);
	    int q = idpcolB->get(x);
            if (p >= noccB && q < noccB) value = AvoB->get(p-noccB,q); 
            else if (p < noccB && q < noccB) value = AooB->get(p-nfrzc,q); 
	    kappaB->set(x, -wogB->get(x)/value);
        }

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
 //outfile->Printf("\n kappa_diag_hess done. \n"); 
	
}// end main
}} // End Namespaces
