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

#include <libqt/qt.h>
#include "defines.h"
#include "dfocc.h"

using namespace psi;
using namespace std;

namespace psi{ namespace dfoccwave{

void DFOCC::kappa_orb_resp()
{ 
//fprintf(outfile,"\n kappa_orb_resp is starting... \n"); fflush(outfile);

    SharedTensor2d K;

if (reference_ == "RESTRICTED") {
    // Set the kappa to -negative of the mo grad
    zvectorA = SharedTensor1d(new Tensor1d("Alpha Z-Vector", noccA * nvirA));
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              int ai = vo_idxAA->get(a,i);
              zvectorA->set(ai, -WorbA->get(a + noccA, i));
         }
    }

    // Build the MO Hessian
    Aorb = SharedTensor2d(new Tensor2d("MO Hessian Matrix", nvirA, noccA, nvirA, noccA));

    // A(ai,bj) = 2 \delta_{ij} f_ab - 2 \delta_{ab} f_ij + 8(ia|jb) - 2(ib|ja)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|OV)", noccA, nvirA, noccA, nvirA));
    timer_on("I/O");
    K->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    #pragma omp parallel for
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              int ai = vo_idxAA->get(a,i);
              int ia = ov_idxAA->get(i,a);
              for (int b = 0; b < nvirA; b++) {
                   int ib = ov_idxAA->get(i,b);
                   for (int j = 0; j < noccA; j++) {
                        int bj = vo_idxAA->get(b,j);
                        int jb = ov_idxAA->get(j,b);
                        int ja = ov_idxAA->get(j,a);
                        double value = (8.0 * K->get(ia,jb)) - (2.0 * K->get(ib,ja)); 
                        if (i == j) value += 2.0 * FockA->get(a + noccA, b + noccA);
                        if (a == b) value -= 2.0 * FockA->get(i, j);
                        Aorb->set(ai, bj, value);
                   }
              }
         }
    }
    K.reset();

    // A(ai,bj) += -2(ij|ab)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|VV)", noccA, noccA, nvirA, nvirA));
    timer_on("I/O");
    K->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    #pragma omp parallel for
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              int ai = vo_idxAA->get(a,i);
              for (int b = 0; b < nvirA; b++) {
                   for (int j = 0; j < noccA; j++) {
                        int bj = vo_idxAA->get(b,j);
                        int ij = oo_idxAA->get(i,j);
                        int ab = vv_idxAA->get(a,b);
                        Aorb->add(ai, bj, -2.0*K->get(ij,ab));
                   }
              }
         }
    }
    K.reset();
    if (print_ > 3) Aorb->print();

    /*
    // Level shifting
    if (level_shift == "TRUE") {
        #pragma omp parallel for
        for (int a = 0; a < nvirA; a++) {
             for (int i = 0; i < noccA; i++) {
                  int ai = vo_idxAA->get(a,i);
                  Aorb->add(ai, ai, lshift_parameter);
             }
        }
    }
    */

    // Solve the orb-resp equations
    pcg_conver = 0;// here 0 means successfull
    if (lineq == "CDGESV") Aorb->cdgesv(zvectorA, pcg_conver);
    else if (lineq == "FLIN") {
         double det = 0.0;      
         Aorb->lineq_flin(zvectorA, &det);
         if (fabs(det) < DIIS_MIN_DET) { 
             fprintf(outfile, "Warning!!! MO Hessian matrix is near-singular\n");
             fprintf(outfile, "Determinant is %6.3E\n", det);
             fflush(outfile);
             pcg_conver = 1;// here 1 means unsuccessful
         }
    }
    else if (lineq == "POPLE") Aorb->lineq_pople(zvectorA, 6, cutoff);
    Aorb.reset();

    // Build kappa for VO block
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
        for (int i = 0; i < naoccA; i++) {
              for (int j = 0; j < nfrzc; j++) {
                   double value = 2.0 * msd_oo_scale * (FockA->get(i + nfrzc, i + nfrzc) - FockA->get(j,j));
                   AooA->set(i, j, value);
              }
        }
     
      // Compute OO-Block orb rot params
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
    if (pcg_conver != 0) {
        // VO Block
        for (int a = 0; a < nvirA; a++) {
              for (int i = 0; i < noccA; i++) {
                   double value = 2.0 * (FockA->get(a + noccA, a + noccA) - FockA->get(i,i));
                   AvoA->set(a, i, value);
              }
        }
        // Build kappa again
        for (int x = 0; x < nidpA; x++) {
	    int p = idprowA->get(x);
	    int q = idpcolA->get(x);
            double value = 0.0;
            if (p >= noccA && q < noccA) value = AvoA->get(p-noccA,q); 
            else if (p < noccA && q < noccA) value = AooA->get(p-nfrzc,q); 
	    kappaA->set(x, -wogA->get(x)/value);
        }
       fprintf(outfile,"\tWarning!!! MO Hessian matrix is near-singular, switching to an approximately diagonal Hartree-Fock Hessian. \n");
       fflush(outfile);
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
    //nidp_tot = nidpA + nidpB;
    nidp_tot = (nvirA*noccA) + (nvirB*noccB);

    // Alpha-Alpha spin cae
    AorbAA = SharedTensor2d(new Tensor2d("MO Hessian Matrix <VO|VO>", nvirA, noccA, nvirA, noccA));

    // A(ai,bj) = 2 \delta_{ij} f_ab - 2 \delta_{ab} f_ij + 4(ia|jb) - 2(ib|ja)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|OV)", noccA, nvirA, noccA, nvirA));
    timer_on("I/O");
    K->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    #pragma omp parallel for
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              int ai = vo_idxAA->get(a,i);
              int ia = ov_idxAA->get(i,a);
              for (int b = 0; b < nvirA; b++) {
                   int ib = ov_idxAA->get(i,b);
                   for (int j = 0; j < noccA; j++) {
                        int bj = vo_idxAA->get(b,j);
                        int jb = ov_idxAA->get(j,b);
                        int ja = ov_idxAA->get(j,a);
                        double value = (4.0 * K->get(ia,jb)) - (2.0 * K->get(ib,ja)); 
                        if (i == j) value += 2.0 * FockA->get(a + noccA, b + noccA);
                        if (a == b) value -= 2.0 * FockA->get(i, j);
                        AorbAA->set(ai, bj, value);
                   }
              }
         }
    }
    K.reset();

    // A(ai,bj) += -2(ij|ab)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|VV)", noccA, noccA, nvirA, nvirA));
    timer_on("I/O");
    K->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    #pragma omp parallel for
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              int ai = vo_idxAA->get(a,i);
              for (int b = 0; b < nvirA; b++) {
                   for (int j = 0; j < noccA; j++) {
                        int bj = vo_idxAA->get(b,j);
                        int ij = oo_idxAA->get(i,j);
                        int ab = vv_idxAA->get(a,b);
                        AorbAA->add(ai, bj, -2.0*K->get(ij,ab));
                   }
              }
         }
    }
    K.reset();
    if (print_ > 3) AorbAA->print();

    // Build the UHF MO Hessian matrix
    Aorb = SharedTensor2d(new Tensor2d("UHF MO Hessian Matrix", nidp_tot, nidp_tot));
    // AAAA part 
    for (int x=0; x<nvoAA; x++) { 
      for (int y=0; y<nvoAA; y++) { 
	Aorb->set(x,y,AorbAA->get(x,y));
      }
    }
    AorbAA.reset();

    // Beta-Beta spin case
    AorbBB = SharedTensor2d(new Tensor2d("MO Hessian Matrix <vo|vo>", nvirB, noccB, nvirB, noccB));

    // A(ai,bj) = 2 \delta_{ij} f_ab - 2 \delta_{ab} f_ij + 4(ia|jb) - 2(ib|ja)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (ov|ov)", noccB, nvirB, noccB, nvirB));
    timer_on("I/O");
    K->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    #pragma omp parallel for
    for (int a = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++) {
              int ai = vo_idxBB->get(a,i);
              int ia = ov_idxBB->get(i,a);
              for (int b = 0; b < nvirB; b++) {
                   int ib = ov_idxBB->get(i,b);
                   for (int j = 0; j < noccB; j++) {
                        int bj = vo_idxBB->get(b,j);
                        int jb = ov_idxBB->get(j,b);
                        int ja = ov_idxBB->get(j,a);
                        double value = (4.0 * K->get(ia,jb)) - (2.0 * K->get(ib,ja));
                        if (i == j) value += 2.0 * FockB->get(a + noccB, b + noccB);
                        if (a == b) value -= 2.0 * FockB->get(i, j);
                        AorbBB->set(ai, bj, value);
                   }
              }
         }
    }
    K.reset();

    // A(ai,bj) += -2(ij|ab)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (oo|vv)", noccB, noccB, nvirB, nvirB));
    timer_on("I/O");
    K->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    #pragma omp parallel for
    for (int a = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++) {
              int ai = vo_idxBB->get(a,i);
              for (int b = 0; b < nvirB; b++) {
                   for (int j = 0; j < noccB; j++) {
                        int bj = vo_idxBB->get(b,j);
                        int ij = oo_idxBB->get(i,j);
                        int ab = vv_idxBB->get(a,b);
                        AorbBB->add(ai, bj, -2.0*K->get(ij,ab));
                   }
              }
         }
    }
    K.reset();
    if (print_ > 3) AorbBB->print();

    // Build the UHF MO Hessian matrix
    // BBBB part 
    for (int x=0; x<nvoBB; x++) { 
      for (int y=0; y<nvoBB; y++) { 
	Aorb->set(x+nvoAA, y+nvoAA, AorbBB->get(x,y));
      }
    }
    AorbBB.reset();

    // Alpha-Beta spin cae
    AorbAB = SharedTensor2d(new Tensor2d("MO Hessian Matrix <VO|vo>", nvirA, noccA, nvirB, noccB));

    // A(AI,bj) = 4(IA|jb)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|ov)", noccA, nvirA, noccB, nvirB));
    timer_on("I/O");
    K->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    #pragma omp parallel for
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              int ai = vo_idxAA->get(a,i);
              int ia = ov_idxAA->get(i,a);
              for (int b = 0; b < nvirB; b++) {
                   for (int j = 0; j < noccB; j++) {
                        int bj = vo_idxBB->get(b,j);
                        int jb = ov_idxBB->get(j,b);
                        AorbAB->set(ai, bj, 4.0*K->get(ia,jb));
                   }
              }
         }
    }
    K.reset();
    if (print_ > 3) AorbAB->print();

    // Build the UHF MO Hessian matrix
    // AABB part 
    for (int x=0; x<nvoAA; x++) { 
      for (int y=0; y<nvoBB; y++) { 
	Aorb->set(x, y+nvoAA, AorbAB->get(x,y));
      }
    }
    
    // BBAA part 
    for (int x=0; x<nvoBB; x++) { 
      for (int y=0; y<nvoAA; y++) { 
	Aorb->set(x+nvoAA, y, AorbAB->get(y,x));
      }
    }
    AorbAB.reset();
    if (print_ > 3) Aorb->print();

    /*
    // Level shifting
    if (level_shift == "TRUE") {
        #pragma omp parallel for
        for (int a = 0; a < nidp_tot; a++) {
             for (int i = 0; i < nidp_tot; i++) {
                  Aorb->add(a, i, lshift_parameter);
             }
        }
    }
    */

    // Build total zvector
    zvector = SharedTensor1d(new Tensor1d("UHF Z-Vector", nidp_tot));
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              int ai = vo_idxAA->get(a,i);
              zvector->set(ai, -WorbA->get(a + noccA, i));
         }
    }
    for (int a = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++) {
              int ai = vo_idxBB->get(a,i);
              zvector->set(ai + nvoAA, -WorbB->get(a + noccB, i));
         }
    }

    // Solve the orb-resp equations
    pcg_conver = 0;// here 0 means successfull
    if (lineq == "CDGESV") Aorb->cdgesv(zvector, pcg_conver);
    else if (lineq == "FLIN") {
         double det = 0.0;      
         Aorb->lineq_flin(zvector, &det);
         if (fabs(det) < DIIS_MIN_DET) { 
             fprintf(outfile, "Warning!!! MO Hessian matrix is near-singular\n");
             fprintf(outfile, "Determinant is %6.3E\n", det);
             fflush(outfile);
             pcg_conver = 1;// here 1 means unsuccessful
         }
    }
    else if (lineq == "POPLE") Aorb->lineq_pople(zvector, 6, cutoff);
    Aorb.reset();

    // Build kappa for VO block
    for (int x = 0; x < nidpA; x++) {
         int p = idprowA->get(x);
	 int q = idpcolA->get(x);
         if (p >= noccA && q < noccA) {
             int ai = vo_idxAA->get(p-noccA,q);
	     kappaA->set(x, zvector->get(ai));
         }
    }

    // Build kappa for vo block
    for (int x = 0; x < nidpB; x++) {
         int p = idprowB->get(x);
	 int q = idpcolB->get(x);
         if (p >= noccB && q < noccB) {
             int ai = vo_idxBB->get(p-noccB,q);
	     kappaB->set(x, zvector->get(ai + nvoAA));
         }
    }
    zvector.reset();

    if (nfrzc > 0) {
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

      // Compute OO-Block orb rot params
      for (int x = 0; x < nidpA; x++) {
	   int p = idprowA->get(x);
	   int q = idpcolA->get(x);
           if (p < noccA && q < noccA) {
               double value = AooA->get(p-nfrzc,q); 
	       kappaA->set(x, -wogA->get(x)/value);
           }
      }
      // Compute oo-Block orb rot params
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
    if (pcg_conver != 0) {
         // VO Block
         for (int a = 0; a < nvirA; a++) {
              for (int i = 0; i < noccA; i++) {
                   double value = 2.0 * (FockA->get(a + noccA, a + noccA) - FockA->get(i,i));
                   AvoA->set(a, i, value);
              }
         }

         // vo Block
         for (int a = 0; a < nvirB; a++) {
              for (int i = 0; i < noccB; i++) {
                   double value = 2.0 * (FockB->get(a + noccB, a + noccB) - FockB->get(i,i));
                   AvoB->set(a, i, value);
              }
         }

	// alpha
        for(int x = 0; x < nidpA; x++) {
	    int p = idprowA->get(x);
	    int q = idpcolA->get(x);
            double value = 0.0;
            if (p >= noccA && q < noccA) value = AvoA->get(p-noccA,q); 
            else if (p < noccA && q < noccA) value = AooA->get(p-nfrzc,q); 
	    kappaA->set(x, -wogA->get(x)/value);
        }
	
	// beta
        for(int x = 0; x < nidpB; x++) {
	    int p = idprowB->get(x);
	    int q = idpcolB->get(x);
            double value = 0.0;
            if (p >= noccB && q < noccB) value = AvoB->get(p-noccB,q); 
            else if (p < noccB && q < noccB) value = AooB->get(p-nfrzc,q); 
	    kappaB->set(x, -wogB->get(x)/value);
        }

       fprintf(outfile,"\tWarning!!! MO Hessian matrix is near-singular, switching to MSD. \n");
       fflush(outfile);
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
 //fprintf(outfile,"\n kappa_orb_resp done. \n"); fflush(outfile);
}// end kappa_orb_resp
}} // End Namespaces


