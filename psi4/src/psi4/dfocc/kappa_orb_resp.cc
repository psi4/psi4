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

#include "psi4/libqt/qt.h"
#include "defines.h"
#include "dfocc.h"

using namespace psi;
using namespace std;

namespace psi{ namespace dfoccwave{

void DFOCC::kappa_orb_resp()
{
//outfile->Printf("\n kappa_orb_resp is starting... \n");

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
    timer_on("MO Hessian");
    Aorb = SharedTensor2d(new Tensor2d("MO Hessian Matrix", nvirA, noccA, nvirA, noccA));

    // A(ai,bj) = 2 \delta_{ij} f_ab
    #pragma omp parallel for
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              int ai = vo_idxAA->get(a,i);
              for (int b = 0; b < nvirA; b++) {
                   int bi = vo_idxAA->get(b,i);
                   double value = 2.0 * FockA->get(a + noccA, b + noccA);
                   Aorb->add(ai, bi, value);
              }
         }
    }

    // A(ai,bj) -= 2 \delta_{ab} f_ij
    #pragma omp parallel for
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              int ai = vo_idxAA->get(a,i);
              for (int j = 0; j < noccA; j++) {
                   int aj = vo_idxAA->get(a,j);
                   double value = -2.0 * FockA->get(i, j);
                   Aorb->add(ai, aj, value);
              }
         }
    }

    // A(ai,bj) += 8(ai|bj) - 2(aj|bi)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (VO|VO)", nvirA, noccA, nvirA, noccA));
    tei_vovo_chem_ref_directAA(K);
    Aorb->sort(1432, K, -2.0, 1.0);
    Aorb->axpy(K, 8.0);
    K.reset();

    // A(ai,bj) += -2(ij|ab)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|VV)", noccA, noccA, nvirA, nvirA));
    tei_oovv_chem_ref_directAA(K);
    Aorb->sort(3142, K, -2.0, 1.0);
    K.reset();
    if (print_ > 3) Aorb->print();
    timer_off("MO Hessian");

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
    pcg_conver = 1;// means successfull
    timer_on("Orb Resp Solver");
    if (lineq == "CDGESV") Aorb->cdgesv(zvectorA, pcg_conver);
    else if (lineq == "FLIN") {
         double det = 0.0;
         Aorb->lineq_flin(zvectorA, &det);
         if (fabs(det) < DIIS_MIN_DET) {
             outfile->Printf( "Warning!!! MO Hessian matrix is near-singular\n");
             outfile->Printf( "Determinant is %6.3E\n", det);

             pcg_conver = 0;// means unsuccessful
         }
    }
    else if (lineq == "POPLE") Aorb->lineq_pople(zvectorA, 6, cutoff);
    Aorb.reset();
    timer_off("Orb Resp Solver");

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
      //approx_diag_hf_mohess_oo();
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
        // VO Block
        #pragma omp parallel for
        for (int a = 0; a < nvirA; a++) {
              for (int i = 0; i < noccA; i++) {
                   double value = 2.0 * (FockA->get(a + noccA, a + noccA) - FockA->get(i,i));
                   AvoA->set(a, i, value);
              }
        }
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
       outfile->Printf("\tWarning!!! MO Hessian matrix is near-singular, switching to an approximately diagonal Hartree-Fock Hessian. \n");

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

    // A(ai,bj) = 2 \delta_{ij} f_ab
    #pragma omp parallel for
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              int ai = vo_idxAA->get(a,i);
              for (int b = 0; b < nvirA; b++) {
                   int bi = vo_idxAA->get(b,i);
                   double value = 2.0 * FockA->get(a + noccA, b + noccA);
                   AorbAA->add(ai, bi, value);
              }
         }
    }

    // A(ai,bj) += - 2 \delta_{ab} f_ij
    #pragma omp parallel for
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              int ai = vo_idxAA->get(a,i);
              for (int j = 0; j < noccA; j++) {
                   int aj = vo_idxAA->get(a,j);
                   double value = -2.0 * FockA->get(i, j);
                   AorbAA->add(ai, aj, value);
              }
         }
    }

    // A(ai,bj) += 4(ai|bj) - 2(aj|bi)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (VO|VO)", nvirA, noccA, nvirA, noccA));
    tei_vovo_chem_ref_directAA(K);
    AorbAA->sort(1432, K, -2.0, 1.0);
    AorbAA->axpy(K, 4.0);
    K.reset();

    // A(ai,bj) += -2(ij|ab)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|VV)", noccA, noccA, nvirA, nvirA));
    tei_oovv_chem_ref_directAA(K);
    AorbAA->sort(3142, K, -2.0, 1.0);
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

    // A(ai,bj) = 2 \delta_{ij} f_ab
    #pragma omp parallel for
    for (int a = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++) {
              int ai = vo_idxBB->get(a,i);
              for (int b = 0; b < nvirB; b++) {
                   int bi = vo_idxBB->get(b,i);
                   double value = 2.0 * FockB->get(a + noccB, b + noccB);
                   AorbBB->add(ai, bi, value);
              }
         }
    }

    // A(ai,bj) -= 2 \delta_{ab} f_ij
    #pragma omp parallel for
    for (int a = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++) {
              int ai = vo_idxBB->get(a,i);
              int ia = ov_idxBB->get(i,a);
              for (int j = 0; j < noccB; j++) {
                   int aj = vo_idxBB->get(a,j);
                   double value = -2.0 * FockB->get(i, j);
                   AorbBB->add(ai, aj, value);
              }
         }
    }

    // A(ai,bj) += 4(ai|bj) - 2(aj|bi)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (vo|vo)", nvirB, noccB, nvirB, noccB));
    tei_vovo_chem_ref_directBB(K);
    AorbBB->sort(1432, K, -2.0, 1.0);
    AorbBB->axpy(K, 4.0);
    K.reset();

    // A(ai,bj) += -2(ij|ab)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (oo|vv)", noccB, noccB, nvirB, nvirB));
    tei_oovv_chem_ref_directBB(K);
    AorbBB->sort(3142, K, -2.0, 1.0);
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

    // A(AI,bj) = 4(AI|bj)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (VO|vo)", nvirA, noccA, nvirB, noccB));
    tei_vovo_chem_ref_directAB(K);
    AorbAB->axpy(K, 4.0);
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
    #pragma omp parallel for
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              int ai = vo_idxAA->get(a,i);
              zvector->set(ai, -WorbA->get(a + noccA, i));
         }
    }
    #pragma omp parallel for
    for (int a = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++) {
              int ai = vo_idxBB->get(a,i);
              zvector->set(ai + nvoAA, -WorbB->get(a + noccB, i));
         }
    }

    // Solve the orb-resp equations
    pcg_conver = 1;// means successfull
    if (lineq == "CDGESV") Aorb->cdgesv(zvector, pcg_conver);
    else if (lineq == "FLIN") {
         double det = 0.0;
         Aorb->lineq_flin(zvector, &det);
         if (fabs(det) < DIIS_MIN_DET) {
             outfile->Printf( "Warning!!! MO Hessian matrix is near-singular\n");
             outfile->Printf( "Determinant is %6.3E\n", det);

             pcg_conver = 0;// means unsuccessful
         }
    }
    else if (lineq == "POPLE") Aorb->lineq_pople(zvector, 6, cutoff);
    Aorb.reset();

    // Build kappa for VO block
    #pragma omp parallel for
    for (int x = 0; x < nidpA; x++) {
         int p = idprowA->get(x);
	 int q = idpcolA->get(x);
         if (p >= noccA && q < noccA) {
             int ai = vo_idxAA->get(p-noccA,q);
	     kappaA->set(x, zvector->get(ai));
         }
    }

    // Build kappa for vo block
    #pragma omp parallel for
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
      // Compute OO-Block orb rot params
      //approx_diag_hf_mohess_oo();
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
         // VO Block
         #pragma omp parallel for
         for (int a = 0; a < nvirA; a++) {
              for (int i = 0; i < noccA; i++) {
                   double value = 2.0 * (FockA->get(a + noccA, a + noccA) - FockA->get(i,i));
                   AvoA->set(a, i, value);
              }
         }

         // vo Block
         #pragma omp parallel for
         for (int a = 0; a < nvirB; a++) {
              for (int i = 0; i < noccB; i++) {
                   double value = 2.0 * (FockB->get(a + noccB, a + noccB) - FockB->get(i,i));
                   AvoB->set(a, i, value);
              }
         }

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

       outfile->Printf("\tWarning!!! MO Hessian matrix is near-singular, switching to MSD. \n");

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
 //outfile->Printf("\n kappa_orb_resp done. \n");
}// end kappa_orb_resp
}} // End Namespaces
