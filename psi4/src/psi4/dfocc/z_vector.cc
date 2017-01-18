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

void DFOCC::z_vector()
{
//outfile->Printf("\n z_vector is starting... \n");

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
    build_rhf_mohess(Aorb);
    timer_off("MO Hessian");

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
        outfile->Printf("\tWarning!!! PCG did NOT converged in %2d iterations. \n", itr_pcg);

    } // end if pcg_conver = 0

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    nidp_tot = nidpA + nidpB;
    //nidp_tot = (nvirA*noccA) + (nvirB*noccB);
    Aorb = SharedTensor2d(new Tensor2d("UHF MO Hessian Matrix", nidp_tot, nidp_tot));
    build_uhf_mohess(Aorb);

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
              zvector->set(ai + nidpA, -WorbB->get(a + noccB, i));
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

    // Build zvector for VO block
    zvectorA = SharedTensor1d(new Tensor1d("Alpha Z-Vector", noccA * nvirA));
    for (int a = 0, ai = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++, ai++) {
	      zvectorA->set(ai, zvector->get(ai));
         }
    }

    // Build zvector for vo block
    zvectorB = SharedTensor1d(new Tensor1d("Beta Z-Vector", noccB * nvirB));
    for (int a = 0, ai = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++, ai++) {
	      zvectorB->set(ai, zvector->get(ai + nidpA));
         }
    }
    zvector.reset();
    //zvectorA->print();
    //zvectorB->print();

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
 //outfile->Printf("\n z_vector done. \n");
}// end z_vector

//=======================================================
//          RHF MO HESSIAN
//=======================================================
void DFOCC::build_rhf_mohess(SharedTensor2d& Aorb_)
{
    SharedTensor2d K;

    // A(ai,bj) = 2 \delta_{ij} f_ab
    #pragma omp parallel for
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              int ai = vo_idxAA->get(a,i);
              for (int b = 0; b < nvirA; b++) {
                   int bi = vo_idxAA->get(b,i);
                   double value = 2.0 * FockA->get(a + noccA, b + noccA);
                   Aorb_->add(ai, bi, value);
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
                   Aorb_->add(ai, aj, value);
              }
         }
    }

    // A(ai,bj) += 8(ai|bj) - 2(aj|bi)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (VO|VO)", nvirA, noccA, nvirA, noccA));
    tei_vovo_chem_ref_directAA(K);
    Aorb_->sort(1432, K, -2.0, 1.0);
    Aorb_->axpy(K, 8.0);
    K.reset();

    // A(ai,bj) += -2(ij|ab)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|VV)", noccA, noccA, nvirA, nvirA));
    tei_oovv_chem_ref_directAA(K);
    Aorb_->sort(3142, K, -2.0, 1.0);
    K.reset();
    if (print_ > 3) Aorb_->print();

    /*
    // Level shifting
    if (level_shift == "TRUE") {
        #pragma omp parallel for
        for (int a = 0; a < nvirA; a++) {
             for (int i = 0; i < noccA; i++) {
                  int ai = vo_idxAA->get(a,i);
                  Aorb_->add(ai, ai, lshift_parameter);
             }
        }
    }
    */

}// end build_rhf_mohess

//=======================================================
//          UHF MO HESSIAN
//=======================================================
void DFOCC::build_uhf_mohess(SharedTensor2d& Aorb_)
{
    SharedTensor2d K;

    // Alpha-Alpha spin cae
    AorbAA = SharedTensor2d(new Tensor2d("MO Hessian Matrix <VO|VO>", nvirA, noccA, nvirA, noccA));

    // A(ai,bj) = 2 \delta_{ij} f_ab => A(ai,bi) = 2 f_ab
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

    // A(ai,bj) += - 2 \delta_{ab} f_ij => A(ai,aj) += -2 f_ij
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
    // AAAA part
    for (int x=0; x<nidpA; x++) {
      for (int y=0; y<nidpA; y++) {
	Aorb_->set(x,y,AorbAA->get(x,y));
      }
    }
    AorbAA.reset();

    // Beta-Beta spin case
    AorbBB = SharedTensor2d(new Tensor2d("MO Hessian Matrix <vo|vo>", nvirB, noccB, nvirB, noccB));

    // A(ai,bj) = 2 \delta_{ij} f_ab => A(ai,bi) = 2 f_ab
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

    // A(ai,bj) += - 2 \delta_{ab} f_ij => A(ai,aj) += -2 f_ij
    #pragma omp parallel for
    for (int a = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++) {
              int ai = vo_idxBB->get(a,i);
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
    for (int x=0; x<nidpB; x++) {
      for (int y=0; y<nidpB; y++) {
	Aorb_->set(x+nidpA, y+nidpA, AorbBB->get(x,y));
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
    for (int x=0; x<nidpA; x++) {
      for (int y=0; y<nidpB; y++) {
	Aorb_->set(x, y+nidpA, AorbAB->get(x,y));
      }
    }

    // BBAA part
    for (int x=0; x<nidpB; x++) {
      for (int y=0; y<nidpA; y++) {
	Aorb_->set(x+nidpA, y, AorbAB->get(y,x));
      }
    }
    AorbAB.reset();
    if (print_ > 3) Aorb_->print();

    /*
    // Level shifting
    if (level_shift == "TRUE") {
        #pragma omp parallel for
        for (int a = 0; a < nidp_tot; a++) {
             for (int i = 0; i < nidp_tot; i++) {
                  Aorb_->add(a, i, lshift_parameter);
             }
        }
    }
    */

}// end build_uhf_mohess

}} // End Namespaces
