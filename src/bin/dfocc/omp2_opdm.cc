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

/** Standard library includes */
#include <libqt/qt.h>
#include "defines.h"
#include "dfocc.h"

using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace dfoccwave{
  
void DFOCC::omp2_opdm()
{   

    timer_on("opdm");
if (reference_ == "RESTRICTED") {
    // Tensors 
    SharedTensor2d t2;
    SharedTensor2d l2;
    t2 = SharedTensor2d(new Tensor2d("T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    l2 = SharedTensor2d(new Tensor2d("T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    t2->read(psio_, PSIF_DFOCC_AMPS);
    l2->read(psio_, PSIF_DFOCC_AMPS);

    // G_ij = 2\sum_{m,e,f} t_im^ef t_jm^ef    
    GijA->contract442(1, 1, t2, l2, 2.0, 0.0);

    // G_ij += -\sum_{m,e,f} t_im^ef t_mj^ef    
    GijA->contract442(1, 2, t2, l2, -1.0, 1.0);

    // G_ab = -2\sum_{m,n,e} t_mn^ea t_mn^eb     
    GabA->contract442(4, 4, t2, l2, -2.0, 0.0);

    // G_ab += \sum_{m,n,e} t_mn^ea t_mn^be    
    GabA->contract442(4, 3, t2, l2, 1.0, 1.0);
    t2.reset();
    l2.reset();

    // Build G1c_oo and G1c_vv
    G1c_oo->set_act_oo(nfrzc, naoccA, GijA);
    G1c_oo->scale(-2.0);
    G1c_vv->set_act_vv(GabA);
    G1c_vv->scale(-2.0);

    //G1c_oo->print();
    //G1c_vv->print();

    // Build G1c
    G1c->set_oo(G1c_oo);
    G1c->set_vv(noccA, G1c_vv);

    // Build G1
    G1->copy(G1c);
    for (int i = 0; i < noccA; i++) G1->add(i, i, 2.0); 

  if(print_ > 2) {
    G1->print();
    double trace = G1->trace();
    fprintf(outfile,"\t trace: %12.12f \n", trace);
    fflush(outfile);
  }

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // tensors
    SharedTensor2d t2;
    SharedTensor2d l2;

    // G_IJ = 1/2 \sum_{M,E,F} t_IM^EF t_JM^EF
    t2 = SharedTensor2d(new Tensor2d("T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    l2 = SharedTensor2d(new Tensor2d("T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    t2->read(psio_, PSIF_DFOCC_AMPS);
    l2->read(psio_, PSIF_DFOCC_AMPS);
    GijA->contract442(1, 1, t2, l2, 0.5, 0.0);
    t2.reset();
    l2.reset();

    // G_IJ = \sum_{m,E,f} t_Im^Ef t_Jm^Ef
    t2 = SharedTensor2d(new Tensor2d("T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    l2 = SharedTensor2d(new Tensor2d("T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    t2->read(psio_, PSIF_DFOCC_AMPS);
    l2->read(psio_, PSIF_DFOCC_AMPS);
    GijA->contract442(1, 1, t2, l2, 1.0, 1.0);
    t2.reset();
    l2.reset();

    // G_ij = 1/2 \sum_{m,e,f} t_im^ef t_jm^ef
    t2 = SharedTensor2d(new Tensor2d("T2_1 <ij|ab>", naoccB, naoccB, navirB, navirB));
    l2 = SharedTensor2d(new Tensor2d("T2_1 <ij|ab>", naoccB, naoccB, navirB, navirB));
    t2->read(psio_, PSIF_DFOCC_AMPS);
    l2->read(psio_, PSIF_DFOCC_AMPS);
    GijB->contract442(1, 1, t2, l2, 0.5, 0.0);
    t2.reset();
    l2.reset();

    // G_ij = \sum_{M,e,F} t_Mi^Fe t_Mj^Fe
    t2 = SharedTensor2d(new Tensor2d("T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    l2 = SharedTensor2d(new Tensor2d("T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    t2->read(psio_, PSIF_DFOCC_AMPS);
    l2->read(psio_, PSIF_DFOCC_AMPS);
    GijB->contract442(2, 2, t2, l2, 1.0, 1.0);
    t2.reset();
    l2.reset();

    // G_AB = -1/2\sum_{M,N,F} t_MN^FA t_MN^FB
    t2 = SharedTensor2d(new Tensor2d("T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    l2 = SharedTensor2d(new Tensor2d("T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    t2->read(psio_, PSIF_DFOCC_AMPS);
    l2->read(psio_, PSIF_DFOCC_AMPS);
    GabA->contract442(4, 4, t2, l2, -0.5, 0.0);
    t2.reset();
    l2.reset();

    // G_AB += -\sum_{M,n,f} t_Mn^Af t_Mn^Bf
    t2 = SharedTensor2d(new Tensor2d("T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    l2 = SharedTensor2d(new Tensor2d("T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    t2->read(psio_, PSIF_DFOCC_AMPS);
    l2->read(psio_, PSIF_DFOCC_AMPS);
    GabA->contract442(3, 3, t2, l2, -1.0, 1.0);
    t2.reset();
    l2.reset();

    // G_ab = -1/2\sum_{m,n,f} t_mn^fa t_mn^fb
    t2 = SharedTensor2d(new Tensor2d("T2_1 <ij|ab>", naoccB, naoccB, navirB, navirB));
    l2 = SharedTensor2d(new Tensor2d("T2_1 <ij|ab>", naoccB, naoccB, navirB, navirB));
    t2->read(psio_, PSIF_DFOCC_AMPS);
    l2->read(psio_, PSIF_DFOCC_AMPS);
    GabB->contract442(4, 4, t2, l2, -0.5, 0.0);
    t2.reset();
    l2.reset();

    // G_ab += -\sum_{m,N,F} t_Mn^Fa t_Mn^Fb
    t2 = SharedTensor2d(new Tensor2d("T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    l2 = SharedTensor2d(new Tensor2d("T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    t2->read(psio_, PSIF_DFOCC_AMPS);
    l2->read(psio_, PSIF_DFOCC_AMPS);
    GabB->contract442(4, 4, t2, l2, -1.0, 1.0);
    t2.reset();
    l2.reset();
    //fprintf(outfile,"\tI am here.\n"); fflush(outfile);

    // Build G1c_oo and G1c_vv
    G1c_ooA->set_act_oo(nfrzc, naoccA, GijA);
    G1c_ooB->set_act_oo(nfrzc, naoccB, GijB);
    G1c_ooA->scale(-1.0);
    G1c_ooB->scale(-1.0);
    G1c_vvA->set_act_vv(GabA);
    G1c_vvB->set_act_vv(GabB);
    G1c_vvA->scale(-1.0);
    G1c_vvB->scale(-1.0);

    // Build G1c
    G1cA->set_oo(G1c_ooA);
    G1cA->set_vv(noccA, G1c_vvA);
    G1cB->set_oo(G1c_ooB);
    G1cB->set_vv(noccB, G1c_vvB);

    // Build G1
    G1A->copy(G1cA);
    G1B->copy(G1cB);
    for (int i = 0; i < noccA; i++) G1A->add(i, i, 1.0); 
    for (int i = 0; i < noccB; i++) G1B->add(i, i, 1.0); 

    // print
  if(print_ > 2) {
    G1A->print();
    G1B->print();
    double trace = G1A->trace();
    fprintf(outfile,"\t Alpha trace: %12.12f \n", trace);
    trace = G1B->trace();
    fprintf(outfile,"\t Beta trace: %12.12f \n", trace);
    fflush(outfile);
  }

}// else if (reference_ == "UNRESTRICTED")
    timer_off("opdm");
} // end omp2_opdm

}} // End Namespaces


