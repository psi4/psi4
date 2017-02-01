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

/** Standard library includes */
#include "psi4/libqt/qt.h"
#include "defines.h"
#include "dfocc.h"


using namespace psi;
using namespace std;


namespace psi{ namespace dfoccwave{

void DFOCC::omp2_opdm()
{

    SharedTensor2d T, U;
    timer_on("opdm");
if (reference_ == "RESTRICTED") {
    // Tensors
    T = SharedTensor2d(new Tensor2d("T2_1 (ia|jb)", naoccA, navirA, naoccA, navirA));
    U = SharedTensor2d(new Tensor2d("U2_1 (ia|jb)", naoccA, navirA, naoccA, navirA));
    if (orb_opt_ == "FALSE" && mp2_amp_type_ == "DIRECT") {
        u2_rmp2_direct(T, U);
    }
    else {
        T->read_symm(psio_, PSIF_DFOCC_AMPS);
        U->read_symm(psio_, PSIF_DFOCC_AMPS);
    }

    // G_ij = \sum_{m,e,f} T'(ie,mf) U'(je,mf)
    GijA->contract442(1, 1, T, U, 1.0, 0.0);

    // G_ab = \sum_{m,n,e} U'(ma,ne) T'(mb,ne)
    GabA->contract442(2, 2, U, T, -1.0, 0.0);
    T.reset();
    U.reset();

    // Build G1c_oo and G1c_vv
    G1c_oo->set_act_oo(nfrzc, naoccA, GijA);
    G1c_oo->scale(-2.0);
    G1c_vv->set_act_vv(GabA);
    G1c_vv->scale(-2.0);

    // Build G1c
    G1c->set_oo(G1c_oo);
    G1c->set_vv(noccA, G1c_vv);

    // Build G1
    G1->copy(G1c);
    for (int i = 0; i < noccA; i++) G1->add(i, i, 2.0);

  if(print_ > 2) {
    G1->print();
    double trace = G1->trace();
    outfile->Printf("\t trace: %12.12f \n", trace);

  }

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // tensors
    SharedTensor2d t2, l2, l1A, l1B;

    // G_IJ = 1/2 \sum_{M,E,F} t_IM^EF t_JM^EF
    t2 = SharedTensor2d(new Tensor2d("T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    l2 = SharedTensor2d(new Tensor2d("T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    if (orb_opt_ == "FALSE" && mp2_amp_type_ == "DIRECT") t2AA_ump2_direct(t2);
    else t2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    l2->copy(t2);
    GijA->contract442(1, 1, t2, l2, 0.5, 0.0);
    // G_AB = -1/2\sum_{M,N,F} t_MN^FA t_MN^FB
    GabA->contract442(4, 4, t2, l2, -0.5, 0.0);
    t2.reset();
    l2.reset();

    // G_IJ = \sum_{m,E,f} t_Im^Ef t_Jm^Ef
    t2 = SharedTensor2d(new Tensor2d("T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    l2 = SharedTensor2d(new Tensor2d("T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    if (orb_opt_ == "FALSE" && mp2_amp_type_ == "DIRECT") t2AB_ump2_direct(t2);
    else t2->read(psio_, PSIF_DFOCC_AMPS);
    l2->copy(t2);
    GijA->contract442(1, 1, t2, l2, 1.0, 1.0);
    // G_ij = \sum_{M,e,F} t_Mi^Fe t_Mj^Fe
    GijB->contract442(2, 2, t2, l2, 1.0, 0.0);
    // G_AB += -\sum_{M,n,f} t_Mn^Af t_Mn^Bf
    GabA->contract442(3, 3, t2, l2, -1.0, 1.0);
    // G_ab += -\sum_{m,N,F} t_Mn^Fa t_Mn^Fb
    GabB->contract442(4, 4, t2, l2, -1.0, 0.0);
    t2.reset();
    l2.reset();

    // G_ij = 1/2 \sum_{m,e,f} t_im^ef t_jm^ef
    t2 = SharedTensor2d(new Tensor2d("T2_1 <ij|ab>", naoccB, naoccB, navirB, navirB));
    l2 = SharedTensor2d(new Tensor2d("T2_1 <ij|ab>", naoccB, naoccB, navirB, navirB));
    if (orb_opt_ == "FALSE" && mp2_amp_type_ == "DIRECT") t2BB_ump2_direct(t2);
    else t2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    l2->copy(t2);
    GijB->contract442(1, 1, t2, l2, 0.5, 1.0);
    // G_ab = -1/2\sum_{m,n,f} t_mn^fa t_mn^fb
    GabB->contract442(4, 4, t2, l2, -0.5, 1.0);
    t2.reset();
    l2.reset();

    //outfile->Printf("\tI am here.\n");

   if (reference == "ROHF" && orb_opt_ == "FALSE") {
       // G_ia = t_i^a
       GiaA->copy(t1A);
       GiaB->copy(t1B);
       GaiA = GiaA->transpose();
       GaiB = GiaB->transpose();

       // G_ij -= \sum_{e} T_i^e L_e^i
       GijA->gemm(false, true, t1A, t1A, -1.0, 1.0);
       GijB->gemm(false, true, t1B, t1B, -1.0, 1.0);

       // G_ab += \sum_{m} L_a^m T_m^b
       GabA->gemm(true, false, t1A, t1A, 1.0, 1.0);
       GabB->gemm(true, false, t1B, t1B, 1.0, 1.0);
   }

    // Build G1c_oo and G1c_vv
    G1c_ooA->set_act_oo(nfrzc, naoccA, GijA);
    G1c_ooB->set_act_oo(nfrzc, naoccB, GijB);
    G1c_ooA->scale(-1.0);
    G1c_ooB->scale(-1.0);
    G1c_vvA->set_act_vv(GabA);
    G1c_vvB->set_act_vv(GabB);
    G1c_vvA->scale(-1.0);
    G1c_vvB->scale(-1.0);

    if (reference == "ROHF" && orb_opt_ == "FALSE") {
        G1c_ovA->set_act_ov(nfrzc, GiaA);
        G1c_ovB->set_act_ov(nfrzc, GiaB);
        G1c_voA->set_act_vo(nfrzc, GaiA);
        G1c_voB->set_act_vo(nfrzc, GaiB);
    }

    // Build G1c
    G1cA->set_oo(G1c_ooA);
    G1cA->set_vv(noccA, G1c_vvA);
    G1cB->set_oo(G1c_ooB);
    G1cB->set_vv(noccB, G1c_vvB);

    if (reference == "ROHF" && orb_opt_ == "FALSE") {
        G1cA->set_ov(G1c_ovA);
        G1cB->set_ov(G1c_ovB);
        G1cA->set_vo(G1c_voA);
        G1cB->set_vo(G1c_voB);
    }

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
    outfile->Printf("\t Alpha trace: %12.12f \n", trace);
    trace = G1B->trace();
    outfile->Printf("\t Beta trace: %12.12f \n", trace);

  }

}// else if (reference_ == "UNRESTRICTED")
    timer_off("opdm");
} // end omp2_opdm

}} // End Namespaces
