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

void DFOCC::effective_mograd()
{
    //outfile->Printf("\n effective_mograd is starting... \n");
    outfile->Printf("\tForming effective orbital gradient...\n");


    if (reference_ == "RESTRICTED") WvoA->form_vo(WorbA);
    else if (reference_ == "UNRESTRICTED") {
        WvoA->form_vo(WorbA);
        WvoB->form_vo(WorbB);
    }

    if (freeze_core_ == "TRUE") {
        z_vector_fc();
        fc_grad_terms();
    }
    //outfile->Printf("\n effective_mograd done. \n");
}// end effective_mograd

//=======================================================
//      Z-Vector: ACO-FC Block
//=======================================================
void DFOCC::z_vector_fc()
{

if (reference_ == "RESTRICTED") {
    // Build Zoo
    ZklA = SharedTensor2d(new Tensor2d("Zvector <I|FC>", naoccA, nfrzc));
    #pragma omp parallel for
    for (int k = 0; k < naoccA; k++) {
         for (int l = 0; l < nfrzc; l++) {
              double value = FockA->get(k + nfrzc, k + nfrzc) - FockA->get(l,l);
              ZklA->set(k, l, -WorbA->get(k + nfrzc, l) / (2.0*value));
         }
    }
    ZlkA = SharedTensor2d(new Tensor2d("Zvector <FC|I>", nfrzc, naoccA));
    ZlkA = ZklA->transpose();

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // Build Zoo
    // Alpha
    ZklA = SharedTensor2d(new Tensor2d("Zvector <I|FC>", naoccA, nfrzc));
    #pragma omp parallel for
    for (int k = 0; k < naoccA; k++) {
         for (int l = 0; l < nfrzc; l++) {
              double value = FockA->get(k + nfrzc, k + nfrzc) - FockA->get(l,l);
              ZklA->set(k, l, -WorbA->get(k + nfrzc, l) / (2.0*value));
         }
    }
    ZlkA = SharedTensor2d(new Tensor2d("Zvector <FC|I>", nfrzc, naoccA));
    ZlkA = ZklA->transpose();

    // Beta
    ZklB = SharedTensor2d(new Tensor2d("Zvector <i|FC>", naoccB, nfrzc));
    #pragma omp parallel for
    for (int k = 0; k < naoccB; k++) {
         for (int l = 0; l < nfrzc; l++) {
              double value = FockB->get(k + nfrzc, k + nfrzc) - FockB->get(l,l);
              ZklB->set(k, l, -WorbB->get(k + nfrzc, l) / (2.0*value));
         }
    }
    ZlkB = SharedTensor2d(new Tensor2d("Zvector <FC|i>", nfrzc, naoccB));
    ZlkB = ZklB->transpose();

}// end if (reference_ == "UNRESTRICTED")

}// end z_vector_fc

//=======================================================
//      FC GRAD TERMS
//=======================================================
void DFOCC::fc_grad_terms()
{

    SharedTensor2d K, L, IvoA, IvoB, G, Gsep, Z, Z2;
    timer_on("fc_grad_terms");
if (reference_ == "RESTRICTED") {
    //=========================
    // OPDM
    //=========================
    G1->add_aocc_fc(ZklA, 2.0, 1.0);
    G1->add_fc_aocc(ZlkA, 2.0, 1.0);

    //=========================
    // Seprable TPDM
    //=========================
    // Z_Q' = 4 \sum_{kl} b_{kl}^{Q} Z_kl
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|KL)", nQ_ref, naoccA, nfrzc));
    L->form_b_kl(K);
    K.reset();
    SharedTensor1d Zq = SharedTensor1d(new Tensor1d("DF_BASIS_SCF Zp_Q", nQ_ref));
    Zq->gemv(false, L, ZklA, 4.0, 0.0);
    L.reset();

    // GFM OO Block
    // F_ij += 2 \sum_{Q} b_ij^Q Z_Q'
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    GFoo->gemv(true, K, Zq, 2.0, 1.0);
    K.reset();

    // GFM VO Block
    // F_ai += 2 \sum_{Q} b_ai^Q Z_Q'
    // W_ai += 2 \sum_{Q} b_ai^Q Z_Q'
    IvoA = SharedTensor2d(new Tensor2d("MO-basis I <V|O>", nvirA, noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VO)", nQ_ref, nvirA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L->swap_3index_col(K);
    K.reset();
    IvoA->gemv(true, L, Zq, 2.0, 0.0);
    L.reset();

    // TPDM
    // G_kl^Q += 2 Z_kl J_Q
    // G_lk^Q += 2 Z_kl J_Q
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OO)", nQ_ref, noccA, noccA));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    #pragma omp parallel for
    for (int Q = 0; Q < nQ_ref; Q++) {
         for (int k = 0; k < naoccA; k++) {
              for (int l = 0; l < nfrzc; l++) {
                   int kl = l + ( (k + nfrzc) * noccA);
                   int lk = k + nfrzc + (l*noccA);
                   double value = 2.0 * ZklA->get(k,l) * Jc->get(Q);
                   Gsep->add(Q, kl, value);
                   Gsep->add(Q, lk, value);
              }
         }
    }

    //  G_ij^Q += 2 Z_Q \delta_{ij}
    #pragma omp parallel for
    for (int Q = 0; Q < nQ_ref; Q++) {
         for (int i = 0; i < noccA; i++) {
              int ii = oo_idxAA->get(i,i);
              Gsep->add(Q, ii, 2.0 * Zq->get(Q));
         }
    }
    Zq.reset();

    // Z_li^Q = 2 * \sum_{k} Z_lk b_ki^Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|KI)", nQ_ref, naoccA, noccA));
    L->form_b_ki(K);
    K.reset();
    Z = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|LI)", nQ_ref, nfrzc, noccA));
    Z->contract233(false, false, nfrzc, noccA, ZlkA, L, 2.0, 0.0);
    L.reset();

    // G_il^Q -= Z_li^Q
    // G_li^Q -= Z_li^Q
    #pragma omp parallel for
    for (int Q = 0; Q < nQ_ref; Q++) {
         for (int i = 0; i < noccA; i++) {
              for (int l = 0; l < nfrzc; l++) {
                   int il = l + (i*noccA);
                   int li = i + (l*noccA);
                   double value = Z->get(Q, li);
                   Gsep->subtract(Q, il, value);
                   Gsep->subtract(Q, li, value);
              }
         }
    }

    // GFM OO Block
    // F_ij -= \sum_{Q} \sum_{l} b_li^Q Z_lj^Q'
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|LI)", nQ_ref, nfrzc, noccA));
    L->form_b_li(K);
    K.reset();
    GFoo->contract(true, false, noccA, noccA, nQ_ref * nfrzc, L, Z, -1.0, 1.0);
    L.reset();

    // GFM VO Block
    // F_ai -= \sum_{Q} \sum_{l} b_la^Q Z_li^Q'
    // W_ai -= \sum_{Q} \sum_{l} b_la^Q Z_li^Q'
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|LA)", nQ_ref, nfrzc, nvirA));
    L->form_b_la(K);
    K.reset();
    IvoA->contract(true, false, nvirA, noccA, nQ_ref * nfrzc, L, Z, -1.0, 1.0);
    L.reset();
    Z.reset();

    // Z_ki^Q = 2 * \sum_{l} Z_kl b_li^Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|LI)", nQ_ref, nfrzc, noccA));
    L->form_b_li(K);
    K.reset();
    Z = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|KI)", nQ_ref, naoccA, noccA));
    Z->contract233(false, false, naoccA, noccA, ZklA, L, 2.0, 0.0);
    L.reset();

    // G_ki^Q -= Z_ki^Q
    // G_ik^Q -= Z_ki^Q
    #pragma omp parallel for
    for (int Q = 0; Q < nQ_ref; Q++) {
         for (int i = 0; i < noccA; i++) {
              for (int k = 0; k < naoccA; k++) {
                   int ik = k + nfrzc + (i*noccA);
                   int ki = i + ( (k+nfrzc) * noccA);
                   int ki2 = i + (k*noccA);
                   double value = Z->get(Q, ki2);
                   Gsep->subtract(Q, ik, value);
                   Gsep->subtract(Q, ki, value);
              }
         }
    }
    Gsep->write(psio_, PSIF_DFOCC_DENS);
    Gsep.reset();

    // GFM OO Block
    // F_ij -= \sum_{Q} \sum_{k} b_ki^Q Z_kj^Q'
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|KI)", nQ_ref, naoccA, noccA));
    L->form_b_ki(K);
    K.reset();
    GFoo->contract(true, false, noccA, noccA, nQ_ref * naoccA, L, Z, -1.0, 1.0);
    L.reset();

    // GFM VO Block
    // F_ai -= \sum_{Q} \sum_{k} b_ka^Q Z_ki^Q'
    // W_ai -= \sum_{Q} \sum_{k} b_ka^Q Z_ki^Q'
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|KA)", nQ_ref, naoccA, nvirA));
    L->form_b_ka(K);
    K.reset();
    IvoA->contract(true, false, nvirA, noccA, nQ_ref * naoccA, L, Z, -1.0, 1.0);
    GFvo->add(IvoA);
    WvoA->add(IvoA);
    IvoA.reset();
    L.reset();
    Z.reset();

    //=========================
    // GFM: AOCC-FC Terms
    //=========================
    // F_kl += 2.0 * z_kl f_kk
    // F_lk += 2.0 * z_kl f_ll
    #pragma omp parallel for
    for (int k = 0; k < naoccA; k++) {
         for (int l = 0; l < nfrzc; l++) {
              GFoo->add(k + nfrzc, l, 2.0 * ZklA->get(k,l) * FockA->get(k + nfrzc, k + nfrzc));
              GFoo->add(l, k + nfrzc, 2.0 * ZklA->get(k,l) * FockA->get(l, l));
	 }
    }

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    //=========================
    // OPDM
    //=========================
    G1A->add_aocc_fc(ZklA, 1.0, 1.0);
    G1A->add_fc_aocc(ZlkA, 1.0, 1.0);
    G1B->add_aocc_fc(ZklB, 1.0, 1.0);
    G1B->add_fc_aocc(ZlkB, 1.0, 1.0);

    //=========================
    // Seprable TPDM : Alpha
    //=========================
    // Z_Q' = 2 \sum_{KL} b_{KL}^{Q} Z_KL
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|KL)", nQ_ref, naoccA, nfrzc));
    L->form_b_kl(K);
    K.reset();
    SharedTensor1d Zq = SharedTensor1d(new Tensor1d("DF_BASIS_SCF Zp_Q", nQ_ref));
    Zq->gemv(false, L, ZklA, 2.0, 0.0);
    L.reset();

    // Z_Q' += 2 \sum_{kl} b_{kl}^{Q} Z_kl
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB, noccB));
    K->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|kl)", nQ_ref, naoccB, nfrzc));
    L->form_b_kl(K);
    K.reset();
    Zq->gemv(false, L, ZklB, 2.0, 1.0);
    L.reset();

    // GFM OO Block
    // F_IJ += \sum_{Q} b_IJ^Q Z_Q'
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    GFooA->gemv(true, K, Zq, 1.0, 1.0);
    K.reset();

    // GFM VO Block
    // F_AI += \sum_{Q} b_AI^Q Z_Q'
    // W_AI += \sum_{Q} b_AI^Q Z_Q'
    IvoA = SharedTensor2d(new Tensor2d("MO-basis I <V|O>", nvirA, noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VO)", nQ_ref, nvirA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L->swap_3index_col(K);
    K.reset();
    IvoA->gemv(true, L, Zq, 1.0, 0.0);
    L.reset();

    // TPDM
    // G_KL^Q += Z_KL J_Q
    // G_LK^Q += Z_KL J_Q
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OO)", nQ_ref, noccA, noccA));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    #pragma omp parallel for
    for (int Q = 0; Q < nQ_ref; Q++) {
         for (int k = 0; k < naoccA; k++) {
              for (int l = 0; l < nfrzc; l++) {
                   int kl = l + ( (k + nfrzc) * noccA);
                   int lk = k + nfrzc + (l*noccA);
                   double value = ZklA->get(k,l) * Jc->get(Q);
                   Gsep->add(Q, kl, value);
                   Gsep->add(Q, lk, value);
              }
         }
    }

    //  G_IJ^Q += Z_Q \delta_{IJ}
    #pragma omp parallel for
    for (int Q = 0; Q < nQ_ref; Q++) {
         for (int i = 0; i < noccA; i++) {
              int ii = oo_idxAA->get(i,i);
              Gsep->add(Q, ii, Zq->get(Q));
         }
    }
    //Zq.reset();

    // Z_LI^Q = \sum_{K} Z_LK b_KI^Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|KI)", nQ_ref, naoccA, noccA));
    L->form_b_ki(K);
    K.reset();
    Z = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|LI)", nQ_ref, nfrzc, noccA));
    Z->contract233(false, false, nfrzc, noccA, ZlkA, L, 1.0, 0.0);
    L.reset();

    // G_IL^Q -= Z_LI^Q
    // G_LI^Q -= Z_LI^Q
    #pragma omp parallel for
    for (int Q = 0; Q < nQ_ref; Q++) {
         for (int i = 0; i < noccA; i++) {
              for (int l = 0; l < nfrzc; l++) {
                   int il = l + (i*noccA);
                   int li = i + (l*noccA);
                   double value = Z->get(Q, li);
                   Gsep->subtract(Q, il, value);
                   Gsep->subtract(Q, li, value);
              }
         }
    }

    // GFM OO Block
    // F_IJ -= \sum_{Q} \sum_{L} b_LI^Q Z_LJ^Q'
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|LI)", nQ_ref, nfrzc, noccA));
    L->form_b_li(K);
    K.reset();
    GFooA->contract(true, false, noccA, noccA, nQ_ref * nfrzc, L, Z, -1.0, 1.0);
    L.reset();

    // GFM VO Block
    // F_AI -= \sum_{Q} \sum_{L} b_LA^Q Z_LI^Q'
    // W_AI -= \sum_{Q} \sum_{L} b_LA^Q Z_LI^Q'
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|LA)", nQ_ref, nfrzc, nvirA));
    L->form_b_la(K);
    K.reset();
    IvoA->contract(true, false, nvirA, noccA, nQ_ref * nfrzc, L, Z, -1.0, 1.0);
    L.reset();
    Z.reset();

    // Z_KI^Q = \sum_{L} Z_KL b_LI^Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|LI)", nQ_ref, nfrzc, noccA));
    L->form_b_li(K);
    K.reset();
    Z = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|KI)", nQ_ref, naoccA, noccA));
    Z->contract233(false, false, naoccA, noccA, ZklA, L, 1.0, 0.0);
    L.reset();

    // G_KI^Q -= Z_KI^Q
    // G_IK^Q -= Z_KI^Q
    #pragma omp parallel for
    for (int Q = 0; Q < nQ_ref; Q++) {
         for (int i = 0; i < noccA; i++) {
              for (int k = 0; k < naoccA; k++) {
                   int ik = k + nfrzc + (i*noccA);
                   int ki = i + ( (k+nfrzc) * noccA);
                   int ki2 = i + (k*noccA);
                   double value = Z->get(Q, ki2);
                   Gsep->subtract(Q, ik, value);
                   Gsep->subtract(Q, ki, value);
              }
         }
    }
    Gsep->write(psio_, PSIF_DFOCC_DENS);
    Gsep.reset();

    // GFM OO Block
    // F_IJ -= \sum_{Q} \sum_{K} b_KI^Q Z_KJ^Q'
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|KI)", nQ_ref, naoccA, noccA));
    L->form_b_ki(K);
    K.reset();
    GFooA->contract(true, false, noccA, noccA, nQ_ref * naoccA, L, Z, -1.0, 1.0);
    L.reset();

    // GFM VO Block
    // F_AI -= \sum_{Q} \sum_{K} b_KA^Q Z_KI^Q'
    // W_AI -= \sum_{Q} \sum_{k} b_KA^Q Z_KI^Q'
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|KA)", nQ_ref, naoccA, nvirA));
    L->form_b_ka(K);
    K.reset();
    IvoA->contract(true, false, nvirA, noccA, nQ_ref * naoccA, L, Z, -1.0, 1.0);
    L.reset();
    Z.reset();
    GFvoA->add(IvoA);
    WvoA->add(2.0, IvoA);
    IvoA.reset();

    //=========================
    // GFM: AOCC-FC Terms
    //=========================
    // F_kl += z_kl f_kk
    // F_lk += z_kl f_ll
    #pragma omp parallel for
    for (int k = 0; k < naoccA; k++) {
         for (int l = 0; l < nfrzc; l++) {
              GFooA->add(k + nfrzc, l, ZklA->get(k,l) * FockA->get(k + nfrzc, k + nfrzc));
              GFooA->add(l, k + nfrzc, ZklA->get(k,l) * FockA->get(l, l));
	 }
    }

    //=========================
    // Seprable TPDM : Beta
    //=========================
    // GFM oo Block
    // F_ij += \sum_{Q} b_ij^Q Z_Q'
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB, noccB));
    K->read(psio_, PSIF_DFOCC_INTS);
    GFooB->gemv(true, K, Zq, 1.0, 1.0);
    K.reset();

    // GFM vo Block
    // F_ai += \sum_{Q} b_ai^Q Z_Q'
    // W_ai += \sum_{Q} b_ai^Q Z_Q'
    IvoB = SharedTensor2d(new Tensor2d("MO-basis I <v|o>", nvirB, noccB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB, nvirB));
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vo)", nQ_ref, nvirB, noccB));
    K->read(psio_, PSIF_DFOCC_INTS);
    L->swap_3index_col(K);
    K.reset();
    IvoB->gemv(true, L, Zq, 1.0, 0.0);
    L.reset();

    // TPDM
    // G_kl^Q += Z_kl J_Q
    // G_lk^Q += Z_kl J_Q
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|oo)", nQ_ref, noccB, noccB));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    #pragma omp parallel for
    for (int Q = 0; Q < nQ_ref; Q++) {
         for (int k = 0; k < naoccB; k++) {
              for (int l = 0; l < nfrzc; l++) {
                   int kl = l + ( (k + nfrzc) * noccB);
                   int lk = k + nfrzc + (l*noccB);
                   double value = ZklB->get(k,l) * Jc->get(Q);
                   Gsep->add(Q, kl, value);
                   Gsep->add(Q, lk, value);
              }
         }
    }

    //  G_ij^Q += Z_Q \delta_{ij}
    #pragma omp parallel for
    for (int Q = 0; Q < nQ_ref; Q++) {
         for (int i = 0; i < noccB; i++) {
              int ii = oo_idxBB->get(i,i);
              Gsep->add(Q, ii, Zq->get(Q));
         }
    }
    Zq.reset();

    // Z_li^Q = \sum_{k} Z_lk b_ki^Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB, noccB));
    K->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ki)", nQ_ref, naoccB, noccB));
    L->form_b_ki(K);
    K.reset();
    Z = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|li)", nQ_ref, nfrzc, noccB));
    Z->contract233(false, false, nfrzc, noccB, ZlkB, L, 1.0, 0.0);
    L.reset();

    // G_il^Q -= Z_li^Q
    // G_li^Q -= Z_li^Q
    #pragma omp parallel for
    for (int Q = 0; Q < nQ_ref; Q++) {
         for (int i = 0; i < noccB; i++) {
              for (int l = 0; l < nfrzc; l++) {
                   int il = l + (i*noccB);
                   int li = i + (l*noccB);
                   double value = Z->get(Q, li);
                   Gsep->subtract(Q, il, value);
                   Gsep->subtract(Q, li, value);
              }
         }
    }

    // GFM oo Block
    // F_ij -= \sum_{Q} \sum_{l} b_li^Q Z_lj^Q'
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB, noccB));
    K->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|li)", nQ_ref, nfrzc, noccB));
    L->form_b_li(K);
    K.reset();
    GFooB->contract(true, false, noccB, noccB, nQ_ref * nfrzc, L, Z, -1.0, 1.0);
    L.reset();

    // GFM vo Block
    // F_ai -= \sum_{Q} \sum_{l} b_la^Q Z_li^Q'
    // W_ai -= \sum_{Q} \sum_{l} b_la^Q Z_li^Q'
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB, nvirB));
    K->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|la)", nQ_ref, nfrzc, nvirB));
    L->form_b_la(K);
    K.reset();
    IvoB->contract(true, false, nvirB, noccB, nQ_ref * nfrzc, L, Z, -1.0, 1.0);
    L.reset();
    Z.reset();

    // Z_ki^Q = \sum_{l} Z_kl b_li^Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB, noccB));
    K->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|li)", nQ_ref, nfrzc, noccB));
    L->form_b_li(K);
    K.reset();
    Z = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|ki)", nQ_ref, naoccB, noccB));
    Z->contract233(false, false, naoccB, noccB, ZklB, L, 1.0, 0.0);
    L.reset();

    // G_ki^Q -= Z_ki^Q
    // G_ik^Q -= Z_ki^Q
    #pragma omp parallel for
    for (int Q = 0; Q < nQ_ref; Q++) {
         for (int i = 0; i < noccB; i++) {
              for (int k = 0; k < naoccB; k++) {
                   int ik = k + nfrzc + (i*noccB);
                   int ki = i + ( (k+nfrzc) * noccB);
                   int ki2 = i + (k*noccB);
                   double value = Z->get(Q, ki2);
                   Gsep->subtract(Q, ik, value);
                   Gsep->subtract(Q, ki, value);
              }
         }
    }
    Gsep->write(psio_, PSIF_DFOCC_DENS);
    Gsep.reset();

    // GFM oo Block
    // F_ij -= \sum_{Q} \sum_{k} b_ki^Q Z_kj^Q'
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB, noccB));
    K->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ki)", nQ_ref, naoccB, noccB));
    L->form_b_ki(K);
    K.reset();
    GFooB->contract(true, false, noccB, noccB, nQ_ref * naoccB, L, Z, -1.0, 1.0);
    L.reset();

    // GFM VO Block
    // F_ai -= \sum_{Q} \sum_{k} b_ka^Q Z_ki^Q'
    // W_ai -= \sum_{Q} \sum_{k} b_ka^Q Z_ki^Q'
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB, nvirB));
    K->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ka)", nQ_ref, naoccB, nvirB));
    L->form_b_ka(K);
    K.reset();
    IvoB->contract(true, false, nvirB, noccB, nQ_ref * naoccB, L, Z, -1.0, 1.0);
    L.reset();
    Z.reset();
    GFvoB->add(IvoB);
    WvoB->add(2.0, IvoB);
    IvoB.reset();

    //=========================
    // GFM: AOCC-FC Terms
    //=========================
    // F_kl += z_kl f_kk
    // F_lk += z_kl f_ll
    #pragma omp parallel for
    for (int k = 0; k < naoccB; k++) {
         for (int l = 0; l < nfrzc; l++) {
              GFooB->add(k + nfrzc, l, ZklB->get(k,l) * FockB->get(k + nfrzc, k + nfrzc));
              GFooB->add(l, k + nfrzc, ZklB->get(k,l) * FockB->get(l, l));
         }
    }

}// end if (reference_ == "UNRESTRICTED")
    timer_off("fc_grad_terms");

}// end fc_grad_terms


}} // End Namespaces
