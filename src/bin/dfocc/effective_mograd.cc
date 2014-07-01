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

void DFOCC::effective_mograd()
{ 
//fprintf(outfile,"\n effective_mograd is starting... \n"); fflush(outfile);
    fprintf(outfile,"\tForming effective orbital gradient...\n");
    fflush(outfile);

    SharedTensor2d K, L;

if (reference_ == "RESTRICTED") {
    // Build Wvo 
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              WvoA->set(a, i, WorbA->get(a + noccA, i));       
         }
    }
 

}// end if (reference_ == "RESTRICTED") 

else if (reference_ == "UNRESTRICTED") {
    // alpha
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              WvoA->set(a, i, WorbA->get(a + noccA, i));       
         }
    }

    // beta
    for (int a = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++) {
              WvoB->set(a, i, WorbB->get(a + noccB, i));       
         }
    }

}// end if (reference_ == "UNRESTRICTED") 

    if (freeze_core_ == "TRUE") {
        z_vector_fc();
        fc_grad_terms();
    }

 //fprintf(outfile,"\n effective_mograd done. \n"); fflush(outfile);
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

    SharedTensor2d K, L, G, Gsep, Z, Z2;
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
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VO)", nQ_ref, nvirA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L->swap_3index_col(K);
    K.reset();
    GFvo->gemv(true, L, Zq, 2.0, 1.0);
    // W_ai += 2 \sum_{Q} b_ai^Q Z_Q'
    WvoA->gemv(true, L, Zq, 2.0, 1.0);
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
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|LA)", nQ_ref, nfrzc, nvirA));
    L->form_b_la(K);
    K.reset();
    GFvo->contract(true, false, nvirA, noccA, nQ_ref * nfrzc, L, Z, -1.0, 1.0);
    // W_ai -= 2\sum_{Q} \sum_{l} b_la^Q Z_li^Q'
    //WvoA->contract(true, false, nvirA, noccA, nQ_ref * nfrzc, L, Z, -2.0, 1.0);
    // W_ai -= \sum_{Q} \sum_{l} b_la^Q Z_li^Q'
    WvoA->contract(true, false, nvirA, noccA, nQ_ref * nfrzc, L, Z, -1.0, 1.0);
    L.reset();
    Z.reset();

    // Z_ki^Q = 2 * \sum_{l} Z_kl b_li^Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|IL)", nQ_ref, noccA, nfrzc));
    L->form_b_il(K);
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
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|KA)", nQ_ref, naoccA, nvirA));
    L->form_b_ka(K);
    K.reset();
    GFvo->contract(true, false, nvirA, noccA, nQ_ref * naoccA, L, Z, -1.0, 1.0);
    // W_ai -= \sum_{Q} \sum_{k} b_ka^Q Z_ki^Q'
    WvoA->contract(true, false, nvirA, noccA, nQ_ref * naoccA, L, Z, -1.0, 1.0);
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
}// end if (reference_ == "UNRESTRICTED") 
    timer_off("fc_grad_terms");

}// end fc_grad_terms


}} // End Namespaces


