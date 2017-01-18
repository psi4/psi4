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

void DFOCC::prepare4grad()
{
    if (wfn_type_ == "DF-OMP2") separable_tpdm();
    else sep_tpdm_cc();
    //if (wfn_type_ == "DF-CCSD" || wfn_type_ == "CD-CCSD" || wfn_type_ == "DF-CCD" || wfn_type_ == "CD-CCD") ccl_energy2();
    idp2();

    outfile->Printf("\tComputing the generalized Fock matrix (GFM)...\n");
    if (wfn_type_ == "DF-OMP2") {
        gfock_vo();
        gfock_ov();
        gfock_oo();
        gfock_vv();
    }
    else {
        gfock_cc_vo();
        gfock_cc_ov();
        gfock_cc_oo();
        gfock_cc_vv();
    }

    outfile->Printf("\tComputing the orbital gradient...\n");
    mograd();
    effective_mograd();
    timer_on("Z-vector");
    z_vector_pcg();
    timer_off("Z-vector");
    effective_pdm_gfm();

}// end prepare4grad

//=======================================================
//          Effective PDMs
//=======================================================
void DFOCC::effective_pdms()
{
    outfile->Printf("\tComputing the relaxed response density matrices...\n");


    SharedTensor2d K, L, G, Gsep, Z, Z2;
    timer_on("effective_pdms");
if (reference_ == "RESTRICTED") {
    //=========================
    // OPDM
    //=========================
    G1->add_vo(ZvoA, 2.0, 1.0);
    G1->add_ov(ZovA, 2.0, 1.0);

    //=========================
    // Seprable TPDM
    //=========================
    // Intermediates
    // Z_Q = 4 \sum_{ai} b_{ai}^{Q} Z_ai
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VO)", nQ_ref, nvirA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L->swap_3index_col(K);
    K.reset();
    SharedTensor1d Zq = SharedTensor1d(new Tensor1d("DF_BASIS_SCF Z_Q", nQ_ref));
    Zq->gemv(false, L, ZvoA, 4.0, 0.0);

    // Z_ij^Q = 4 * \sum_{a} Z_ia b_aj^Q
    Z = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|OO)", nQ_ref, noccA, noccA));
    Z->contract233(false, false, noccA, noccA, ZovA, L, 4.0, 0.0);
    L.reset();
    Z->symmetrize3(Z);

    // OO Block
    //  G_ij^Q += 2 Z_Q \delta_{ij} - 1/2 (Z_ij^Q + Z_ji^Q)
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OO)", nQ_ref, noccA, noccA));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    Gsep->subtract(Z);
    Z.reset();
    #pragma omp parallel for
    for (int Q = 0; Q < nQ_ref; Q++) {
         for (int i = 0; i < noccA; i++) {
              int ii = oo_idxAA->get(i,i);
              Gsep->add(Q, ii, 2.0 * Zq->get(Q));
         }
    }
    Zq.reset();
    Gsep->write(psio_, PSIF_DFOCC_DENS);
    Gsep.reset();

    // OV and VO Blocks
    // Z_ai^Q = 2 * \sum_{m} Z_am b_mi^Q
    Z = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|VO)", nQ_ref, nvirA, noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    Z->contract233(false, false, nvirA, noccA, ZvoA, K, 2.0, 0.0);
    K.reset();

    // G_ai^Q += 2 * J_Q Zai - Z_ai^Q
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VO)", nQ_ref, nvirA, noccA));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    Gsep->subtract(Z);
    Gsep->dirprd123(Jc, ZvoA, 2.0, 1.0);
    Gsep->write(psio_, PSIF_DFOCC_DENS);
    Gsep.reset();

    // G_ia^Q += 2 * J_Q Zia - Z_ia^Q : where Z_ia^Q = Z_ai^Q
    Z2 = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|OV)", nQ_ref, noccA, nvirA));
    Z2->swap_3index_col(Z);
    Z.reset();
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OV)", nQ_ref, noccA, nvirA));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    Gsep->subtract(Z2);
    Z2.reset();
    Gsep->dirprd123(Jc, ZovA, 2.0, 1.0);
    Gsep->write(psio_, PSIF_DFOCC_DENS);
    Gsep.reset();

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    //=========================
    // OPDM
    //=========================
    // Alpha
    G1A->add_vo(ZvoA, 1.0, 1.0);
    G1A->add_ov(ZovA, 1.0, 1.0);

    // Beta
    G1B->add_vo(ZvoB, 1.0, 1.0);
    G1B->add_ov(ZovB, 1.0, 1.0);

    //=========================
    // Seprable TPDM
    //=========================
    // Intermediates
    // Z_Q = 2 \sum_{AI} b_{AI}^{Q} Z_AI
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VO)", nQ_ref, nvirA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L->swap_3index_col(K);
    K.reset();
    SharedTensor1d Zq = SharedTensor1d(new Tensor1d("DF_BASIS_SCF Z_Q", nQ_ref));
    Zq->gemv(false, L, ZvoA, 2.0, 0.0);

    // Z_IJ^Q = 2 * \sum_{A} Z_IA b_AJ^Q
    Z = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|OO)", nQ_ref, noccA, noccA));
    Z->contract233(false, false, noccA, noccA, ZovA, L, 2.0, 0.0);
    L.reset();
    Z->symmetrize3(Z);
    Z->write(psio_, PSIF_DFOCC_DENS);
    Z.reset();

    // Z_Q += 2 \sum_{ai} b_{ai}^{Q} Z_ai
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB, nvirB));
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vo)", nQ_ref, nvirB, noccB));
    K->read(psio_, PSIF_DFOCC_INTS);
    L->swap_3index_col(K);
    K.reset();
    Zq->gemv(false, L, ZvoB, 2.0, 1.0);

    // Z_ij^Q = 2 * \sum_{a} Z_ia b_aj^Q
    Z = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|oo)", nQ_ref, noccB, noccB));
    Z->contract233(false, false, noccB, noccB, ZovB, L, 2.0, 0.0);
    L.reset();
    Z->symmetrize3(Z);
    Z->write(psio_, PSIF_DFOCC_DENS);
    Z.reset();

    // OO Block
    // G_IJ^Q += Z_Q \delta_{IJ} - 1/2 (Z_IJ^Q + Z_JI^Q)
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OO)", nQ_ref, noccA, noccA));
    Z = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|OO)", nQ_ref, noccA, noccA));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    Z->read(psio_, PSIF_DFOCC_DENS);
    Gsep->subtract(Z);
    Z.reset();
    #pragma omp parallel for
    for (int Q = 0; Q < nQ_ref; Q++) {
         for (int i = 0; i < noccA; i++) {
              int ii = oo_idxAA->get(i,i);
              Gsep->add(Q, ii, Zq->get(Q));
         }
    }
    Gsep->write(psio_, PSIF_DFOCC_DENS);
    Gsep.reset();

    // oo Block
    //  G_ij^Q += Z_Q \delta_{ij} - 1/2 (Z_ij^Q + Z_ji^Q)
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|oo)", nQ_ref, noccB, noccB));
    Z = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|oo)", nQ_ref, noccB, noccB));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    Z->read(psio_, PSIF_DFOCC_DENS);
    Gsep->subtract(Z);
    Z.reset();
    #pragma omp parallel for
    for (int Q = 0; Q < nQ_ref; Q++) {
         for (int i = 0; i < noccB; i++) {
              int ii = oo_idxBB->get(i,i);
              Gsep->add(Q, ii, Zq->get(Q));
         }
    }
    Zq.reset();
    Gsep->write(psio_, PSIF_DFOCC_DENS);
    Gsep.reset();

    // OV and VO Blocks
    // Z_AI^Q = \sum_{M} Z_AM b_MI^Q
    Z = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|VO)", nQ_ref, nvirA, noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    Z->contract233(false, false, nvirA, noccA, ZvoA, K, 1.0, 0.0);
    K.reset();

    // G_AI^Q += J_Q Z_AI - Z_AI^Q
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VO)", nQ_ref, nvirA, noccA));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    Gsep->subtract(Z);
    Gsep->dirprd123(Jc, ZvoA, 1.0, 1.0);
    Gsep->write(psio_, PSIF_DFOCC_DENS);
    Gsep.reset();

    // G_IA^Q += J_Q Z_IA - Z_IA^Q : where Z_IA^Q = Z_AI^Q
    Z2 = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|OV)", nQ_ref, noccA, nvirA));
    Z2->swap_3index_col(Z);
    Z.reset();
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OV)", nQ_ref, noccA, nvirA));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    Gsep->subtract(Z2);
    Z2.reset();
    Gsep->dirprd123(Jc, ZovA, 1.0, 1.0);
    Gsep->write(psio_, PSIF_DFOCC_DENS);
    Gsep.reset();

    // ov and vo Blocks
    // Z_ai^Q = \sum_{m} Z_am b_mi^Q
    Z = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|vo)", nQ_ref, nvirB, noccB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB, noccB));
    K->read(psio_, PSIF_DFOCC_INTS);
    Z->contract233(false, false, nvirB, noccB, ZvoB, K, 1.0, 0.0);
    K.reset();

    // G_ai^Q += J_Q Zai - Z_ai^Q
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|vo)", nQ_ref, nvirB, noccB));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    Gsep->subtract(Z);
    Gsep->dirprd123(Jc, ZvoB, 1.0, 1.0);
    Gsep->write(psio_, PSIF_DFOCC_DENS);
    Gsep.reset();

    // G_ia^Q += J_Q Zia - Z_ia^Q : where Z_ia^Q = Z_ai^Q
    Z2 = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|ov)", nQ_ref, noccB, nvirB));
    Z2->swap_3index_col(Z);
    Z.reset();
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|ov)", nQ_ref, noccB, nvirB));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    Gsep->subtract(Z2);
    Z2.reset();
    Gsep->dirprd123(Jc, ZovB, 1.0, 1.0);
    Gsep->write(psio_, PSIF_DFOCC_DENS);
    Gsep.reset();

}// else if (reference_ == "UNRESTRICTED")
    timer_off("effective_pdms");
} // end effective_pdms

//=======================================================
//          Effective GFM
//=======================================================
void DFOCC::effective_gfm()
{
    outfile->Printf("\tComputing the relaxed GFM...\n");


    SharedTensor2d K, L;
    timer_on("effective_gfm");
if (reference_ == "RESTRICTED") {
    // OO Block
    // F_ij += 8 * \sum_{em} (ij|em) z_em
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|OV)", noccA, noccA, noccA, nvirA));
    tei_ooov_chem_ref_directAA(K);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|VO)", noccA, noccA, nvirA, noccA));
    L->sort(1243, K, 1.0, 0.0);
    GFoo->gemv(false, L, ZvoA, 8.0, 1.0);
    L.reset();

    // F_ij += -2 * \sum_{em} <ij|em> z_em
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OO|VO>", noccA, noccA, nvirA, noccA));
    // <ij|em> = (ie|jm) = (jm|ie)
    L->sort(3142, K, 1.0, 0.0);
    K.reset();
    GFoo->gemv(false, L, ZvoA, -2.0, 1.0);
    L.reset();

    // F_ij += -2 * \sum_{em} <ij|me> z_me
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OO|OV>", noccA, noccA, noccA, nvirA));
    tei_ooov_phys_ref_directAA(K);
    GFoo->gemv(false, K, ZovA, -2.0, 1.0);
    K.reset();

    // Set global GF
    GF->set_oo(GFoo);

    // VO Block
    // F_ai += 8 * \sum_{em} (ai|em) z_em
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|OV)", noccA, nvirA, noccA, nvirA));
    tei_ovov_chem_ref_directAA(K);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (VO|VO)", nvirA, noccA, nvirA, noccA));
    L->sort(2143, K, 1.0, 0.0);
    GFvo->gemv(false, L, ZvoA, 8.0, 1.0);
    L.reset();

    // F_ai += -2 * \sum_{em} <ai|me> z_me: (ie|ma) => <ai|me>
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <VO|OV>", nvirA, noccA, noccA, nvirA));
    L->sort(4132, K, 1.0, 0.0);
    K.reset();
    GFvo->gemv(false, L, ZovA, -2.0, 1.0);
    L.reset();

    // F_ai += -2 * \sum_{em} <ai|em> z_em
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|VV)", noccA, noccA, nvirA, nvirA));
    tei_oovv_chem_ref_directAA(K);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <VO|VO>", nvirA, noccA, nvirA, noccA));
    L->sort(3142, K, 1.0, 0.0);
    K.reset();
    GFvo->gemv(false, L, ZvoA, -2.0, 1.0);
    L.reset();

    // Set global GF
    GF->set_vo(GFvo);

    // OV Block and remaining part of the VO block
    // F_ia += 2.0 * z_ai f_ii
    // F_ai += 2.0 * z_ai f_aa
    #pragma omp parallel for
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              GF->add(i, a + noccA, 2.0 * ZvoA->get(a, i) * FockA->get(i,i));
              GF->add(a + noccA, i, 2.0 * ZvoA->get(a, i) * FockA->get(a + noccA, a + noccA));
	 }
    }

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // OO Block
    // F_IJ += 2 * \sum_{EM} (IJ|EM) Z_EM
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|OV)", noccA, noccA, noccA, nvirA));
    tei_ooov_chem_ref_directAA(K);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|VO)", noccA, noccA, nvirA, noccA));
    L->sort(1243, K, 1.0, 0.0);
    GFooA->gemv(false, L, ZvoA, 2.0, 1.0);
    L.reset();

    // F_IJ -=  \sum_{EM} <IJ|EM> Z_EM. Note: <IJ|EM> = (IE|JM) = (JM|IE)
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OO|VO>", noccA, noccA, nvirA, noccA));
    L->sort(3142, K, 1.0, 0.0);
    K.reset();
    GFooA->gemv(false, L, ZvoA, -1.0, 1.0);
    L.reset();

    // F_IJ -= \sum_{EM} <IJ|ME> Z_ME
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OO|OV>", noccA, noccA, noccA, nvirA));
    tei_ooov_phys_ref_directAA(K);
    GFooA->gemv(false, K, ZovA, -1.0, 1.0);
    K.reset();

    // F_IJ += 2 * \sum_{em} (IJ|em) Z_em =  2 * \sum_{em} (IJ|me) Z_me
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|ov)", noccA, noccA, noccB, nvirB));
    tei_ooov_chem_ref_directAB(K);
    GFooA->gemv(false, K, ZovB, 2.0, 1.0);
    K.reset();

    // Set global GF
    GFA->set_oo(GFooA);

    // oo Block
    // F_ij += 2 * \sum_{em} (ij|em) Z_em
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (oo|oV)", noccB, noccB, noccB, nvirB));
    tei_ooov_chem_ref_directBB(K);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (oo|vo)", noccB, noccB, nvirB, noccB));
    L->sort(1243, K, 1.0, 0.0);
    GFooB->gemv(false, L, ZvoB, 2.0, 1.0);
    L.reset();

    // F_ij -=  \sum_{em} <ij|em> Z_em. Note: <ij|em> = (ie|jm) = (jm|ie)
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <oo|vo>", noccB, noccB, nvirB, noccB));
    L->sort(3142, K, 1.0, 0.0);
    K.reset();
    GFooB->gemv(false, L, ZvoB, -1.0, 1.0);
    L.reset();

    // F_ij -= \sum_{em} <ij|me> Z_me
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <oo|ov>", noccB, noccB, noccB, nvirB));
    tei_ooov_phys_ref_directBB(K);
    GFooB->gemv(false, K, ZovB, -1.0, 1.0);
    K.reset();

    // F_ij += 2 * \sum_{EM} (ij|EM) Z_EM = 2 * \sum_{ME} (ij|ME) Z_ME
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|oo)", noccA, nvirA, noccB, noccB));
    tei_ovoo_chem_ref_directAB(K);
    GFooB->gemv(true, K, ZovA, 2.0, 1.0);
    K.reset();

    // Set global GF
    GFB->set_oo(GFooB);

    // VO Block
    // F_AI += 2 * \sum_{EM} (AI|EM) Z_EM
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|OV)", noccA, nvirA, noccA, nvirA));
    tei_ovov_chem_ref_directAA(K);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (VO|VO)", nvirA, noccA, nvirA, noccA));
    L->sort(2143, K, 1.0, 0.0);
    GFvoA->gemv(false, L, ZvoA, 2.0, 1.0);
    L.reset();

    // F_AI -= \sum_{EM} <AI|ME> Z_ME: (IE|MA) => <AI|ME>
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <VO|OV>", nvirA, noccA, noccA, nvirA));
    L->sort(4132, K, 1.0, 0.0);
    K.reset();
    GFvoA->gemv(false, L, ZovA, -1.0, 1.0);
    L.reset();

    // F_AI -= \sum_{EM} <AI|EM> Z_EM
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|VV)", noccA, noccA, nvirA, nvirA));
    tei_oovv_chem_ref_directAA(K);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <VO|VO>", nvirA, noccA, nvirA, noccA));
    L->sort(3142, K, 1.0, 0.0);
    K.reset();
    GFvoA->gemv(false, L, ZvoA, -1.0, 1.0);
    L.reset();

    // F_AI += 2 * \sum_{em} (AI|em) Z_em
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|ov)", noccA, nvirA, noccB, nvirB));
    tei_ovov_chem_ref_directAB(K);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (VO|vo)", nvirA, noccA, nvirB, noccB));
    L->sort(2143, K, 1.0, 0.0);
    K.reset();
    GFvoA->gemv(false, L, ZvoB, 2.0, 1.0);
    // F_ai += 2 * \sum_{EM} (ai|EM) Z_EM
    GFvoB->gemv(true, L, ZvoA, 2.0, 1.0);
    L.reset();

    // Set global GF
    GFA->set_vo(GFvoA);

    // vo Block
    // F_ai += 2 * \sum_{em} (ai|em) Z_em
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (ov|ov)", noccB, nvirB, noccB, nvirB));
    tei_ovov_chem_ref_directBB(K);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (vo|vo)", nvirB, noccB, nvirB, noccB));
    L->sort(2143, K, 1.0, 0.0);
    GFvoB->gemv(false, L, ZvoB, 2.0, 1.0);
    L.reset();

    // F_ai -= \sum_{em} <ai|me> Z_me: (ie|ma) => <ai|me>
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <vo|ov>", nvirB, noccB, noccB, nvirB));
    L->sort(4132, K, 1.0, 0.0);
    K.reset();
    GFvoB->gemv(false, L, ZovB, -1.0, 1.0);
    L.reset();

    // F_ai -= \sum_{em} <ai|em> Z_em
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (oo|vv)", noccB, noccB, nvirB, nvirB));
    tei_oovv_chem_ref_directBB(K);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <vo|vo>", nvirB, noccB, nvirB, noccB));
    L->sort(3142, K, 1.0, 0.0);
    K.reset();
    GFvoB->gemv(false, L, ZvoB, -1.0, 1.0);
    L.reset();

    /*
    // F_ai += 2 * \sum_{EM} (ai|EM) Z_EM
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|ov)", noccA, nvirA, noccB, nvirB));
    tei_ovov_chem_ref_directAB(K);
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (VO|vo)", nvirA, noccA, nvirB, noccB));
    L->sort(2143, K, 1.0, 0.0);
    K.reset();
    GFvoB->gemv(true, L, ZvoA, 2.0, 1.0);
    L.reset();
    */

    // Set global GF
    GFB->set_vo(GFvoB);

    // OV Block and remaining part of the VO block
    // F_IA += Z_AI f_II
    // F_AI += Z_AI f_AA
    #pragma omp parallel for
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              GFA->add(i, a + noccA, ZvoA->get(a, i) * FockA->get(i,i));
              GFA->add(a + noccA, i, ZvoA->get(a, i) * FockA->get(a + noccA, a + noccA));
	 }
    }

    // ov Block and remaining part of the vo block
    // F_ia += z_ai f_ii
    // F_ai += z_ai f_aa
    #pragma omp parallel for
    for (int a = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++) {
              GFB->add(i, a + noccB, ZvoB->get(a, i) * FockB->get(i,i));
              GFB->add(a + noccB, i, ZvoB->get(a, i) * FockB->get(a + noccB, a + noccB));
	 }
    }

}// else if (reference_ == "UNRESTRICTED")
    timer_off("effective_gfm");
} // end effective_gfm

//=======================================================
//          Effective PDMs & GFM
//=======================================================
void DFOCC::effective_pdm_gfm()
{
    outfile->Printf("\tComputing the relaxed PDMs and GFM...\n");


    SharedTensor2d K, L, G, Gsep, Z, Z2;
    timer_on("Gamma^eff");
if (reference_ == "RESTRICTED") {
    //=========================
    // OPDM
    //=========================
    G1->add_vo(ZvoA, 2.0, 1.0);
    G1->add_ov(ZovA, 2.0, 1.0);

    //=========================
    // Seprable TPDM & GFM
    //=========================
    // Intermediates
    // Z_Q = 4 \sum_{ai} b_{ai}^{Q} Z_ai
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VO)", nQ_ref, nvirA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L->swap_3index_col(K);
    K.reset();
    SharedTensor1d Zq = SharedTensor1d(new Tensor1d("DF_BASIS_SCF Z_Q", nQ_ref));
    Zq->gemv(false, L, ZvoA, 4.0, 0.0);

    // Z_ij^Q = 4 * \sum_{a} Z_ia b_aj^Q
    Z = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|OO)", nQ_ref, noccA, noccA));
    Z->contract233(false, false, noccA, noccA, ZovA, L, 4.0, 0.0);

    // GFM VO Block
    // F_ai += 2 \sum_{Q} b_ai^Q Z_Q
    GFvo->gemv(true, L, Zq, 2.0, 1.0);
    L.reset();

    // F_ai -= 1/2 \sum_{Q}\sum_{m} b_ma^Q Z_mi^Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS);
    GFvo->contract(true, false, nvirA, noccA, nQ_ref * noccA, K, Z, -0.5, 1.0);
    K.reset();

    // GFM OO Block
    // F_ij += 2 \sum_{Q} b_ij^Q Z_Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    GFoo->gemv(true, K, Zq, 2.0, 1.0);

    // F_ij -= 1/2 \sum_{Q}\sum_{m} b_mi^Q Z_mj^Q
    GFoo->contract(true, false, noccA, noccA, nQ_ref * noccA, K, Z, -0.5, 1.0);
    K.reset();

    // Symmetrize Z_ij^Q for TPDM
    Z->symmetrize3(Z);

    // TPDM OO Block
    //  G_ij^Q += 2 Z_Q \delta_{ij} - 1/2 (Z_ij^Q + Z_ji^Q)
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OO)", nQ_ref, noccA, noccA));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    Gsep->subtract(Z);
    Z.reset();
    #pragma omp parallel for
    for (int Q = 0; Q < nQ_ref; Q++) {
         for (int i = 0; i < noccA; i++) {
              int ii = oo_idxAA->get(i,i);
              Gsep->add(Q, ii, 2.0 * Zq->get(Q));
         }
    }
    Zq.reset();
    Gsep->write(psio_, PSIF_DFOCC_DENS);
    Gsep.reset();

    // Z intermediate VO Block
    // Z_ai^Q = 2 * \sum_{m} Z_am b_mi^Q
    Z = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|VO)", nQ_ref, nvirA, noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    Z->contract233(false, false, nvirA, noccA, ZvoA, K, 2.0, 0.0);
    K.reset();

    // GFM VO Block
    // F_ai -= \sum_{Q} \sum_{e} b_ea^Q Z_ei^Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    GFvo->contract(true, false, nvirA, noccA, nQ_ref * nvirA, K, Z, -1.0, 1.0);
    K.reset();

    // GFM OO Block
    // F_ij -= \sum_{Q} \sum_{e} b_ei^Q Z_ej^Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VO)", nQ_ref, nvirA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L->swap_3index_col(K);
    K.reset();
    GFoo->contract(true, false, noccA, noccA, nQ_ref * nvirA, L, Z, -1.0, 1.0);
    L.reset();

    // TPDM OV and VO Blocks
    // G_ai^Q += 2 * J_Q Zai - Z_ai^Q
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VO)", nQ_ref, nvirA, noccA));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    Gsep->subtract(Z);
    Gsep->dirprd123(Jc, ZvoA, 2.0, 1.0);
    Gsep->write(psio_, PSIF_DFOCC_DENS);
    Gsep.reset();

    // G_ia^Q += 2 * J_Q Zia - Z_ia^Q : where Z_ia^Q = Z_ai^Q
    Z2 = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|OV)", nQ_ref, noccA, nvirA));
    Z2->swap_3index_col(Z);
    Z.reset();
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OV)", nQ_ref, noccA, nvirA));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    Gsep->subtract(Z2);
    Z2.reset();
    Gsep->dirprd123(Jc, ZovA, 2.0, 1.0);
    Gsep->write(psio_, PSIF_DFOCC_DENS);
    Gsep.reset();

    // Final GFM terms
    // Set global GF
    GF->set_oo(GFoo);
    GF->set_vo(GFvo);

    // OV Block and remaining part of the VO block
    // F_ia += 2.0 * z_ai f_ii
    // F_ai += 2.0 * z_ai f_aa
    #pragma omp parallel for
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              GF->add(i, a + noccA, 2.0 * ZvoA->get(a, i) * FockA->get(i,i));
              GF->add(a + noccA, i, 2.0 * ZvoA->get(a, i) * FockA->get(a + noccA, a + noccA));
	 }
    }

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    //=========================
    // OPDM
    //=========================
    // Alpha
    G1A->add_vo(ZvoA, 1.0, 1.0);
    G1A->add_ov(ZovA, 1.0, 1.0);

    // Beta
    G1B->add_vo(ZvoB, 1.0, 1.0);
    G1B->add_ov(ZovB, 1.0, 1.0);

    //=========================
    // Seprable TPDM
    //=========================
    // Intermediates
    // Z_Q = 2 \sum_{AI} b_{AI}^{Q} Z_AI
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VO)", nQ_ref, nvirA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L->swap_3index_col(K);
    K.reset();
    SharedTensor1d Zq = SharedTensor1d(new Tensor1d("DF_BASIS_SCF Z_Q", nQ_ref));
    Zq->gemv(false, L, ZvoA, 2.0, 0.0);

    // Z_IJ^Q = 2 * \sum_{A} Z_IA b_AJ^Q
    Z = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|OO)", nQ_ref, noccA, noccA));
    Z->contract233(false, false, noccA, noccA, ZovA, L, 2.0, 0.0);
    L.reset();

    // GFM VO Block
    // F_AI -= 1/2 \sum_{Q}\sum_{M} b_MA^Q Z_MI^Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS);
    GFvoA->contract(true, false, nvirA, noccA, nQ_ref * noccA, K, Z, -0.5, 1.0);
    K.reset();

    // GFM OO Block
    // F_IJ -= 1/2 \sum_{Q}\sum_{M} b_MI^Q Z_MJ^Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    GFooA->contract(true, false, noccA, noccA, nQ_ref * noccA, K, Z, -0.5, 1.0);
    K.reset();

    // Symmetrize Z_IJ^Q
    Z->symmetrize3(Z);
    Z->write(psio_, PSIF_DFOCC_DENS);
    Z.reset();

    // Z_Q += 2 \sum_{ai} b_{ai}^{Q} Z_ai
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB, nvirB));
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vo)", nQ_ref, nvirB, noccB));
    K->read(psio_, PSIF_DFOCC_INTS);
    L->swap_3index_col(K);
    K.reset();
    Zq->gemv(false, L, ZvoB, 2.0, 1.0);

    // Z_ij^Q = 2 * \sum_{a} Z_ia b_aj^Q
    Z = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|oo)", nQ_ref, noccB, noccB));
    Z->contract233(false, false, noccB, noccB, ZovB, L, 2.0, 0.0);

    // GFM VO Block
    // F_ai += \sum_{Q} b_ai^Q Z_Q
    GFvoB->gemv(true, L, Zq, 1.0, 1.0);
    L.reset();

    // F_AI += \sum_{Q} b_AI^Q Z_Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VO)", nQ_ref, nvirA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L->swap_3index_col(K);
    K.reset();
    GFvoA->gemv(true, L, Zq, 1.0, 1.0);
    L.reset();

    // F_ai -= 1/2 \sum_{Q}\sum_{m} b_ma^Q Z_mi^Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB, nvirB));
    K->read(psio_, PSIF_DFOCC_INTS);
    GFvoB->contract(true, false, nvirB, noccB, nQ_ref * noccB, K, Z, -0.5, 1.0);
    K.reset();

    // GFM OO Block
    // F_IJ += \sum_{Q} b_IJ^Q Z_Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    GFooA->gemv(true, K, Zq, 1.0, 1.0);
    K.reset();

    // F_ij += \sum_{Q} b_ij^Q Z_Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB, noccB));
    K->read(psio_, PSIF_DFOCC_INTS);
    GFooB->gemv(true, K, Zq, 1.0, 1.0);

    // F_ij -= 1/2 \sum_{Q}\sum_{m} b_mi^Q Z_mj^Q
    GFooB->contract(true, false, noccB, noccB, nQ_ref * noccB, K, Z, -0.5, 1.0);
    K.reset();

    // Symmetrize Z_ij^Q
    Z->symmetrize3(Z);
    Z->write(psio_, PSIF_DFOCC_DENS);
    Z.reset();

    // TPDM OO Block
    // G_IJ^Q += Z_Q \delta_{IJ} - 1/2 (Z_IJ^Q + Z_JI^Q)
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OO)", nQ_ref, noccA, noccA));
    Z = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|OO)", nQ_ref, noccA, noccA));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    Z->read(psio_, PSIF_DFOCC_DENS);
    Gsep->subtract(Z);
    Z.reset();
    #pragma omp parallel for
    for (int Q = 0; Q < nQ_ref; Q++) {
         for (int i = 0; i < noccA; i++) {
              int ii = oo_idxAA->get(i,i);
              Gsep->add(Q, ii, Zq->get(Q));
         }
    }
    Gsep->write(psio_, PSIF_DFOCC_DENS);
    Gsep.reset();

    // TPDM oo Block
    //  G_ij^Q += Z_Q \delta_{ij} - 1/2 (Z_ij^Q + Z_ji^Q)
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|oo)", nQ_ref, noccB, noccB));
    Z = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|oo)", nQ_ref, noccB, noccB));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    Z->read(psio_, PSIF_DFOCC_DENS);
    Gsep->subtract(Z);
    Z.reset();
    #pragma omp parallel for
    for (int Q = 0; Q < nQ_ref; Q++) {
         for (int i = 0; i < noccB; i++) {
              int ii = oo_idxBB->get(i,i);
              Gsep->add(Q, ii, Zq->get(Q));
         }
    }
    Zq.reset();
    Gsep->write(psio_, PSIF_DFOCC_DENS);
    Gsep.reset();

    // Z Int VO block
    // Z_AI^Q = \sum_{M} Z_AM b_MI^Q
    Z = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|VO)", nQ_ref, nvirA, noccA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    Z->contract233(false, false, nvirA, noccA, ZvoA, K, 1.0, 0.0);
    K.reset();

    // GFM VO Block
    // F_AI -= \sum_{Q} \sum_{E} b_EA^Q Z_EI^Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    GFvoA->contract(true, false, nvirA, noccA, nQ_ref * nvirA, K, Z, -1.0, 1.0);
    K.reset();

    // GFM OO Block
    // F_IJ -= \sum_{Q} \sum_{E} b_EI^Q Z_EJ^Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA, nvirA));
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VO)", nQ_ref, nvirA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    L->swap_3index_col(K);
    K.reset();
    GFooA->contract(true, false, noccA, noccA, nQ_ref * nvirA, L, Z, -1.0, 1.0);
    L.reset();

    // TPDM OV and VO Blocks
    // G_AI^Q += J_Q Z_AI - Z_AI^Q
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|VO)", nQ_ref, nvirA, noccA));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    Gsep->subtract(Z);
    Gsep->dirprd123(Jc, ZvoA, 1.0, 1.0);
    Gsep->write(psio_, PSIF_DFOCC_DENS);
    Gsep.reset();

    // G_IA^Q += J_Q Z_IA - Z_IA^Q : where Z_IA^Q = Z_AI^Q
    Z2 = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|OV)", nQ_ref, noccA, nvirA));
    Z2->swap_3index_col(Z);
    Z.reset();
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|OV)", nQ_ref, noccA, nvirA));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    Gsep->subtract(Z2);
    Z2.reset();
    Gsep->dirprd123(Jc, ZovA, 1.0, 1.0);
    Gsep->write(psio_, PSIF_DFOCC_DENS);
    Gsep.reset();

    // Z Int vo block
    // Z_ai^Q = \sum_{m} Z_am b_mi^Q
    Z = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|vo)", nQ_ref, nvirB, noccB));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB, noccB));
    K->read(psio_, PSIF_DFOCC_INTS);
    Z->contract233(false, false, nvirB, noccB, ZvoB, K, 1.0, 0.0);
    K.reset();

    // GFM vo
    // F_ai -= \sum_{Q} \sum_{e} b_ea^Q Z_ei^Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vv)", nQ_ref, nvirB, nvirB));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    GFvoB->contract(true, false, nvirB, noccB, nQ_ref * nvirB, K, Z, -1.0, 1.0);
    K.reset();

    // GFM oo
    // F_ij -= \sum_{Q} \sum_{e} b_ei^Q Z_ej^Q
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB, nvirB));
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vo)", nQ_ref, nvirB, noccB));
    K->read(psio_, PSIF_DFOCC_INTS);
    L->swap_3index_col(K);
    K.reset();
    GFooB->contract(true, false, noccB, noccB, nQ_ref * nvirB, L, Z, -1.0, 1.0);
    L.reset();

    // TPDM ov and vo Blocks
    // G_ai^Q += J_Q Zai - Z_ai^Q
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|vo)", nQ_ref, nvirB, noccB));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    Gsep->subtract(Z);
    Gsep->dirprd123(Jc, ZvoB, 1.0, 1.0);
    Gsep->write(psio_, PSIF_DFOCC_DENS);
    Gsep.reset();

    // G_ia^Q += J_Q Zia - Z_ia^Q : where Z_ia^Q = Z_ai^Q
    Z2 = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Z (Q|ov)", nQ_ref, noccB, nvirB));
    Z2->swap_3index_col(Z);
    Z.reset();
    Gsep = SharedTensor2d(new Tensor2d("3-Index Separable TPDM (Q|ov)", nQ_ref, noccB, nvirB));
    Gsep->read(psio_, PSIF_DFOCC_DENS);
    Gsep->subtract(Z2);
    Z2.reset();
    Gsep->dirprd123(Jc, ZovB, 1.0, 1.0);
    Gsep->write(psio_, PSIF_DFOCC_DENS);
    Gsep.reset();

    // Final GFM terms
    // Set global GF
    GFA->set_oo(GFooA);
    GFB->set_oo(GFooB);
    GFA->set_vo(GFvoA);
    GFB->set_vo(GFvoB);

    // OV Block and remaining part of the VO block
    // F_IA += Z_AI f_II
    // F_AI += Z_AI f_AA
    #pragma omp parallel for
    for (int a = 0; a < nvirA; a++) {
         for (int i = 0; i < noccA; i++) {
              GFA->add(i, a + noccA, ZvoA->get(a, i) * FockA->get(i,i));
              GFA->add(a + noccA, i, ZvoA->get(a, i) * FockA->get(a + noccA, a + noccA));
	 }
    }

    // ov Block and remaining part of the vo block
    // F_ia += z_ai f_ii
    // F_ai += z_ai f_aa
    #pragma omp parallel for
    for (int a = 0; a < nvirB; a++) {
         for (int i = 0; i < noccB; i++) {
              GFB->add(i, a + noccB, ZvoB->get(a, i) * FockB->get(i,i));
              GFB->add(a + noccB, i, ZvoB->get(a, i) * FockB->get(a + noccB, a + noccB));
	 }
    }

}// else if (reference_ == "UNRESTRICTED")
    timer_off("Gamma^eff");
} // end effective_pdm_gfm

}} // End Namespaces
