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

void DFOCC::mp2_direct()
{
    SharedTensor2d K, L, M, T, U;
    timer_on("MP2");
if (reference_ == "RESTRICTED") {
    // Build amplitudes in Mulliken order
    T = SharedTensor2d(new Tensor2d("T2_1 (ia|jb)", naoccA, navirA, naoccA, navirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(K);
    T->copy(K);
    T->apply_denom_chem(nfrzc, noccA, FockA);

    // form U(ia,jb)
    U = SharedTensor2d(new Tensor2d("U2_1 (ia|jb)", naoccA, navirA, naoccA, navirA));
    U->sort(1432, T, 1.0, 0.0);
    U->scale(-1.0);
    U->axpy(T, 2.0);
    T.reset();
    Ecorr = U->vector_dot(K);
    U.reset();
    K.reset();
    Emp2 = Eref + Ecorr;

    /*
    // Form U(Q,ia) = \sum_{jb} b_jb^Q u_jb^ia
    T = SharedTensor2d(new Tensor2d("U2_1 (Q|IA)", nQ, naoccA * navirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA * navirA));
    T->gemm(false, false, K, U, 1.0, 0.0);
    U.reset();
    Ecorr = T->vector_dot(K);
    T.reset();
    K.reset();
    Emp2 = Eref + Ecorr;
    */

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // T2AA
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(L);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ||AB>", naoccA, naoccA, navirA, navirA));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    t2_1AA = SharedTensor2d(new Tensor2d("T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    t2_1AA->copy(K);
    t2_1AA->apply_denom(nfrzc, noccA, FockA);
    Emp2AA = 0.25 * t2_1AA->vector_dot(K);
    K.reset();
    t2_1AA.reset();
    Escsmp2AA = ss_scale * Emp2AA;
    Escsnmp2AA = 1.76 * Emp2AA;

    // T2BB
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB));
    tei_iajb_chem_directBB(L);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|ab>", naoccB, naoccB, navirB, navirB));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij||ab>", naoccB, naoccB, navirB, navirB));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    t2_1BB = SharedTensor2d(new Tensor2d("T2_1 <ij|ab>", naoccB, naoccB, navirB, navirB));
    t2_1BB->copy(K);
    t2_1BB->apply_denom(nfrzc, noccB, FockB);
    Emp2BB = 0.25 * t2_1BB->vector_dot(K);
    K.reset();
    t2_1BB.reset();
    Escsmp2BB = ss_scale * Emp2BB;
    Escsnmp2BB = 1.76 * Emp2BB;

    // T2AB
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|jb)", naoccA, navirA, naoccB, navirB));
    tei_iajb_chem_directAB(L);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    t2_1AB = SharedTensor2d(new Tensor2d("T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    t2_1AB->copy(K);
    t2_1AB->apply_denom_os(nfrzc, noccA, noccB, FockA, FockB);
    Emp2AB = t2_1AB->vector_dot(K);
    K.reset();
    t2_1AB.reset();
    Escsmp2AB = os_scale * Emp2AB;
    if (mo_optimized == 0) Esosmp2AB = sos_scale * Emp2AB;
    else if (mo_optimized == 1) Esosmp2AB = sos_scale2 * Emp2AB;

    Emp2_t1 = 0.0;
    if (reference == "ROHF" && orb_opt_ == "FALSE") {
        if (wfn_type_ == "DF-OMP2" || wfn_type_ == "CD-OMP2") {
        //Singles-contribution
        //Alpha
        for(int i = 0 ; i < naoccA; ++i){
            for(int a = 0 ; a < navirA; ++a){
                Emp2_t1 += t1A->get(i, a) * FockA->get(a + noccA, i + nfrzc);
            }
        }

        // Beta
        for(int i = 0 ; i < naoccB; ++i){
            for(int a = 0 ; a < navirB; ++a){
                Emp2_t1 += t1B->get(i, a) * FockB->get(a + noccB, i + nfrzc);
            }
        }
        }
    }// end if (reference == "ROHF")

    Ecorr = Emp2AA + Emp2AB + Emp2BB + Emp2_t1;
    Emp2 = Eref + Ecorr;
    Escsmp2 = Eref + Escsmp2AA + Escsmp2AB + Escsmp2BB;
    Esosmp2 = Eref + Esosmp2AB;
    Escsnmp2 = Eref + Escsnmp2AA + Escsnmp2BB;
}// else if (reference_ == "UNRESTRICTED")
    timer_off("MP2");
} // end mp2_direct

//======================================================================
//    MP2 for CCSD
//======================================================================
void DFOCC::ccsd_mp2()
{
    SharedTensor2d K, L, M, T, U, Tau;
    timer_on("MP2");
if (reference_ == "RESTRICTED") {
    // Build amplitudes in Mulliken order
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    t2->copy(K);
    t2->apply_denom_chem(nfrzc, noccA, FockA);

    // Form U(ia,jb) = 2*T(ia,jb) - T (ib,ja)
    U = SharedTensor2d(new Tensor2d("U2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_u2_amps(U,t2);
    Ecorr = U->vector_dot(K);
    U.reset();
    K.reset();
    Emp2 = Eref + Ecorr;
    Eccsd = Emp2;

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // T2AA
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(L);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ||AB>", naoccA, naoccA, navirA, navirA));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    t2_1AA = SharedTensor2d(new Tensor2d("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    t2_1AA->copy(K);
    t2_1AA->apply_denom(nfrzc, noccA, FockA);
    t2_1AA->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Emp2AA = 0.25 * t2_1AA->vector_dot(K);
    K.reset();
    t2_1AA.reset();
    Escsmp2AA = ss_scale * Emp2AA;
    Escsnmp2AA = 1.76 * Emp2AA;

    // T2BB
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB));
    tei_iajb_chem_directBB(L);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|ab>", naoccB, naoccB, navirB, navirB));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij||ab>", naoccB, naoccB, navirB, navirB));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    t2_1BB = SharedTensor2d(new Tensor2d("T2 <ij|ab>", naoccB, naoccB, navirB, navirB));
    t2_1BB->copy(K);
    t2_1BB->apply_denom(nfrzc, noccB, FockB);
    t2_1BB->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Emp2BB = 0.25 * t2_1BB->vector_dot(K);
    K.reset();
    t2_1BB.reset();
    Escsmp2BB = ss_scale * Emp2BB;
    Escsnmp2BB = 1.76 * Emp2BB;

    // T2AB
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|jb)", naoccA, navirA, naoccB, navirB));
    tei_iajb_chem_directAB(L);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    t2_1AB = SharedTensor2d(new Tensor2d("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    t2_1AB->copy(K);
    t2_1AB->apply_denom_os(nfrzc, noccA, noccB, FockA, FockB);
    t2_1AB->write(psio_, PSIF_DFOCC_AMPS);
    Emp2AB = t2_1AB->vector_dot(K);
    K.reset();
    t2_1AB.reset();
    Escsmp2AB = os_scale * Emp2AB;
    if (mo_optimized == 0) Esosmp2AB = sos_scale * Emp2AB;
    else if (mo_optimized == 1) Esosmp2AB = sos_scale2 * Emp2AB;

    Emp2_t1 = 0.0;
    if (reference == "ROHF" && orb_opt_ == "FALSE") {
        //Singles-contribution
        //Alpha
        for(int i = 0 ; i < naoccA; ++i){
            for(int a = 0 ; a < navirA; ++a){
                Emp2_t1 += t1A->get(i, a) * FockA->get(a + noccA, i + nfrzc);
            }
        }

        // Beta
        for(int i = 0 ; i < naoccB; ++i){
            for(int a = 0 ; a < navirB; ++a){
                Emp2_t1 += t1B->get(i, a) * FockB->get(a + noccB, i + nfrzc);
            }
        }
    }// end if (reference == "ROHF")

    Ecorr = Emp2AA + Emp2AB + Emp2BB + Emp2_t1;
    Emp2 = Eref + Ecorr;
    Escsmp2 = Eref + Escsmp2AA + Escsmp2AB + Escsmp2BB;
    Esosmp2 = Eref + Esosmp2AB;
    Escsnmp2 = Eref + Escsnmp2AA + Escsnmp2BB;
    Eccsd = Emp2;
}// else if (reference_ == "UNRESTRICTED")
    timer_off("MP2");
}// end of ccsd_mp2

//======================================================================
//    MP2 for CCSD
//======================================================================
void DFOCC::ccsd_mp2_low()
{
    SharedTensor2d K, L, M, T, U, Tau;
    timer_on("MP2");
if (reference_ == "RESTRICTED") {
    // Build amplitudes in Mulliken order
    T = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(K);
    T->copy(K);
    T->apply_denom_chem(nfrzc, noccA, FockA);
    T->write_symm(psio_, PSIF_DFOCC_AMPS);

    // Form T'(ib,ja) = T(ia,jb)
    U = SharedTensor2d(new Tensor2d("T2p (IA|JB)", naoccA, navirA, naoccA, navirA));
    U->sort(1432, T, 1.0, 0.0);
    U->write_symm(psio_, PSIF_DFOCC_AMPS);
    U.reset();

    // Form Tau(ia,jb) = T(ia,jb) + t(ia) * t(jb)
    U = SharedTensor2d(new Tensor2d("Tau (IA|JB)", naoccA, navirA, naoccA, navirA));
    U->dirprd224(t1A, t1A);
    U->add(T);
    U->write_symm(psio_, PSIF_DFOCC_AMPS);
    U.reset();

    // Form U(ia,jb) = 2*T(ia,jb) - T (ib,ja)
    U = SharedTensor2d(new Tensor2d("U2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    U->sort(1432, T, 1.0, 0.0);
    U->scale(-1.0);
    U->axpy(T, 2.0);
    T.reset();
    U->write_symm(psio_, PSIF_DFOCC_AMPS);
    Ecorr = U->vector_dot(K);
    U.reset();
    K.reset();
    Emp2 = Eref + Ecorr;
    Eccsd = Emp2;

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // T2AA
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(L);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ||AB>", naoccA, naoccA, navirA, navirA));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    t2_1AA = SharedTensor2d(new Tensor2d("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    t2_1AA->copy(K);
    t2_1AA->apply_denom(nfrzc, noccA, FockA);
    t2_1AA->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Emp2AA = 0.25 * t2_1AA->vector_dot(K);
    K.reset();
    t2_1AA.reset();
    Escsmp2AA = ss_scale * Emp2AA;
    Escsnmp2AA = 1.76 * Emp2AA;

    // T2BB
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB));
    tei_iajb_chem_directBB(L);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|ab>", naoccB, naoccB, navirB, navirB));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij||ab>", naoccB, naoccB, navirB, navirB));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    t2_1BB = SharedTensor2d(new Tensor2d("T2 <ij|ab>", naoccB, naoccB, navirB, navirB));
    t2_1BB->copy(K);
    t2_1BB->apply_denom(nfrzc, noccB, FockB);
    t2_1BB->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Emp2BB = 0.25 * t2_1BB->vector_dot(K);
    K.reset();
    t2_1BB.reset();
    Escsmp2BB = ss_scale * Emp2BB;
    Escsnmp2BB = 1.76 * Emp2BB;

    // T2AB
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|jb)", naoccA, navirA, naoccB, navirB));
    tei_iajb_chem_directAB(L);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    t2_1AB = SharedTensor2d(new Tensor2d("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    t2_1AB->copy(K);
    t2_1AB->apply_denom_os(nfrzc, noccA, noccB, FockA, FockB);
    t2_1AB->write(psio_, PSIF_DFOCC_AMPS);
    Emp2AB = t2_1AB->vector_dot(K);
    K.reset();
    t2_1AB.reset();
    Escsmp2AB = os_scale * Emp2AB;
    if (mo_optimized == 0) Esosmp2AB = sos_scale * Emp2AB;
    else if (mo_optimized == 1) Esosmp2AB = sos_scale2 * Emp2AB;

    Emp2_t1 = 0.0;
    if (reference == "ROHF" && orb_opt_ == "FALSE") {
        //Singles-contribution
        //Alpha
        for(int i = 0 ; i < naoccA; ++i){
            for(int a = 0 ; a < navirA; ++a){
                Emp2_t1 += t1A->get(i, a) * FockA->get(a + noccA, i + nfrzc);
            }
        }

        // Beta
        for(int i = 0 ; i < naoccB; ++i){
            for(int a = 0 ; a < navirB; ++a){
                Emp2_t1 += t1B->get(i, a) * FockB->get(a + noccB, i + nfrzc);
            }
        }
    }// end if (reference == "ROHF")

    Ecorr = Emp2AA + Emp2AB + Emp2BB + Emp2_t1;
    Emp2 = Eref + Ecorr;
    Escsmp2 = Eref + Escsmp2AA + Escsmp2AB + Escsmp2BB;
    Esosmp2 = Eref + Esosmp2AB;
    Escsnmp2 = Eref + Escsnmp2AA + Escsnmp2BB;
    Eccsd = Emp2;
}// else if (reference_ == "UNRESTRICTED")
    timer_off("MP2");
}// end of ccsd_mp2_low

//======================================================================
//    MP2 for CCD
//======================================================================
void DFOCC::ccd_mp2()
{
    SharedTensor2d K, L, M, T, U, Tau;
    timer_on("MP2");
if (reference_ == "RESTRICTED") {
    // Build amplitudes in Mulliken order
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    t2->copy(K);
    t2->apply_denom_chem(nfrzc, noccA, FockA);

    // Form U(ia,jb) = 2*T(ia,jb) - T (ib,ja)
    U = SharedTensor2d(new Tensor2d("U2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    ccsd_u2_amps(U,t2);
    Ecorr = U->vector_dot(K);
    U.reset();
    K.reset();
    Emp2 = Eref + Ecorr;

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // T2AA
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(L);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ||AB>", naoccA, naoccA, navirA, navirA));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    t2_1AA = SharedTensor2d(new Tensor2d("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    t2_1AA->copy(K);
    t2_1AA->apply_denom(nfrzc, noccA, FockA);
    t2_1AA->write_anti_symm(psio_, PSIF_DFOCC_AMPS);

    // Form T2(IA|JB)
    t2p_1 = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    t2p_1->sort(1324, t2_1AA, 1.0, 0.0);
    t2p_1->write_symm(psio_, PSIF_DFOCC_AMPS);
    //t2p_1->write(psio_, PSIF_DFOCC_AMPS);
    t2p_1.reset();

    // Energy
    Emp2AA = 0.25 * t2_1AA->vector_dot(K);
    K.reset();
    t2_1AA.reset();
    Escsmp2AA = ss_scale * Emp2AA;
    Escsnmp2AA = 1.76 * Emp2AA;

    // T2BB
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB));
    tei_iajb_chem_directBB(L);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|ab>", naoccB, naoccB, navirB, navirB));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij||ab>", naoccB, naoccB, navirB, navirB));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    t2_1BB = SharedTensor2d(new Tensor2d("T2 <ij|ab>", naoccB, naoccB, navirB, navirB));
    t2_1BB->copy(K);
    t2_1BB->apply_denom(nfrzc, noccB, FockB);
    t2_1BB->write_anti_symm(psio_, PSIF_DFOCC_AMPS);

    // Form T2(ia|jb)
    t2p_1 = SharedTensor2d(new Tensor2d("T2 (ia|jb)", naoccB, navirB, naoccB, navirB));
    t2p_1->sort(1324, t2_1BB, 1.0, 0.0);
    t2p_1->write_symm(psio_, PSIF_DFOCC_AMPS);
    t2p_1.reset();

    // Energy
    Emp2BB = 0.25 * t2_1BB->vector_dot(K);
    K.reset();
    t2_1BB.reset();
    Escsmp2BB = ss_scale * Emp2BB;
    Escsnmp2BB = 1.76 * Emp2BB;

    // T2AB
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|jb)", naoccA, navirA, naoccB, navirB));
    tei_iajb_chem_directAB(L);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    t2_1AB = SharedTensor2d(new Tensor2d("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    t2_1AB->copy(K);
    t2_1AB->apply_denom_os(nfrzc, noccA, noccB, FockA, FockB);
    t2_1AB->write(psio_, PSIF_DFOCC_AMPS);

    // Form T2(IA|jb)
    t2p_1 = SharedTensor2d(new Tensor2d("T2 (IA|jb)", naoccA, navirA, naoccB, navirB));
    t2p_1->sort(1324, t2_1AB, 1.0, 0.0);
    t2p_1->write(psio_, PSIF_DFOCC_AMPS);
    t2p_1.reset();

    // Energy
    Emp2AB = t2_1AB->vector_dot(K);
    K.reset();
    t2_1AB.reset();
    Escsmp2AB = os_scale * Emp2AB;
    if (mo_optimized == 0) Esosmp2AB = sos_scale * Emp2AB;
    else if (mo_optimized == 1) Esosmp2AB = sos_scale2 * Emp2AB;

    Emp2_t1 = 0.0;
    if (reference == "ROHF" && orb_opt_ == "FALSE") {
        //Singles-contribution
        //Alpha
        for(int i = 0 ; i < naoccA; ++i){
            for(int a = 0 ; a < navirA; ++a){
                Emp2_t1 += t1A->get(i, a) * FockA->get(a + noccA, i + nfrzc);
            }
        }

        // Beta
        for(int i = 0 ; i < naoccB; ++i){
            for(int a = 0 ; a < navirB; ++a){
                Emp2_t1 += t1B->get(i, a) * FockB->get(a + noccB, i + nfrzc);
            }
        }
    }// end if (reference == "ROHF")

    Ecorr = Emp2AA + Emp2AB + Emp2BB + Emp2_t1;
    Emp2 = Eref + Ecorr;
    Escsmp2 = Eref + Escsmp2AA + Escsmp2AB + Escsmp2BB;
    Esosmp2 = Eref + Esosmp2AB;
    Escsnmp2 = Eref + Escsnmp2AA + Escsnmp2BB;
}// else if (reference_ == "UNRESTRICTED")
    Eccd = Emp2;
    Elccd = Emp2;
    timer_off("MP2");
}// end of ccd_mp2

//======================================================================
//    MP2 for CCD
//======================================================================
void DFOCC::ccd_mp2_low()
{
    SharedTensor2d K, L, M, T, U, Tau;
    timer_on("MP2");
if (reference_ == "RESTRICTED") {
    // Build amplitudes in Mulliken order
    T = SharedTensor2d(new Tensor2d("T2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(K);
    T->copy(K);
    T->apply_denom_chem(nfrzc, noccA, FockA);
    T->write_symm(psio_, PSIF_DFOCC_AMPS);

    // Form T'(ib,ja) = T(ia,jb)
    U = SharedTensor2d(new Tensor2d("T2p (IA|JB)", naoccA, navirA, naoccA, navirA));
    U->sort(1432, T, 1.0, 0.0);
    U->write_symm(psio_, PSIF_DFOCC_AMPS);
    U.reset();

    // Form U(ia,jb) = 2*T(ia,jb) - T (ib,ja)
    U = SharedTensor2d(new Tensor2d("U2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    U->sort(1432, T, 1.0, 0.0);
    U->scale(-1.0);
    U->axpy(T, 2.0);
    T.reset();
    U->write_symm(psio_, PSIF_DFOCC_AMPS);
    Ecorr = U->vector_dot(K);
    U.reset();
    K.reset();
    Emp2 = Eref + Ecorr;

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // T2AA
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(L);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ||AB>", naoccA, naoccA, navirA, navirA));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    t2_1AA = SharedTensor2d(new Tensor2d("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    t2_1AA->copy(K);
    t2_1AA->apply_denom(nfrzc, noccA, FockA);
    t2_1AA->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Emp2AA = 0.25 * t2_1AA->vector_dot(K);
    K.reset();
    t2_1AA.reset();
    Escsmp2AA = ss_scale * Emp2AA;
    Escsnmp2AA = 1.76 * Emp2AA;

    // T2BB
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB));
    tei_iajb_chem_directBB(L);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|ab>", naoccB, naoccB, navirB, navirB));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij||ab>", naoccB, naoccB, navirB, navirB));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    t2_1BB = SharedTensor2d(new Tensor2d("T2 <ij|ab>", naoccB, naoccB, navirB, navirB));
    t2_1BB->copy(K);
    t2_1BB->apply_denom(nfrzc, noccB, FockB);
    t2_1BB->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Emp2BB = 0.25 * t2_1BB->vector_dot(K);
    K.reset();
    t2_1BB.reset();
    Escsmp2BB = ss_scale * Emp2BB;
    Escsnmp2BB = 1.76 * Emp2BB;

    // T2AB
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|jb)", naoccA, navirA, naoccB, navirB));
    tei_iajb_chem_directAB(L);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    t2_1AB = SharedTensor2d(new Tensor2d("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    t2_1AB->copy(K);
    t2_1AB->apply_denom_os(nfrzc, noccA, noccB, FockA, FockB);
    t2_1AB->write(psio_, PSIF_DFOCC_AMPS);
    Emp2AB = t2_1AB->vector_dot(K);
    K.reset();
    t2_1AB.reset();
    Escsmp2AB = os_scale * Emp2AB;
    if (mo_optimized == 0) Esosmp2AB = sos_scale * Emp2AB;
    else if (mo_optimized == 1) Esosmp2AB = sos_scale2 * Emp2AB;

    Emp2_t1 = 0.0;
    if (reference == "ROHF" && orb_opt_ == "FALSE") {
        //Singles-contribution
        //Alpha
        for(int i = 0 ; i < naoccA; ++i){
            for(int a = 0 ; a < navirA; ++a){
                Emp2_t1 += t1A->get(i, a) * FockA->get(a + noccA, i + nfrzc);
            }
        }

        // Beta
        for(int i = 0 ; i < naoccB; ++i){
            for(int a = 0 ; a < navirB; ++a){
                Emp2_t1 += t1B->get(i, a) * FockB->get(a + noccB, i + nfrzc);
            }
        }
    }// end if (reference == "ROHF")

    Ecorr = Emp2AA + Emp2AB + Emp2BB + Emp2_t1;
    Emp2 = Eref + Ecorr;
    Escsmp2 = Eref + Escsmp2AA + Escsmp2AB + Escsmp2BB;
    Esosmp2 = Eref + Esosmp2AB;
    Escsnmp2 = Eref + Escsnmp2AA + Escsnmp2BB;
}// else if (reference_ == "UNRESTRICTED")
    Eccd = Emp2;
    Elccd = Emp2;
    timer_off("MP2");
}// end of ccd_mp2_low

}} // End Namespaces
