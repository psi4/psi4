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

using namespace psi;
using namespace std;


namespace psi{ namespace dfoccwave{
  
void DFOCC::t2_1st_sc()
{   
    SharedTensor2d K, L, M, T, U;
    timer_on("1st-order T2");
if (reference_ == "RESTRICTED") {
    // Build amplitudes in Mulliken order 
    T = SharedTensor2d(new Tensor2d("T2_1(ia,jb)", naoccA, navirA, naoccA, navirA));
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    if (conv_tei_type == "DISK") K->read(psio_, PSIF_DFOCC_INTS);
    else tei_iajb_chem_directAA(K);
    T->copy(K);
    if (regularization == "FALSE") T->apply_denom_chem(nfrzc, noccA, FockA);
    else if (regularization == "TRUE") T->reg_denom_chem(nfrzc, noccA, FockA, reg_param);
    T->write(psio_, PSIF_DFOCC_AMPS);

    /*
    // Sort amplitudes to Dirac order
    t2_1 = SharedTensor2d(new Tensor2d("T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    t2_1->sort(1324, t2p_1, 1.0, 0.0);
    t2_1->write(psio_, PSIF_DFOCC_AMPS);
    t2_1.reset();
    */
 
    // form U(ia,jb)
    U = SharedTensor2d(new Tensor2d("2*T2_1(ia,jb) - T2_1(ib,ja)", naoccA, navirA, naoccA, navirA));
    U->sort(1432, T, 1.0, 0.0);
    U->scale(-1.0);
    U->axpy(T, 2.0);
    T.reset();
    Ecorr = U->vector_dot(K);
    U->write(psio_, PSIF_DFOCC_AMPS);
    U.reset();
    K.reset();
    Emp2 = Eref + Ecorr;

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // T2AA
    if (conv_tei_type == "DISK") {
        K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ||AB>", naoccA, naoccA, navirA, navirA));
        K->read(psio_, PSIF_DFOCC_INTS);
    }
    else {
        L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
        tei_iajb_chem_directAA(L);
        M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA));
        M->sort(1324, L, 1.0, 0.0);
        L.reset();
        K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ||AB>", naoccA, naoccA, navirA, navirA));
        tei_pqrs_anti_symm_direct(K, M);
        M.reset();
    }
    t2_1AA = SharedTensor2d(new Tensor2d("T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    t2_1AA->copy(K);
    if (regularization == "FALSE") t2_1AA->apply_denom(nfrzc, noccA, FockA);
    else if (regularization == "TRUE") t2_1AA->reg_denom(nfrzc, noccA, FockA, reg_param);
    t2_1AA->write(psio_, PSIF_DFOCC_AMPS);

    // Enegy AA 
    Emp2AA = 0.25 * t2_1AA->vector_dot(K);
    K.reset();
    Escsmp2AA = ss_scale * Emp2AA;
    Escsnmp2AA = 1.76 * Emp2AA;

    // sort AA
    t2p_1 = SharedTensor2d(new Tensor2d("T2_1(IA,JB)", naoccA, navirA, naoccA, navirA));
    t2p_1->sort(1324, t2_1AA, 1.0, 0.0);
    t2_1AA.reset();
    t2p_1->write(psio_, PSIF_DFOCC_AMPS);
    t2p_1.reset();

    // T2BB
    if (conv_tei_type == "DISK") {
        K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij||ab>", naoccB, naoccB, navirB, navirB));
        K->read(psio_, PSIF_DFOCC_INTS);
    }
    else {
        L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB));
        tei_iajb_chem_directBB(L);
        M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|ab>", naoccB, naoccB, navirB, navirB));
        M->sort(1324, L, 1.0, 0.0);
        L.reset();
        K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij||ab>", naoccB, naoccB, navirB, navirB));
        tei_pqrs_anti_symm_direct(K, M);
        M.reset();
    }
    t2_1BB = SharedTensor2d(new Tensor2d("T2_1 <ij|ab>", naoccB, naoccB, navirB, navirB));
    t2_1BB->copy(K);
    if (regularization == "FALSE") t2_1BB->apply_denom(nfrzc, noccB, FockB);
    else if (regularization == "TRUE") t2_1BB->reg_denom(nfrzc, noccB, FockB, reg_param);
    t2_1BB->write(psio_, PSIF_DFOCC_AMPS);

    // Energy BB
    Emp2BB = 0.25 * t2_1BB->vector_dot(K);
    K.reset();
    Escsmp2BB = ss_scale * Emp2BB;
    Escsnmp2BB = 1.76 * Emp2BB;

    // sort BB
    t2p_1 = SharedTensor2d(new Tensor2d("T2_1(ia,jb)", naoccB, navirB, naoccB, navirB));
    t2p_1->sort(1324, t2_1BB, 1.0, 0.0);
    t2_1BB.reset();
    t2p_1->write(psio_, PSIF_DFOCC_AMPS);
    t2p_1.reset();

    // T2AB
    if (conv_tei_type == "DISK") {
        K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ij|Ab>", naoccA, naoccB, navirA, navirB));
        K->read(psio_, PSIF_DFOCC_INTS);
    }
    else {
        L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|jb)", naoccA, navirA, naoccB, navirB));
        tei_iajb_chem_directAB(L);
        K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ij|Ab>", naoccA, naoccB, navirA, navirB));
        K->sort(1324, L, 1.0, 0.0);
        L.reset();
    }
    t2_1AB = SharedTensor2d(new Tensor2d("T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    t2_1AB->copy(K);
    if (regularization == "FALSE") t2_1AB->apply_denom_os(nfrzc, noccA, noccB, FockA, FockB);
    else if (regularization == "TRUE") t2_1AB->reg_denom_os(nfrzc, noccA, noccB, FockA, FockB, reg_param);
    t2_1AB->write(psio_, PSIF_DFOCC_AMPS);

    // Energy AB
    Emp2AB = t2_1AB->vector_dot(K);
    K.reset();
    Escsmp2AB = os_scale * Emp2AB;
    if (mo_optimized == 0) Esosmp2AB = sos_scale * Emp2AB;
    else if (mo_optimized == 1) Esosmp2AB = sos_scale2 * Emp2AB;

    // sort AB
    t2p_1 = SharedTensor2d(new Tensor2d("T2_1(IA,jb)", naoccA, navirA, naoccB, navirB));
    t2p_1->sort(1324, t2_1AB, 1.0, 0.0);
    t2_1AB.reset();
    t2p_1->write(psio_, PSIF_DFOCC_AMPS);
    t2p_1.reset();

    Ecorr = Emp2AA + Emp2AB + Emp2BB;
    if (reference == "ROHF") Ecorr += Emp2_t1;
    Emp2 = Eref + Ecorr;
    Escsmp2 = Eref + Escsmp2AA + Escsmp2AB + Escsmp2BB;
    Esosmp2 = Eref + Esosmp2AB;
    Escsnmp2 = Eref + Escsnmp2AA + Escsnmp2BB;
}// else if (reference_ == "UNRESTRICTED")
    timer_off("1st-order T2");
} // end t2_1st_sc

}} // End Namespaces


