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

void DFOCC::ref_energy()
{
     double Ehf;
     Ehf=0.0;

 if (reference_ == "RESTRICTED") {
      for (int i=0; i < noccA;i++) {
	Ehf+=HmoA->get(i,i) + FockA->get(i,i);
      }
    Eref = Ehf + Enuc;
 }// end rhf

 else if (reference_ == "UNRESTRICTED") {

     // alpha contribution
      for (int i=0; i<noccA;i++) {
           Ehf+=HmoA->get(i,i) + FockA->get(i,i);
      }

    // beta contribution
      for (int i=0; i<noccB;i++) {
           Ehf+=HmoB->get(i,i) + FockB->get(i,i);
      }

    Eref = (0.5 * Ehf) + Enuc;
 }// end uhf

} // end of ref_energy

//=======================================================
//          MP2 Energy
//=======================================================
void DFOCC::mp2_energy()
{
    SharedTensor2d K, L, M;
    timer_on("MP2 Energy");
if (reference_ == "RESTRICTED") {
    Ecorr = 0.0;
    JiajbAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(JiajbAA);
    u2p_1 = SharedTensor2d(new Tensor2d("U2_1 (ia|jb)", naoccA, navirA, naoccA, navirA));
    u2p_1->read_symm(psio_, PSIF_DFOCC_AMPS);
    Ecorr = u2p_1->vector_dot(JiajbAA);
    JiajbAA.reset();
    u2p_1.reset();
    Emp2 = Eref + Ecorr;
}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    Emp2AA = 0.0;
    Emp2BB = 0.0;
    Emp2AB = 0.0;

    // AA part
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(L);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ||AB>", naoccA, naoccA, navirA, navirA));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    t2_1AA = SharedTensor2d(new Tensor2d("T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    t2_1AA->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Emp2AA = 0.25 * t2_1AA->vector_dot(K);
    K.reset();
    t2_1AA.reset();
    Escsmp2AA = ss_scale * Emp2AA;
    Escsnmp2AA = 1.76 * Emp2AA;

    // BB part
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB));
    tei_iajb_chem_directBB(L);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|ab>", naoccB, naoccB, navirB, navirB));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij||ab>", naoccB, naoccB, navirB, navirB));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    t2_1BB = SharedTensor2d(new Tensor2d("T2_1 <ij|ab>", naoccB, naoccB, navirB, navirB));
    t2_1BB->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Emp2BB = 0.25 * t2_1BB->vector_dot(K);
    K.reset();
    t2_1BB.reset();
    Escsmp2BB = ss_scale * Emp2BB;
    Escsnmp2BB = 1.76 * Emp2BB;

    // AB part
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|jb)", naoccA, navirA, naoccB, navirB));
    tei_iajb_chem_directAB(L);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    t2_1AB = SharedTensor2d(new Tensor2d("T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    t2_1AB->read(psio_, PSIF_DFOCC_AMPS);
    Emp2AB = t2_1AB->vector_dot(K);
    K.reset();
    t2_1AB.reset();
    Escsmp2AB = os_scale * Emp2AB;
    if (mo_optimized == 0) Esosmp2AB = sos_scale * Emp2AB;
    else if (mo_optimized == 1) Esosmp2AB = sos_scale2 * Emp2AB;


    //Singles-contribution
    Emp2_t1 = 0.0;
    if (reference == "ROHF" && orb_opt_ == "FALSE") {
        if (wfn_type_ == "DF-OMP2" || wfn_type_ == "CD-OMP2") {
        //Alpha
        Emp2_t1 = 0.0;
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
    timer_off("MP2 Energy");
} // end mp2_energy

//=======================================================
//          SCS-MP2 Energy
//=======================================================
void DFOCC::scs_mp2_energy()
{

    SharedTensor2d K, L, M;
    timer_on("MP2 Energy");
if (reference_ == "RESTRICTED") {
    JiajbAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(JiajbAA);

    // Same spin part
    SharedTensor2d temp = SharedTensor2d(new Tensor2d("U2_1 (ia|jb)", naoccA, navirA, naoccA, navirA));
    temp->read(psio_, PSIF_DFOCC_AMPS);
    Ecorr = 0.5*temp->vector_dot(JiajbAA);
    temp.reset();
    Emp2AA = Ecorr;
    Emp2BB = Emp2AA;
    Escsmp2AA = ss_scale * Emp2AA;
    Escsnmp2AA = 1.76 * Emp2AA;
    Escsmp2BB = ss_scale * Emp2BB;
    Escsnmp2BB = 1.76 * Emp2BB;

    // Opposit spin part
    u2p_1 = SharedTensor2d(new Tensor2d("U2_1 (ia|jb)", naoccA, navirA, naoccA, navirA));
    u2p_1->read(psio_, PSIF_DFOCC_AMPS);
    Ecorr = u2p_1->vector_dot(JiajbAA);
    JiajbAA.reset();
    u2p_1.reset();
    Emp2AB = Ecorr - Emp2AA - Emp2BB;
    Escsmp2AB = os_scale * Emp2AB;
    if (mo_optimized == 0) Esosmp2AB = sos_scale * Emp2AB;
    else if (mo_optimized == 1) Esosmp2AB = sos_scale2 * Emp2AB;
}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // AA part
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    tei_iajb_chem_directAA(L);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ||AB>", naoccA, naoccA, navirA, navirA));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    t2_1AA = SharedTensor2d(new Tensor2d("T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    t2_1AA->read(psio_, PSIF_DFOCC_AMPS);
    Emp2AA = 0.25 * t2_1AA->vector_dot(K);
    K.reset();
    t2_1AA.reset();
    Escsmp2AA = ss_scale * Emp2AA;
    Escsnmp2AA = 1.76 * Emp2AA;

    // BB part
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB));
    tei_iajb_chem_directBB(L);
    M = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij|ab>", naoccB, naoccB, navirB, navirB));
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij||ab>", naoccB, naoccB, navirB, navirB));
    tei_pqrs_anti_symm_direct(K, M);
    M.reset();
    t2_1BB = SharedTensor2d(new Tensor2d("T2_1 <ij|ab>", naoccB, naoccB, navirB, navirB));
    t2_1BB->read(psio_, PSIF_DFOCC_AMPS);
    Emp2BB = 0.25 * t2_1BB->vector_dot(K);
    K.reset();
    t2_1BB.reset();
    Escsmp2BB = ss_scale * Emp2BB;
    Escsnmp2BB = 1.76 * Emp2BB;

    // AB part
    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|jb)", naoccA, navirA, naoccB, navirB));
    tei_iajb_chem_directAB(L);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    t2_1AB = SharedTensor2d(new Tensor2d("T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    t2_1AB->read(psio_, PSIF_DFOCC_AMPS);
    Emp2AB = t2_1AB->vector_dot(K);
    K.reset();
    t2_1AB.reset();
    Escsmp2AB = os_scale * Emp2AB;
    if (mo_optimized == 0) Esosmp2AB = sos_scale * Emp2AB;
    else if (mo_optimized == 1) Esosmp2AB = sos_scale2 * Emp2AB;


    Emp2_t1 = 0.0;
    if (reference == "ROHF" && orb_opt_ == "FALSE" && wfn_type_ == "DF-OMP2") {
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
}// else if (reference_ == "UNRESTRICTED")

    Ecorr = Emp2AA + Emp2AB + Emp2BB + Emp2_t1;
    Emp2 = Eref + Ecorr;
    Escsmp2 = Eref + Escsmp2AA + Escsmp2AB + Escsmp2BB;
    Esosmp2 = Eref + Esosmp2AB;
    Escsnmp2 = Eref + Escsnmp2AA + Escsnmp2BB;
    //outfile->Printf("\tDF-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
    //outfile->Printf("\tDF-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
    //
    timer_off("MP2 Energy");
} // end mp2_energy

//=======================================================
//          CCSD Energy
//=======================================================
void DFOCC::ccsd_energy()
{
} // end ccsd_energy

}} // End Namespaces
