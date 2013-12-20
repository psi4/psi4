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

    timer_on("MP2 Energy");
if (reference_ == "RESTRICTED") {
    JiajbAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA));
    timer_on("I/O");
    JiajbAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");

    // Same spin part
    SharedTensor2d temp = SharedTensor2d(new Tensor2d("T2_1(ia,jb) - T2_1(ib,ja)", naoccA, navirA, naoccA, navirA));
    timer_on("I/O");
    temp->read(psio_, PSIF_DFOCC_AMPS);
    timer_off("I/O");
    Ecorr = 0.5*temp->vector_dot(JiajbAA); 
    temp.reset();
    Emp2AA = Ecorr;
    Emp2BB = Emp2AA;
    Escsmp2AA = ss_scale * Emp2AA;
    Escsnmp2AA = 1.76 * Emp2AA;
    Escsmp2BB = ss_scale * Emp2BB;
    Escsnmp2BB = 1.76 * Emp2BB;

    // Opposit spin part
    u2p_1 = SharedTensor2d(new Tensor2d("2*T2_1(ia,jb) - T2_1(ib,ja)", naoccA, navirA, naoccA, navirA));
    timer_on("I/O");
    u2p_1->read(psio_, PSIF_DFOCC_AMPS);
    timer_off("I/O");
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
    t2_1AA = SharedTensor2d(new Tensor2d("T2_1 <IJ|AB>", naoccA, naoccA, navirA, navirA));
    AIijabAA = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <IJ||AB>", naoccA, naoccA, navirA, navirA));
    timer_on("I/O");
    t2_1AA->read(psio_, PSIF_DFOCC_AMPS);
    AIijabAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    Emp2AA = 0.25 * t2_1AA->vector_dot(AIijabAA);
    AIijabAA.reset();
    t2_1AA.reset();
    Escsmp2AA = ss_scale * Emp2AA;
    Escsnmp2AA = 1.76 * Emp2AA;

    // BB part
    t2_1BB = SharedTensor2d(new Tensor2d("T2_1 <ij|ab>", naoccB, naoccB, navirB, navirB));
    AIijabBB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <ij||ab>", naoccB, naoccB, navirB, navirB));
    timer_on("I/O");
    t2_1BB->read(psio_, PSIF_DFOCC_AMPS);
    AIijabBB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    Emp2BB = 0.25 * t2_1BB->vector_dot(AIijabBB);
    AIijabBB.reset();
    t2_1BB.reset();
    Escsmp2BB = ss_scale * Emp2BB;
    Escsnmp2BB = 1.76 * Emp2BB;

    // AB part
    t2_1AB = SharedTensor2d(new Tensor2d("T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    IijabAB = SharedTensor2d(new Tensor2d("DF_BASIS_CC MO Ints <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    timer_on("I/O");
    t2_1AB->read(psio_, PSIF_DFOCC_AMPS);
    IijabAB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    Emp2AB = t2_1AB->vector_dot(IijabAB);
    IijabAB.reset();
    t2_1AB.reset();
    Escsmp2AB = os_scale * Emp2AB;
    if (mo_optimized == 0) Esosmp2AB = sos_scale * Emp2AB;
    else if (mo_optimized == 1) Esosmp2AB = sos_scale2 * Emp2AB;


    if (reference == "ROHF" && orb_opt_ == "FALSE" && wfn_type_ == "DF-OMP2") {
        //Singles-contribution
        Emp2_t1 = 0.0;
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
    //fprintf(outfile,"\tDF-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
    //fprintf(outfile,"\tDF-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
    //fflush(outfile);
    timer_off("MP2 Energy");
} // end mp2_energy

//=======================================================
//          MP3 Energy
//=======================================================          

}} // End Namespaces


