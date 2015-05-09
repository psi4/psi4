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
  
void DFOCC::mp2l_energy()
{   

    timer_on("CCL Energy");

    EcorrL = 0.0;
    double Eoei = 0.0;
    double Eoo = 0.0;
    double Eov = 0.0;
    double Evv = 0.0;
    Emp2L_old = Emp2L;

if (reference_ == "RESTRICTED") {
    // DE = \sum_{p,q} G_pq f_pq
    EcorrL += G1c->vector_dot(FockA);
    Eoei = EcorrL;

    // DE += \sum_{Q} \sum_{i,a} G_ia^Q b_ia^Q
    G2c_ov = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA * nvirA));
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA * nvirA));
    G2c_ov->read(psio_, PSIF_DFOCC_DENS);
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += G2c_ov->vector_dot(bQovA);
    G2c_ov.reset();
    bQovA.reset();
    Eov = EcorrL - Eoei;

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // DE = \sum_{p,q} G_pq f_pq
    EcorrL += G1cA->vector_dot(FockA);
    EcorrL += G1cB->vector_dot(FockB);
    Eoei = EcorrL;

    // DE += \sum_{Q} \sum_{I,A} G_IA^Q b_IA^Q
    G2c_ovA = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA * nvirA));
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA * nvirA));
    G2c_ovA->read(psio_, PSIF_DFOCC_DENS);
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += G2c_ovA->vector_dot(bQovA);
    G2c_ovA.reset();
    bQovA.reset();

    // DE += \sum_{Q} \sum_{i,a} G_ia^Q b_ia^Q
    G2c_ovB = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|ov)", nQ, noccB * nvirB));
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ov)", nQ, noccB * nvirB));
    G2c_ovB->read(psio_, PSIF_DFOCC_DENS);
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += G2c_ovB->vector_dot(bQovB);
    G2c_ovB.reset();
    bQovB.reset();
    Eov = EcorrL - Eoei;

}// else if (reference_ == "UNRESTRICTED")

    Emp2L = Eref + EcorrL;
    DE = Emp2L - Emp2L_old;
    /*
    outfile->Printf("\tDF-MP2L One-Electron Energy (a.u.) : %20.14f\n", Eoei);
    outfile->Printf("\tDF-MP2L OV Energy (a.u.)           : %20.14f\n", Eov);
    outfile->Printf("\tDF-MP2L Correlation Energy (a.u.)  : %20.14f\n", EcorrL);
    outfile->Printf("\tDF-MP2L Total Energy (a.u.)        : %20.14f\n", Emp2L);
    
    */
    timer_off("CCL Energy");
} // end mp2_energy

//=======================================================
//          CCSD Energy
//=======================================================          
void DFOCC::ccsdl_energy()
{   

    timer_on("CCL Energy");

    SharedTensor2d G, K;

    EcorrL = 0.0;
    double Eoei = 0.0;
    double Eoo = 0.0;
    double Eov = 0.0;
    double Evv = 0.0;

if (reference_ == "RESTRICTED") {
    // DE = \sum_{p,q} G_pq f_pq
    EcorrL += G1c->vector_dot(FockA);
    Eoei = EcorrL;

    // DE += 1/2 \sum_{Q} \sum_{i,j} G_ij^Q b_ij^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OO)", nQ, noccA, noccA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA, noccA));
    K->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += 0.5 * G->vector_dot(K);
    G.reset();
    K.reset();
    Eoo = EcorrL - Eoei;

    // DE += \sum_{Q} \sum_{i,a} G_ia^Q b_ia^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += G->vector_dot(K);
    G.reset();
    K.reset();
    Eov = EcorrL - Eoei - Eoo;

    // DE += 1/2 \sum_{Q} \sum_{a,b} G_ab^Q b_ab^Q
    G = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|VV)", nQ, nvirA, nvirA));
    G->read(psio_, PSIF_DFOCC_DENS, true, true);
    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA));
    K->read(psio_, PSIF_DFOCC_INTS, true, true);
    EcorrL += 0.5 * G->vector_dot(K);
    G.reset();
    K.reset();
    Evv = EcorrL - Eoei - Eoo - Eov;

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // DE = \sum_{p,q} G_pq f_pq
    EcorrL += G1cA->vector_dot(FockA);
    EcorrL += G1cB->vector_dot(FockB);
    Eoei = EcorrL;

    // DE += \sum_{Q} \sum_{I,A} G_IA^Q b_IA^Q
    G2c_ovA = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|OV)", nQ, noccA * nvirA));
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA * nvirA));
    G2c_ovA->read(psio_, PSIF_DFOCC_DENS);
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += G2c_ovA->vector_dot(bQovA);
    G2c_ovA.reset();
    bQovA.reset();

    // DE += \sum_{Q} \sum_{i,a} G_ia^Q b_ia^Q
    G2c_ovB = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|ov)", nQ, noccB * nvirB));
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ov)", nQ, noccB * nvirB));
    G2c_ovB->read(psio_, PSIF_DFOCC_DENS);
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    EcorrL += G2c_ovB->vector_dot(bQovB);
    G2c_ovB.reset();
    bQovB.reset();
    Eov = EcorrL - Eoei;

}// else if (reference_ == "UNRESTRICTED")

    EccsdL = Eref + EcorrL;
    outfile->Printf("\tDF-CCSDL One-Electron Energy (a.u.): %20.14f\n", Eoei);
    outfile->Printf("\tDF-CCSDL OO Energy (a.u.)          : %20.14f\n", Eoo);
    outfile->Printf("\tDF-CCSDL OV Energy (a.u.)          : %20.14f\n", Eov);
    outfile->Printf("\tDF-CCSDL VV Energy (a.u.)          : %20.14f\n", Evv);
    outfile->Printf("\tDF-CCSDL Correlation Energy (a.u.) : %20.14f\n", EcorrL);
    outfile->Printf("\tDF-CCSDL Total Energy (a.u.)       : %20.14f\n", EccsdL);
    
    timer_off("CCL Energy");
} // end ccsdl_energy


}} // End Namespaces


