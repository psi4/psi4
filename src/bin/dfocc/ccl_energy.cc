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
    timer_on("I/O");
    G2c_ov->read(psio_, PSIF_DFOCC_DENS);
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
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
    timer_on("I/O");
    G2c_ovA->read(psio_, PSIF_DFOCC_DENS);
    bQovA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    EcorrL += G2c_ovA->vector_dot(bQovA);
    G2c_ovA.reset();
    bQovA.reset();

    // DE += \sum_{Q} \sum_{i,a} G_ia^Q b_ia^Q
    G2c_ovB = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|ov)", nQ, noccB * nvirB));
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ov)", nQ, noccB * nvirB));
    timer_on("I/O");
    G2c_ovB->read(psio_, PSIF_DFOCC_DENS);
    bQovB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    EcorrL += G2c_ovB->vector_dot(bQovB);
    G2c_ovB.reset();
    bQovB.reset();
    Eov = EcorrL - Eoei;

}// else if (reference_ == "UNRESTRICTED")

    Emp2L = Eref + EcorrL;
    DE = Emp2L - Emp2L_old;
    /*
    fprintf(outfile,"\tDF-MP2L One-Electron Energy (a.u.) : %20.14f\n", Eoei);
    fprintf(outfile,"\tDF-MP2L OV Energy (a.u.)           : %20.14f\n", Eov);
    fprintf(outfile,"\tDF-MP2L Correlation Energy (a.u.)  : %20.14f\n", EcorrL);
    fprintf(outfile,"\tDF-MP2L Total Energy (a.u.)        : %20.14f\n", Emp2L);
    fflush(outfile);
    */
    timer_off("CCL Energy");
} // end mp2_energy

//=======================================================
//          MP3 Energy
//=======================================================          

}} // End Namespaces


