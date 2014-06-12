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
  
void DFOCC::s2_response()
{   

    fprintf(outfile,"\tComputing <S**2>_response...\n");  
    fflush(outfile);
    timer_on("s2_response");

    //=========================
    // Read AO basis SO 
    //=========================
    //Sso = SharedTensor2d(new Tensor2d("SO-basis Overlap Ints", nso_, nso_));
    SharedTensor2d SmoAB = SharedTensor2d(new Tensor2d("MO-basis Alpha-Beta Overlap Ints", nmo_, nmo_));
    SharedTensor2d SmoBA = SharedTensor2d(new Tensor2d("MO-basis Beta-Alpha Overlap Ints", nmo_, nmo_));
    SharedTensor2d temp = SharedTensor2d(new Tensor2d("Temp", nso_, nmo_));

    // AB
    temp->gemm(false, false, Sso, CmoB, 1.0, 0.0);
    SmoAB->gemm(true, false, CmoA, temp, 1.0, 0.0);

    // BA
    temp->zero();
    temp->gemm(false, false, Sso, CmoA, 1.0, 0.0);
    SmoBA->gemm(true, false, CmoB, temp, 1.0, 0.0);
    temp.reset();

    //=========================
    // Comput <S2>_ref 
    //=========================
    // <S2>_ref = 1/4 (Na - Nb)^2 + 1/2 (Na + Nb) - \sum_{Ij} S_Ij * S_jI
    s2_ref = 0.0;
    s2_ref = 0.25 * (noccA - noccB) * (noccA - noccB);
    s2_ref += 0.5 * (noccA + noccB);
    for (int i = 0; i < noccA; i++) {
         for (int j = 0; j < noccB; j++) {
              s2_ref -= SmoAB->get(i,j) * SmoBA->get(j,i);
         }
    }

    //=========================
    // Comput <S2>_resp
    //=========================
    // <S2> = <S2>_ref - 1/2 \sum_{Ij} \sum_{Ab} S_Ib S_jA t_Ij^Ab
    s2_resp = 0.0;
    s2_proj = 0.0;
    s2_resp += s2_ref;
    s2_proj += s2_ref;

    // Form overlap blocks
    // S_Ia
    SharedTensor2d SovAB = SharedTensor2d(new Tensor2d("S <I|a>", naoccA, navirB));
    for (int i = 0; i < naoccA; i++) {
         for (int a = 0; a < navirB; a++) {
              SovAB->set(i, a, SmoAB->get(i + nfrzc, a + noccB));
         }
    }

    // S_iA
    SharedTensor2d SovBA = SharedTensor2d(new Tensor2d("S <i|A>", naoccB, navirA));
    for (int i = 0; i < naoccB; i++) {
         for (int a = 0; a < navirA; a++) {
              SovBA->set(i, a, SmoBA->get(i + nfrzc, a + noccA));
         }
    }

    // Compute amplitude contribution
    SharedTensor2d T = SharedTensor2d(new Tensor2d("T2_1 <Ij|Ab>", naoccA, naoccB, navirA, navirB));
    SharedTensor2d T2 = SharedTensor2d(new Tensor2d("T2_1 (Ib,jA)", naoccA, navirB, naoccB, navirA));
    T->read(psio_, PSIF_DFOCC_AMPS);
    T2->sort(1423, T, 1.0, 0.0);
    T.reset();

    // X_Ib = -1/2 \sum_{jA} T2(Ib,jA) S_jA. NoTe: For projected value there is no 1/2 factor
    SharedTensor2d X = SharedTensor2d(new Tensor2d("X <O|v>", noccA, nvirB));
    //X->gemv(false, T2, SovBA, -0.5, 0.0);
    X->gemv(false, T2, SovBA, -1.0, 0.0);
    T2.reset();

    // <S2> += \sum_{Ib} X_Ib S_Ib
    double value = 0.0;
    value = X->vector_dot(SovAB);
    X.reset();
    s2_proj += value;
    s2_resp += 0.5 * value;

    double value2 = 1.0 + (4.0*s2_resp);
    double approx_s = 0.5 * (sqrt(value2) - 1.0);
    double approx_s2 = sqrt(value2);
 
    // Print
    fprintf(outfile,"\t<S**2>_reference (a.u.): %12.8f\n", s2_ref);
    fprintf(outfile,"\t<S**2>_response (a.u.) : %12.8f\n", s2_resp);
    fprintf(outfile,"\t<S**2>_projected (a.u.): %12.8f\n", s2_proj);
    fprintf(outfile,"\tApproximate spin q.n.  : %12.8f\n", approx_s);
    fprintf(outfile,"\tApproximate spin mult. : %12.8f\n", approx_s2);
    fflush(outfile);

    timer_off("s2_response");
} // s2_response

}} // End Namespaces


