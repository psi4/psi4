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
  
void DFOCC::ccsd_opdm()
{   

    SharedTensor2d T, U;
    timer_on("opdm");
//if (reference_ == "RESTRICTED") {

    // G1_ij = -(G_ij + G_ji)
    T = SharedTensor2d(new Tensor2d("G Intermediate <I|J>", naoccA, naoccA));
    T->symmetrize(GtijA);
    T->scale(-2.0);
    G1c_oo->set_act_oo(nfrzc, naoccA, T);
    T.reset();

    //  G1_ab = -(G_ab + G_ba)
    T = SharedTensor2d(new Tensor2d("G Intermediate <A|B>", navirA, navirA));
    T->symmetrize(GtabA);
    T->scale(-2.0);
    G1c_vv->set_act_vv(T);
    T.reset();

    // G1_ia = t_i^a + l_i^a
    T = SharedTensor2d(new Tensor2d("Corr OPDM <I|A>", naoccA, navirA));
    T->copy(t1A);
    T->add(l1A);

    // G1_ia += \sum(me) [U(ia,me) - t_m^a * t_i^e] l_m^e
    U = SharedTensor2d(new Tensor2d("U2 (IA|JB)", naoccA, navirA, naoccA, navirA));
    U->read_symm(psio_, PSIF_DFOCC_AMPS);
    #pragma omp parallel for
    for(int i = 0 ; i < naoccA; ++i){
        for(int j = 0 ; j < naoccA; ++j){
            for(int a = 0 ; a < navirA; ++a){
                int ia = ia_idxAA->get(i,a);
                for(int b = 0 ; b < navirA; ++b){
                    int jb = ia_idxAA->get(j,b);
                    int bj = ai_idxAA->get(b,j);
                    U->subtract(ia, jb, t1A->get(i,b) * t1A->get(j,a) );
                }
            }
        }
    }
    T->gemv(false, U, l1A, 1.0, 1.0);
    U.reset();

    // G1_ia -= \sum(m) t_m^a G_im
    T->gemm(false, false, GijA, t1A, -1.0, 1.0);

    // G1_ia += \sum(e) t_i^e G_ea
    T->gemm(false, false, t1A, GabA, 1.0, 1.0);
    G1c_ov->set_act_ov(nfrzc, T);
    T.reset();

    // Build G1_ai
    G1c_vo = G1c_ov->transpose();

    // Build G1c
    G1c->set_oo(G1c_oo);
    G1c->set_ov(G1c_ov);
    G1c->set_vo(G1c_vo);
    G1c->set_vv(noccA, G1c_vv);

    // Build G1
    G1->copy(G1c);
    for (int i = 0; i < noccA; i++) G1->add(i, i, 2.0); 

  if(print_ > 2) {
    G1->print();
    double trace = G1->trace();
    outfile->Printf("\t trace: %12.12f \n", trace);
    
  }

//}// end if (reference_ == "RESTRICTED")

//else if (reference_ == "UNRESTRICTED") {
//}// else if (reference_ == "UNRESTRICTED")
    timer_off("opdm");
} // end ccsd_opdm

}} // End Namespaces


