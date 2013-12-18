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
  
void DFOCC::fock()
{   

    timer_on("Fock");
if (reference_ == "RESTRICTED") {
    // OEI contribution
    FockA->copy(HmoA);

    // F_IJ
    JooooAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|OO)", noccA, noccA, noccA, noccA));
    timer_on("I/O");
    JooooAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    for(int i = 0 ; i < noccA; ++i){
        for(int j = 0 ; j < noccA; ++j){
            int ij = oo_idxAA->get(i,j);
            double sum = 0.0 ;
            for(int m = 0 ; m < noccA; ++m){
                int im = oo_idxAA->get(i,m);
                int jm = oo_idxAA->get(j,m);
                int mm = oo_idxAA->get(m,m);
                sum += (2.0*JooooAA->get(ij,mm)) - JooooAA->get(im,jm); 
            }
            FockA->add(i,j,sum);
        }
    }
    JooooAA.reset();

    // F_IA
    JooovAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|OV)", noccA, noccA, noccA, nvirA));
    timer_on("I/O");
    JooovAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    for(int i = 0 ; i < noccA; ++i){
        for(int a = 0 ; a < nvirA; ++a){
            int ia = ov_idxAA->get(i,a);
            double sum = 0.0 ;
            for(int m = 0 ; m < noccA; ++m){
                int im = oo_idxAA->get(i,m);
                int ma = ov_idxAA->get(m,a);
                int mm = oo_idxAA->get(m,m);
                sum += (2.0*JooovAA->get(mm,ia)) - JooovAA->get(im,ma); 
            }
            FockA->add(i, a + noccA, sum);
            FockA->add(a + noccA, i,sum);
        }
    }
    JooovAA.reset();

    // F_AB
    JoovvAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|VV)", noccA, noccA, nvirA, nvirA));
    JovovAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|OV)", noccA, nvirA, noccA, nvirA));
    timer_on("I/O");
    JoovvAA->read(psio_, PSIF_DFOCC_INTS);
    JovovAA->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    for(int a = 0 ; a < nvirA; ++a){
        for(int b = 0 ; b < nvirA; ++b){
            int ab = vv_idxAA->get(a,b);
            double sum = 0.0 ;
            for(int m = 0 ; m < noccA; ++m){
                int ma = ov_idxAA->get(m,a);
                int mb = ov_idxAA->get(m,b);
                int mm = oo_idxAA->get(m,m);
                sum += (2.0*JoovvAA->get(mm,ab)) - JovovAA->get(ma,mb); 
            }
            FockA->add(a + noccA, b + noccA, sum);
        }
    }
    JoovvAA.reset();
    JovovAA.reset();

    // Fock Blocks
    FooA->form_oo(FockA);
    FvoA->form_vo(FockA);
    FvvA->form_vv(noccA, FockA);

    if (print_ > 2) FockA->print();
}// end if (reference_ == "RESTRICTED")



else if (reference_ == "UNRESTRICTED") {
    // OEI contribution
    FockA->copy(HmoA);
    FockB->copy(HmoB);

    // F_IJ = h_IJ + \sum_{M} <IM||JM> + \sum_{m} <Im|Jm>
    AIooooAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OO||OO>", noccA, noccA, noccA, noccA));
    IooooAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <Oo|Oo>", noccA, noccB, noccA, noccB));
    timer_on("I/O");
    AIooooAA->read(psio_, PSIF_DFOCC_INTS);
    IooooAB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    for(int i = 0 ; i < noccA; ++i){
        for(int j = 0 ; j < noccA; ++j){
            double sum1 = 0.0 ;
            for(int m = 0 ; m < noccA; ++m){
                int im = oo_idxAA->get(i,m);
                int jm = oo_idxAA->get(j,m);
                sum1 += AIooooAA->get(im,jm); 
            }
            double sum2 = 0.0 ;
            for(int m = 0 ; m < noccB; ++m){
                int im = oo_idxAB->get(i,m);
                int jm = oo_idxAB->get(j,m);
                sum2 += IooooAB->get(im,jm); 
            }
            FockA->add(i,j,sum1+sum2);
        }
    }
    AIooooAA.reset();
    IooooAB.reset();

    // F_ij = h_ij + \sum_{m} <im||jm> + \sum_{M} <Mi|Mj>
    AIooooBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <oo||oo>", noccB, noccB, noccB, noccB));
    IooooAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <Oo|Oo>", noccA, noccB, noccA, noccB));
    timer_on("I/O");
    AIooooBB->read(psio_, PSIF_DFOCC_INTS);
    IooooAB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    for(int i = 0 ; i < noccB; ++i){
        for(int j = 0 ; j < noccB; ++j){
            double sum1 = 0.0 ;
            for(int m = 0 ; m < noccB; ++m){
                int im = oo_idxBB->get(i,m);
                int jm = oo_idxBB->get(j,m);
                sum1 += AIooooBB->get(im,jm); 
            }
            double sum2 = 0.0 ;
            for(int m = 0 ; m < noccA; ++m){
                int mi = oo_idxAB->get(m,i);
                int mj = oo_idxAB->get(m,j);
                sum2 += IooooAB->get(mi,mj); 
            }
            FockB->add(i,j,sum1+sum2);
        }
    }
    AIooooBB.reset();
    IooooAB.reset();

    // F_IA = h_IA + \sum_{M} <MI||MA> + \sum_{m} <Im|Am>
    AIooovAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OO||OV>", noccA, noccA, noccA, nvirA));
    IoovoAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <Oo|Vo>", noccA, noccB, nvirA, noccB));
    timer_on("I/O");
    AIooovAA->read(psio_, PSIF_DFOCC_INTS);
    IoovoAB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    for(int i = 0 ; i < noccA; ++i){
        for(int a = 0 ; a < nvirA; ++a){
            double sum1 = 0.0 ;
            for(int m = 0 ; m < noccA; ++m){
                int mi = oo_idxAA->get(m,i);
                int ma = ov_idxAA->get(m,a);
                sum1 += AIooovAA->get(mi,ma); 
            }
            double sum2 = 0.0 ;
            for(int m = 0 ; m < noccB; ++m){
                int im = oo_idxAB->get(i,m);
                int am = vo_idxAB->get(a,m);
                sum2 += IoovoAB->get(im,am); 
            }
            FockA->add(i, a + noccA, sum1+sum2);
            FockA->add(a + noccA, i,sum1+sum2);
        }
    }
    AIooovAA.reset();
    IoovoAB.reset();

    // F_ia = h_ia + \sum_{m} <mi||ma> + \sum_{M} <Mi|Ma>
    AIooovBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <oo||ov>", noccB, noccB, noccB, nvirB));
    IooovAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <Oo|Ov>", noccA, noccB, noccA, nvirB));
    timer_on("I/O");
    AIooovBB->read(psio_, PSIF_DFOCC_INTS);
    IooovAB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    for(int i = 0 ; i < noccB; ++i){
        for(int a = 0 ; a < nvirB; ++a){
            double sum1 = 0.0 ;
            for(int m = 0 ; m < noccB; ++m){
                int mi = oo_idxBB->get(m,i);
                int ma = ov_idxBB->get(m,a);
                sum1 += AIooovBB->get(mi,ma); 
            }
            double sum2 = 0.0 ;
            for(int m = 0 ; m < noccA; ++m){
                int mi = oo_idxAB->get(m,i);
                int ma = ov_idxAB->get(m,a);
                sum2 += IooovAB->get(mi,ma); 
            }
            FockB->add(i, a + noccB, sum1+sum2);
            FockB->add(a + noccB, i, sum1+sum2);
        }
    }
    AIooovBB.reset();
    IooovAB.reset();

    // F_AB = h_AB + \sum_{M} <MA||MB> + \sum_{m} <Am|Bm>
    AIovovAA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <OV||OV>", noccA, nvirA, noccA, nvirA));
    IvovoAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <Vo|Vo>", nvirA, noccB, nvirA, noccB));
    timer_on("I/O");
    AIovovAA->read(psio_, PSIF_DFOCC_INTS);
    IvovoAB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    for(int a = 0 ; a < nvirA; ++a){
        for(int b = 0 ; b < nvirA; ++b){
            double sum1 = 0.0 ;
            for(int m = 0 ; m < noccA; ++m){
                int ma = ov_idxAA->get(m,a);
                int mb = ov_idxAA->get(m,b);
                sum1 += AIovovAA->get(ma,mb); 
            }
            double sum2 = 0.0 ;
            for(int m = 0 ; m < noccB; ++m){
                int am = vo_idxAB->get(a,m);
                int bm = vo_idxAB->get(b,m);
                sum2 += IvovoAB->get(am,bm); 
            }
            FockA->add(a + noccA, b + noccA, sum1+sum2);
        }
    }
    AIovovAA.reset();
    IvovoAB.reset();

    // F_ab = h_ab + \sum_{m} <ma||mb> + \sum_{M} <Ma|Mb>
    AIovovBB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <ov||ov>", noccB, nvirB, noccB, nvirB));
    IovovAB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints <Ov|Ov>", noccA, nvirB, noccA, nvirB));
    timer_on("I/O");
    AIovovBB->read(psio_, PSIF_DFOCC_INTS);
    IovovAB->read(psio_, PSIF_DFOCC_INTS);
    timer_off("I/O");
    for(int a = 0 ; a < nvirB; ++a){
        for(int b = 0 ; b < nvirB; ++b){
            double sum1 = 0.0 ;
            for(int m = 0 ; m < noccB; ++m){
                int ma = ov_idxBB->get(m,a);
                int mb = ov_idxBB->get(m,b);
                sum1 += AIovovBB->get(ma,mb); 
            }
            double sum2 = 0.0 ;
            for(int m = 0 ; m < noccA; ++m){
                int ma = ov_idxAB->get(m,a);
                int mb = ov_idxAB->get(m,b);
                sum2 += IovovAB->get(ma,mb); 
            }
            FockB->add(a + noccB, b + noccB, sum1+sum2);
        }
    }
    AIovovBB.reset();
    IvovoAB.reset();

    // Fock Blocks
    FooA->form_oo(FockA);
    FooB->form_oo(FockB);
    FvoA->form_vo(FockA);
    FvoB->form_vo(FockB);
    FvvA->form_vv(noccA, FockA);
    FvvB->form_vv(noccB, FockB);

    if (print_ > 2) {
        FockA->print();
        FockB->print();
    }

}// end else if (reference_ == "UNRESTRICTED") 
    timer_off("Fock");
} // end fock

}} // End Namespaces


