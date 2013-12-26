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

    SharedTensor2d K, L, M;
    timer_on("Fock");
if (reference_ == "RESTRICTED") {
    // OEI contribution
    FockA->copy(HmoA);

    // F_IJ
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|OO)", noccA, noccA, noccA, noccA));
    if (conv_tei_type == "DISK") K->read(psio_, PSIF_DFOCC_INTS);
    else tei_oooo_chem_ref_directAA(K);
    #pragma omp parallel for
    for(int i = 0 ; i < noccA; ++i){
        for(int j = 0 ; j < noccA; ++j){
            int ij = oo_idxAA->get(i,j);
            double sum = 0.0 ;
            for(int m = 0 ; m < noccA; ++m){
                int im = oo_idxAA->get(i,m);
                int jm = oo_idxAA->get(j,m);
                int mm = oo_idxAA->get(m,m);
                sum += (2.0*K->get(ij,mm)) - K->get(im,jm); 
            }
            FockA->add(i,j,sum);
        }
    }
    K.reset();

    // F_IA
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|OV)", noccA, noccA, noccA, nvirA));
    if (conv_tei_type == "DISK") K->read(psio_, PSIF_DFOCC_INTS);
    else tei_ooov_chem_ref_directAA(K);
    #pragma omp parallel for
    for(int i = 0 ; i < noccA; ++i){
        for(int a = 0 ; a < nvirA; ++a){
            int ia = ov_idxAA->get(i,a);
            double sum = 0.0 ;
            for(int m = 0 ; m < noccA; ++m){
                int im = oo_idxAA->get(i,m);
                int ma = ov_idxAA->get(m,a);
                int mm = oo_idxAA->get(m,m);
                sum += (2.0*K->get(mm,ia)) - K->get(im,ma); 
            }
            FockA->add(i, a + noccA, sum);
            FockA->add(a + noccA, i,sum);
        }
    }
    K.reset();

    // F_AB
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|VV)", noccA, noccA, nvirA, nvirA));
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|OV)", noccA, nvirA, noccA, nvirA));
    if (conv_tei_type == "DISK") {
        K->read(psio_, PSIF_DFOCC_INTS);
        L->read(psio_, PSIF_DFOCC_INTS);
    }
    else {
       tei_oovv_chem_ref_directAA(K);
       tei_ovov_chem_ref_directAA(L);
    }
    #pragma omp parallel for
    for(int a = 0 ; a < nvirA; ++a){
        for(int b = 0 ; b < nvirA; ++b){
            int ab = vv_idxAA->get(a,b);
            double sum = 0.0 ;
            for(int m = 0 ; m < noccA; ++m){
                int ma = ov_idxAA->get(m,a);
                int mb = ov_idxAA->get(m,b);
                int mm = oo_idxAA->get(m,m);
                sum += (2.0*K->get(mm,ab)) - L->get(ma,mb); 
            }
            FockA->add(a + noccA, b + noccA, sum);
        }
    }
    K.reset();
    L.reset();

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

    // F_IJ += \sum_{M} <IM||JM> + \sum_{m} <Im|Jm> = \sum_{M} (IJ|MM) - (IM|JM) + \sum_{m} (IJ|mm)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|OO)", noccA, noccA, noccA, noccA));
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|oo)", noccA, noccA, noccB, noccB));
    if (conv_tei_type == "DISK") {
        K->read(psio_, PSIF_DFOCC_INTS);
        L->read(psio_, PSIF_DFOCC_INTS);
    }
    else {
       tei_oooo_chem_ref_directAA(K);
       tei_oooo_chem_ref_directAB(L);
    }
    #pragma omp parallel for
    for(int i = 0 ; i < noccA; ++i){
        for(int j = 0 ; j < noccA; ++j){
            int ij = oo_idxAA->get(i,j);
            double sum1 = 0.0 ;
            for(int m = 0 ; m < noccA; ++m){
                int im = oo_idxAA->get(i,m);
                int jm = oo_idxAA->get(j,m);
                int mm = oo_idxAA->get(m,m);
                sum1 += K->get(ij,mm) - K->get(im,jm); 
            }
            double sum2 = 0.0 ;
            for(int m = 0 ; m < noccB; ++m){
                int mm = oo_idxBB->get(m,m);
                sum2 += L->get(ij,mm); 
            }
            FockA->add(i,j,sum1+sum2);
        }
    }
    K.reset();

    // F_ij = \sum_{m} <im||jm> + \sum_{M} <Mi|Mj> = \sum_{m} (ij|mm) - (im|jm) + \sum_{M} (MM|ij)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (oo|oo)", noccB, noccB, noccB, noccB));
    if (conv_tei_type == "DISK") K->read(psio_, PSIF_DFOCC_INTS);
    else tei_oooo_chem_ref_directBB(K);
    #pragma omp parallel for
    for(int i = 0 ; i < noccB; ++i){
        for(int j = 0 ; j < noccB; ++j){
            int ij = oo_idxBB->get(i,j);
            double sum1 = 0.0 ;
            for(int m = 0 ; m < noccB; ++m){
                int im = oo_idxBB->get(i,m);
                int jm = oo_idxBB->get(j,m);
                int mm = oo_idxBB->get(m,m);
                sum1 += K->get(ij,mm) - K->get(im,jm); 
            }
            double sum2 = 0.0 ;
            for(int m = 0 ; m < noccA; ++m){
                int mm = oo_idxAA->get(m,m);
                sum2 += L->get(mm,ij); 
            }
            FockB->add(i,j,sum1+sum2);
        }
    }
    K.reset();
    L.reset();

    // F_IA = \sum_{M} <MI||MA> + \sum_{m} <Im|Am> = \sum_{M} (MM|IA) - (IM|MA) + \sum_{m} (IA|mm)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|OV)", noccA, noccA, noccA, nvirA));
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|oo)", noccA, nvirA, noccB, noccB));
    if (conv_tei_type == "DISK") {
        K->read(psio_, PSIF_DFOCC_INTS);
        L->read(psio_, PSIF_DFOCC_INTS);
    }
    else {
       tei_ooov_chem_ref_directAA(K);
       tei_ovoo_chem_ref_directAB(L);
    }
    #pragma omp parallel for
    for(int i = 0 ; i < noccA; ++i){
        for(int a = 0 ; a < nvirA; ++a){
            int ia = ov_idxAA->get(i,a);
            double sum1 = 0.0 ;
            for(int m = 0 ; m < noccA; ++m){
                int mi = oo_idxAA->get(m,i);
                int ma = ov_idxAA->get(m,a);
                int mm = oo_idxAA->get(m,m);
                sum1 += K->get(mm,ia) - K->get(mi,ma); 
            }
            double sum2 = 0.0 ;
            for(int m = 0 ; m < noccB; ++m){
                int mm = oo_idxBB->get(m,m);
                sum2 += L->get(ia,mm); 
            }
            FockA->add(i, a + noccA, sum1+sum2);
            FockA->add(a + noccA, i,sum1+sum2);
        }
    }
    K.reset();
    L.reset();

    // F_ia = \sum_{m} <mi||ma> + \sum_{M} <Mi|Ma> = \sum_{m} (mm|ia) - (im|ma) + \sum_{M} (MM|ia) 
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (oo|ov)", noccB, noccB, noccB, nvirB));
    L = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|ov)", noccA, noccA, noccB, nvirB));
    if (conv_tei_type == "DISK") {
        K->read(psio_, PSIF_DFOCC_INTS);
        L->read(psio_, PSIF_DFOCC_INTS);
    }
    else {
       tei_ooov_chem_ref_directBB(K);
       tei_ooov_chem_ref_directAB(L);
    }
    #pragma omp parallel for
    for(int i = 0 ; i < noccB; ++i){
        for(int a = 0 ; a < nvirB; ++a){
            int ia = ov_idxBB->get(i,a);
            double sum1 = 0.0 ;
            for(int m = 0 ; m < noccB; ++m){
                int mi = oo_idxBB->get(m,i);
                int ma = ov_idxBB->get(m,a);
                int mm = oo_idxBB->get(m,m);
                sum1 += K->get(mm,ia) - K->get(mi,ma); 
            }
            double sum2 = 0.0 ;
            for(int m = 0 ; m < noccA; ++m){
                int mm = oo_idxAA->get(m,m);
                sum2 += L->get(mm,ia); 
            }
            FockB->add(i, a + noccB, sum1+sum2);
            FockB->add(a + noccB, i, sum1+sum2);
        }
    }
    K.reset();
    L.reset();

    // F_AB = \sum_{M} <MA||MB> + \sum_{m} <Am|Bm> = \sum_{M} (MM|AB) - (MA|MB) + \sum_{m} (AB|mm)
    // F_AB += \sum_{M} (MM|AB)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|VV)", noccA, noccA, nvirA, nvirA));
    if (conv_tei_type == "DISK") K->read(psio_, PSIF_DFOCC_INTS);
    else tei_oovv_chem_ref_directAA(K);
    #pragma omp parallel for
    for(int a = 0 ; a < nvirA; ++a){
        for(int b = 0 ; b < nvirA; ++b){
            int ab = vv_idxAA->get(a,b);
            double sum1 = 0.0 ;
            for(int m = 0 ; m < noccA; ++m){
                int mm = oo_idxAA->get(m,m);
                sum1 += K->get(mm,ab); 
            }
            FockA->add(a + noccA, b + noccA, sum1);
        }
    }
    K.reset();

    // F_AB -= \sum_{M} (MA|MB)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OV|OV)", noccA, nvirA, noccA, nvirA));
    if (conv_tei_type == "DISK") K->read(psio_, PSIF_DFOCC_INTS);
    else tei_ovov_chem_ref_directAA(K);
    #pragma omp parallel for
    for(int a = 0 ; a < nvirA; ++a){
        for(int b = 0 ; b < nvirA; ++b){
            double sum1 = 0.0 ;
            for(int m = 0 ; m < noccA; ++m){
                int ma = ov_idxAA->get(m,a);
                int mb = ov_idxAA->get(m,b);
                sum1 -= K->get(ma,mb); 
            }
            FockA->add(a + noccA, b + noccA, sum1);
        }
    }
    K.reset();

    // F_AB += \sum_{m} (AB|mm)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (VV|oo)", nvirA, nvirA, noccB, noccB));
    if (conv_tei_type == "DISK") K->read(psio_, PSIF_DFOCC_INTS);
    else tei_vvoo_chem_ref_directAB(K);
    #pragma omp parallel for
    for(int a = 0 ; a < nvirA; ++a){
        for(int b = 0 ; b < nvirA; ++b){
            int ab = vv_idxAA->get(a,b);
            double sum1 = 0.0 ;
            for(int m = 0 ; m < noccB; ++m){
                int mm = oo_idxBB->get(m,m);
                sum1 += K->get(ab,mm); 
            }
            FockA->add(a + noccA, b + noccA, sum1);
        }
    }
    K.reset();

    // F_ab = \sum_{m} <ma||mb> + \sum_{M} <Ma|Mb> = \sum_{m} (mm|ab) - (ma|mb) + \sum_{M} (MM|ab)
    // F_ab += \sum_{m} (mm|ab)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (oo|vv)", noccB, noccB, nvirB, nvirB));
    if (conv_tei_type == "DISK") K->read(psio_, PSIF_DFOCC_INTS);
    else tei_oovv_chem_ref_directBB(K);
    #pragma omp parallel for
    for(int a = 0 ; a < nvirB; ++a){
        for(int b = 0 ; b < nvirB; ++b){
            int ab = vv_idxBB->get(a,b);
            double sum1 = 0.0 ;
            for(int m = 0 ; m < noccB; ++m){
                int mm = oo_idxBB->get(m,m);
                sum1 += K->get(mm,ab); 
            }
            FockB->add(a + noccB, b + noccB, sum1);
        }
    }
    K.reset();

    // F_ab -= \sum_{m} (ma|mb)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (ov|ov)", noccB, nvirB, noccB, nvirB));
    if (conv_tei_type == "DISK") K->read(psio_, PSIF_DFOCC_INTS);
    else tei_ovov_chem_ref_directBB(K);
    #pragma omp parallel for
    for(int a = 0 ; a < nvirB; ++a){
        for(int b = 0 ; b < nvirB; ++b){
            double sum1 = 0.0 ;
            for(int m = 0 ; m < noccB; ++m){
                int ma = ov_idxBB->get(m,a);
                int mb = ov_idxBB->get(m,b);
                sum1 -= K->get(ma,mb); 
            }
            FockB->add(a + noccB, b + noccB, sum1);
        }
    }
    K.reset();

    // F_ab += \sum_{M} (MM|ab)
    K = SharedTensor2d(new Tensor2d("DF_BASIS_SCF MO Ints (OO|vv)", noccA, noccA, nvirB, nvirB));
    if (conv_tei_type == "DISK") K->read(psio_, PSIF_DFOCC_INTS);
    else tei_oovv_chem_ref_directAB(K);
    #pragma omp parallel for
    for(int a = 0 ; a < nvirB; ++a){
        for(int b = 0 ; b < nvirB; ++b){
            int ab = vv_idxBB->get(a,b);
            double sum1 = 0.0 ;
            for(int m = 0 ; m < noccA; ++m){
                int mm = oo_idxAA->get(m,m);
                sum1 += K->get(mm,ab); 
            }
            FockB->add(a + noccB, b + noccB, sum1);
        }
    }
    K.reset();

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


