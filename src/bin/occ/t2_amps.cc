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

#include <libqt/qt.h>
#include <libtrans/integraltransform.h>

#include "defines.h"
#include "occwave.h"

using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace occwave{

void OCCWave::t2_amps()
{   
     //fprintf(outfile,"\n t2_amps is starting... \n"); fflush(outfile);

//===========================================================================================
//========================= RHF =============================================================
//===========================================================================================
if (reference_ == "RESTRICTED") {

     dpdbuf4 K, T, Tau, Tnew, D, R, Tp, W, TAA, TAB, TBB, Ttemp;
     dpdfile2 Fo,Fv;
     int nElements;

     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);

    // Build T(IA,JB)    
    // T_IJ^AB(2) = \sum_{M,E} Tau_IM^AE(1) W_MBEJ(1) => T(IA,JB)(2) = \sum_{M,E} Tau'(IA,ME) (ME|JB)
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2 (IA|JB)");
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Tau (OV|OV)");
    dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    dpd_->contract444(&Tp, &K, &T, 0, 0, 1.0, 0.0);
    // T_IJ^AB(2) += \sum_{M,E} Tau_JM^BE(1) W_MAEI(1) => T(IA,JB)(2) = \sum_{M,E} Tau'(JB,ME) (ME|IA)
    dpd_->contract444(&K, &Tp, &T, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&K);
    dpd_->buf4_close(&Tp);
    
    // T_IJ^AB(2) -= \sum_{m,e} T_im^ae(1) W_mbje(1) => T(IA,JB)(2) -= \sum_{m,e} T'(ia,me) <me|jb>
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2 (OV|OV)");
    dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>");
    dpd_->contract444(&Tp, &K, &T, 0, 0, -1.0, 1.0);
    // T_IJ^AB(2) -= \sum_{m,e} T_jm^be(1) W_maie(1) => T(IA,JB)(2) -= \sum_{m,e} T'(jb,me) <me|ia>
    dpd_->contract444(&K, &Tp, &T, 0, 0, -1.0, 1.0);
    dpd_->buf4_close(&K);
    dpd_->buf4_close(&Tp);
    
    // T(IA,JB) => T_IJ^AB(2)
    dpd_->buf4_sort(&T, PSIF_OCC_DPD , prqs, ID("[O,O]"), ID("[V,V]"), "T2 <IJ|AB>");
    dpd_->buf4_close(&T);
    
    
    // Build T(JA,IB)    
    // T_IJ^AB(2) = -\sum_{M,E} T_MJ^AE(1) W_MBIE(1) => T(JA,IB)(2) = -\sum_{M,E} T"(JA,ME) <ME|IB>
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2 (JA|IB)");
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2pp (OV|OV)");
    dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>");
    dpd_->contract444(&Tp, &K, &T, 0, 0, -1.0, 0.0);
    // T_IJ^AB(2) = -\sum_{M,E} T_IM^EB(1) W_MAJE(1) => T(JA,IB)(2) = -\sum_{M,E} T"(ME,IB) <ME|JA>
    dpd_->contract444(&K, &Tp, &T, 1, 1, -1.0, 1.0);
    dpd_->buf4_close(&K);
    dpd_->buf4_close(&Tp);
    
    // T(JA,IB) => T_IJ^AB(2)
    dpd_->buf4_sort(&T, PSIF_OCC_DPD , rpqs, ID("[O,O]"), ID("[V,V]"), "T2 (IJ|AB)");
    dpd_->buf4_close(&T);    
    

    // Build T2AAnew
    // T_IJ^AB = <IJ||AB>
    dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    dpd_->buf4_copy(&K, PSIF_OCC_DPD, "T2new <OO|VV>");
    dpd_->buf4_close(&K);

    // T_IJ^AB += T(IA,JB)
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2new <OO|VV>");
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <IJ|AB>");
    dpd_->buf4_axpy(&Tp, &T, 1.0); // 1.0*Tp + T -> T
    dpd_->buf4_close(&Tp); 
    
    // T_IJ^AB += T(JA,IB)
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 (IJ|AB)");
    dpd_->buf4_axpy(&Tp, &T, 1.0); // 1.0*Tp + T -> T
    dpd_->buf4_close(&Tp);


    // T_IJ^AB(2) += \sum_{M,N} T_MN^AB(1) W_MNIJ(1) = \sum_{M,N} T_MN^AB(1) <MN|IJ>
    dpd_->buf4_init(&TAA, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
    dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "MO Ints <OO|OO>");
    dpd_->contract444(&K, &TAA, &T, 1, 1, 1.0, 1.0);
    dpd_->buf4_close(&K);
    dpd_->buf4_close(&T);
   
    // T_IJ^AB(2) += 1/2 \sum_{E,F} T_IJ^EF(1) <EF||AB> = \sum_{E,F} T_IJ^EF(1) <AB|EF>
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[V,V]"), ID("[O,O]"),
                  ID("[V,V]"), ID("[O,O]"), 0, "Z2 <VV|OO>");
    dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV|VV>");
    dpd_->contract444(&K, &TAA, &T, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&K);
    dpd_->buf4_close(&TAA);
    dpd_->buf4_sort(&T, PSIF_OCC_DPD , rspq, ID("[O,O]"), ID("[V,V]"), "Z2 <OO|VV>");
    dpd_->buf4_close(&T);
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2new <OO|VV>");
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Z2 <OO|VV>");
    dpd_->buf4_axpy(&Tp, &T, 1.0); // 1.0*Tp + T -> T
    dpd_->buf4_close(&T);
    dpd_->buf4_close(&Tp);


    // initalize Tnew and Told
    dpd_->buf4_init(&Tnew, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2new <OO|VV>");
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
    

    // T_Ij^Ab = \sum_{e} T_Ij^Ae * F_be 
    dpd_->file2_init(&Fv, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");  
    dpd_->contract424(&T, &Fv, &Tnew, 3, 1, 0, 1.0, 1.0); 

    // T_Ij^Ab = \sum_{E} T_Ij^Eb * F_AE
    dpd_->contract244(&Fv, &T, &Tnew, 1, 2, 1, 1.0, 1.0); 
    dpd_->file2_close(&Fv);
    
    // T_Ij^Ab = -\sum_{m} T_Im^Ab * F_mj 
    dpd_->file2_init(&Fo, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    dpd_->contract424(&T, &Fo, &Tnew, 1, 0, 1, -1.0, 1.0);

    // T_Ij^Ab = - \sum_{M} T_Mj^Ab * F_MI
    dpd_->contract244(&Fo, &T, &Tnew, 0, 0, 0, -1.0, 1.0);
    dpd_->file2_close(&Fo);
    dpd_->buf4_close(&T);
    
    // T_IJ^AB = T_IJ^AB / D_IJ^AB
    dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "D <OO|VV>");
    dpd_->buf4_dirprd(&D, &Tnew);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&Tnew);
    

    // Check convergence?
    dpd_->buf4_init(&Tnew, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2new <OO|VV>");
    dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "RT2 <OO|VV>");
    dpd_->buf4_init(&R, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "RT2 <OO|VV>");
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
    dpd_->buf4_axpy(&T, &R, -1.0); // -1.0*T + R -> R
    dpd_->buf4_close(&T);
    
    nElements = 0;
    for(int h = 0; h < nirrep_; h++) nElements += R.params->coltot[h] * R.params->rowtot[h];
    rms_t2 = 0.0;
    rms_t2 = dpd_->buf4_dot_self(&R);
    dpd_->buf4_close(&R);
    rms_t2 = sqrt(rms_t2) / nElements;
    
    // Reset
    dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "T2 <OO|VV>");
    if (print_ > 3) dpd_->buf4_print(&Tnew, outfile, 1);
    dpd_->buf4_close(&Tnew);

    // DIIS
    if (orb_opt_ == "FALSE" || mo_optimized == 1) {
    dpd_->buf4_init(&R, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "RT2 <OO|VV>");
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
    t2DiisManager->add_entry(2, &R, &T);
    if (t2DiisManager->subspace_size() >= cc_mindiis_) t2DiisManager->extrapolate(1, &T);
    dpd_->buf4_close(&R);
    dpd_->buf4_close(&T);
    }

     // Build Tau(ij,ab) = 2*T(ij,ab) - T(ji,ab)
     // Build TAA(ij,ab) = T(ij,ab) - T(ji,ab)
     dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
     dpd_->buf4_copy(&T, PSIF_OCC_DPD, "Tau <OO|VV>");
     dpd_->buf4_copy(&T, PSIF_OCC_DPD, "T2AA <OO|VV>");
     dpd_->buf4_sort(&T, PSIF_OCC_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "T2jiab <OO|VV>");
     dpd_->buf4_close(&T);
     dpd_->buf4_init(&Tau, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tau <OO|VV>");
     dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2AA <OO|VV>");
     dpd_->buf4_init(&Ttemp, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2jiab <OO|VV>");
     dpd_->buf4_scm(&Tau, 2.0);
     dpd_->buf4_axpy(&Ttemp, &Tau, -1.0); // -1.0*Ttemp + Tau -> Tau
     dpd_->buf4_axpy(&Ttemp, &Tp, -1.0); // -1.0*Ttemp + Tp -> Tp
     dpd_->buf4_close(&Ttemp);
     dpd_->buf4_close(&Tp);
     dpd_->buf4_close(&Tau);
 

    // Build amplitudes in chemist notation
    // T_IJ^AB => T'(IA,JB), T"(JA,IB)
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
    dpd_->buf4_sort(&T, PSIF_OCC_DPD , prqs, ID("[O,V]"), ID("[O,V]"), "T2 (OV|OV)");
    dpd_->buf4_sort(&T, PSIF_OCC_DPD , qrps, ID("[O,V]"), ID("[O,V]"), "T2pp (OV|OV)");
    dpd_->buf4_close(&T);

    // Tau(IJ,AB) => Tau'(IA,JB), Tau"(JA,IB)
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tau <OO|VV>");
    dpd_->buf4_sort(&T, PSIF_OCC_DPD , prqs, ID("[O,V]"), ID("[O,V]"), "Tau (OV|OV)");
    dpd_->buf4_sort(&T, PSIF_OCC_DPD , qrps, ID("[O,V]"), ID("[O,V]"), "Taupp (OV|OV)");
    dpd_->buf4_close(&T);
 

    psio_->close(PSIF_LIBTRANS_DPD, 1);
    psio_->close(PSIF_OCC_DPD, 1);

}// end if (reference_ == "RESTRICTED") 


//===========================================================================================
//========================= UHF =============================================================
//===========================================================================================
else if (reference_ == "UNRESTRICTED") {
     
/********************************************************************************************/
/************************** Build W intermediates *******************************************/
/********************************************************************************************/     
     timer_on("W int");
     w_int();
     timer_off("W int");
     
     dpdbuf4 K, T, Tnew, D, R, Tp, W, TAA, TAB, TBB;
     dpdfile2 Fo,Fv;
     int nElements;
     
     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
     
/********************************************************************************************/
/************************** Alpha-Alpha spin case *******************************************/
/********************************************************************************************/        
    // Build T(IA,JB)    
    // T_IJ^AB(2) = \sum_{M,E} T_IM^AE(1) W_MBEJ(1) => T(IA,JB)(2) = \sum_{M,E} T(IA,ME) W(ME,JB)
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2 (IA|JB)");
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2 (OV|OV)");
    dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W (OV|OV)");
    dpd_->contract444(&Tp, &W, &T, 0, 0, 1.0, 0.0);
    // T_IJ^AB(2) += \sum_{M,E} T_JM^BE(1) W_MAEI(1) => T(IA,JB)(2) = \sum_{M,E} T(JB,ME) W(ME,IA)
    dpd_->contract444(&W, &Tp, &T, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&W);
    dpd_->buf4_close(&Tp);
    
    // T_IJ^AB(2) += \sum_{m,e} T_Im^Ae(1) W_JeBm(1) => T(IA,JB)(2) += \sum_{m,e} T(IA,me) W(JB,me)
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "T2 (OV|ov)");
    dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "W (OV|ov)");
    dpd_->contract444(&Tp, &W, &T, 0, 0, 1.0, 1.0);
    // T_IJ^AB(2) += \sum_{m,e} T_Jm^Be(1) W_IeAm(1) => T(IA,JB)(2) += \sum_{m,e} T(JB,me) W(IA,me)
    dpd_->contract444(&W, &Tp, &T, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&W);
    dpd_->buf4_close(&Tp);
    
    // T(IA,JB) => T_IJ^AB(2)
    dpd_->buf4_sort(&T, PSIF_OCC_DPD , prqs, ID("[O,O]"), ID("[V,V]"), "T2 <IJ|AB>");
    dpd_->buf4_close(&T);
    
    
    
    // Build T(JA,IB)    
    // T_IJ^AB(2) = -\sum_{M,E} T_JM^AE(1) W_MBEI(1) => T(JA,IB)(2) = -\sum_{M,E} T(JA,ME) W(ME,IB)
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2 (JA|IB)");
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2 (OV|OV)");
    dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W (OV|OV)");
    dpd_->contract444(&Tp, &W, &T, 0, 0, -1.0, 0.0);
    // T_IJ^AB(2) = -\sum_{M,E} T_IM^BE(1) W_MAEJ(1) => T(JA,IB)(2) = -\sum_{M,E} T(IB,ME) W(ME,JA)
    dpd_->contract444(&W, &Tp, &T, 0, 0, -1.0, 1.0);
    dpd_->buf4_close(&W);
    dpd_->buf4_close(&Tp);
    
    // T_IJ^AB(2) = -\sum_{m,e} T_Jm^Ae(1) W_IeBm(1) => T(JA,IB)(2) -= \sum_{m,e} T(JA,me) W(IB,me)
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "T2 (OV|ov)");
    dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "W (OV|ov)");
    dpd_->contract444(&Tp, &W, &T, 0, 0, -1.0, 1.0);
    // T_IJ^AB(2) = -\sum_{m,e} T_Im^Be(1) W_JeAm(1) => T(JA,IB)(2) -= \sum_{m,e} T(IB,me) W(JA,me)
    dpd_->contract444(&W, &Tp, &T, 0, 0, -1.0, 1.0);
    dpd_->buf4_close(&W);
    dpd_->buf4_close(&Tp);
    
    // T(JA,IB) => T_IJ^AB(2)
    dpd_->buf4_sort(&T, PSIF_OCC_DPD , rpqs, ID("[O,O]"), ID("[V,V]"), "T2 (IJ|AB)");
    dpd_->buf4_close(&T);    
    
     
     
    // Build T2AAnew
    // T_IJ^AB = <IJ||AB>
    dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
    dpd_->buf4_copy(&K, PSIF_OCC_DPD, "T2new <OO|VV>");
    dpd_->buf4_close(&K);

    // T_IJ^AB += T(IA,JB)
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2new <OO|VV>");
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <IJ|AB>");
    dpd_->buf4_axpy(&Tp, &T, 1.0); // 1.0*Tp + T -> T
    dpd_->buf4_close(&Tp); 
    
    // T_IJ^AB += T(JA,IB)
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 (IJ|AB)");
    dpd_->buf4_axpy(&Tp, &T, 1.0); // 1.0*Tp + T -> T
    dpd_->buf4_close(&Tp);


    // T_IJ^AB(2) += 1/2 \sum_{M,N} T_MN^AB(1) W_MNIJ(1) = 1/2 \sum_{M,N} T_MN^AB(1) <MN||IJ>
    dpd_->buf4_init(&TAA, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
    //dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[O,O]"),
    //              ID("[O,O]"), ID("[O,O]"), 0, "W_1 <OO|OO>");
    dpd_->buf4_init(&W, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "MO Ints <OO||OO>");
    dpd_->contract444(&W, &TAA, &T, 1, 1, 0.5, 1.0);
    dpd_->buf4_close(&W);
    dpd_->buf4_close(&T);
   

    // T_IJ^AB(2) += 1/2 \sum_{E,F} T_IJ^EF(1) <EF||AB> = \sum_{E,F} T_IJ^EF(1) <AB|EF>
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[V,V]"), ID("[O,O]"),
                  ID("[V,V]"), ID("[O,O]"), 0, "Z2 <VV|OO>");
    dpd_->buf4_init(&W, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV|VV>");
    dpd_->contract444(&W, &TAA, &T, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&W);
    dpd_->buf4_close(&TAA);
    dpd_->buf4_sort(&T, PSIF_OCC_DPD , rspq, ID("[O,O]"), ID("[V,V]"), "Z2 <OO|VV>");
    dpd_->buf4_close(&T);
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2new <OO|VV>");
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Z2 <OO|VV>");
    dpd_->buf4_axpy(&Tp, &T, 1.0); // 1.0*Tp + T -> T
    dpd_->buf4_close(&T);
    dpd_->buf4_close(&Tp);


    // initalize Tnew and Told
    dpd_->buf4_init(&Tnew, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2new <OO|VV>");
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
    
    
    // T_IJ^AB = \sum_{E} T_IJ^AE * F_EB + \sum_{E} T_IJ^EB * F_AE
    dpd_->file2_init(&Fv, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");  
    dpd_->contract424(&T, &Fv, &Tnew, 3, 1, 0, 1.0, 1.0); 
    dpd_->contract244(&Fv, &T, &Tnew, 1, 2, 1, 1.0, 1.0); 
    dpd_->file2_close(&Fv);
    
    // T_IJ^AB = -\sum_{M} T_IM^AB * F_MJ - \sum_{M} T_MJ^AB * F_MI
    dpd_->file2_init(&Fo, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    dpd_->contract424(&T, &Fo, &Tnew, 1, 0, 1, -1.0, 1.0);
    dpd_->contract244(&Fo, &T, &Tnew, 0, 0, 0, -1.0, 1.0);
    dpd_->file2_close(&Fo);
    dpd_->buf4_close(&T);
    
  
    // T_IJ^AB = T_IJ^AB / D_IJ^AB
    dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "D <OO|VV>");
    dpd_->buf4_dirprd(&D, &Tnew);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&Tnew);
    
/********************************************************************************************/
/************************** Beta-Beta spin case *********************************************/
/********************************************************************************************/ 
    // Build T(ia,jb)    
    // T_ij^ab(2) = \sum_{m,e} T_im^ae(1) W_mbej(1) => T(ia,jb)(2) = \sum_{m,e} T(ia,me) W(me,jb)
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "T2 (ia|jb)");
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "T2 (ov|ov)");
    dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "W (ov|ov)");
    dpd_->contract444(&Tp, &W, &T, 0, 0, 1.0, 0.0);
    // T_ij^ab(2) += \sum_{m,e} T_jm^be(1) W_maei(1) => T(ia,jb)(2) = \sum_{m,e} T(jb,me) W(me,ia)
    dpd_->contract444(&W, &Tp, &T, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&W);
    dpd_->buf4_close(&Tp);
    
    // T_ij^ab(2) += \sum_{M,E} T_Mi^Ea(1) W_MbEj(1) => T(ia,jb)(2) += \sum_{M,E} T(ME,ia) W(ME,jb)
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "T2 (OV|ov)");
    dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "W (OV|ov)");
    dpd_->contract444(&Tp, &W, &T, 1, 1, 1.0, 1.0);
    // T_ij^ab(2) += \sum_{M,E} T_Mj^Eb(1) W_MaEi(1) => T(ia,jb)(2) += \sum_{M,E} T(ME,jb) W(ME,ia)
    dpd_->contract444(&W, &Tp, &T, 1, 1, 1.0, 1.0);
    dpd_->buf4_close(&W);
    dpd_->buf4_close(&Tp);
    
    // T(ia,jb) => T_ij^ab(2)
    dpd_->buf4_sort(&T, PSIF_OCC_DPD , prqs, ID("[o,o]"), ID("[v,v]"), "T2 <ij|ab>");
    dpd_->buf4_close(&T);
    
    
    
    // Build T(ja,ib)    
    // T_ij^ab(2) = -\sum_{m,e} T_jm^ae(1) W_mbei(1) => T(ja,ib)(2) = -\sum_{m,e} T(ja,me) W(me,ib)
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "T2 (ja|ib)");
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "T2 (ov|ov)");
    dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "W (ov|ov)");
    dpd_->contract444(&Tp, &W, &T, 0, 0, -1.0, 0.0);
    // T_ij^ab(2) = -\sum_{m,e} T_im^be(1) W_maej(1) => T(ja,ib)(2) = -\sum_{m,e} T(ib,me) W(me,ja)
    dpd_->contract444(&W, &Tp, &T, 0, 0, -1.0, 1.0);
    dpd_->buf4_close(&W);
    dpd_->buf4_close(&Tp);
    
    // T_ij^ab(2) = -\sum_{M,E} T_Mj^Ea(1) W_MbEi(1) => T(ja,ib)(2) -= \sum_{M,E} T(ME,ja) W(ME,ib)
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "T2 (OV|ov)");
    dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "W (OV|ov)");
    dpd_->contract444(&Tp, &W, &T, 1, 1, -1.0, 1.0);
    // T_ij^ab(2) = -\sum_{M,E} T_Mi^Eb(1) W_MaEj(1) => T(ja,ib)(2) -= \sum_{M,E} T(ME,ib) W(ME,ja)
    dpd_->contract444(&W, &Tp, &T, 1, 1, -1.0, 1.0);
    dpd_->buf4_close(&W);
    dpd_->buf4_close(&Tp);
    
    // T(ja,ib) => T_ij^ab(2)
    dpd_->buf4_sort(&T, PSIF_OCC_DPD , rpqs, ID("[o,o]"), ID("[v,v]"), "T2 (ij|ab)");
    dpd_->buf4_close(&T);    
    
     
    // Build T2BBnew
    // T_ij^ab = <ij|ab>
    dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
    dpd_->buf4_copy(&K, PSIF_OCC_DPD, "T2new <oo|vv>");
    dpd_->buf4_close(&K);


    // T_ij^ab += T(ia,jb)
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2new <oo|vv>");
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <ij|ab>");
    dpd_->buf4_axpy(&Tp, &T, 1.0); // 1.0*Tp + T -> T
    dpd_->buf4_close(&Tp); 
    
    // T_ij^ab += T(ja,ib)
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 (ij|ab)");
    dpd_->buf4_axpy(&Tp, &T, 1.0); // 1.0*Tp + T -> T
    dpd_->buf4_close(&Tp); 
 
    
    // T_ij^ab(2) += 1/2 \sum_{m,n} T_mn^ab(1) W_mnij(1)
    dpd_->buf4_init(&TBB, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
    //dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[o,o]"),
    //              ID("[o,o]"), ID("[o,o]"), 0, "W_1 <oo|oo>");
    dpd_->buf4_init(&W, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "MO Ints <oo||oo>");
    dpd_->contract444(&W, &TBB, &T, 1, 1, 0.5, 1.0);
    dpd_->buf4_close(&W);
    dpd_->buf4_close(&T);
   
   
    // T_ij^ab(2) += 1/2 \sum_{e,f} T_ij^ef(1) <ef||ab> = \sum_{e,f} T_ij^ef(1) <ab|ef>
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[v,v]"), ID("[o,o]"),
                  ID("[v,v]"), ID("[o,o]"), 0, "Z2 <vv|oo>");
    dpd_->buf4_init(&W, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "MO Ints <vv|vv>");
    dpd_->contract444(&W, &TBB, &T, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&W);
    dpd_->buf4_close(&TBB);
    dpd_->buf4_sort(&T, PSIF_OCC_DPD , rspq, ID("[o,o]"), ID("[v,v]"), "Z2 <oo|vv>");
    dpd_->buf4_close(&T);
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2new <oo|vv>");
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Z2 <oo|vv>");
    dpd_->buf4_axpy(&Tp, &T, 1.0); // 1.0*Tp + T -> T
    dpd_->buf4_close(&T);
    dpd_->buf4_close(&Tp);


    // initalize Tnew and Told
    dpd_->buf4_init(&Tnew, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2new <oo|vv>");
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
    
    
    // T_ij^ab = \sum_{e} T_ij^ae * F_eb + \sum_{e} T_ij^eb * F_ae
    dpd_->file2_init(&Fv, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "F <v|v>");  
    dpd_->contract424(&T, &Fv, &Tnew, 3, 1, 0, 1.0, 1.0); 
    dpd_->contract244(&Fv, &T, &Tnew, 1, 2, 1, 1.0, 1.0); 
    dpd_->file2_close(&Fv);
    
    // T_ij^ab = -\sum_{m} T_im^ab * F_mj - \sum_{m} T_mj^ab * F_mi
    dpd_->file2_init(&Fo, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "F <o|o>");
    dpd_->contract424(&T, &Fo, &Tnew, 1, 0, 1, -1.0, 1.0);
    dpd_->contract244(&Fo, &T, &Tnew, 0, 0, 0, -1.0, 1.0);
    dpd_->file2_close(&Fo);
    dpd_->buf4_close(&T);
    
  
    // T_ij^ab = T_ij^ab / D_ij^ab
    dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "D <oo|vv>");
    dpd_->buf4_dirprd(&D, &Tnew);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&Tnew);
    
/********************************************************************************************/
/************************** Alpha-Beta spin case ********************************************/
/********************************************************************************************/
    // Build T(IA,jb)    
    // T_Ij^Ab(2) = \sum_{M,E} T_IM^AE(1) W_MbEj(1) => T(IA,jb)(2) = \sum_{M,E} T(IA,ME) W(jb,ME)
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "T2 (IA|jb)");
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2 (OV|OV)");
    dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[O,V]"),
                  ID("[o,v]"), ID("[O,V]"), 0, "W (ov|OV)");
    dpd_->contract444(&Tp, &W, &T, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&Tp);
    dpd_->buf4_close(&W);

    
    // T_Ij^Ab(2) += \sum_{m,e} T_jm^be(1) W_IeAm(1) => T(IA,jb)(2) = \sum_{m,e} W(IA,me) T(jb,me)  
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "T2 (ov|ov)");
    dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "W (OV|ov)");
    dpd_->contract444(&W, &Tp, &T, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&W);
    dpd_->buf4_close(&Tp);
    
    // T_Ij^Ab(2) += \sum_{m,e} T_Im^Ae(1) W_mbej(1) => T(IA,jb)(2) += \sum_{m,e} T(IA,me) W(me,jb)
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "T2 (OV|ov)");
    dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "W (ov|ov)");
    dpd_->contract444(&Tp, &W, &T, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&W);
    
    // T_Ij^Ab(2) += \sum_{M,E} T_Mj^Eb(1) W_MAEI(1) => T(IA,jb)(2) += \sum_{M,E} W(ME,IA) T(ME,jb)  
    dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W (OV|OV)");
    dpd_->contract444(&W, &Tp, &T, 1, 1, 1.0, 1.0);
    dpd_->buf4_close(&W);
    dpd_->buf4_close(&Tp);
    
    // T(IA,jb) => T_Ij^Ab(2)
    dpd_->buf4_sort(&T, PSIF_OCC_DPD , prqs, ID("[O,o]"), ID("[V,v]"), "T2 <Ij|Ab>");
    dpd_->buf4_close(&T);
    
    
    // Build T(jA,Ib)    
    // T_Ij^Ab(2) = \sum_{M,e} T_Mj^Ae(1) W_MbeI(1) => T(jA,Ib)(2) = \sum_{M,e} T(jA,Me) W(Me,Ib)
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,V]"), ID("[O,v]"),
                  ID("[o,V]"), ID("[O,v]"), 0, "T2 (jA|Ib)");
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[o,V]"), ID("[O,v]"),
                  ID("[o,V]"), ID("[O,v]"), 0, "T2 (oV|Ov)");
    dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "W (Ov|Ov)");
    dpd_->contract444(&Tp, &W, &T, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&W);
    
    // T_Ij^Ab(2) = +\sum_{m,E} T_Im^Eb(1) W_mAEj(1) => T(jA,Ib)(2) = +\sum_{m,E} W(mE,jA) T(mE,Ib) 
    dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "W (oV|oV)");
    dpd_->contract444(&W, &Tp, &T, 1, 1, 1.0, 1.0); 
    dpd_->buf4_close(&W);
    dpd_->buf4_close(&Tp);
    
    
    // T(jA,Ib) => T_Ij^Ab(2)
    dpd_->buf4_sort(&T, PSIF_OCC_DPD , rpqs, ID("[O,o]"), ID("[V,v]"), "T2 (Ij|Ab)");
    dpd_->buf4_close(&T);    
    
     
    // Build T2ABnew
    // T_Ij^Ab = <Ij|Ab>
    dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    dpd_->buf4_copy(&K, PSIF_OCC_DPD, "T2new <Oo|Vv>");
    dpd_->buf4_close(&K);


    // T_Ij^Ab += T(IA,jb)
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2new <Oo|Vv>");
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2 <Ij|Ab>");
    dpd_->buf4_axpy(&Tp, &T, 1.0); // 1.0*Tp + T -> T
    dpd_->buf4_close(&Tp); 
    
    // T_Ij^Ab += T(jA,Ib)
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2 (Ij|Ab)");
    dpd_->buf4_axpy(&Tp, &T, 1.0); // 1.0*Tp + T -> T
    dpd_->buf4_close(&Tp); 
    

    // T_Ij^Ab(2) += \sum_{M,n} T_Mn^Ab(1) W_MnIj(1) = \sum_{M,n} W(Mn,Ij) T(Mn,Ab)
    dpd_->buf4_init(&TAB, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
    //dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[O,o]"),
    //              ID("[O,o]"), ID("[O,o]"), 0, "W_1 <Oo|Oo>");
    dpd_->buf4_init(&W, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "MO Ints <Oo|Oo>");
    dpd_->contract444(&W, &TAB, &T, 1, 1, 1.0, 1.0);
    dpd_->buf4_close(&W);
    dpd_->buf4_close(&T);
    
    
    // T_Ij^Ab(2) +=  \sum_{E,f} T_Ij^Ef(1) W_AbEf(1) =  \sum_{E,f} T(Ij,Ef) W(Ab,Ef)
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[V,v]"), ID("[O,o]"),
                  ID("[V,v]"), ID("[O,o]"), 0, "Z2 <Vv|Oo>");
    dpd_->buf4_init(&W, PSIF_LIBTRANS_DPD, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "MO Ints <Vv|Vv>");
    dpd_->contract444(&W, &TAB, &T, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&W);
    dpd_->buf4_close(&TAB);
    dpd_->buf4_sort(&T, PSIF_OCC_DPD , rspq, ID("[O,o]"), ID("[V,v]"), "Z2 <Oo|Vv>");
    dpd_->buf4_close(&T);
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2new <Oo|Vv>");
    dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Z2 <Oo|Vv>");
    dpd_->buf4_axpy(&Tp, &T, 1.0); // 1.0*Tp + T -> T
    dpd_->buf4_close(&T);
    dpd_->buf4_close(&Tp);
 

    // initalize Tnew and Told
    dpd_->buf4_init(&Tnew, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2new <Oo|Vv>");
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
    
    
    // T_Ij^Ab = \sum_{e} T_Ij^Ae * F_be + \sum_{E} T_Ij^Eb * F_AE
    dpd_->file2_init(&Fv, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "F <v|v>");  
    dpd_->contract424(&T, &Fv, &Tnew, 3, 1, 0, 1.0, 1.0); 
    dpd_->file2_close(&Fv);
    dpd_->file2_init(&Fv, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");  
    dpd_->contract244(&Fv, &T, &Tnew, 1, 2, 1, 1.0, 1.0); 
    dpd_->file2_close(&Fv);
    
    // T_Ij^Ab = -\sum_{m} T_Im^Ab * F_mj - \sum_{M} T_Mj^Ab * F_MI
    dpd_->file2_init(&Fo, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "F <o|o>");
    dpd_->contract424(&T, &Fo, &Tnew, 1, 0, 1, -1.0, 1.0);
    dpd_->file2_close(&Fo);
    dpd_->file2_init(&Fo, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    dpd_->contract244(&Fo, &T, &Tnew, 0, 0, 0, -1.0, 1.0);
    dpd_->file2_close(&Fo);
    dpd_->buf4_close(&T);
  
    // T_Ij^Ab = T_Ij^Ab / D_Ij^Ab
    dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "D <Oo|Vv>");
    dpd_->buf4_dirprd(&D, &Tnew);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&Tnew);
    
/********************************************************************************************/
/************************** Compute amplitude residual to Check Convergence *****************/
/********************************************************************************************/    
    // Alpha-Alpha spin case
    dpd_->buf4_init(&Tnew, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2new <OO|VV>");
    dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "RT2 <OO|VV>");
    dpd_->buf4_init(&R, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "RT2 <OO|VV>");
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
    dpd_->buf4_axpy(&T, &R, -1.0); // -1.0*T + R -> R
    dpd_->buf4_close(&T);
    
    nElements = 0;
    for(int h = 0; h < nirrep_; h++) nElements += R.params->coltot[h] * R.params->rowtot[h];
    rms_t2AA = 0.0;
    rms_t2AA = dpd_->buf4_dot_self(&R);
    dpd_->buf4_close(&R);
    rms_t2AA = sqrt(rms_t2AA) / nElements;
    
    // Reset
    dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "T2 <OO|VV>");
    if (print_ > 3) dpd_->buf4_print(&Tnew, outfile, 1);
    dpd_->buf4_close(&Tnew);
    
    
    // Beta-Beta spin case
    dpd_->buf4_init(&Tnew, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2new <oo|vv>");
    dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "RT2 <oo|vv>");
    dpd_->buf4_init(&R, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "RT2 <oo|vv>");
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
    dpd_->buf4_axpy(&T, &R, -1.0); // -1.0*T + R -> R
    dpd_->buf4_close(&T);
    
    nElements = 0;
    for(int h = 0; h < nirrep_; h++) nElements += R.params->coltot[h] * R.params->rowtot[h];
    rms_t2BB = 0.0;
    rms_t2BB = dpd_->buf4_dot_self(&R);
    dpd_->buf4_close(&R);
    rms_t2BB = sqrt(rms_t2BB) / nElements;
    
    
    // Reset
    dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "T2 <oo|vv>");
    if (print_ > 3) dpd_->buf4_print(&Tnew, outfile, 1);
    dpd_->buf4_close(&Tnew);
    
    
    // Alpha-Beta spin case
    dpd_->buf4_init(&Tnew, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2new <Oo|Vv>");
    dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "RT2 <Oo|Vv>");
    dpd_->buf4_init(&R, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "RT2 <Oo|Vv>");
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
    dpd_->buf4_axpy(&T, &R, -1.0); // -1.0*T + R -> R
    dpd_->buf4_close(&T);
    
    nElements = 0;
    for(int h = 0; h < nirrep_; h++) nElements += R.params->coltot[h] * R.params->rowtot[h];
    rms_t2AB = 0.0;
    rms_t2AB = dpd_->buf4_dot_self(&R);
    dpd_->buf4_close(&R);
    rms_t2AB = sqrt(rms_t2AA) / nElements;
    
    // Reset
    dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "T2 <Oo|Vv>");
    if (print_ > 3) dpd_->buf4_print(&Tnew, outfile, 1);
    dpd_->buf4_close(&Tnew);
    
    // DIIS
    if (orb_opt_ == "FALSE" || mo_optimized == 1) {
    dpdbuf4 Raa, Rbb, Rab, Taa, Tbb, Tab;
    dpd_->buf4_init(&Raa, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "RT2 <OO|VV>");
    dpd_->buf4_init(&Taa, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
    dpd_->buf4_init(&Rbb, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "RT2 <oo|vv>");
    dpd_->buf4_init(&Tbb, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
    dpd_->buf4_init(&Rab, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "RT2 <Oo|Vv>");
    dpd_->buf4_init(&Tab, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
    t2DiisManager->add_entry(6, &Raa, &Rbb, &Rab, &Taa, &Tbb, &Tab);
    if (t2DiisManager->subspace_size() >= cc_mindiis_) t2DiisManager->extrapolate(3, &Taa, &Tbb, &Tab);
    dpd_->buf4_close(&Raa);
    dpd_->buf4_close(&Rbb);
    dpd_->buf4_close(&Rab);
    dpd_->buf4_close(&Taa);
    dpd_->buf4_close(&Tbb);
    dpd_->buf4_close(&Tab);
    }

    // Build amplitudes in chemist notation
    // T_IJ^AB => T(IA,JB)
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
    dpd_->buf4_sort(&T, PSIF_OCC_DPD , prqs, ID("[O,V]"), ID("[O,V]"), "T2 (OV|OV)");
    dpd_->buf4_close(&T);
    
    // T_ij^ab => T(ia,jb)
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
    dpd_->buf4_sort(&T, PSIF_OCC_DPD , prqs, ID("[o,v]"), ID("[o,v]"), "T2 (ov|ov)");
    dpd_->buf4_close(&T);
    
    
    // T_Ij^Ab => T(IA,jb), T(jA,Ib)
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
    dpd_->buf4_sort(&T, PSIF_OCC_DPD , prqs, ID("[O,V]"), ID("[o,v]"), "T2 (OV|ov)");
    dpd_->buf4_sort(&T, PSIF_OCC_DPD , qrps, ID("[o,V]"), ID("[O,v]"), "T2 (oV|Ov)");
    dpd_->buf4_close(&T);    
   
    // T(IA,jb) => T(jb,IA)   
    dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "T2 (OV|ov)");
    dpd_->buf4_sort(&T, PSIF_OCC_DPD , rspq, ID("[o,v]"), ID("[O,V]"), "T2 (ov|OV)");
    dpd_->buf4_close(&T);

    // close files
    psio_->close(PSIF_LIBTRANS_DPD, 1);
    psio_->close(PSIF_OCC_DPD, 1);
    
}// end if (reference_ == "UNRESTRICTED") 
 //fprintf(outfile,"\n t2_amps done. \n"); fflush(outfile);

} // end t2_amps
}} // End Namespaces

