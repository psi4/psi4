/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

#include <libtrans/integraltransform.h>

#include "defines.h"
#include "occwave.h"

using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace occwave{
  
void OCCWave::omp3_t2_1st_general()
{   
     //outfile->Printf("\n omp3_t2_1st_general is starting... \n"); 

//===========================================================================================
//========================= RHF =============================================================
//===========================================================================================
if (reference_ == "RESTRICTED") {
     dpdbuf4 K, T, Tnew, D, R, Tau, Ttemp, Tss;
     dpdfile2 Fo,Fv;
     int nElements;
     
     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
     
     // T_ij^ab = <ij|ab>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    global_dpd_->buf4_copy(&K, PSIF_OCC_DPD, "T2_1new <OO|VV>");
    global_dpd_->buf4_close(&K);
    
    
    // initalize Tnew and Told
    global_dpd_->buf4_init(&Tnew, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1new <OO|VV>");
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    
    // T_ij^ab = \sum_{e} T_ij^ae * F_be + \sum_{e} T_ij^eb * F_ae
    global_dpd_->file2_init(&Fv, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");  
    global_dpd_->contract424(&T, &Fv, &Tnew, 3, 1, 0, 1.0, 1.0); 
    global_dpd_->contract244(&Fv, &T, &Tnew, 1, 2, 1, 1.0, 1.0); 
    //dpd_contract244(&Fv, &T, &Tnew, 1, 0, 1, 1.0, 1.0); 
    global_dpd_->file2_close(&Fv);
    
    // T_ij^ab = -\sum_{m} T_im^ab * F_mj - \sum_{m} T_mj^ab * F_mi
    global_dpd_->file2_init(&Fo, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    global_dpd_->contract424(&T, &Fo, &Tnew, 1, 0, 1, -1.0, 1.0);
    global_dpd_->contract244(&Fo, &T, &Tnew, 0, 0, 0, -1.0, 1.0);
    //dpd_contract244(&Fo, &T, &Tnew, 0, 2, 0, -1.0, 1.0);
    global_dpd_->file2_close(&Fo);
    
  
    // T_ij^ab = T_ij^ab / D_ij^ab
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "D <OO|VV>");
    global_dpd_->buf4_init(&Tnew, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1new <OO|VV>");
    global_dpd_->buf4_dirprd(&D, &Tnew);
    global_dpd_->buf4_close(&D);
    
     // Build Tau(ij,ab) = 2*T(ij,ab) - T(ji,ab)
     // Build TAA(ij,ab) = T(ij,ab) - T(ji,ab)
     global_dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "Tau_1 <OO|VV>");
     global_dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "T2_1AA <OO|VV>");
     global_dpd_->buf4_sort(&Tnew, PSIF_OCC_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "T2_1jiab <OO|VV>");
     global_dpd_->buf4_init(&Tau, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tau_1 <OO|VV>");
     global_dpd_->buf4_init(&Tss, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1AA <OO|VV>");
     global_dpd_->buf4_init(&Ttemp, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1jiab <OO|VV>");
     global_dpd_->buf4_scm(&Tau, 2.0);
     global_dpd_->buf4_axpy(&Ttemp, &Tau, -1.0); // -1.0*Ttemp + Tau -> Tau
     global_dpd_->buf4_axpy(&Ttemp, &Tss, -1.0); // -1.0*Ttemp + Tss -> Tss
     global_dpd_->buf4_close(&Ttemp);
     global_dpd_->buf4_close(&Tau);
     global_dpd_->buf4_close(&Tss);
    
    // Compute amplitude residual to Check Convergence
    global_dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "RT2_1 <OO|VV>");
    global_dpd_->buf4_init(&R, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "RT2_1 <OO|VV>");
    global_dpd_->buf4_axpy(&T, &R, -1.0); // -1.0*T + R -> R
    global_dpd_->buf4_close(&T);
    
    nElements = 0;
    for(int h = 0; h < nirrep_; h++) nElements += R.params->coltot[h] * R.params->rowtot[h];
    rms_t2 = 0.0;
    rms_t2 = global_dpd_->buf4_dot_self(&R);
    global_dpd_->buf4_close(&R);
    rms_t2 = sqrt(rms_t2 / nElements);
    
    // Reset
    global_dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "T2_1 <OO|VV>");
    if (print_ > 2) global_dpd_->buf4_print(&Tnew, "outfile", 1);
    global_dpd_->buf4_close(&Tnew);

    // Build amplitudes in chemist notation
    // T_IJ^AB => T'(IA,JB), T"(JA,IB)
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    global_dpd_->buf4_sort(&T, PSIF_OCC_DPD , prqs, ID("[O,V]"), ID("[O,V]"), "T2_1 (OV|OV)");
    global_dpd_->buf4_sort(&T, PSIF_OCC_DPD , qrps, ID("[O,V]"), ID("[O,V]"), "T2_1pp (OV|OV)");
    global_dpd_->buf4_close(&T);
  
    // Tau(IJ,AB) => Tau'(IA,JB), Tau"(JA,IB)
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tau_1 <OO|VV>");
    global_dpd_->buf4_sort(&T, PSIF_OCC_DPD , prqs, ID("[O,V]"), ID("[O,V]"), "Tau_1 (OV|OV)");
    global_dpd_->buf4_sort(&T, PSIF_OCC_DPD , qrps, ID("[O,V]"), ID("[O,V]"), "Tau_1pp (OV|OV)");
    global_dpd_->buf4_close(&T);
    
    /* 
    // Tau'(IA,JB) = 2T_ij^ab - T_ji^ab = 2T' - T"
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1 (OV|OV)");
    dpd_buf4_copy(&T, PSIF_OCC_DPD, "Tau_1 (OV|OV)");
    dpd_buf4_close(&T);
    dpd_buf4_init(&Tau, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Tau_1 (OV|OV)");
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1pp (OV|OV)");
    dpd_buf4_scm(&Tau, 2.0);
    dpd_buf4_axpy(&T, &Tau, -1.0); // -1.0*T + Tau -> Tau
    dpd_buf4_close(&Tau);
    dpd_buf4_close(&T);
     
    // Tau"(JA,IB) = 2T_ij^ab - T_ji^ab = 2T" - T'
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1pp (OV|OV)");
    dpd_buf4_copy(&T, PSIF_OCC_DPD, "Tau_1pp (OV|OV)");
    dpd_buf4_close(&T);
    dpd_buf4_init(&Tau, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Tau_1pp (OV|OV)");
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1 (OV|OV)");
    dpd_buf4_scm(&Tau, 2.0);
    dpd_buf4_axpy(&T, &Tau, -1.0); // -1.0*T + Tau -> Tau
    dpd_buf4_close(&Tau);
    dpd_buf4_close(&T);

    // Tau'(IA,JB) = 2T_ij^ab - T_ji^ab = 2T'(ia,jb) - T"(ja,ib)
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1pp (OV|OV)");
    dpd_buf4_sort(&T, PSIF_OCC_DPD , rqps, ID("[O,V]"), ID("[O,V]"), "Tau_1 (OV|OV)");
    dpd_buf4_close(&T);
    dpd_buf4_init(&Tau, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Tau_1 (OV|OV)");
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1 (OV|OV)");
    dpd_buf4_axpy(&T, &Tau, 2.0); // 2.0*T + Tau -> Tau
    dpd_buf4_close(&Tau);
    dpd_buf4_close(&T);

    // Tau_tilde'(IA,JB) = 2T_ij^ab - T_ji^ab = 2T'(ia,jb) - T"(ja,ib)
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1pp (OV|OV)");
    dpd_buf4_sort(&T, PSIF_OCC_DPD , rqps, ID("[O,V]"), ID("[O,V]"), "Tau_tilde_1 (OV|OV)");
    dpd_buf4_close(&T);
    dpd_buf4_init(&Tau, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Tau_tilde_1 (OV|OV)");
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1 (OV|OV)");
    dpd_buf4_axpy(&T, &Tau, 2.0); // 2.0*T + Tau -> Tau
    dpd_buf4_close(&Tau);
    dpd_buf4_close(&T);

    // Tau_tilde"(JA,IB) = 2T_ij^ab - T_ji^ab = 2T"(ja,ib) - T'(ja,ib)
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1 (OV|OV)");
    dpd_buf4_sort(&T, PSIF_OCC_DPD , rqps, ID("[O,V]"), ID("[O,V]"), "Tau_tilde_1pp (OV|OV)");
    dpd_buf4_close(&T);
    dpd_buf4_init(&Tau, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Tau_tilde_1pp (OV|OV)");
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1pp (OV|OV)");
    dpd_buf4_axpy(&T, &Tau, 2.0); // 2.0*T + Tau -> Tau
    dpd_buf4_close(&Tau);
    dpd_buf4_close(&T);
    */

    psio_->close(PSIF_LIBTRANS_DPD, 1);
    psio_->close(PSIF_OCC_DPD, 1);
 
}// end if (reference_ == "RESTRICTED") 


//===========================================================================================
//========================= UHF =============================================================
//===========================================================================================
else if (reference_ == "UNRESTRICTED") {
     dpdbuf4 K, T, Tnew, D, R;
     dpdfile2 Fo,Fv;
     int nElements;
     
     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
     
     
    // Build new T2AA 
    // T_IJ^AB = <IJ||AB>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
    global_dpd_->buf4_copy(&K, PSIF_OCC_DPD, "T2_1new <OO|VV>");
    global_dpd_->buf4_close(&K);
    
    
    // initalize Tnew and Told
    global_dpd_->buf4_init(&Tnew, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1new <OO|VV>");
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    
    
    // T_IJ^AB = \sum_{E} T_IJ^AE * F_EB + \sum_{E} T_IJ^EB * F_AE
    global_dpd_->file2_init(&Fv, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");  
    global_dpd_->contract424(&T, &Fv, &Tnew, 3, 1, 0, 1.0, 1.0); 
    global_dpd_->contract244(&Fv, &T, &Tnew, 1, 2, 1, 1.0, 1.0); 
    global_dpd_->file2_close(&Fv);
    
    // T_IJ^AB = -\sum_{M} T_IM^AB * F_MJ - \sum_{M} T_MJ^AB * F_MI
    global_dpd_->file2_init(&Fo, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    global_dpd_->contract424(&T, &Fo, &Tnew, 1, 0, 1, -1.0, 1.0);
    global_dpd_->contract244(&Fo, &T, &Tnew, 0, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&Fo);
    global_dpd_->buf4_close(&T);
    
  
    // T_IJ^AB = T_IJ^AB / D_IJ^AB
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "D <OO|VV>");
    global_dpd_->buf4_dirprd(&D, &Tnew);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Tnew);
    
    
    // Build new T2BB 
    // T_ij^ab = <ij||ab>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
    global_dpd_->buf4_copy(&K, PSIF_OCC_DPD, "T2_1new <oo|vv>");
    global_dpd_->buf4_close(&K);
    
    
    // initalize Tnew and Told
    global_dpd_->buf4_init(&Tnew, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1new <oo|vv>");
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
    
    
    // T_ij^ab = \sum_{e} T_ij^ae * F_eb + \sum_{e} T_ij^eb * F_ae
    global_dpd_->file2_init(&Fv, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "F <v|v>");  
    global_dpd_->contract424(&T, &Fv, &Tnew, 3, 1, 0, 1.0, 1.0); 
    global_dpd_->contract244(&Fv, &T, &Tnew, 1, 2, 1, 1.0, 1.0); 
    global_dpd_->file2_close(&Fv);
    
    // T_ij^ab = -\sum_{m} T_im^ab * F_mj - \sum_{m} T_mj^ab * F_mi
    global_dpd_->file2_init(&Fo, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "F <o|o>");
    global_dpd_->contract424(&T, &Fo, &Tnew, 1, 0, 1, -1.0, 1.0);
    global_dpd_->contract244(&Fo, &T, &Tnew, 0, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&Fo);
    global_dpd_->buf4_close(&T);
    
  
    // T_ij^ab = T_ij^ab / D_ij^ab
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "D <oo|vv>");
    global_dpd_->buf4_dirprd(&D, &Tnew);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Tnew);
    
    
    // Build new T2AB
    // T_Ij^Ab = <Ij||Ab>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    global_dpd_->buf4_copy(&K, PSIF_OCC_DPD, "T2_1new <Oo|Vv>");
    global_dpd_->buf4_close(&K);
    
    
    // initalize Tnew and Told
    global_dpd_->buf4_init(&Tnew, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1new <Oo|Vv>");
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
    
    
    // T_Ij^Ab = \sum_{e} T_Ij^Ae * F_be + \sum_{E} T_Ij^Eb * F_AE
    global_dpd_->file2_init(&Fv, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "F <v|v>");  
    global_dpd_->contract424(&T, &Fv, &Tnew, 3, 1, 0, 1.0, 1.0); 
    global_dpd_->file2_close(&Fv);
    global_dpd_->file2_init(&Fv, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");  
    global_dpd_->contract244(&Fv, &T, &Tnew, 1, 2, 1, 1.0, 1.0); 
    global_dpd_->file2_close(&Fv);
    
    // T_Ij^Ab = -\sum_{m} T_Im^Ab * F_mj - \sum_{M} T_Mj^Ab * F_MI
    global_dpd_->file2_init(&Fo, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "F <o|o>");
    global_dpd_->contract424(&T, &Fo, &Tnew, 1, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&Fo);
    global_dpd_->file2_init(&Fo, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    global_dpd_->contract244(&Fo, &T, &Tnew, 0, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&Fo);
    global_dpd_->buf4_close(&T);
    
  
    // T_Ij^Ab = T_Ij^Ab / D_Ij^Ab
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "D <Oo|Vv>");
    global_dpd_->buf4_dirprd(&D, &Tnew);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Tnew);
    
    
    // Compute amplitude residual to Check Convergence
    // Alpha-Alpha spin case
    global_dpd_->buf4_init(&Tnew, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1new <OO|VV>");
    global_dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "RT2_1 <OO|VV>");
    global_dpd_->buf4_init(&R, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "RT2_1 <OO|VV>");
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    global_dpd_->buf4_axpy(&T, &R, -1.0); // -1.0*T + R -> R
    global_dpd_->buf4_close(&T);
    
    nElements = 0;
    for(int h = 0; h < nirrep_; h++) nElements += R.params->coltot[h] * R.params->rowtot[h];
    rms_t2AA = 0.0;
    rms_t2AA = global_dpd_->buf4_dot_self(&R);
    global_dpd_->buf4_close(&R);
    rms_t2AA = sqrt(rms_t2AA) / nElements;
    
    // Reset
    global_dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "T2_1 <OO|VV>");
    if (print_ > 1) global_dpd_->buf4_print(&Tnew, "outfile", 1);
    global_dpd_->buf4_close(&Tnew);
    
    
    // Beta-Beta spin case
    global_dpd_->buf4_init(&Tnew, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1new <oo|vv>");
    global_dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "RT2_1 <oo|vv>");
    global_dpd_->buf4_init(&R, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "RT2_1 <oo|vv>");
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
    global_dpd_->buf4_axpy(&T, &R, -1.0); // -1.0*T + R -> R
    global_dpd_->buf4_close(&T);
    
    nElements = 0;
    for(int h = 0; h < nirrep_; h++) nElements += R.params->coltot[h] * R.params->rowtot[h];
    rms_t2BB = 0.0;
    rms_t2BB = global_dpd_->buf4_dot_self(&R);
    global_dpd_->buf4_close(&R);
    rms_t2BB = sqrt(rms_t2BB) / nElements;
    
    // Reset
    global_dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "T2_1 <oo|vv>");
    if (print_ > 1) global_dpd_->buf4_print(&Tnew, "outfile", 1);
    global_dpd_->buf4_close(&Tnew);
    
    
    // Alpha-Beta spin case
    global_dpd_->buf4_init(&Tnew, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1new <Oo|Vv>");
    global_dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "RT2_1 <Oo|Vv>");
    global_dpd_->buf4_init(&R, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "RT2_1 <Oo|Vv>");
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
    global_dpd_->buf4_axpy(&T, &R, -1.0); // -1.0*T + R -> R
    global_dpd_->buf4_close(&T);
    
    nElements = 0;
    for(int h = 0; h < nirrep_; h++) nElements += R.params->coltot[h] * R.params->rowtot[h];
    rms_t2AB = 0.0;
    rms_t2AB = global_dpd_->buf4_dot_self(&R);
    global_dpd_->buf4_close(&R);
    rms_t2AB = sqrt(rms_t2AA) / nElements;
    
    // Reset
    global_dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "T2_1 <Oo|Vv>");
    if (print_ > 1) global_dpd_->buf4_print(&Tnew, "outfile", 1);
    global_dpd_->buf4_close(&Tnew);
    
    
    // Build amplitudes in chemist notation
    // T_IJ^AB => T(IA,JB)
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    global_dpd_->buf4_sort(&T, PSIF_OCC_DPD , prqs, ID("[O,V]"), ID("[O,V]"), "T2_1 (OV|OV)");
    global_dpd_->buf4_close(&T);
    
    // T_ij^ab => T(ia,jb)
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
    global_dpd_->buf4_sort(&T, PSIF_OCC_DPD , prqs, ID("[o,v]"), ID("[o,v]"), "T2_1 (ov|ov)");
    global_dpd_->buf4_close(&T);
    
    // T_Ij^Ab => T(IA,jb), T(jA,Ib)
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
    global_dpd_->buf4_sort(&T, PSIF_OCC_DPD , prqs, ID("[O,V]"), ID("[o,v]"), "T2_1 (OV|ov)");
    global_dpd_->buf4_sort(&T, PSIF_OCC_DPD , qrps, ID("[o,V]"), ID("[O,v]"), "T2_1 (oV|Ov)");
    global_dpd_->buf4_close(&T);

    // T(IA,jb) => T(jb,IA)   
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "T2_1 (OV|ov)");
    global_dpd_->buf4_sort(&T, PSIF_OCC_DPD , rspq, ID("[o,v]"), ID("[O,V]"), "T2_1 (ov|OV)");
    global_dpd_->buf4_close(&T);

     /* 
    // Build Lambda amplitudes 
    // T_IJ^AB => L_AB^IJ
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    dpd_buf4_sort(&T, PSIF_OCC_DPD , rspq, ID("[V,V]"), ID("[O,O]"), "L2_1 <VV|OO>");
    dpd_buf4_close(&T);
    
    // T_ij^ab => L_ab^ij
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
    dpd_buf4_sort(&T, PSIF_OCC_DPD , rspq, ID("[v,v]"), ID("[o,o]"), "L2_1 <vv|oo>");
    dpd_buf4_close(&T);
     
    // T_Ij^Ab => L_Ab^Ij
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
    dpd_buf4_sort(&T, PSIF_OCC_DPD , rspq, ID("[V,v]"), ID("[O,o]"), "L2_1 <Vv|Oo>");
    dpd_buf4_close(&T);    
    */
    
    // close files
    psio_->close(PSIF_LIBTRANS_DPD, 1);
    psio_->close(PSIF_OCC_DPD, 1);
}// end if (reference_ == "UNRESTRICTED") 
    
 //outfile->Printf("\n omp3_t2_1st_general done. \n"); 

} // end omp3_t2_1st_general
}} // End Namespaces