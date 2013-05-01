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

#include <libtrans/integraltransform.h>

#include "defines.h"
#include "occwave.h"

using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace occwave{
  
void OCCWave::ocepa_t2_1st_sc()
{   

//===========================================================================================
//========================= RHF =============================================================
//===========================================================================================
if (reference_ == "RESTRICTED") {
     dpdbuf4 K, T, D, Tau, Ttemp, Tss;
     
     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
     
    // T_ij^ab = <ij|ab>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    dpd_buf4_copy(&K, PSIF_OCC_DPD, "T2 <OO|VV>");
    dpd_buf4_close(&K);
    
    // T_ij^ab = T_ij^ab / D_ij^ab
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "D <OO|VV>");
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
    dpd_buf4_dirprd(&D, &T);
    dpd_buf4_close(&D);
    
     // Build Tau(ij,ab) = 2*T(ij,ab) - T(ji,ab)
     // Build TAA(ij,ab) = T(ij,ab) - T(ji,ab)
     dpd_buf4_copy(&T, PSIF_OCC_DPD, "Tau <OO|VV>");
     dpd_buf4_copy(&T, PSIF_OCC_DPD, "T2AA <OO|VV>");
     dpd_buf4_sort(&T, PSIF_OCC_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "T2jiab <OO|VV>");
     dpd_buf4_init(&Tau, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tau <OO|VV>");
     dpd_buf4_init(&Tss, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2AA <OO|VV>");
     dpd_buf4_init(&Ttemp, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2jiab <OO|VV>");
     dpd_buf4_scm(&Tau, 2.0);
     dpd_buf4_axpy(&Ttemp, &Tau, -1.0); // -1.0*Ttemp + Tau -> Tau
     dpd_buf4_axpy(&Ttemp, &Tss, -1.0); // -1.0*Ttemp + Tss -> Tss
     dpd_buf4_close(&Ttemp);
     dpd_buf4_close(&Tau);
     dpd_buf4_close(&Tss);
     
     if (print_ > 4) dpd_buf4_print(&T, outfile, 1);
     dpd_buf4_close(&T);


    // Build amplitudes in chemist notation
    // T_IJ^AB => T'(IA,JB), T"(JA,IB)
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
    dpd_buf4_sort(&T, PSIF_OCC_DPD , prqs, ID("[O,V]"), ID("[O,V]"), "T2 (OV|OV)");
    dpd_buf4_sort(&T, PSIF_OCC_DPD , qrps, ID("[O,V]"), ID("[O,V]"), "T2pp (OV|OV)");
    dpd_buf4_close(&T);

    // Tau(IJ,AB) => Tau'(IA,JB), Tau"(JA,IB)
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tau <OO|VV>");
    dpd_buf4_sort(&T, PSIF_OCC_DPD , prqs, ID("[O,V]"), ID("[O,V]"), "Tau (OV|OV)");
    dpd_buf4_sort(&T, PSIF_OCC_DPD , qrps, ID("[O,V]"), ID("[O,V]"), "Taupp (OV|OV)");
    dpd_buf4_close(&T);
    

     psio_->close(PSIF_LIBTRANS_DPD, 1);
     psio_->close(PSIF_OCC_DPD, 1);

}// end if (reference_ == "RESTRICTED") 


//===========================================================================================
//========================= UHF =============================================================
//===========================================================================================
else if (reference_ == "UNRESTRICTED") {


     dpdbuf4 K, T, D;
     
     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
     
     // Build T2AA
     // T_IJ^AB = <IJ||AB>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
    dpd_buf4_copy(&K, PSIF_OCC_DPD, "T2 <OO|VV>");
    dpd_buf4_close(&K);
    
    
    // T_IJ^AB = T_IJ^AB / D_IJ^AB
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "D <OO|VV>");
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
    dpd_buf4_dirprd(&D, &T);
    dpd_buf4_close(&D);
    if (print_ > 1) dpd_buf4_print(&T, outfile, 1);
    dpd_buf4_close(&T);
    
    
    // Build T2BB
    // T_ij^ab = <ij|ab>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
    dpd_buf4_copy(&K, PSIF_OCC_DPD, "T2 <oo|vv>");
    dpd_buf4_close(&K);
    
    
    // T_ij^ab = T_ij^ab / D_ij^ab
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "D <oo|vv>");
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
    dpd_buf4_dirprd(&D, &T);
    dpd_buf4_close(&D);
    if (print_ > 1) dpd_buf4_print(&T, outfile, 1);
    dpd_buf4_close(&T);
    
    
    // Build T2AB
    // T_Ij^Ab = <Ij|Ab>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    dpd_buf4_copy(&K, PSIF_OCC_DPD, "T2 <Oo|Vv>");
    dpd_buf4_close(&K);
    
    
    // T_Ij^Ab = T_Ij^Ab / D_Ij^Ab
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "D <Oo|Vv>");
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
    dpd_buf4_dirprd(&D, &T);
    dpd_buf4_close(&D);
    if (print_ > 1) dpd_buf4_print(&T, outfile, 1);
    dpd_buf4_close(&T);
    
    
    
    // Build amplitudes in chemist notation
    // T_IJ^AB => T(IA,JB)
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
    dpd_buf4_sort(&T, PSIF_OCC_DPD , prqs, ID("[O,V]"), ID("[O,V]"), "T2 (OV|OV)");
    dpd_buf4_close(&T);
    
    // T_ij^ab => T(ia,jb)
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
    dpd_buf4_sort(&T, PSIF_OCC_DPD , prqs, ID("[o,v]"), ID("[o,v]"), "T2 (ov|ov)");
    dpd_buf4_close(&T);
    
    
    // T_Ij^Ab => T(IA,jb), T(jA,Ib)
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
    dpd_buf4_sort(&T, PSIF_OCC_DPD , prqs, ID("[O,V]"), ID("[o,v]"), "T2 (OV|ov)");
    dpd_buf4_sort(&T, PSIF_OCC_DPD , qrps, ID("[o,V]"), ID("[O,v]"), "T2 (oV|Ov)");
    dpd_buf4_close(&T);    
   
    // T(IA,jb) => T(jb,IA)   
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "T2 (OV|ov)");
    dpd_buf4_sort(&T, PSIF_OCC_DPD , rspq, ID("[o,v]"), ID("[O,V]"), "T2 (ov|OV)");
    dpd_buf4_close(&T);

    psio_->close(PSIF_LIBTRANS_DPD, 1);
    psio_->close(PSIF_OCC_DPD, 1);
}// end if (reference_ == "UNRESTRICTED") 

} // end t2_1st_sc
}} // End Namespaces




