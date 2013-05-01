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
  
void OCCWave::v_2nd_order()
{   

 if (reference_ == "RESTRICTED") {

     dpdbuf4 T, Tau, V;     
     
     psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD); 
     psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);

    // Build V(IJ,KL)  
    // V_IJKL(2) = 2\sum_{E,F} T_IJ^EF(1) (2T_KL^EF(1) - T_LK^EF(1))
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    dpd_buf4_init(&Tau, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tau_1 <OO|VV>");
    dpd_buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "V <OO|OO>");
    dpd_contract444(&T, &Tau, &V, 0, 0, 2.0, 0.0);
    dpd_buf4_close(&V);
    dpd_buf4_close(&T);
    dpd_buf4_close(&Tau);

    //Print 
    if (print_ > 3) {
      dpd_buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "V <OO|OO>");
      dpd_buf4_print(&V, outfile, 1);
      dpd_buf4_close(&V);
    }

    // Build V(IA,JB)  
    // V_IAJB(2) => V(IB,JA) = \sum_{M,E} (2T_MI^BE(1) - T_IM^BE(1)) T_JM^EA(1) 
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1pp (OV|OV)");
    dpd_buf4_init(&Tau, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Tau_1pp (OV|OV)");
    dpd_buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "V <IB|JA>");
    dpd_contract444(&Tau, &T, &V, 0, 0, 1.0, 0.0); 
    dpd_buf4_close(&T);
    dpd_buf4_close(&Tau);

    // V_IAJB(2) => V(IB,JA) += \sum_{M,E} (2T_IM^BE(1) - T_MI^BE(1)) T_JM^AE(1) 
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1 (OV|OV)");
    dpd_buf4_init(&Tau, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Tau_1 (OV|OV)");
    dpd_contract444(&Tau, &T, &V, 0, 0, 1.0, 1.0); 
    dpd_buf4_close(&T);
    dpd_buf4_close(&Tau);

    // V(IB,JA) => V(IA,JB)
    dpd_buf4_sort(&V, PSIF_OCC_DENSITY , psrq, ID("[O,V]"), ID("[O,V]"), "V <OV|OV>");
    dpd_buf4_close(&V);


    // Build V(IA,BJ)  
    // V_IABJ(2) => V(IB,JA) = \sum_{M,E} (2T_IM^BE(1) - T_MI^BE(1)) (2T_JM^AE(1) - T_JM^EA(1)) 
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Tau_1 (OV|OV)");
    dpd_buf4_init(&Tau, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Tau_1 (OV|OV)");
    dpd_buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "V (IB|JA)");
    dpd_contract444(&Tau, &T, &V, 0, 0, 1.0, 0.0); 
    dpd_buf4_close(&T);
    dpd_buf4_close(&Tau);

    // V(IB,JA) => V(IA,BJ)
    dpd_buf4_sort(&V, PSIF_OCC_DENSITY , psqr, ID("[O,V]"), ID("[V,O]"), "V <OV|VO>");
    dpd_buf4_close(&V);

     psio_->close(PSIF_OCC_DENSITY, 1);
     psio_->close(PSIF_OCC_DPD, 1);

 }// end if (reference_ == "RESTRICTED") 

 else if (reference_ == "UNRESTRICTED") {

     dpdbuf4 T, L, V, V1, V2;     
     
     psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD); 
     psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);
     
    // Build V(IJ,KL)  
    // V_IJKL(2) = 1/2 \sum_{E,F} T_IJ^EF(1) L_EF^KL(1) = 1/2 \sum_{E,F} T_IJ^EF(1) T_KL^EF(1)
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    dpd_buf4_init(&L, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    dpd_buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "V <OO|OO>");
    dpd_contract444(&T, &L, &V, 0, 0, 0.5, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    dpd_buf4_close(&V);

    // Build V(ij,kl)
    // V_ijkl(2) = 1/2 \sum_{e,f} T_ij^ef(1) L_ef^kl(1) = 1/2 \sum_{e,f} T_ij^ef(1) T_kl^ef(1)
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
    dpd_buf4_init(&L, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
    dpd_buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "V <oo|oo>");    
    dpd_contract444(&T, &L, &V, 0, 0, 0.5, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    dpd_buf4_close(&V);
    
    
    // Build V(Ij,Kl)
    // V_IjKl(2) = \sum_{E,f} T_Ij^Ef(1) L_EF^KL(1) = \sum_{E,f} T_Ij^Ef(1) T_Kl^Ef(1)
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
    dpd_buf4_init(&L, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
    dpd_buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "V <Oo|Oo>");    
    dpd_contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    dpd_buf4_close(&V);
    
    
    //Print 
    if (print_ > 3) {
      dpd_buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "V <OO|OO>");
      dpd_buf4_print(&V, outfile, 1);
      dpd_buf4_close(&V);
    
      dpd_buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "V <oo|oo>");
      dpd_buf4_print(&V, outfile, 1);
      dpd_buf4_close(&V);
    
      dpd_buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "V <Oo|Oo>");
      dpd_buf4_print(&V, outfile, 1);
      dpd_buf4_close(&V); 
    }
    
    // Build V(IA,JB)  
    // V_IAJB(2) => V(IB,JA) = 1/2 \sum_{M,E} T_IM^BE(1) T_JM^AE(1) = 1/2 \sum_{M,E} T(IB,ME) T(JA,ME) 
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1 (OV|OV)");
    dpd_buf4_init(&L, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1 (OV|OV)");
    dpd_buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "V <IB|JA>");
    dpd_contract444(&T, &L, &V, 0, 0, 0.5, 0.0); 
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
   
    // V_IAJB(2) => V(IB,JA) += 1/2 \sum_{m,e} T_Im^Be(1) T_Jm^Ae(1) = 1/2 \sum_{m,e} T(IB,me) T(JA,me) 
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "T2_1 (OV|ov)");
    dpd_buf4_init(&L, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "T2_1 (OV|ov)");
    dpd_contract444(&T, &L, &V, 0, 0, 0.5, 1.0); 
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);

    // V(IB,JA) => V(IA,JB)
    dpd_buf4_sort(&V, PSIF_OCC_DENSITY , psrq, ID("[O,V]"), ID("[O,V]"), "V <OV|OV>");
    dpd_buf4_close(&V);
    
    
    // Build V(ia,jb)  
    // V_iajb(2) => V(ib,ja) = 1/2 \sum_{m,e} T_im^be(1) T_jm^ae(1) = 1/2 \sum_{m,e} T(ib,me) T(ja,me) 
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "T2_1 (ov|ov)");
    dpd_buf4_init(&L, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "T2_1 (ov|ov)");
    dpd_buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "V <ib|ja>");
    dpd_contract444(&T, &L, &V, 0, 0, 0.5, 0.0); 
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
   
    // V_iajb(2) => V(ib,ja) += 1/2 \sum_{M,E} T_Mi^Eb(1) T_Mj^Ea(1) = 1/2 \sum_{M,E} T(ME,ib) T(ME,ja) 
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "T2_1 (OV|ov)");
    dpd_buf4_init(&L, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "T2_1 (OV|ov)");
    dpd_contract444(&T, &L, &V, 1, 1, 0.5, 1.0); 
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);

    // V(ib,ja) => V(ia,jb)
    dpd_buf4_sort(&V, PSIF_OCC_DENSITY , psrq, ID("[o,v]"), ID("[o,v]"), "V <ov|ov>");
    dpd_buf4_close(&V);



    // Build V(Ia,Jb) & V(iA,jB)    
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,V]"), ID("[O,v]"),
                  ID("[o,V]"), ID("[O,v]"), 0, "T2_1 (oV|Ov)");
    dpd_buf4_init(&L, PSIF_OCC_DPD, 0, ID("[o,V]"), ID("[O,v]"),
                  ID("[o,V]"), ID("[O,v]"), 0, "T2_1 (oV|Ov)");
    dpd_buf4_init(&V1, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "V <Ib|Ja>");
    dpd_buf4_init(&V2, PSIF_OCC_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "V <iB|jA>");
    // V_IaJb(2) => V(Ib,Ja) = 1/2 \sum_{m,E} T_Im^Eb(1) T_Jm^Ea(1) = 1/2 \sum_{m,E} T(mE,IB) T(mE,Ja) 
    dpd_contract444(&T, &L, &V1, 1, 1, 0.5, 0.0); 
    // V_iAjB(2) => V(iB,jA) = 1/2 \sum_{M,e} T_Mi^Be(1) T_Mj^Ae(1) = 1/2 \sum_{M,e} T(iB,Me) T(jA,Me) 
    dpd_contract444(&T, &L, &V2, 0, 0, 0.5, 0.0); 
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    
    // V(Ib,Ja) => V(Ia,Jb) & V(iB,jA) => V(iA,jB)
    dpd_buf4_sort(&V1, PSIF_OCC_DENSITY , psrq, ID("[O,v]"), ID("[O,v]"), "V <Ov|Ov>");
    dpd_buf4_close(&V1);
    dpd_buf4_sort(&V2, PSIF_OCC_DENSITY , psrq, ID("[o,V]"), ID("[o,V]"), "V <oV|oV>");
    dpd_buf4_close(&V2);
    
    
    
    // Build V(Ia,jB)  
    /* old
    // V_IajB(2) => V(IB,ja) = 1/2 \sum_{M,E} T_IM^BE(1) T_Mj^Ea(1) = 1/2 \sum_{M,E} T(IB,ME) T(ME,ja) 
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1 (OV|OV)");
    dpd_buf4_init(&L, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "T2_1 (OV|ov)");
    dpd_buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "V <IB|ja>");
    dpd_contract444(&T, &L, &V, 0, 1, 0.5, 0.0); 
    dpd_buf4_close(&T);
    */
   
    // Build V(Ia,jB)  
    // V_IajB(2) => V(IB,ja) = 1/2 \sum_{M,E} T_IM^BE(1) T_Mj^Ea(1) = 1/2 \sum_{M,E} T(IB,ME) T(ME,ja) 
    //                       = 1/2 \sum_{M,E} T(IB,ME) T(ja,ME)
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1 (OV|OV)");
    dpd_buf4_init(&L, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[O,V]"),
                  ID("[o,v]"), ID("[O,V]"), 0, "T2_1 (ov|OV)");
    dpd_buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "V <IB|ja>");
    dpd_contract444(&T, &L, &V, 0, 0, 0.5, 0.0); 
    dpd_buf4_close(&T);
    //

    // V_IajB(2) => V(IB,ja) += 1/2 \sum_{m,e} T_Im^Be(1) T_jm^ae(1) = 1/2 \sum_{m,e} T(IB,me) T(ja,me) 
    //                        = 1/2 \sum_{m,e} T(me,IB) T(me,ja)
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "T2_1 (ov|ov)");
    dpd_contract444(&L, &T, &V, 1, 1, 0.5, 1.0); 
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    
    // V(IB,ja) => V(Ia,jB)
    dpd_buf4_sort(&V, PSIF_OCC_DENSITY , psrq, ID("[O,v]"), ID("[o,V]"), "V <Ov|oV>");
    dpd_buf4_close(&V);

    // Note: V(iA,Jb) = V(Jb,iA)  

     psio_->close(PSIF_OCC_DENSITY, 1);
     psio_->close(PSIF_OCC_DPD, 1);

 }// end if (reference_ == "UNRESTRICTED") 

} // end V_2nd_order
}} // End Namespaces

