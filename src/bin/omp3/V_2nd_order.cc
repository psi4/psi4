/** Standard library includes */
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>

/** Required PSI3 includes */ 
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <libtrans/mospace.h>
#include <libtrans/integraltransform.h>

/** Required libmints includes */
#include <libmints/mints.h>
#include <libmints/factory.h>
#include <libmints/wavefunction.h>

#include "defines.h"
#include "omp3wave.h"

using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace omp3wave{
  
void OMP3Wave::V_2nd_order()
{   
     dpdbuf4 T, L, V, V1, V2;     
     
     psio_->open(PSIF_OMP3_DPD, PSIO_OPEN_OLD); 
     psio_->open(PSIF_OMP3_DENSITY, PSIO_OPEN_OLD);
     
/********************************************************************************************/
/************************** OOOO-block ******************************************************/
/********************************************************************************************/     
    // Build V(IJ,KL)  
    // V_IJKL(2) = 1/2 \sum_{E,F} T_IJ^EF(1) L_EF^KL(1) = 1/2 \sum_{E,F} T_IJ^EF(1) T_KL^EF(1)
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    dpd_buf4_init(&L, PSIF_OMP3_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    dpd_buf4_init(&V, PSIF_OMP3_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "V_2 <OO|OO>");
    dpd_contract444(&T, &L, &V, 0, 0, 0.5, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    dpd_buf4_close(&V);
    
    
    // Build V(ij,kl)
    // V_ijkl(2) = 1/2 \sum_{e,f} T_ij^ef(1) L_ef^kl(1) = 1/2 \sum_{e,f} T_ij^ef(1) T_kl^ef(1)
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
    dpd_buf4_init(&L, PSIF_OMP3_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
    dpd_buf4_init(&V, PSIF_OMP3_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "V_2 <oo|oo>");    
    dpd_contract444(&T, &L, &V, 0, 0, 0.5, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    dpd_buf4_close(&V);
    
    
    // Build V(Ij,Kl)
    // V_IjKl(2) = \sum_{E,f} T_Ij^Ef(1) L_EF^KL(1) = \sum_{E,f} T_Ij^Ef(1) T_Kl^Ef(1)
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
    dpd_buf4_init(&L, PSIF_OMP3_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
    dpd_buf4_init(&V, PSIF_OMP3_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "V_2 <Oo|Oo>");    
    dpd_contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    dpd_buf4_close(&V);
  
    
    
    /*
    //Print 
    if (print_ > 3) {
      dpd_buf4_init(&V, PSIF_OMP3_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "V_2 <OO|OO>");
      dpd_buf4_print(&V, outfile, 1);
      dpd_buf4_close(&V);
    
      dpd_buf4_init(&V, PSIF_OMP3_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "V_2 <oo|oo>");
      dpd_buf4_print(&V, outfile, 1);
      dpd_buf4_close(&V);
    
      dpd_buf4_init(&V, PSIF_OMP3_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "V_2 <Oo|Oo>");
      dpd_buf4_print(&V, outfile, 1);
      dpd_buf4_close(&V); 
    }
    */
    
/********************************************************************************************/
/************************** OVOV-block ******************************************************/
/********************************************************************************************/    
    // Build V(IA,JB)  
    // V_IAJB(2) => V(IB,JA) = 1/2 \sum_{M,E} T_IM^BE(1) T_JM^AE(1) = 1/2 \sum_{M,E} T(IB,ME) T(JA,ME) 
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1 (OV|OV)");
    dpd_buf4_init(&L, PSIF_OMP3_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1 (OV|OV)");
    dpd_buf4_init(&V, PSIF_OMP3_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "V_2 <IB|JA>");
    dpd_contract444(&T, &L, &V, 0, 0, 0.5, 0.0); 
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    
    // V_IAJB(2) => V(IB,JA) += 1/2 \sum_{m,e} T_Im^Be(1) T_Jm^Ae(1) = 1/2 \sum_{m,e} T(IB,me) T(JA,me) 
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "T2_1 (OV|ov)");
    dpd_buf4_init(&L, PSIF_OMP3_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "T2_1 (OV|ov)");
    dpd_contract444(&T, &L, &V, 0, 0, 0.5, 1.0); 
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    
    // V(IB,JA) => V(IA,JB)
    dpd_buf4_sort(&V, PSIF_OMP3_DENSITY , psrq, ID("[O,V]"), ID("[O,V]"), "V_2 <OV|OV>");
    dpd_buf4_close(&V);
    
    
    
    
    // Build V(ia,jb)  
    // V_iajb(2) => V(ib,ja) = 1/2 \sum_{m,e} T_im^be(1) T_jm^ae(1) = 1/2 \sum_{m,e} T(ib,me) T(ja,me) 
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "T2_1 (ov|ov)");
    dpd_buf4_init(&L, PSIF_OMP3_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "T2_1 (ov|ov)");
    dpd_buf4_init(&V, PSIF_OMP3_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "V_2 <ib|ja>");
    dpd_contract444(&T, &L, &V, 0, 0, 0.5, 0.0); 
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    
    // V_iajb(2) => V(ib,ja) += 1/2 \sum_{M,E} T_Mi^Eb(1) T_Mj^Ea(1) = 1/2 \sum_{M,E} T(ME,ib) T(ME,ja) 
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "T2_1 (OV|ov)");
    dpd_buf4_init(&L, PSIF_OMP3_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "T2_1 (OV|ov)");
    dpd_contract444(&T, &L, &V, 1, 1, 0.5, 1.0); 
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    
    // V(ib,ja) => V(ia,jb)
    dpd_buf4_sort(&V, PSIF_OMP3_DENSITY , psrq, ID("[o,v]"), ID("[o,v]"), "V_2 <ov|ov>");
    dpd_buf4_close(&V);



    // Build V(Ia,Jb) & V(iA,jB)    
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[o,V]"), ID("[O,v]"),
                  ID("[o,V]"), ID("[O,v]"), 0, "T2_1 (oV|Ov)");
    dpd_buf4_init(&L, PSIF_OMP3_DPD, 0, ID("[o,V]"), ID("[O,v]"),
                  ID("[o,V]"), ID("[O,v]"), 0, "T2_1 (oV|Ov)");
    dpd_buf4_init(&V1, PSIF_OMP3_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "V_2 <Ib|Ja>");
    dpd_buf4_init(&V2, PSIF_OMP3_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "V_2 <iB|jA>");
    // V_IaJb(2) => V(Ib,Ja) = 1/2 \sum_{m,E} T_Im^Eb(1) T_Jm^Ea(1) = 1/2 \sum_{m,E} T(mE,IB) T(mE,Ja) 
    dpd_contract444(&T, &L, &V1, 1, 1, 0.5, 0.0); 
    // V_iAjB(2) => V(iB,jA) = 1/2 \sum_{M,e} T_Mi^Be(1) T_Mj^Ae(1) = 1/2 \sum_{M,e} T(iB,Me) T(jA,Me) 
    dpd_contract444(&T, &L, &V2, 0, 0, 0.5, 0.0); 
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    
    // V(Ib,Ja) => V(Ia,Jb) & V(iB,jA) => V(iA,jB)
    dpd_buf4_sort(&V1, PSIF_OMP3_DENSITY , psrq, ID("[O,v]"), ID("[O,v]"), "V_2 <Ov|Ov>");
    dpd_buf4_close(&V1);
    dpd_buf4_sort(&V2, PSIF_OMP3_DENSITY , psrq, ID("[o,V]"), ID("[o,V]"), "V_2 <oV|oV>");
    dpd_buf4_close(&V2);
    
    
    
    // Build V(Ia,jB)  
    // V_IajB(2) => V(IB,ja) = 1/2 \sum_{M,E} T_IM^BE(1) T_Mj^Ea(1) = 1/2 \sum_{M,E} T(IB,ME) T(ME,ja) 
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1 (OV|OV)");
    dpd_buf4_init(&L, PSIF_OMP3_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "T2_1 (OV|ov)");
    dpd_buf4_init(&V, PSIF_OMP3_DENSITY, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "V_2 <IB|ja>");
    dpd_contract444(&T, &L, &V, 0, 1, 0.5, 0.0); 
    dpd_buf4_close(&T);
    
    // V_IajB(2) => V(IB,ja) += 1/2 \sum_{m,e} T_Im^Be(1) T_jm^ae(1) = 1/2 \sum_{m,e} T(IB,me) T(ja,me) 
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "T2_1 (ov|ov)");
    dpd_contract444(&L, &T, &V, 0, 0, 0.5, 1.0); 
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    
    // V(IB,ja) => V(Ia,jB)
    dpd_buf4_sort(&V, PSIF_OMP3_DENSITY , psrq, ID("[O,v]"), ID("[o,V]"), "V_2 <Ov|oV>");
    dpd_buf4_close(&V);


    // Note: V(iA,Jb) = V(Jb,iA)  


     psio_->close(PSIF_OMP3_DENSITY, 1);
     psio_->close(PSIF_OMP3_DPD, 1);

} // end V_2nd_order
}} // End Namespaces

