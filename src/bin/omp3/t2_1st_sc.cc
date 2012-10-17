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
  
void OMP3Wave::t2_1st_sc()
{   

//===========================================================================================
//========================= RHF =============================================================
//===========================================================================================
if (reference == "RHF") {
     dpdbuf4 K, T, D, Tau, Ttemp, Tss;
     
     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OMP3_DPD, PSIO_OPEN_OLD);
     
    // T_ij^ab = <ij|ab>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    dpd_buf4_copy(&K, PSIF_OMP3_DPD, "T2_1 <OO|VV>");
    dpd_buf4_close(&K);
    
    // T_ij^ab = T_ij^ab / D_ij^ab
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "D <OO|VV>");
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    dpd_buf4_dirprd(&D, &T);
    dpd_buf4_close(&D);
    
     // Build Tau(ij,ab) = 2*T(ij,ab) - T(ji,ab)
     // Build TAA(ij,ab) = T(ij,ab) - T(ji,ab)
     dpd_buf4_copy(&T, PSIF_OMP3_DPD, "Tau_1 <OO|VV>");
     dpd_buf4_copy(&T, PSIF_OMP3_DPD, "T2_1AA <OO|VV>");
     dpd_buf4_sort(&T, PSIF_OMP3_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "T2_1jiab <OO|VV>");
     dpd_buf4_init(&Tau, PSIF_OMP3_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tau_1 <OO|VV>");
     dpd_buf4_init(&Tss, PSIF_OMP3_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1AA <OO|VV>");
     dpd_buf4_init(&Ttemp, PSIF_OMP3_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1jiab <OO|VV>");
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
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    dpd_buf4_sort(&T, PSIF_OMP3_DPD , prqs, ID("[O,V]"), ID("[O,V]"), "T2_1 (OV|OV)");
    dpd_buf4_sort(&T, PSIF_OMP3_DPD , qrps, ID("[O,V]"), ID("[O,V]"), "T2_1pp (OV|OV)");
    dpd_buf4_close(&T);

    // Tau(IJ,AB) => Tau'(IA,JB), Tau"(JA,IB)
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tau_1 <OO|VV>");
    dpd_buf4_sort(&T, PSIF_OMP3_DPD , prqs, ID("[O,V]"), ID("[O,V]"), "Tau_1 (OV|OV)");
    dpd_buf4_sort(&T, PSIF_OMP3_DPD , qrps, ID("[O,V]"), ID("[O,V]"), "Tau_1pp (OV|OV)");
    dpd_buf4_close(&T);
    
    /*   
    // Tau'(IA,JB) = 2T_ij^ab - T_ji^ab = 2T' - T"
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1 (OV|OV)");
    dpd_buf4_copy(&T, PSIF_OMP3_DPD, "Tau_1 (OV|OV)");
    dpd_buf4_close(&T);
    dpd_buf4_init(&Tau, PSIF_OMP3_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Tau_1 (OV|OV)");
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1pp (OV|OV)");
    dpd_buf4_scm(&Tau, 2.0);
    dpd_buf4_axpy(&T, &Tau, -1.0); // -1.0*T + Tau -> Tau
    dpd_buf4_close(&Tau);
    dpd_buf4_close(&T);
     
    // Tau"(JA,IB) = 2T_ij^ab - T_ji^ab = 2T" - T'
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1pp (OV|OV)");
    dpd_buf4_copy(&T, PSIF_OMP3_DPD, "Tau_1pp (OV|OV)");
    dpd_buf4_close(&T);
    dpd_buf4_init(&Tau, PSIF_OMP3_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Tau_1pp (OV|OV)");
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1 (OV|OV)");
    dpd_buf4_scm(&Tau, 2.0);
    dpd_buf4_axpy(&T, &Tau, -1.0); // -1.0*T + Tau -> Tau
    dpd_buf4_close(&Tau);
    dpd_buf4_close(&T);

    // Tau_tilde'(IA,JB) = 2T_ij^ab - T_ji^ab = 2T'(ia,jb) - T"(ja,ib)
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1pp (OV|OV)");
    dpd_buf4_sort(&T, PSIF_OMP3_DPD , rqps, ID("[O,V]"), ID("[O,V]"), "Tau_tilde_1 (OV|OV)");
    dpd_buf4_close(&T);
    dpd_buf4_init(&Tau, PSIF_OMP3_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Tau_tilde_1 (OV|OV)");
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1 (OV|OV)");
    dpd_buf4_axpy(&T, &Tau, 2.0); // 2.0*T + Tau -> Tau
    dpd_buf4_close(&Tau);
    dpd_buf4_close(&T);

    // Tau_tilde"(JA,IB) = 2T_ij^ab - T_ji^ab = 2T"(ja,ib) - T'(ja,ib)
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1 (OV|OV)");
    dpd_buf4_sort(&T, PSIF_OMP3_DPD , rqps, ID("[O,V]"), ID("[O,V]"), "Tau_tilde_1pp (OV|OV)");
    dpd_buf4_close(&T);
    dpd_buf4_init(&Tau, PSIF_OMP3_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Tau_tilde_1pp (OV|OV)");
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2_1pp (OV|OV)");
    dpd_buf4_axpy(&T, &Tau, 2.0); // 2.0*T + Tau -> Tau
    dpd_buf4_close(&Tau);
    dpd_buf4_close(&T);
    */

     psio_->close(PSIF_LIBTRANS_DPD, 1);
     psio_->close(PSIF_OMP3_DPD, 1);

}// end if (reference == "RHF") 


//===========================================================================================
//========================= UHF =============================================================
//===========================================================================================
else if (reference == "UHF") {


     dpdbuf4 K, T, D;
     
     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OMP3_DPD, PSIO_OPEN_OLD);
     
     // Build T2AA
     // T_IJ^AB = <IJ||AB>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
    dpd_buf4_copy(&K, PSIF_OMP3_DPD, "T2_1 <OO|VV>");
    dpd_buf4_close(&K);
    
    
    // T_IJ^AB = T_IJ^AB / D_IJ^AB
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "D <OO|VV>");
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    dpd_buf4_dirprd(&D, &T);
    dpd_buf4_close(&D);
    if (print_ > 1) dpd_buf4_print(&T, outfile, 1);
    dpd_buf4_close(&T);
    
    
    // Build T2BB
    // T_ij^ab = <ij|ab>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
    dpd_buf4_copy(&K, PSIF_OMP3_DPD, "T2_1 <oo|vv>");
    dpd_buf4_close(&K);
    
    
    // T_ij^ab = T_ij^ab / D_ij^ab
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "D <oo|vv>");
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
    dpd_buf4_dirprd(&D, &T);
    dpd_buf4_close(&D);
    if (print_ > 1) dpd_buf4_print(&T, outfile, 1);
    dpd_buf4_close(&T);
    
    
    // Build T2AB
    // T_Ij^Ab = <Ij|Ab>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    dpd_buf4_copy(&K, PSIF_OMP3_DPD, "T2_1 <Oo|Vv>");
    dpd_buf4_close(&K);
    
    
    // T_Ij^Ab = T_Ij^Ab / D_Ij^Ab
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "D <Oo|Vv>");
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
    dpd_buf4_dirprd(&D, &T);
    dpd_buf4_close(&D);
    if (print_ > 1) dpd_buf4_print(&T, outfile, 1);
    dpd_buf4_close(&T);
    
    
    
    // Build amplitudes in chemist notation
    // T_IJ^AB => T(IA,JB)
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    dpd_buf4_sort(&T, PSIF_OMP3_DPD , prqs, ID("[O,V]"), ID("[O,V]"), "T2_1 (OV|OV)");
    dpd_buf4_close(&T);
    
    // T_ij^ab => T(ia,jb)
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
    dpd_buf4_sort(&T, PSIF_OMP3_DPD , prqs, ID("[o,v]"), ID("[o,v]"), "T2_1 (ov|ov)");
    dpd_buf4_close(&T);
    
    /*
    // T_Ij^Ab => T(IA,jb), T(Ib,jA)
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
    dpd_buf4_sort(&T, PSIF_OMP3_DPD , prqs, ID("[O,V]"), ID("[o,v]"), "T2_1 (OV|ov)");
    dpd_buf4_sort(&T, PSIF_OMP3_DPD , psqr, ID("[O,v]"), ID("[o,V]"), "T2_1 (Ov|oV)");
    dpd_buf4_close(&T);
    */
    
    // T_Ij^Ab => T(IA,jb), T(jA,Ib)
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
    dpd_buf4_sort(&T, PSIF_OMP3_DPD , prqs, ID("[O,V]"), ID("[o,v]"), "T2_1 (OV|ov)");
    dpd_buf4_sort(&T, PSIF_OMP3_DPD , qrps, ID("[o,V]"), ID("[O,v]"), "T2_1 (oV|Ov)");
    dpd_buf4_close(&T);    
   
    // T(IA,jb) => T(jb,IA)   
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "T2_1 (OV|ov)");
    dpd_buf4_sort(&T, PSIF_OMP3_DPD , rspq, ID("[o,v]"), ID("[O,V]"), "T2_1 (ov|OV)");
    dpd_buf4_close(&T);
     
    /* 
    // Build Lambda amplitudes 
    // T_IJ^AB => L_AB^IJ
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    dpd_buf4_sort(&T, PSIF_OMP3_DPD , rspq, ID("[V,V]"), ID("[O,O]"), "L2_1 <VV|OO>");
    dpd_buf4_close(&T);
    
    // T_ij^ab => L_ab^ij
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
    dpd_buf4_sort(&T, PSIF_OMP3_DPD , rspq, ID("[v,v]"), ID("[o,o]"), "L2_1 <vv|oo>");
    dpd_buf4_close(&T);
     
    // T_Ij^Ab => L_Ab^Ij
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
    dpd_buf4_sort(&T, PSIF_OMP3_DPD , rspq, ID("[V,v]"), ID("[O,o]"), "L2_1 <Vv|Oo>");
    dpd_buf4_close(&T);    
    */

    psio_->close(PSIF_LIBTRANS_DPD, 1);
    psio_->close(PSIF_OMP3_DPD, 1);
}// end if (reference == "UHF") 

} // end t2_1st_sc
}} // End Namespaces




