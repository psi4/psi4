#include <libtrans/integraltransform.h>

#include "defines.h"
#include "occwave.h"

using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace occwave{
  
void OCCWave::w_int()
{   

     dpdbuf4 K, W;
     
     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);

//===========================================================================================
//========================= RHF =============================================================
//===========================================================================================
if (reference_ == "RESTRICTED") {
    // W_mbej => W(me,jb) = <mb|ej> = (me|jb)
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    dpd_buf4_copy(&K, PSIF_OCC_DPD, "W (OV|OV)");
    dpd_buf4_close(&K);
    
    // W_mbje => W'(me,jb) = <mb|je> = <me|jb>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>");
    dpd_buf4_copy(&K, PSIF_OCC_DPD, "W <OV|OV>");
    dpd_buf4_close(&K);
    
}// end if (reference_ == "RESTRICTED") 


//===========================================================================================
//========================= UHF =============================================================
//===========================================================================================
else if (reference_ == "UNRESTRICTED") {
    /*
    // W_MNIJ = <MN||IJ>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "MO Ints <OO||OO>");
    dpd_buf4_copy(&K, PSIF_OCC_DPD, "W_1 <OO|OO>");
    dpd_buf4_close(&K);    
  
    // W_mnij = <mn||ij>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "MO Ints <oo||oo>");
    dpd_buf4_copy(&K, PSIF_OCC_DPD, "W_1 <oo|oo>");
    dpd_buf4_close(&K);    

    // W_MnIj = <Mn|Ij>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "MO Ints <Oo|Oo>");
    dpd_buf4_copy(&K, PSIF_OCC_DPD, "W_1 <Oo|Oo>");
    dpd_buf4_close(&K);
    
    //Print 
    if (print_ > 3) {
      dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "W_1 <OO|OO>");
      dpd_buf4_print(&W, outfile, 1);
      dpd_buf4_close(&W);
    
      dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "W_1 <oo|oo>");
      dpd_buf4_print(&W, outfile, 1);
      dpd_buf4_close(&W);
    
      dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "W_1 <Oo|Oo>");
      dpd_buf4_print(&W, outfile, 1);
      dpd_buf4_close(&W); 
    }
    
    
    
    // W_ABEF = <AB||EF>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV||VV>");
    dpd_buf4_copy(&K, PSIF_OCC_DPD, "W_1 <VV|VV>");
    dpd_buf4_close(&K);
    
    // W_abef = <ab||ef>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "MO Ints <vv||vv>");
    dpd_buf4_copy(&K, PSIF_OCC_DPD, "W_1 <vv|vv>");
    dpd_buf4_close(&K);
    
    // W_AbEf = <Ab|Ef>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "MO Ints <Vv|Vv>");
    dpd_buf4_copy(&K, PSIF_OCC_DPD, "W_1 <Vv|Vv>");
    dpd_buf4_close(&K);
    
    //Print 
    if (print_ > 3) {
      dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "W_1 <VV|VV>");
      dpd_buf4_print(&W, outfile, 1);
      dpd_buf4_close(&W);
    
      dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "W_1 <vv|vv>");
      dpd_buf4_print(&W, outfile, 1);
      dpd_buf4_close(&W);
    
      dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "W_1 <Vv|Vv>");
      dpd_buf4_print(&W, outfile, 1);
      dpd_buf4_close(&W); 
    }
    */
    
    
    // W_MBEJ => W(ME,JB) = <MB||EJ> = (ME|JB) - <ME|JB>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    dpd_buf4_copy(&K, PSIF_OCC_DPD, "W (OV|OV)");
    dpd_buf4_close(&K);
    dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W (OV|OV)");
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>");
    dpd_buf4_axpy(&K, &W, -1.0); // -1.0*K + W -> W
    dpd_buf4_close(&K);
    dpd_buf4_close(&W);    
    
    
    // W_mbej => W(me,jb) = <mb||ej> = (me|jb) - <me|jb>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov|ov)");
    dpd_buf4_copy(&K, PSIF_OCC_DPD, "W (ov|ov)");
    dpd_buf4_close(&K);
    dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "W (ov|ov)");
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov|ov>");
    dpd_buf4_axpy(&K, &W, -1.0); // -1.0*K + W -> W
    dpd_buf4_close(&K);
    dpd_buf4_close(&W);
    
    
    // W_MbEj => W(ME,jb) = <Mb||Ej> = (ME|jb)
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
    dpd_buf4_copy(&K, PSIF_OCC_DPD, "W (OV|ov)");
    dpd_buf4_close(&K);
    
    
    // W_MbeJ => W(Me,Jb) = <Mb||eJ> = -(MJ|be)
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[v,v]"),
                  ID("[O>=O]+"), ID("[v>=v]+"), 0, "MO Ints (OO|vv)");
    dpd_buf4_sort(&K, PSIF_OCC_DPD , psqr, ID("[O,v]"), ID("[O,v]"), "W (Ov|Ov)");
    dpd_buf4_close(&K);
    dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "W (Ov|Ov)");
    dpd_buf4_scm(&W, -1.0);
    dpd_buf4_close(&W);
    
    
    // W_mBEj => W(mE,jB) = <mB||Ej> = -(BE|mj)
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[o,o]"),
                  ID("[V>=V]+"), ID("[o>=o]+"), 0, "MO Ints (VV|oo)");
    dpd_buf4_sort(&K, PSIF_OCC_DPD , rqsp, ID("[o,V]"), ID("[o,V]"), "W (oV|oV)");
    dpd_buf4_close(&K);
    dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "W (oV|oV)");
    dpd_buf4_scm(&W, -1.0);
    dpd_buf4_close(&W);
    
    
    // it is unnecessary for ocepa, but i will create it so that can use DPD with OOC
    // W_mBeJ => W(me,JB) = <mB||eJ> = (JB|me) = W(JB,me) = W_JeBm
    dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "W (OV|ov)");
    dpd_buf4_sort(&W, PSIF_OCC_DPD , rspq, ID("[o,v]"), ID("[O,V]"), "W (ov|OV)");
    dpd_buf4_close(&W);
    
    //Print 
    if (print_ > 3) {
      dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W (OV|OV)");
      dpd_buf4_print(&W, outfile, 1);
      dpd_buf4_close(&W);

      dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "W (ov|ov)");
      dpd_buf4_print(&W, outfile, 1);
      dpd_buf4_close(&W);

      dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "W (OV|ov)");
      dpd_buf4_print(&W, outfile, 1);
      dpd_buf4_close(&W);

      dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "W (Ov|Ov)");
      dpd_buf4_print(&W, outfile, 1);
      dpd_buf4_close(&W);

      dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "W (oV|oV)");
      dpd_buf4_print(&W, outfile, 1);
      dpd_buf4_close(&W);
    }
    

}// end if (reference_ == "UNRESTRICTED") 

     psio_->close(PSIF_LIBTRANS_DPD, 1);
     psio_->close(PSIF_OCC_DPD, 1);

} // end W_int
}} // End Namespaces

