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
  
void OMP3Wave::W_1st_order()
{   
     dpdbuf4 K, W;
     
     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OMP3_DPD, PSIO_OPEN_OLD);
     
    /*
    // W_MNIJ = <MN||IJ>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "MO Ints <OO||OO>");
    dpd_buf4_copy(&K, PSIF_OMP3_DPD, "W_1 <OO|OO>");
    dpd_buf4_close(&K);    
  
    // W_mnij = <mn||ij>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "MO Ints <oo||oo>");
    dpd_buf4_copy(&K, PSIF_OMP3_DPD, "W_1 <oo|oo>");
    dpd_buf4_close(&K);    

    // W_MnIj = <Mn|Ij>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "MO Ints <Oo|Oo>");
    dpd_buf4_copy(&K, PSIF_OMP3_DPD, "W_1 <Oo|Oo>");
    dpd_buf4_close(&K);
    
    //Print 
    if (print_ > 3) {
      dpd_buf4_init(&W, PSIF_OMP3_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "W_1 <OO|OO>");
      dpd_buf4_print(&W, outfile, 1);
      dpd_buf4_close(&W);
    
      dpd_buf4_init(&W, PSIF_OMP3_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "W_1 <oo|oo>");
      dpd_buf4_print(&W, outfile, 1);
      dpd_buf4_close(&W);
    
      dpd_buf4_init(&W, PSIF_OMP3_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "W_1 <Oo|Oo>");
      dpd_buf4_print(&W, outfile, 1);
      dpd_buf4_close(&W); 
    }
    
    
    
    // W_ABEF = <AB||EF>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV||VV>");
    dpd_buf4_copy(&K, PSIF_OMP3_DPD, "W_1 <VV|VV>");
    dpd_buf4_close(&K);
    
    // W_abef = <ab||ef>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "MO Ints <vv||vv>");
    dpd_buf4_copy(&K, PSIF_OMP3_DPD, "W_1 <vv|vv>");
    dpd_buf4_close(&K);
    
    // W_AbEf = <Ab|Ef>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "MO Ints <Vv|Vv>");
    dpd_buf4_copy(&K, PSIF_OMP3_DPD, "W_1 <Vv|Vv>");
    dpd_buf4_close(&K);
    
    //Print 
    if (print_ > 3) {
      dpd_buf4_init(&W, PSIF_OMP3_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "W_1 <VV|VV>");
      dpd_buf4_print(&W, outfile, 1);
      dpd_buf4_close(&W);
    
      dpd_buf4_init(&W, PSIF_OMP3_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "W_1 <vv|vv>");
      dpd_buf4_print(&W, outfile, 1);
      dpd_buf4_close(&W);
    
      dpd_buf4_init(&W, PSIF_OMP3_DPD, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "W_1 <Vv|Vv>");
      dpd_buf4_print(&W, outfile, 1);
      dpd_buf4_close(&W); 
    }
    */
    
    
    // W_MBEJ => W(ME,JB) = <MB||EJ> = (ME|JB) - <ME|JB>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    dpd_buf4_copy(&K, PSIF_OMP3_DPD, "W_1 (OV|OV)");
    dpd_buf4_close(&K);
    dpd_buf4_init(&W, PSIF_OMP3_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W_1 (OV|OV)");
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>");
    dpd_buf4_axpy(&K, &W, -1.0); // -1.0*K + W -> W
    dpd_buf4_close(&K);
    dpd_buf4_close(&W);    
    
    //Print 
    if (print_ > 3) {
      dpd_buf4_init(&W, PSIF_OMP3_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W_1 (OV|OV)");
      dpd_buf4_print(&W, outfile, 1);
      dpd_buf4_close(&W);
    }
    
    // W_mbej => W(me,jb) = <mb||ej> = (me|jb) - <me|jb>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov|ov)");
    dpd_buf4_copy(&K, PSIF_OMP3_DPD, "W_1 (ov|ov)");
    dpd_buf4_close(&K);
    dpd_buf4_init(&W, PSIF_OMP3_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "W_1 (ov|ov)");
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov|ov>");
    dpd_buf4_axpy(&K, &W, -1.0); // -1.0*K + W -> W
    dpd_buf4_close(&K);
    dpd_buf4_close(&W);
    
    //Print 
    if (print_ > 3) {
      dpd_buf4_init(&W, PSIF_OMP3_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "W_1 (ov|ov)");
      dpd_buf4_print(&W, outfile, 1);
      dpd_buf4_close(&W);
    }
    
    // W_MbEj => W(ME,jb) = <Mb||Ej> = (ME|jb)
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
    dpd_buf4_copy(&K, PSIF_OMP3_DPD, "W_1 (OV|ov)");
    dpd_buf4_close(&K);
    
    //Print 
    if (print_ > 3) {
      dpd_buf4_init(&W, PSIF_OMP3_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "W_1 (OV|ov)");
      dpd_buf4_print(&W, outfile, 1);
      dpd_buf4_close(&W);
    }
    
    
    // W_MbeJ => W(Me,Jb) = <Mb||eJ> = -(MJ|be)
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[v,v]"),
                  ID("[O>=O]+"), ID("[v>=v]+"), 0, "MO Ints (OO|vv)");
    dpd_buf4_sort(&K, PSIF_OMP3_DPD , psqr, ID("[O,v]"), ID("[O,v]"), "W_1 (Ov|Ov)");
    dpd_buf4_close(&K);
    dpd_buf4_init(&W, PSIF_OMP3_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "W_1 (Ov|Ov)");
    dpd_buf4_scm(&W, -1.0);
    dpd_buf4_close(&W);
    
    //Print 
    if (print_ > 3) {
      dpd_buf4_init(&W, PSIF_OMP3_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "W_1 (Ov|Ov)");
      dpd_buf4_print(&W, outfile, 1);
      dpd_buf4_close(&W);
    }
      
    
    // W_mBEj => W(mE,jB) = <mB||Ej> = -(BE|mj)
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[o,o]"),
                  ID("[V>=V]+"), ID("[o>=o]+"), 0, "MO Ints (VV|oo)");
    dpd_buf4_sort(&K, PSIF_OMP3_DPD , rqsp, ID("[o,V]"), ID("[o,V]"), "W_1 (oV|oV)");
    dpd_buf4_close(&K);
    dpd_buf4_init(&W, PSIF_OMP3_DPD, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "W_1 (oV|oV)");
    dpd_buf4_scm(&W, -1.0);
    dpd_buf4_close(&W);
    
    //Print 
    if (print_ > 3) {
      dpd_buf4_init(&W, PSIF_OMP3_DPD, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "W_1 (oV|oV)");
      dpd_buf4_print(&W, outfile, 1);
      dpd_buf4_close(&W);
    }
    
    // it is unnecessary for omp3, but i will create it so that can use DPD with OOC
    // W_mBeJ => W(me,JB) = <mB||eJ> = (JB|me) = W(JB,me) = W_JeBm
    dpd_buf4_init(&W, PSIF_OMP3_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "W_1 (OV|ov)");
    dpd_buf4_sort(&W, PSIF_OMP3_DPD , rspq, ID("[o,v]"), ID("[O,V]"), "W_1 (ov|OV)");
    dpd_buf4_close(&W);
    
    
     psio_->close(PSIF_LIBTRANS_DPD, 1);
     psio_->close(PSIF_OMP3_DPD, 1);

} // end W_1st_order
}} // End Namespaces

