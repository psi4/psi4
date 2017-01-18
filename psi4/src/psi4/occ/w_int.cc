/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
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

#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "defines.h"
#include "occwave.h"


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
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    global_dpd_->buf4_copy(&K, PSIF_OCC_DPD, "W (OV|OV)");
    global_dpd_->buf4_close(&K);

    // W_mbje => W'(me,jb) = <mb|je> = <me|jb>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>");
    global_dpd_->buf4_copy(&K, PSIF_OCC_DPD, "W <OV|OV>");
    global_dpd_->buf4_close(&K);

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
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    global_dpd_->buf4_copy(&K, PSIF_OCC_DPD, "W (OV|OV)");
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W (OV|OV)");
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>");
    global_dpd_->buf4_axpy(&K, &W, -1.0); // -1.0*K + W -> W
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&W);


    // W_mbej => W(me,jb) = <mb||ej> = (me|jb) - <me|jb>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov|ov)");
    global_dpd_->buf4_copy(&K, PSIF_OCC_DPD, "W (ov|ov)");
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "W (ov|ov)");
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov|ov>");
    global_dpd_->buf4_axpy(&K, &W, -1.0); // -1.0*K + W -> W
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&W);


    // W_MbEj => W(ME,jb) = <Mb||Ej> = (ME|jb)
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
    global_dpd_->buf4_copy(&K, PSIF_OCC_DPD, "W (OV|ov)");
    global_dpd_->buf4_close(&K);


    // W_MbeJ => W(Me,Jb) = <Mb||eJ> = -(MJ|be)
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[v,v]"),
                  ID("[O>=O]+"), ID("[v>=v]+"), 0, "MO Ints (OO|vv)");
    global_dpd_->buf4_sort(&K, PSIF_OCC_DPD , psqr, ID("[O,v]"), ID("[O,v]"), "W (Ov|Ov)");
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "W (Ov|Ov)");
    global_dpd_->buf4_scm(&W, -1.0);
    global_dpd_->buf4_close(&W);


    // W_mBEj => W(mE,jB) = <mB||Ej> = -(BE|mj)
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[o,o]"),
                  ID("[V>=V]+"), ID("[o>=o]+"), 0, "MO Ints (VV|oo)");
    global_dpd_->buf4_sort(&K, PSIF_OCC_DPD , rqsp, ID("[o,V]"), ID("[o,V]"), "W (oV|oV)");
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "W (oV|oV)");
    global_dpd_->buf4_scm(&W, -1.0);
    global_dpd_->buf4_close(&W);


    // it is unnecessary for ocepa, but i will create it so that can use DPD with OOC
    // W_mBeJ => W(me,JB) = <mB||eJ> = (JB|me) = W(JB,me) = W_JeBm
    global_dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "W (OV|ov)");
    global_dpd_->buf4_sort(&W, PSIF_OCC_DPD , rspq, ID("[o,v]"), ID("[O,V]"), "W (ov|OV)");
    global_dpd_->buf4_close(&W);

    //Print
    if (print_ > 3) {
      global_dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "W (OV|OV)");
      global_dpd_->buf4_print(&W, "outfile", 1);
      global_dpd_->buf4_close(&W);

      global_dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "W (ov|ov)");
      global_dpd_->buf4_print(&W, "outfile", 1);
      global_dpd_->buf4_close(&W);

      global_dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "W (OV|ov)");
      global_dpd_->buf4_print(&W, "outfile", 1);
      global_dpd_->buf4_close(&W);

      global_dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "W (Ov|Ov)");
      global_dpd_->buf4_print(&W, "outfile", 1);
      global_dpd_->buf4_close(&W);

      global_dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "W (oV|oV)");
      global_dpd_->buf4_print(&W, "outfile", 1);
      global_dpd_->buf4_close(&W);
    }


}// end if (reference_ == "UNRESTRICTED")

     psio_->close(PSIF_LIBTRANS_DPD, 1);
     psio_->close(PSIF_OCC_DPD, 1);

} // end W_int
}} // End Namespaces
