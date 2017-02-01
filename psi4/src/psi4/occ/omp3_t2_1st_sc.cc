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


using namespace std;


namespace psi{ namespace occwave{

void OCCWave::omp3_t2_1st_sc()
{

//===========================================================================================
//========================= RHF =============================================================
//===========================================================================================
if (reference_ == "RESTRICTED") {
     dpdbuf4 K, T, D, Tau, Ttemp, Tss;

     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);

    // T_ij^ab = <ij|ab>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    global_dpd_->buf4_copy(&K, PSIF_OCC_DPD, "T2_1 <OO|VV>");
    global_dpd_->buf4_close(&K);

    // T_ij^ab = T_ij^ab / D_ij^ab
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "D <OO|VV>");
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    global_dpd_->buf4_dirprd(&D, &T);
    global_dpd_->buf4_close(&D);

     // Build Tau(ij,ab) = 2*T(ij,ab) - T(ji,ab)
     // Build TAA(ij,ab) = T(ij,ab) - T(ji,ab)
     global_dpd_->buf4_copy(&T, PSIF_OCC_DPD, "Tau_1 <OO|VV>");
     global_dpd_->buf4_copy(&T, PSIF_OCC_DPD, "T2_1AA <OO|VV>");
     global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "T2_1jiab <OO|VV>");
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

     if (print_ > 4) global_dpd_->buf4_print(&T, "outfile", 1);
     global_dpd_->buf4_close(&T);


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


     dpdbuf4 K, T, D;

     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);

     // Build T2AA
     // T_IJ^AB = <IJ||AB>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
    global_dpd_->buf4_copy(&K, PSIF_OCC_DPD, "T2_1 <OO|VV>");
    //if (print_ > 1) dpd_buf4_print(&K, outfile, 1);
    global_dpd_->buf4_close(&K);


    // T_IJ^AB = T_IJ^AB / D_IJ^AB
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "D <OO|VV>");
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    global_dpd_->buf4_dirprd(&D, &T);
    global_dpd_->buf4_close(&D);
    if (print_ > 1) global_dpd_->buf4_print(&T, "outfile", 1);
    global_dpd_->buf4_close(&T);


    // Build T2BB
    // T_ij^ab = <ij|ab>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
    global_dpd_->buf4_copy(&K, PSIF_OCC_DPD, "T2_1 <oo|vv>");
    global_dpd_->buf4_close(&K);


    // T_ij^ab = T_ij^ab / D_ij^ab
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "D <oo|vv>");
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
    global_dpd_->buf4_dirprd(&D, &T);
    global_dpd_->buf4_close(&D);
    if (print_ > 1) global_dpd_->buf4_print(&T, "outfile", 1);
    global_dpd_->buf4_close(&T);


    // Build T2AB
    // T_Ij^Ab = <Ij|Ab>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    global_dpd_->buf4_copy(&K, PSIF_OCC_DPD, "T2_1 <Oo|Vv>");
    global_dpd_->buf4_close(&K);


    // T_Ij^Ab = T_Ij^Ab / D_Ij^Ab
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "D <Oo|Vv>");
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
    global_dpd_->buf4_dirprd(&D, &T);
    global_dpd_->buf4_close(&D);
    if (print_ > 1) global_dpd_->buf4_print(&T, "outfile", 1);
    global_dpd_->buf4_close(&T);



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

    /*
    // T_Ij^Ab => T(IA,jb), T(Ib,jA)
    dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
    dpd_buf4_sort(&T, PSIF_OCC_DPD , prqs, ID("[O,V]"), ID("[o,v]"), "T2_1 (OV|ov)");
    dpd_buf4_sort(&T, PSIF_OCC_DPD , psqr, ID("[O,v]"), ID("[o,V]"), "T2_1 (Ov|oV)");
    dpd_buf4_close(&T);
    */

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

    psio_->close(PSIF_LIBTRANS_DPD, 1);
    psio_->close(PSIF_OCC_DPD, 1);
}// end if (reference_ == "UNRESTRICTED")

} // end t2_1st_sc
}} // End Namespaces


