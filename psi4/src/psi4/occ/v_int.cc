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

void OCCWave::v_int()
{

 if (reference_ == "RESTRICTED") {

     dpdbuf4 T, Tau, V;

     psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);

    // Build V(IJ,KL)
    // V_IJKL(2) = 2\sum_{E,F} T_IJ^EF(1) (2T_KL^EF(1) - T_LK^EF(1))
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
    global_dpd_->buf4_init(&Tau, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tau <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "V <OO|OO>");
    global_dpd_->contract444(&T, &Tau, &V, 0, 0, 2.0, 0.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&Tau);

    //Print
    if (print_ > 3) {
      global_dpd_->buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "V <OO|OO>");
      global_dpd_->buf4_print(&V, "outfile", 1);
      global_dpd_->buf4_close(&V);
    }

    // Build V(IA,JB)
    // V_IAJB(2) => V(IB,JA) = \sum_{M,E} (2T_MI^BE(1) - T_IM^BE(1)) T_JM^EA(1)
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2pp (OV|OV)");
    global_dpd_->buf4_init(&Tau, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Taupp (OV|OV)");
    global_dpd_->buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "V <IB|JA>");
    global_dpd_->contract444(&Tau, &T, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&Tau);

    // V_IAJB(2) => V(IB,JA) += \sum_{M,E} (2T_IM^BE(1) - T_MI^BE(1)) T_JM^AE(1)
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2 (OV|OV)");
    global_dpd_->buf4_init(&Tau, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Tau (OV|OV)");
    global_dpd_->contract444(&Tau, &T, &V, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&Tau);

    // V(IB,JA) => V(IA,JB)
    global_dpd_->buf4_sort(&V, PSIF_OCC_DENSITY , psrq, ID("[O,V]"), ID("[O,V]"), "V <OV|OV>");
    global_dpd_->buf4_close(&V);


    // Build V(IA,BJ)
    // V_IABJ(2) => V(IB,JA) = \sum_{M,E} (2T_IM^BE(1) - T_MI^BE(1)) (2T_JM^AE(1) - T_JM^EA(1))
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Tau (OV|OV)");
    global_dpd_->buf4_init(&Tau, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Tau (OV|OV)");
    global_dpd_->buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "V (IB|JA)");
    global_dpd_->contract444(&Tau, &T, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&Tau);

    // V(IB,JA) => V(IA,BJ)
    global_dpd_->buf4_sort(&V, PSIF_OCC_DENSITY , psqr, ID("[O,V]"), ID("[V,O]"), "V <OV|VO>");
    global_dpd_->buf4_close(&V);

     psio_->close(PSIF_OCC_DENSITY, 1);
     psio_->close(PSIF_OCC_DPD, 1);

 }// end if (reference_ == "RESTRICTED")

 else if (reference_ == "UNRESTRICTED") {

     dpdbuf4 T, L, V, V1, V2;

     psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);

    // Build V(IJ,KL)
    // V_IJKL(2) = 1/2 \sum_{E,F} T_IJ^EF(1) L_EF^KL(1) = 1/2 \sum_{E,F} T_IJ^EF(1) T_KL^EF(1)
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
    global_dpd_->buf4_init(&L, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "V <OO|OO>");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 0.5, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    // Build V(ij,kl)
    // V_ijkl(2) = 1/2 \sum_{e,f} T_ij^ef(1) L_ef^kl(1) = 1/2 \sum_{e,f} T_ij^ef(1) T_kl^ef(1)
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
    global_dpd_->buf4_init(&L, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "V <oo|oo>");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 0.5, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);


    // Build V(Ij,Kl)
    // V_IjKl(2) = \sum_{E,f} T_Ij^Ef(1) L_EF^KL(1) = \sum_{E,f} T_Ij^Ef(1) T_Kl^Ef(1)
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
    global_dpd_->buf4_init(&L, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
    global_dpd_->buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "V <Oo|Oo>");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);


    //Print
    if (print_ > 3) {
      global_dpd_->buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "V <OO|OO>");
      global_dpd_->buf4_print(&V, "outfile", 1);
      global_dpd_->buf4_close(&V);

      global_dpd_->buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "V <oo|oo>");
      global_dpd_->buf4_print(&V, "outfile", 1);
      global_dpd_->buf4_close(&V);

      global_dpd_->buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "V <Oo|Oo>");
      global_dpd_->buf4_print(&V, "outfile", 1);
      global_dpd_->buf4_close(&V);
    }

    // Build V(IA,JB)
    // V_IAJB(2) => V(IB,JA) = 1/2 \sum_{M,E} T_IM^BE(1) T_JM^AE(1) = 1/2 \sum_{M,E} T(IB,ME) T(JA,ME)
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2 (OV|OV)");
    global_dpd_->buf4_init(&L, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2 (OV|OV)");
    global_dpd_->buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "V <IB|JA>");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 0.5, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);

    // V_IAJB(2) => V(IB,JA) += 1/2 \sum_{m,e} T_Im^Be(1) T_Jm^Ae(1) = 1/2 \sum_{m,e} T(IB,me) T(JA,me)
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "T2 (OV|ov)");
    global_dpd_->buf4_init(&L, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "T2 (OV|ov)");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 0.5, 1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);

    // V(IB,JA) => V(IA,JB)
    global_dpd_->buf4_sort(&V, PSIF_OCC_DENSITY , psrq, ID("[O,V]"), ID("[O,V]"), "V <OV|OV>");
    global_dpd_->buf4_close(&V);


    // Build V(ia,jb)
    // V_iajb(2) => V(ib,ja) = 1/2 \sum_{m,e} T_im^be(1) T_jm^ae(1) = 1/2 \sum_{m,e} T(ib,me) T(ja,me)
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "T2 (ov|ov)");
    global_dpd_->buf4_init(&L, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "T2 (ov|ov)");
    global_dpd_->buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "V <ib|ja>");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 0.5, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);

    // V_iajb(2) => V(ib,ja) += 1/2 \sum_{M,E} T_Mi^Eb(1) T_Mj^Ea(1) = 1/2 \sum_{M,E} T(ME,ib) T(ME,ja)
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "T2 (OV|ov)");
    global_dpd_->buf4_init(&L, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "T2 (OV|ov)");
    global_dpd_->contract444(&T, &L, &V, 1, 1, 0.5, 1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);

    // V(ib,ja) => V(ia,jb)
    global_dpd_->buf4_sort(&V, PSIF_OCC_DENSITY , psrq, ID("[o,v]"), ID("[o,v]"), "V <ov|ov>");
    global_dpd_->buf4_close(&V);



    // Build V(Ia,Jb) & V(iA,jB)
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,V]"), ID("[O,v]"),
                  ID("[o,V]"), ID("[O,v]"), 0, "T2 (oV|Ov)");
    global_dpd_->buf4_init(&L, PSIF_OCC_DPD, 0, ID("[o,V]"), ID("[O,v]"),
                  ID("[o,V]"), ID("[O,v]"), 0, "T2 (oV|Ov)");
    global_dpd_->buf4_init(&V1, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "V <Ib|Ja>");
    global_dpd_->buf4_init(&V2, PSIF_OCC_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "V <iB|jA>");
    // V_IaJb(2) => V(Ib,Ja) = 1/2 \sum_{m,E} T_Im^Eb(1) T_Jm^Ea(1) = 1/2 \sum_{m,E} T(mE,IB) T(mE,Ja)
    global_dpd_->contract444(&T, &L, &V1, 1, 1, 0.5, 0.0);
    // V_iAjB(2) => V(iB,jA) = 1/2 \sum_{M,e} T_Mi^Be(1) T_Mj^Ae(1) = 1/2 \sum_{M,e} T(iB,Me) T(jA,Me)
    global_dpd_->contract444(&T, &L, &V2, 0, 0, 0.5, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);

    // V(Ib,Ja) => V(Ia,Jb) & V(iB,jA) => V(iA,jB)
    global_dpd_->buf4_sort(&V1, PSIF_OCC_DENSITY , psrq, ID("[O,v]"), ID("[O,v]"), "V <Ov|Ov>");
    global_dpd_->buf4_close(&V1);
    global_dpd_->buf4_sort(&V2, PSIF_OCC_DENSITY , psrq, ID("[o,V]"), ID("[o,V]"), "V <oV|oV>");
    global_dpd_->buf4_close(&V2);


    // Build V(Ia,jB)
    // V_IajB(2) => V(IB,ja) = 1/2 \sum_{M,E} T_IM^BE(1) T_Mj^Ea(1) = 1/2 \sum_{M,E} T(IB,ME) T(ME,ja)
    //                       = 1/2 \sum_{M,E} T(IB,ME) T(ja,ME)
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "T2 (OV|OV)");
    global_dpd_->buf4_init(&L, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[O,V]"),
                  ID("[o,v]"), ID("[O,V]"), 0, "T2 (ov|OV)");
    global_dpd_->buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "V <IB|ja>");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 0.5, 0.0);
    global_dpd_->buf4_close(&T);
    //

    // V_IajB(2) => V(IB,ja) += 1/2 \sum_{m,e} T_Im^Be(1) T_jm^ae(1) = 1/2 \sum_{m,e} T(IB,me) T(ja,me)
    //                        = 1/2 \sum_{m,e} T(me,IB) T(me,ja)
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "T2 (ov|ov)");
    global_dpd_->contract444(&L, &T, &V, 1, 1, 0.5, 1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);

    // V(IB,ja) => V(Ia,jB)
    global_dpd_->buf4_sort(&V, PSIF_OCC_DENSITY , psrq, ID("[O,v]"), ID("[o,V]"), "V <Ov|oV>");
    global_dpd_->buf4_close(&V);

    // Note: V(iA,Jb) = V(Jb,iA)

     psio_->close(PSIF_OCC_DENSITY, 1);
     psio_->close(PSIF_OCC_DPD, 1);

 }// end if (reference_ == "UNRESTRICTED")

} // end V_int
}} // End Namespaces
