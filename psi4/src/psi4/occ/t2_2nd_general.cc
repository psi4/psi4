/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/libqt/qt.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "defines.h"
#include "occwave.h"

#include <cmath>

namespace psi {
namespace occwave {

void OCCWave::t2_2nd_general() {
    double omp2p5_factor = (wfn_type_ == "OMP2.5" ? 0.5 : 1);

    //===========================================================================================
    //========================= RHF =============================================================
    //===========================================================================================
    if (reference_ == "RESTRICTED") {
        dpdbuf4 K, T, Tau, Tnew, D, R, Tp, W, TAA, TAB, TBB, Ttemp;
        dpdfile2 Fo, Fv;
        int nElements;

        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
        psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);

        // Build T(IA,JB)
        // T_IJ^AB(2) = \sum_{M,E} Tau_IM^AE(1) W_MBEJ(1) => T(IA,JB)(2) = \sum_{M,E} Tau'(IA,ME) (ME|JB)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "T2_2 (IA|JB)");
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "Tau_1 (OV|OV)");
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "MO Ints (OV|OV)");
        global_dpd_->contract444(&Tp, &K, &T, 0, 0, 1.0, 0.0);
        // T_IJ^AB(2) += \sum_{M,E} Tau_JM^BE(1) W_MAEI(1) => T(IA,JB)(2) = \sum_{M,E} Tau'(JB,ME) (ME|IA)
        global_dpd_->contract444(&K, &Tp, &T, 0, 0, 1.0, 1.0);
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_close(&Tp);

        // T_IJ^AB(2) -= \sum_{m,e} T_im^ae(1) W_mbje(1) => T(IA,JB)(2) -= \sum_{m,e} T'(ia,me) <me|jb>
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "T2_1 (OV|OV)");
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "MO Ints <OV|OV>");
        global_dpd_->contract444(&Tp, &K, &T, 0, 0, -1.0, 1.0);
        // T_IJ^AB(2) -= \sum_{m,e} T_jm^be(1) W_maie(1) => T(IA,JB)(2) -= \sum_{m,e} T'(jb,me) <me|ia>
        global_dpd_->contract444(&K, &Tp, &T, 0, 0, -1.0, 1.0);
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_close(&Tp);

        // T(IA,JB) => T_IJ^AB(2)
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "T2_2 <IJ|AB>");
        global_dpd_->buf4_close(&T);

        // Build T(JA,IB)
        // T_IJ^AB(2) = -\sum_{M,E} T_MJ^AE(1) W_MBIE(1) => T(JA,IB)(2) = -\sum_{M,E} T"(JA,ME) <ME|IB>
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "T2_2 (JA|IB)");
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "T2_1pp (OV|OV)");
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "MO Ints <OV|OV>");
        global_dpd_->contract444(&Tp, &K, &T, 0, 0, -1.0, 0.0);
        // T_IJ^AB(2) = -\sum_{M,E} T_IM^EB(1) W_MAJE(1) => T(JA,IB)(2) = -\sum_{M,E} T"(ME,IB) <ME|JA>
        global_dpd_->contract444(&K, &Tp, &T, 1, 1, -1.0, 1.0);
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_close(&Tp);

        // T(JA,IB) => T_IJ^AB(2)
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, rpqs, ID("[O,O]"), ID("[V,V]"), "T2_2 (IJ|AB)");
        global_dpd_->buf4_close(&T);

        // Build T2AAnew
        // T_IJ^AB(2) = T(IA,JB)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_2 <IJ|AB>");
        global_dpd_->buf4_copy(&T, PSIF_OCC_DPD, "T2_2new <OO|VV>");
        global_dpd_->buf4_close(&T);

        // T_IJ^AB(2) += T(JA,IB)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_2new <OO|VV>");
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_2 (IJ|AB)");
        global_dpd_->buf4_axpy(&Tp, &T, 1.0);  // 1.0*Tp + T -> T
        global_dpd_->buf4_close(&Tp);

        // T_IJ^AB(2) += \sum_{M,N} T_MN^AB(1) W_MNIJ(1) = \sum_{M,N} T_MN^AB(1) <MN|IJ>
        global_dpd_->buf4_init(&TAA, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_1 <OO|VV>");
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), 0,
                               "MO Ints <OO|OO>");
        global_dpd_->contract444(&K, &TAA, &T, 1, 1, 1.0, 1.0);
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_close(&T);

        // T_IJ^AB(2) += 1/2 \sum_{E,F} T_IJ^EF(1) <EF||AB> = \sum_{E,F} T_IJ^EF(1) <AB|EF>
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), 0,
                               "Z2_2 <VV|OO>");
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"), ID("[V,V]"), ID("[V,V]"), 0,
                               "MO Ints <VV|VV>");
        global_dpd_->contract444(&K, &TAA, &T, 0, 0, 1.0, 0.0);
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_close(&TAA);
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, rspq, ID("[O,O]"), ID("[V,V]"), "Z2_2 <OO|VV>");
        global_dpd_->buf4_close(&T);
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_2new <OO|VV>");
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "Z2_2 <OO|VV>");
        global_dpd_->buf4_axpy(&Tp, &T, 1.0);  // 1.0*Tp + T -> T
        global_dpd_->buf4_close(&T);
        global_dpd_->buf4_close(&Tp);

        // initialize Tnew and Told
        global_dpd_->buf4_init(&Tnew, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_2new <OO|VV>");
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_2 <OO|VV>");

        // T_Ij^Ab = \sum_{e} T_Ij^Ae * F_be
        global_dpd_->file2_init(&Fv, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");
        global_dpd_->contract424(&T, &Fv, &Tnew, 3, 1, 0, 1.0, 1.0);

        // T_Ij^Ab = \sum_{E} T_Ij^Eb * F_AE
        global_dpd_->contract244(&Fv, &T, &Tnew, 1, 2, 1, 1.0, 1.0);
        global_dpd_->file2_close(&Fv);

        // T_Ij^Ab = -\sum_{m} T_Im^Ab * F_mj
        global_dpd_->file2_init(&Fo, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
        global_dpd_->contract424(&T, &Fo, &Tnew, 1, 0, 1, -1.0, 1.0);

        // T_Ij^Ab = - \sum_{M} T_Mj^Ab * F_MI
        global_dpd_->contract244(&Fo, &T, &Tnew, 0, 0, 0, -1.0, 1.0);
        global_dpd_->file2_close(&Fo);
        global_dpd_->buf4_close(&T);

        // T_IJ^AB = T_IJ^AB / D_IJ^AB
        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "D <OO|VV>");
        global_dpd_->buf4_dirprd(&D, &Tnew);
        global_dpd_->buf4_close(&D);
        if (print_ > 2) global_dpd_->buf4_print(&Tnew, "outfile", 1);
        global_dpd_->buf4_close(&Tnew);

        // Check convergence?
        global_dpd_->buf4_init(&Tnew, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_2new <OO|VV>");
        global_dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "RT2_2 <OO|VV>");
        global_dpd_->buf4_init(&R, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "RT2_2 <OO|VV>");
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_2 <OO|VV>");
        global_dpd_->buf4_axpy(&T, &R, -1.0);  // -1.0*T + R -> R
        global_dpd_->buf4_close(&T);

        nElements = 0;
        for (int h = 0; h < nirrep_; h++) nElements += R.params->coltot[h] * R.params->rowtot[h];
        rms_t2 = 0.0;
        rms_t2 = global_dpd_->buf4_dot_self(&R);
        global_dpd_->buf4_close(&R);
        rms_t2 = std::sqrt(rms_t2) / nElements;

        // Reset
        global_dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "T2_2 <OO|VV>");
        if (print_ > 1) global_dpd_->buf4_print(&Tnew, "outfile", 1);
        global_dpd_->buf4_close(&Tnew);


        psio_->close(PSIF_LIBTRANS_DPD, 1);
        psio_->close(PSIF_OCC_DPD, 1);

    }  // end if (reference_ == "RESTRICTED")

    //===========================================================================================
    //========================= UHF =============================================================
    //===========================================================================================
    else if (reference_ == "UNRESTRICTED") {
        /********************************************************************************************/
        /************************** Build W intermediates *******************************************/
        /********************************************************************************************/
        timer_on("W int");
        w_1st_order();
        timer_off("W int");

        dpdbuf4 K, T, Tnew, D, R, Tp, W, TAA, TAB, TBB;
        dpdfile2 Fo, Fv;
        int nElements;

        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
        psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);

        /********************************************************************************************/
        /************************** Alpha-Alpha spin case *******************************************/
        /********************************************************************************************/
        // Build T(IA,JB)
        // T_IJ^AB(2) = \sum_{M,E} T_IM^AE(1) W_MBEJ(1) => T(IA,JB)(2) = \sum_{M,E} T(IA,ME) W(ME,JB)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "T2_2 (IA|JB)");
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "T2_1 (OV|OV)");
        global_dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "W_1 (OV|OV)");
        global_dpd_->contract444(&Tp, &W, &T, 0, 0, 1.0, 0.0);
        // T_IJ^AB(2) += \sum_{M,E} T_JM^BE(1) W_MAEI(1) => T(IA,JB)(2) = \sum_{M,E} T(JB,ME) W(ME,IA)
        global_dpd_->contract444(&W, &Tp, &T, 0, 0, 1.0, 1.0);
        global_dpd_->buf4_close(&W);
        global_dpd_->buf4_close(&Tp);

        // T_IJ^AB(2) += \sum_{m,e} T_Im^Ae(1) W_JeBm(1) => T(IA,JB)(2) += \sum_{m,e} T(IA,me) W(JB,me)
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                               "T2_1 (OV|ov)");
        global_dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                               "W_1 (OV|ov)");
        global_dpd_->contract444(&Tp, &W, &T, 0, 0, 1.0, 1.0);
        // T_IJ^AB(2) += \sum_{m,e} T_Jm^Be(1) W_IeAm(1) => T(IA,JB)(2) += \sum_{m,e} T(JB,me) W(IA,me)
        global_dpd_->contract444(&W, &Tp, &T, 0, 0, 1.0, 1.0);
        global_dpd_->buf4_close(&W);
        global_dpd_->buf4_close(&Tp);

        // T(IA,JB) => T_IJ^AB(2)
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "T2_2 <IJ|AB>");
        global_dpd_->buf4_close(&T);

        // Build T(JA,IB)
        // T_IJ^AB(2) = -\sum_{M,E} T_JM^AE(1) W_MBEI(1) => T(JA,IB)(2) = -\sum_{M,E} T(JA,ME) W(ME,IB)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "T2_2 (JA|IB)");
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "T2_1 (OV|OV)");
        global_dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "W_1 (OV|OV)");
        global_dpd_->contract444(&Tp, &W, &T, 0, 0, -1.0, 0.0);
        // T_IJ^AB(2) = -\sum_{M,E} T_IM^BE(1) W_MAEJ(1) => T(JA,IB)(2) = -\sum_{M,E} T(IB,ME) W(ME,JA)
        global_dpd_->contract444(&W, &Tp, &T, 0, 0, -1.0, 1.0);
        global_dpd_->buf4_close(&W);
        global_dpd_->buf4_close(&Tp);

        // T_IJ^AB(2) = -\sum_{m,e} T_Jm^Ae(1) W_IeBm(1) => T(JA,IB)(2) -= \sum_{m,e} T(JA,me) W(IB,me)
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                               "T2_1 (OV|ov)");
        global_dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                               "W_1 (OV|ov)");
        global_dpd_->contract444(&Tp, &W, &T, 0, 0, -1.0, 1.0);
        // T_IJ^AB(2) = -\sum_{m,e} T_Im^Be(1) W_JeAm(1) => T(JA,IB)(2) -= \sum_{m,e} T(IB,me) W(JA,me)
        global_dpd_->contract444(&W, &Tp, &T, 0, 0, -1.0, 1.0);
        global_dpd_->buf4_close(&W);
        global_dpd_->buf4_close(&Tp);

        // T(JA,IB) => T_IJ^AB(2)
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, rpqs, ID("[O,O]"), ID("[V,V]"), "T2_2 (IJ|AB)");
        global_dpd_->buf4_close(&T);

        // Build T2AAnew
        // T_IJ^AB(2) = T(IA,JB)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_2 <IJ|AB>");
        global_dpd_->buf4_copy(&T, PSIF_OCC_DPD, "T2_2new <OO|VV>");
        global_dpd_->buf4_close(&T);

        // T_IJ^AB(2) += T(JA,IB)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_2new <OO|VV>");
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_2 (IJ|AB)");
        global_dpd_->buf4_axpy(&Tp, &T, 1.0);  // 1.0*Tp + T -> T
        global_dpd_->buf4_close(&Tp);

        // T_IJ^AB(2) += 1/2 \sum_{M,N} T_MN^AB(1) W_MNIJ(1) = 1/2 \sum_{M,N} T_MN^AB(1) <MN||IJ>
        global_dpd_->buf4_init(&TAA, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_1 <OO|VV>");
        // dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[O,O]"),
        //              ID("[O,O]"), ID("[O,O]"), 0, "W_1 <OO|OO>");
        global_dpd_->buf4_init(&W, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), 0,
                               "MO Ints <OO||OO>");
        global_dpd_->contract444(&W, &TAA, &T, 1, 1, 0.5, 1.0);
        global_dpd_->buf4_close(&W);
        global_dpd_->buf4_close(&T);

        /*
        // T_IJ^AB(2) += 1/2 \sum_{E,F} T_IJ^EF(1) W_ABEF(1)
        //dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[V,V]"), ID("[V,V]"),
        //              ID("[V,V]"), ID("[V,V]"), 0, "W_1 <VV|VV>");
        dpd_buf4_init(&W, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                      ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV||VV>");
        dpd_contract444(&TAA, &W, &T, 0, 0, 0.5, 1.0);
        dpd_buf4_close(&W);
        dpd_buf4_close(&T);
        dpd_buf4_close(&TAA);
        */

        /*
        // T_IJ^AB(2) += 1/2 \sum_{E,F} T_IJ^EF(1) <EF||AB> = \sum_{E,F} T_IJ^EF(1) <AB|EF>
        dpd_buf4_init(&W, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                      ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV|VV>");
        dpd_contract444(&TAA, &W, &T, 0, 0, 1.0, 1.0);
        dpd_buf4_close(&W);
        dpd_buf4_close(&T);
        dpd_buf4_close(&TAA);
        */

        // T_IJ^AB(2) += 1/2 \sum_{E,F} T_IJ^EF(1) <EF||AB> = \sum_{E,F} T_IJ^EF(1) <AB|EF>
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), 0,
                               "Z2_2 <VV|OO>");
        global_dpd_->buf4_init(&W, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"), ID("[V,V]"), ID("[V,V]"), 0,
                               "MO Ints <VV|VV>");
        global_dpd_->contract444(&W, &TAA, &T, 0, 0, 1.0, 0.0);
        global_dpd_->buf4_close(&W);
        global_dpd_->buf4_close(&TAA);
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, rspq, ID("[O,O]"), ID("[V,V]"), "Z2_2 <OO|VV>");
        global_dpd_->buf4_close(&T);
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_2new <OO|VV>");
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "Z2_2 <OO|VV>");
        global_dpd_->buf4_axpy(&Tp, &T, 1.0);  // 1.0*Tp + T -> T
        global_dpd_->buf4_close(&T);
        global_dpd_->buf4_close(&Tp);

        // initialize Tnew and Told
        global_dpd_->buf4_init(&Tnew, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_2new <OO|VV>");
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_2 <OO|VV>");

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
        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "D <OO|VV>");
        global_dpd_->buf4_dirprd(&D, &Tnew);
        global_dpd_->buf4_close(&D);
        if (print_ > 2) global_dpd_->buf4_print(&Tnew, "outfile", 1);
        global_dpd_->buf4_close(&Tnew);

        /********************************************************************************************/
        /************************** Beta-Beta spin case *********************************************/
        /********************************************************************************************/
        // Build T(ia,jb)
        // T_ij^ab(2) = \sum_{m,e} T_im^ae(1) W_mbej(1) => T(ia,jb)(2) = \sum_{m,e} T(ia,me) W(me,jb)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                               "T2_2 (ia|jb)");
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                               "T2_1 (ov|ov)");
        global_dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                               "W_1 (ov|ov)");
        global_dpd_->contract444(&Tp, &W, &T, 0, 0, 1.0, 0.0);
        // T_ij^ab(2) += \sum_{m,e} T_jm^be(1) W_maei(1) => T(ia,jb)(2) = \sum_{m,e} T(jb,me) W(me,ia)
        global_dpd_->contract444(&W, &Tp, &T, 0, 0, 1.0, 1.0);
        global_dpd_->buf4_close(&W);
        global_dpd_->buf4_close(&Tp);

        // T_ij^ab(2) += \sum_{M,E} T_Mi^Ea(1) W_MbEj(1) => T(ia,jb)(2) += \sum_{M,E} T(ME,ia) W(ME,jb)
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                               "T2_1 (OV|ov)");
        global_dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                               "W_1 (OV|ov)");
        global_dpd_->contract444(&Tp, &W, &T, 1, 1, 1.0, 1.0);
        // T_ij^ab(2) += \sum_{M,E} T_Mj^Eb(1) W_MaEi(1) => T(ia,jb)(2) += \sum_{M,E} T(ME,jb) W(ME,ia)
        global_dpd_->contract444(&W, &Tp, &T, 1, 1, 1.0, 1.0);
        global_dpd_->buf4_close(&W);
        global_dpd_->buf4_close(&Tp);

        // T(ia,jb) => T_ij^ab(2)
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, prqs, ID("[o,o]"), ID("[v,v]"), "T2_2 <ij|ab>");
        global_dpd_->buf4_close(&T);

        // Build T(ja,ib)
        // T_ij^ab(2) = -\sum_{m,e} T_jm^ae(1) W_mbei(1) => T(ja,ib)(2) = -\sum_{m,e} T(ja,me) W(me,ib)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                               "T2_2 (ja|ib)");
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                               "T2_1 (ov|ov)");
        global_dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                               "W_1 (ov|ov)");
        global_dpd_->contract444(&Tp, &W, &T, 0, 0, -1.0, 0.0);
        // T_ij^ab(2) = -\sum_{m,e} T_im^be(1) W_maej(1) => T(ja,ib)(2) = -\sum_{m,e} T(ib,me) W(me,ja)
        global_dpd_->contract444(&W, &Tp, &T, 0, 0, -1.0, 1.0);
        global_dpd_->buf4_close(&W);
        global_dpd_->buf4_close(&Tp);

        // T_ij^ab(2) = -\sum_{M,E} T_Mj^Ea(1) W_MbEi(1) => T(ja,ib)(2) -= \sum_{M,E} T(ME,ja) W(ME,ib)
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                               "T2_1 (OV|ov)");
        global_dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                               "W_1 (OV|ov)");
        global_dpd_->contract444(&Tp, &W, &T, 1, 1, -1.0, 1.0);
        // T_ij^ab(2) = -\sum_{M,E} T_Mi^Eb(1) W_MaEj(1) => T(ja,ib)(2) -= \sum_{M,E} T(ME,ib) W(ME,ja)
        global_dpd_->contract444(&W, &Tp, &T, 1, 1, -1.0, 1.0);
        global_dpd_->buf4_close(&W);
        global_dpd_->buf4_close(&Tp);

        // T(ja,ib) => T_ij^ab(2)
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, rpqs, ID("[o,o]"), ID("[v,v]"), "T2_2 (ij|ab)");
        global_dpd_->buf4_close(&T);

        // Build T2BBnew
        // T_ij^ab(2) = T(ia,jb)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "T2_2 <ij|ab>");
        global_dpd_->buf4_copy(&T, PSIF_OCC_DPD, "T2_2new <oo|vv>");
        global_dpd_->buf4_close(&T);

        // T_ij^ab(2) += T(ja,ib)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "T2_2new <oo|vv>");
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "T2_2 (ij|ab)");
        global_dpd_->buf4_axpy(&Tp, &T, 1.0);  // 1.0*Tp + T -> T
        global_dpd_->buf4_close(&Tp);

        // T_ij^ab(2) += 1/2 \sum_{m,n} T_mn^ab(1) W_mnij(1)
        global_dpd_->buf4_init(&TBB, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "T2_1 <oo|vv>");
        // dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[o,o]"),
        //              ID("[o,o]"), ID("[o,o]"), 0, "W_1 <oo|oo>");
        global_dpd_->buf4_init(&W, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"), ID("[o,o]"), ID("[o,o]"), 0,
                               "MO Ints <oo||oo>");
        global_dpd_->contract444(&W, &TBB, &T, 1, 1, 0.5, 1.0);
        global_dpd_->buf4_close(&W);
        global_dpd_->buf4_close(&T);

        /*
        // T_ij^ab(2) += 1/2 \sum_{e,f} T_ij^ef(1) W_abef(1)
        //dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[v,v]"), ID("[v,v]"),
        //              ID("[v,v]"), ID("[v,v]"), 0, "W_1 <vv|vv>");
        dpd_buf4_init(&W, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                      ID("[v,v]"), ID("[v,v]"), 0, "MO Ints <vv||vv>");
        dpd_contract444(&TBB, &W, &T, 0, 0, 0.5, 1.0);
        dpd_buf4_close(&W);
        dpd_buf4_close(&T);
        dpd_buf4_close(&TBB);
        */

        /*
        // T_ij^ab(2) += 1/2 \sum_{e,f} T_ij^ef(1) <ef||ab> = \sum_{e,f} T_ij^ef(1) <ab|ef>
        dpd_buf4_init(&W, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                      ID("[v,v]"), ID("[v,v]"), 0, "MO Ints <vv|vv>");
        dpd_contract444(&TBB, &W, &T, 0, 0, 1.0, 1.0);
        dpd_buf4_close(&W);
        dpd_buf4_close(&T);
        dpd_buf4_close(&TBB);
        */

        // T_ij^ab(2) += 1/2 \sum_{e,f} T_ij^ef(1) <ef||ab> = \sum_{e,f} T_ij^ef(1) <ab|ef>
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), 0,
                               "Z2_2 <vv|oo>");
        global_dpd_->buf4_init(&W, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"), ID("[v,v]"), ID("[v,v]"), 0,
                               "MO Ints <vv|vv>");
        global_dpd_->contract444(&W, &TBB, &T, 0, 0, 1.0, 0.0);
        global_dpd_->buf4_close(&W);
        global_dpd_->buf4_close(&TBB);
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, rspq, ID("[o,o]"), ID("[v,v]"), "Z2_2 <oo|vv>");
        global_dpd_->buf4_close(&T);
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "T2_2new <oo|vv>");
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "Z2_2 <oo|vv>");
        global_dpd_->buf4_axpy(&Tp, &T, 1.0);  // 1.0*Tp + T -> T
        global_dpd_->buf4_close(&T);
        global_dpd_->buf4_close(&Tp);

        // initialize Tnew and Told
        global_dpd_->buf4_init(&Tnew, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "T2_2new <oo|vv>");
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "T2_2 <oo|vv>");

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
        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "D <oo|vv>");
        global_dpd_->buf4_dirprd(&D, &Tnew);
        global_dpd_->buf4_close(&D);
        if (print_ > 2) global_dpd_->buf4_print(&Tnew, "outfile", 1);
        global_dpd_->buf4_close(&Tnew);

        /********************************************************************************************/
        /************************** Alpha-Beta spin case ********************************************/
        /********************************************************************************************/
        // Build T(IA,jb)
        /*
        // T_Ij^Ab(2) = \sum_{M,E} T_IM^AE(1) W_MbEj(1) => T(IA,jb)(2) = \sum_{M,E} T(IA,ME) W(ME,jb)
        dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                      ID("[O,V]"), ID("[o,v]"), 0, "T2_2 (IA|jb)");
        dpd_buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "T2_1 (OV|OV)");
        dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                      ID("[O,V]"), ID("[o,v]"), 0, "W_1 (OV|ov)");
        dpd_contract444(&Tp, &W, &T, 0, 1, 1.0, 0.0);
        dpd_buf4_close(&Tp);
        */

        // T_Ij^Ab(2) = \sum_{M,E} T_IM^AE(1) W_MbEj(1) => T(IA,jb)(2) = \sum_{M,E} T(IA,ME) W(jb,ME)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                               "T2_2 (IA|jb)");
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "T2_1 (OV|OV)");
        global_dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), 0,
                               "W_1 (ov|OV)");
        global_dpd_->contract444(&Tp, &W, &T, 0, 0, 1.0, 0.0);
        global_dpd_->buf4_close(&Tp);
        global_dpd_->buf4_close(&W);

        // T_Ij^Ab(2) += \sum_{m,e} T_jm^be(1) W_IeAm(1) => T(IA,jb)(2) = \sum_{m,e} W(IA,me) T(jb,me)
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                               "T2_1 (ov|ov)");
        global_dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                               "W_1 (OV|ov)");
        global_dpd_->contract444(&W, &Tp, &T, 0, 0, 1.0, 1.0);
        global_dpd_->buf4_close(&W);
        global_dpd_->buf4_close(&Tp);

        // T_Ij^Ab(2) += \sum_{m,e} T_Im^Ae(1) W_mbej(1) => T(IA,jb)(2) += \sum_{m,e} T(IA,me) W(me,jb)
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                               "T2_1 (OV|ov)");
        global_dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                               "W_1 (ov|ov)");
        global_dpd_->contract444(&Tp, &W, &T, 0, 0, 1.0, 1.0);
        global_dpd_->buf4_close(&W);

        // T_Ij^Ab(2) += \sum_{M,E} T_Mj^Eb(1) W_MAEI(1) => T(IA,jb)(2) += \sum_{M,E} W(ME,IA) T(ME,jb)
        global_dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "W_1 (OV|OV)");
        global_dpd_->contract444(&W, &Tp, &T, 1, 1, 1.0, 1.0);
        global_dpd_->buf4_close(&W);
        global_dpd_->buf4_close(&Tp);

        // T(IA,jb) => T_Ij^Ab(2)
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, prqs, ID("[O,o]"), ID("[V,v]"), "T2_2 <Ij|Ab>");
        global_dpd_->buf4_close(&T);

        // Build T(jA,Ib)
        // T_Ij^Ab(2) = \sum_{M,e} T_Mj^Ae(1) W_MbeI(1) => T(jA,Ib)(2) = \sum_{M,e} T(jA,Me) W(Me,Ib)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,V]"), ID("[O,v]"), ID("[o,V]"), ID("[O,v]"), 0,
                               "T2_2 (jA|Ib)");
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[o,V]"), ID("[O,v]"), ID("[o,V]"), ID("[O,v]"), 0,
                               "T2_1 (oV|Ov)");
        global_dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), 0,
                               "W_1 (Ov|Ov)");
        global_dpd_->contract444(&Tp, &W, &T, 0, 0, 1.0, 0.0);
        global_dpd_->buf4_close(&W);

        // T_Ij^Ab(2) = +\sum_{m,E} T_Im^Eb(1) W_mAEj(1) => T(jA,Ib)(2) = +\sum_{m,E} W(mE,jA) T(mE,Ib)
        global_dpd_->buf4_init(&W, PSIF_OCC_DPD, 0, ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), 0,
                               "W_1 (oV|oV)");
        global_dpd_->contract444(&W, &Tp, &T, 1, 1, 1.0, 1.0);
        global_dpd_->buf4_close(&W);
        global_dpd_->buf4_close(&Tp);

        // T(jA,Ib) => T_Ij^Ab(2)
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, rpqs, ID("[O,o]"), ID("[V,v]"), "T2_2 (Ij|Ab)");
        global_dpd_->buf4_close(&T);

        // Build T2ABnew
        // T_Ij^Ab(2) = T(IA,jb)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "T2_2 <Ij|Ab>");
        global_dpd_->buf4_copy(&T, PSIF_OCC_DPD, "T2_2new <Oo|Vv>");
        global_dpd_->buf4_close(&T);

        // T_Ij^Ab(2) += T(jA,Ib)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "T2_2new <Oo|Vv>");
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "T2_2 (Ij|Ab)");
        global_dpd_->buf4_axpy(&Tp, &T, 1.0);  // 1.0*Tp + T -> T
        global_dpd_->buf4_close(&Tp);

        // T_Ij^Ab(2) += \sum_{M,n} T_Mn^Ab(1) W_MnIj(1) = \sum_{M,n} W(Mn,Ij) T(Mn,Ab)
        global_dpd_->buf4_init(&TAB, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "T2_1 <Oo|Vv>");
        // dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[O,o]"),
        //              ID("[O,o]"), ID("[O,o]"), 0, "W_1 <Oo|Oo>");
        global_dpd_->buf4_init(&W, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), 0,
                               "MO Ints <Oo|Oo>");
        global_dpd_->contract444(&W, &TAB, &T, 1, 1, 1.0, 1.0);
        global_dpd_->buf4_close(&W);
        global_dpd_->buf4_close(&T);

        /*
       // T_Ij^Ab(2) +=  \sum_{E,f} T_Ij^Ef(1) W_AbEf(1) =  \sum_{E,f} T(Ij,Ef) W(Ab,Ef)
       //dpd_buf4_init(&W, PSIF_OCC_DPD, 0, ID("[V,v]"), ID("[V,v]"),
       //              ID("[V,v]"), ID("[V,v]"), 0, "W_1 <Vv|Vv>");
       dpd_buf4_init(&W, PSIF_LIBTRANS_DPD, 0, ID("[V,v]"), ID("[V,v]"),
                     ID("[V,v]"), ID("[V,v]"), 0, "MO Ints <Vv|Vv>");
       dpd_contract444(&TAB, &W, &T, 0, 0, 1.0, 1.0);
       dpd_buf4_close(&W);
       dpd_buf4_close(&T);
       dpd_buf4_close(&TAB);
       */

        // T_Ij^Ab(2) +=  \sum_{E,f} T_Ij^Ef(1) W_AbEf(1) =  \sum_{E,f} T(Ij,Ef) W(Ab,Ef)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), 0,
                               "Z2_2 <Vv|Oo>");
        global_dpd_->buf4_init(&W, PSIF_LIBTRANS_DPD, 0, ID("[V,v]"), ID("[V,v]"), ID("[V,v]"), ID("[V,v]"), 0,
                               "MO Ints <Vv|Vv>");
        global_dpd_->contract444(&W, &TAB, &T, 0, 0, 1.0, 0.0);
        global_dpd_->buf4_close(&W);
        global_dpd_->buf4_close(&TAB);
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, rspq, ID("[O,o]"), ID("[V,v]"), "Z2_2 <Oo|Vv>");
        global_dpd_->buf4_close(&T);
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "T2_2new <Oo|Vv>");
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "Z2_2 <Oo|Vv>");
        global_dpd_->buf4_axpy(&Tp, &T, 1.0);  // 1.0*Tp + T -> T
        global_dpd_->buf4_close(&T);
        global_dpd_->buf4_close(&Tp);

        // initialize Tnew and Told
        global_dpd_->buf4_init(&Tnew, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "T2_2new <Oo|Vv>");
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "T2_2 <Oo|Vv>");

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
        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "D <Oo|Vv>");
        global_dpd_->buf4_dirprd(&D, &Tnew);
        global_dpd_->buf4_close(&D);
        if (print_ > 2) global_dpd_->buf4_print(&Tnew, "outfile", 1);
        global_dpd_->buf4_close(&Tnew);

        /********************************************************************************************/
        /************************** Compute amplitude residual to Check Convergence *****************/
        /********************************************************************************************/
        // Alpha-Alpha spin case
        global_dpd_->buf4_init(&Tnew, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_2new <OO|VV>");
        global_dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "RT2_2 <OO|VV>");
        global_dpd_->buf4_init(&R, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "RT2_2 <OO|VV>");
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_2 <OO|VV>");
        global_dpd_->buf4_axpy(&T, &R, -1.0);  // -1.0*T + R -> R
        global_dpd_->buf4_close(&T);

        nElements = 0;
        for (int h = 0; h < nirrep_; h++) nElements += R.params->coltot[h] * R.params->rowtot[h];
        rms_t2AA = 0.0;
        rms_t2AA = global_dpd_->buf4_dot_self(&R);
        global_dpd_->buf4_close(&R);
        rms_t2AA = std::sqrt(rms_t2AA) / nElements;

        // Reset
        global_dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "T2_2 <OO|VV>");
        if (print_ > 1) global_dpd_->buf4_print(&Tnew, "outfile", 1);
        global_dpd_->buf4_close(&Tnew);

        // Beta-Beta spin case
        global_dpd_->buf4_init(&Tnew, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "T2_2new <oo|vv>");
        global_dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "RT2_2 <oo|vv>");
        global_dpd_->buf4_init(&R, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "RT2_2 <oo|vv>");
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "T2_2 <oo|vv>");
        global_dpd_->buf4_axpy(&T, &R, -1.0);  // -1.0*T + R -> R
        global_dpd_->buf4_close(&T);

        nElements = 0;
        for (int h = 0; h < nirrep_; h++) nElements += R.params->coltot[h] * R.params->rowtot[h];
        rms_t2BB = 0.0;
        rms_t2BB = global_dpd_->buf4_dot_self(&R);
        global_dpd_->buf4_close(&R);
        rms_t2BB = std::sqrt(rms_t2BB) / nElements;

        // Reset
        global_dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "T2_2 <oo|vv>");
        if (print_ > 1) global_dpd_->buf4_print(&Tnew, "outfile", 1);
        global_dpd_->buf4_close(&Tnew);

        // Alpha-Beta spin case
        global_dpd_->buf4_init(&Tnew, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "T2_2new <Oo|Vv>");
        global_dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "RT2_2 <Oo|Vv>");
        global_dpd_->buf4_init(&R, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "RT2_2 <Oo|Vv>");
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "T2_2 <Oo|Vv>");
        global_dpd_->buf4_axpy(&T, &R, -1.0);  // -1.0*T + R -> R
        global_dpd_->buf4_close(&T);

        nElements = 0;
        for (int h = 0; h < nirrep_; h++) nElements += R.params->coltot[h] * R.params->rowtot[h];
        rms_t2AB = 0.0;
        rms_t2AB = global_dpd_->buf4_dot_self(&R);
        global_dpd_->buf4_close(&R);
        rms_t2AB = std::sqrt(rms_t2AA) / nElements;

        // Reset
        global_dpd_->buf4_copy(&Tnew, PSIF_OCC_DPD, "T2_2 <Oo|Vv>");
        if (print_ > 1) global_dpd_->buf4_print(&Tnew, "outfile", 1);
        global_dpd_->buf4_close(&Tnew);


        // close files
        psio_->close(PSIF_LIBTRANS_DPD, 1);
        psio_->close(PSIF_OCC_DPD, 1);

    }  // end if (reference_ == "UNRESTRICTED")
       // outfile->Printf("\n t2_2nd_general done. \n");

}  // end t2_2nd_general

// For iterative MP methods, construct post-diis quantities.
// These are the sum of amplitudes and resorted amplitudes.
// At present, only OMP2.5 and OMP3 use this method.
void OCCWave::iterative_mp_postdiis_amplitudes() {
    double omp2p5_factor = (wfn_type_ == "OMP2.5" ? 0.5 : 1);
    psio_open(PSIF_OCC_DPD, PSIO_OPEN_OLD);

    if (reference_ == "RESTRICTED") {
        dpdbuf4 T, Tau, Tp, Ttemp;
        /********************************************************************************************/
        /******************************* Chemist T2(1) Amplitudes ***********************************/
        /********************************************************************************************/
        // T_IJ^AB => T'(IA,JB), T"(JA,IB)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_1 <OO|VV>");
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, prqs, ID("[O,V]"), ID("[O,V]"), "T2_1 (OV|OV)");
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, qrps, ID("[O,V]"), ID("[O,V]"), "T2_1pp (OV|OV)");
        global_dpd_->buf4_close(&T);

        // Tau(IJ,AB) => Tau'(IA,JB), Tau"(JA,IB)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "Tau_1 <OO|VV>");
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, prqs, ID("[O,V]"), ID("[O,V]"), "Tau_1 (OV|OV)");
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, qrps, ID("[O,V]"), ID("[O,V]"), "Tau_1pp (OV|OV)");
        global_dpd_->buf4_close(&T);

        /********************************************************************************************/
        /************************** Sum up 1st & 2nd order amplitudes *******************************/
        /********************************************************************************************/
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_1 <OO|VV>");
        global_dpd_->buf4_copy(&T, PSIF_OCC_DPD, "T2 <OO|VV>");
        global_dpd_->buf4_close(&T);
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2 <OO|VV>");
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_2 <OO|VV>");
        global_dpd_->buf4_axpy(&Tp, &T, omp2p5_factor);  // 1.0*Tp + T -> T
        global_dpd_->buf4_close(&T);
        global_dpd_->buf4_close(&Tp);

        // Build Tau(ij,ab) = 2*T(ij,ab) - T(ji,ab)
        // Build TAA(ij,ab) = T(ij,ab) - T(ji,ab)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_2 <OO|VV>");
        global_dpd_->buf4_copy(&T, PSIF_OCC_DPD, "Tau_2 <OO|VV>");
        global_dpd_->buf4_copy(&T, PSIF_OCC_DPD, "T2_2AA <OO|VV>");
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "T2_2jiab <OO|VV>");
        global_dpd_->buf4_close(&T);
        global_dpd_->buf4_init(&Tau, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "Tau_2 <OO|VV>");
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_2AA <OO|VV>");
        global_dpd_->buf4_init(&Ttemp, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_2jiab <OO|VV>");
        global_dpd_->buf4_scm(&Tau, 2.0);
        global_dpd_->buf4_axpy(&Ttemp, &Tau, -1.0);  // -1.0*Ttemp + Tau -> Tau
        global_dpd_->buf4_axpy(&Ttemp, &Tp, -1.0);   // -1.0*Ttemp + Tp -> Tp
        global_dpd_->buf4_close(&Ttemp);
        global_dpd_->buf4_close(&Tp);
        global_dpd_->buf4_close(&Tau);

        // Build Tau(ij,ab) = 2*T(ij,ab) - T(ji,ab)
        // Build TAA(ij,ab) = T(ij,ab) - T(ji,ab)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2 <OO|VV>");
        global_dpd_->buf4_copy(&T, PSIF_OCC_DPD, "Tau <OO|VV>");
        global_dpd_->buf4_copy(&T, PSIF_OCC_DPD, "T2AA <OO|VV>");
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "T2jiab <OO|VV>");
        global_dpd_->buf4_close(&T);
        global_dpd_->buf4_init(&Tau, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "Tau <OO|VV>");
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2AA <OO|VV>");
        global_dpd_->buf4_init(&Ttemp, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2jiab <OO|VV>");
        global_dpd_->buf4_scm(&Tau, 2.0);
        global_dpd_->buf4_axpy(&Ttemp, &Tau, -1.0);  // -1.0*Ttemp + Tau -> Tau
        global_dpd_->buf4_axpy(&Ttemp, &Tp, -1.0);   // -1.0*Ttemp + Tp -> Tp
        global_dpd_->buf4_close(&Ttemp);
        global_dpd_->buf4_close(&Tp);
        global_dpd_->buf4_close(&Tau);
    } else {
        dpdbuf4 T, Tp;
        /********************************************************************************************/
        /******************************* Chemist T2(1) Amplitudes ***********************************/
        /********************************************************************************************/
        // T_IJ^AB => T(IA,JB)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_1 <OO|VV>");
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, prqs, ID("[O,V]"), ID("[O,V]"), "T2_1 (OV|OV)");
        global_dpd_->buf4_close(&T);

        // T_ij^ab => T(ia,jb)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "T2_1 <oo|vv>");
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, prqs, ID("[o,v]"), ID("[o,v]"), "T2_1 (ov|ov)");
        global_dpd_->buf4_close(&T);

        // T_Ij^Ab => T(IA,jb), T(jA,Ib)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "T2_1 <Oo|Vv>");
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, prqs, ID("[O,V]"), ID("[o,v]"), "T2_1 (OV|ov)");
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, qrps, ID("[o,V]"), ID("[O,v]"), "T2_1 (oV|Ov)");
        global_dpd_->buf4_close(&T);

        // T(IA,jb) => T(jb,IA)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                               "T2_1 (OV|ov)");
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, rspq, ID("[o,v]"), ID("[O,V]"), "T2_1 (ov|OV)");
        global_dpd_->buf4_close(&T);

        /********************************************************************************************/
        /************************** Sum up 1st & 2nd order amplitudes *******************************/
        /********************************************************************************************/
        // Build T2 = T2(1) + T2(2)
        // Alpha-Alpha spin case
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_1 <OO|VV>");
        global_dpd_->buf4_copy(&T, PSIF_OCC_DPD, "T2 <OO|VV>");
        global_dpd_->buf4_close(&T);
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2 <OO|VV>");
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_2 <OO|VV>");
        global_dpd_->buf4_axpy(&Tp, &T, omp2p5_factor);  // 1.0*Tp + T -> T
        global_dpd_->buf4_close(&T);
        global_dpd_->buf4_close(&Tp);

        // Beta-Beta spin case
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "T2_1 <oo|vv>");
        global_dpd_->buf4_copy(&T, PSIF_OCC_DPD, "T2 <oo|vv>");
        global_dpd_->buf4_close(&T);
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "T2 <oo|vv>");
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "T2_2 <oo|vv>");
        global_dpd_->buf4_axpy(&Tp, &T, omp2p5_factor);  // 1.0*Tp + T -> T
        global_dpd_->buf4_close(&T);
        global_dpd_->buf4_close(&Tp);

        // Alpha-Beta spin case
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "T2_1 <Oo|Vv>");
        global_dpd_->buf4_copy(&T, PSIF_OCC_DPD, "T2 <Oo|Vv>");
        global_dpd_->buf4_close(&T);
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "T2 <Oo|Vv>");
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "T2_2 <Oo|Vv>");
        global_dpd_->buf4_axpy(&Tp, &T, omp2p5_factor);  // 1.0*Tp + T -> T
        global_dpd_->buf4_close(&T);
        global_dpd_->buf4_close(&Tp);
    }
    psio_->close(PSIF_OCC_DPD, 1);
}
}
}  // End Namespaces
