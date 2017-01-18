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

/* This code includes correlation TPDMs. */

#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "occwave.h"
#include "defines.h"


using namespace std;

namespace psi{ namespace occwave{

void OCCWave::omp2_tpdm_oovv()
{
    dpdbuf4 G, T, Tau;

    psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
    psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);

 if (reference_ == "RESTRICTED") {
    global_dpd_->buf4_init(&Tau, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tau <OO|VV>");
    global_dpd_->buf4_copy(&Tau, PSIF_OCC_DENSITY, "TPDM <OO|VV>");
    global_dpd_->buf4_close(&Tau);

    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
    global_dpd_->buf4_scm(&G, 0.25);
    global_dpd_->buf4_close(&G);

 }// end RHF

 else if (reference_ == "UNRESTRICTED") {
    // Alpha-Alpha spin case
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    global_dpd_->buf4_copy(&T, PSIF_OCC_DENSITY, "TPDM <OO|VV>");
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
    global_dpd_->buf4_scm(&G, 0.25);
    global_dpd_->buf4_close(&G);

    // Beta-Beta spin case
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
    global_dpd_->buf4_copy(&T, PSIF_OCC_DENSITY, "TPDM <oo|vv>");
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "TPDM <oo|vv>");
    global_dpd_->buf4_scm(&G, 0.25);
    global_dpd_->buf4_close(&G);

    // Alpha-Beta spin case
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
    global_dpd_->buf4_copy(&T, PSIF_OCC_DENSITY, "TPDM <Oo|Vv>");
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
    global_dpd_->buf4_scm(&G, 0.25);
    global_dpd_->buf4_close(&G);

 }// end UHF

    psio_->close(PSIF_OCC_DPD, 1);
    psio_->close(PSIF_OCC_DENSITY, 1);

} // end of omp2_tpdm_oovv


/*=======================*/
/*  tpdm_oovv()          */
/*=======================*/
void OCCWave::tpdm_oovv()
{

    dpdbuf4 G, T, V, Tau;

    psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
    psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);

 if (reference_ == "RESTRICTED") {
    // G (IJ,AB) = 1/4 V(IB,AJ)
    global_dpd_->buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[V,O]"),
                  ID("[O,V]"), ID("[V,O]"), 0, "V <OV|VO>");
    global_dpd_->buf4_sort(&V, PSIF_OCC_DENSITY , psrq, ID("[O,O]"), ID("[V,V]"), "TPDM <OO|VV>");
    global_dpd_->buf4_close(&V);


    // For OMP2.5 G(IJ,AB) = 1/8 V(IB,AJ)
    if (wfn_type_ == "OMP2.5") {
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
    global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_close(&G);
    }

    // G (IJ,AB) += 1/4 (2T_IJ^AB - T_JI^AB)
    global_dpd_->buf4_init(&Tau, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tau <OO|VV>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
    global_dpd_->buf4_axpy(&Tau, &G, 1.0); // 1.0*Tau + G -> G
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_scm(&G, 0.25);
    global_dpd_->buf4_close(&G);
 }// end if (reference_ == "RESTRICTED")

 else if (reference_ == "UNRESTRICTED") {
    // Alpha-Alpha spin case
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
    global_dpd_->buf4_copy(&T, PSIF_OCC_DENSITY, "TPDM <OO|VV>");
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
    global_dpd_->buf4_scm(&G, 0.25);
    global_dpd_->buf4_close(&G);

    // Beta-Beta spin case
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
    global_dpd_->buf4_copy(&T, PSIF_OCC_DENSITY, "TPDM <oo|vv>");
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "TPDM <oo|vv>");
    global_dpd_->buf4_scm(&G, 0.25);
    global_dpd_->buf4_close(&G);

    // Alpha-Beta spin case
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
    global_dpd_->buf4_copy(&T, PSIF_OCC_DENSITY, "TPDM <Oo|Vv>");
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
    global_dpd_->buf4_scm(&G, 0.25);
    global_dpd_->buf4_close(&G);

    //Print
    if (print_ > 3) {
      global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
      global_dpd_->buf4_print(&G, "outfile", 1);
      global_dpd_->buf4_close(&G);

      global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "TPDM <oo|vv>");
      global_dpd_->buf4_print(&G, "outfile", 1);
      global_dpd_->buf4_close(&G);

      global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
      global_dpd_->buf4_print(&G, "outfile", 1);
      global_dpd_->buf4_close(&G);
    }

 }// end if (reference_ == "UNRESTRICTED")

    psio_->close(PSIF_OCC_DPD, 1);
    psio_->close(PSIF_OCC_DENSITY, 1);

} // end of twopdm_oovv



/*=======================*/
/*  tpdm_oooo()        */
/*=======================*/
void OCCWave::tpdm_oooo()
{
    dpdbuf4 G, T, V;

    psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);

 if (reference_ == "RESTRICTED") {
    // G(IJ,KL) = 1/8 (V_IJKL +  V_ILKJ)
    global_dpd_->buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "V <OO|OO>");
    global_dpd_->buf4_sort(&V, PSIF_OCC_DENSITY , psrq, ID("[O,O]"), ID("[O,O]"), "TPDM <OO|OO>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
    global_dpd_->buf4_axpy(&V, &G, 1.0); // 1.0*V + G -> G
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_scm(&G, 0.125);
    // For OMP2.5 G(IJ,KL) = 1/16 (V_IJKL + V_ILKJ)
    if (wfn_type_ == "OMP2.5") global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_close(&G);
 }// end if (reference_ == "RESTRICTED")

 else if (reference_ == "UNRESTRICTED") {
    // Alpha-Alpha spin case
    global_dpd_->buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "V <OO|OO>");
    global_dpd_->buf4_copy(&V, PSIF_OCC_DENSITY, "TPDM <OO|OO>");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
    global_dpd_->buf4_scm(&G, 0.25);
    // For OMP2.5 G(IJ,KL) = 1/8 (V_IJKL + V_ILKJ)
    if (wfn_type_ == "OMP2.5") global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_close(&G);

    // Beta-Beta spin case
    global_dpd_->buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "V <oo|oo>");
    global_dpd_->buf4_copy(&V, PSIF_OCC_DENSITY, "TPDM <oo|oo>");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "TPDM <oo|oo>");
    global_dpd_->buf4_scm(&G, 0.25);
    // For OMP2.5 G(IJ,KL) = 1/8 (V_IJKL + V_ILKJ)
    if (wfn_type_ == "OMP2.5") global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_close(&G);

    // Alpha-Beta spin case
    global_dpd_->buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "V <Oo|Oo>");
    global_dpd_->buf4_copy(&V, PSIF_OCC_DENSITY, "TPDM <Oo|Oo>");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "TPDM <Oo|Oo>");
    global_dpd_->buf4_scm(&G, 0.25);
    // For OMP2.5 G(IJ,KL) = 1/8 (V_IJKL + V_ILKJ)
    if (wfn_type_ == "OMP2.5") global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_close(&G);

 }// end if (reference_ == "UNRESTRICTED")

    psio_->close(PSIF_OCC_DENSITY, 1);

} // end of twopdm_oooo


/*=======================*/
/*  omp3_tpdm_vvvv()     */
/*=======================*/
void OCCWave::omp3_tpdm_vvvv()
{
    // NOTE: contract444 can handle only TN and NT type contractions, which means (0,0) and (1,1) type target indices,
    //  with out-of-core algorithm!!!!
    dpdbuf4  T, L, G, V;

    psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
    psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);

 if (reference_ == "RESTRICTED") {
    // NOTE: A VVVV-sort takes too long time, and there is no OOC for psrq type sorting.
    // Hence I will not symmetrize the TPDM VVVV-block.
    // However, in for analytical gradients I need to symmetrize it.

    if (time4grad == 0) {
    // G_ABCD(2) = 1/2\sum_{M,N} T_MN^CD(1) (2T_MN^AB(1) - T_MN^BA(1))
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    global_dpd_->buf4_init(&L, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tau_1 <OO|VV>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
    global_dpd_->contract444(&L, &T, &G, 1, 1, 0.5, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&G);

    // For OMP2.5 G_ABCD(2) = 1/4\sum_{M,N} T_MN^CD(1) (2T_MN^AB(1) - T_MN^BA(1))
    if (wfn_type_ == "OMP2.5") {
        global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                      ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
        global_dpd_->buf4_scm(&G, 0.5);
        global_dpd_->buf4_close(&G);
    }
    }// end if (time4grad == 0)

    else if (time4grad == 1) {
    // Symmetrize the VVVV block
    // buf_axpy assume that half of each buffer can be stored incore, hence it is a problem
    // for 2 VVVV type matrices.
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    global_dpd_->buf4_init(&L, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tau_1 <OO|VV>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <AD|CB>");
    // G_ABCD(2) = 1/4\sum_{M,N} T_MN^CB(1) (2T_MN^AD(1) - T_MN^DA(1))
    global_dpd_->contract444(&L, &T, &G, 1, 1, 0.25, 0.0);
    global_dpd_->buf4_sort(&G, PSIF_OCC_DENSITY, prsq, ID("[V,V]"), ID("[V,V]"), "TPDM <AC|BD>");
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <AC|BD>");
    global_dpd_->buf4_sort(&G, PSIF_OCC_DENSITY, prqs, ID("[V,V]"), ID("[V,V]"), "TPDM <VV|VV>");
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
              ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
    // G_ABCD(2) += 1/4\sum_{M,N} T_MN^CD(1) (2T_MN^AB(1) - T_MN^BA(1))
    global_dpd_->contract444(&L, &T, &G, 1, 1, 0.25, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);

    if (wfn_type_ == "OMP2.5") {
        global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                      ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
        global_dpd_->buf4_scm(&G, 0.5);
        global_dpd_->buf4_close(&G);
    }
    } // end else if (time4grad == 1)

    //Print
    if (print_ > 3) {
      global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
      global_dpd_->buf4_print(&G, "outfile", 1);
      global_dpd_->buf4_close(&G);
    }
 }// end if (reference_ == "RESTRICTED")

 else if (reference_ == "UNRESTRICTED") {
    // Alpha-Alpha spin case
    // G_ABCD(2) = 1/8 \sum_{M,N} T_MN^CD(1) L_AB^MN(1) = 1/8 \sum_{M,N} T_MN^AB(1) T_MN^CD(1)
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    global_dpd_->buf4_init(&L, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
    global_dpd_->contract444(&L, &T, &G, 1, 1, 0.125, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&G);

    // For OMP2.5 G_ABCD(2) = 1/16 \sum_{M,N} T_MN^CD(1) L_AB^MN(1) = 1/16 \sum_{M,N} T_MN^AB(1) T_MN^CD(1)
    if (wfn_type_ == "OMP2.5") {
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
    global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_close(&G);
    }

    // Beta-Beta spin case
    // G_abcd(2) = 1/8 \sum_{m,n} T_mn^cd(1) L_ab^mn(1) = 1/8 \sum_{m,n} T_mn^ab(1) T_mn^cd(1)
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
    global_dpd_->buf4_init(&L, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "TPDM <vv|vv>");
    global_dpd_->contract444(&L, &T, &G, 1, 1, 0.125, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&G);

    // OMP2.5: G_abcd(2) = 1/16 \sum_{m,n} T_mn^cd(1) L_ab^mn(1) = 1/16 \sum_{m,n} T_mn^ab(1) T_mn^cd(1)
    if (wfn_type_ == "OMP2.5") {
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "TPDM <vv|vv>");
    global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_close(&G);
    }


    // Alpha-Beta spin case
    // G_AbCd(2) = 1/4 \sum_{M,n} T_Mn^Cd(1) L_Ab^Mn(1) = 1/4 \sum_{M,n} T_Mn^Ab(1) T_Mn^Cd(1)
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
    global_dpd_->buf4_init(&L, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "TPDM <Vv|Vv>");
    global_dpd_->contract444(&L, &T, &G, 1, 1, 0.25, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&G);

    // OMP2.5: G_AbCd(2) = 1/8 \sum_{M,n} T_Mn^Cd(1) L_Ab^Mn(1) = 1/8 \sum_{M,n} T_Mn^Ab(1) T_Mn^Cd(1)
    if (wfn_type_ == "OMP2.5") {
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "TPDM <Vv|Vv>");
    global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_close(&G);
    }


    //Print
    if (print_ > 3) {
      global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
      global_dpd_->buf4_print(&G, "outfile", 1);
      global_dpd_->buf4_close(&G);

      global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "TPDM <vv|vv>");
      global_dpd_->buf4_print(&G, "outfile", 1);
      global_dpd_->buf4_close(&G);

      global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "TPDM <Vv|Vv>");
      global_dpd_->buf4_print(&G, "outfile", 1);
      global_dpd_->buf4_close(&G);
    }
 }// end if (reference_ == "UNRESTRICTED")

    psio_->close(PSIF_OCC_DENSITY, 1);
    psio_->close(PSIF_OCC_DPD, 1);

} // end of omp3_tpdm_vvvv

/*=======================*/
/*  ocepa_tpdm_vvvv()    */
/*=======================*/
void OCCWave::ocepa_tpdm_vvvv()
{
    //outfile->Printf("\n ocepa_tpdm_vvvv is starting... \n");
    // NOTE: contract444 can handle only TN and NT type contractions, which means (0,0) and (1,1) type target indices,
    //  with out-of-core algorithm!!!!
    dpdbuf4  T, L, G, V;

    psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
    psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);

 if (reference_ == "RESTRICTED") {
    // NOTE: A VVVV-sort takes too long time, and there is no OOC for psrq type sorting.
    // Hence I will not symmetrize the TPDM VVVV-block.
    // However, in for analytical gradients I need to symmetrize it.

    if (time4grad == 0) {
    // G_ABCD(2) = 1/2\sum_{M,N} T_MN^CD(1) (2T_MN^AB(1) - T_MN^BA(1))
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
    global_dpd_->buf4_init(&L, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tau <OO|VV>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
    global_dpd_->contract444(&L, &T, &G, 1, 1, 0.5, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&G);
    }// end if (time4grad == 0)

    else if (time4grad == 1) {
    // Symmetrize the VVVV block
    // buf_axpy assume that half of each buffer can be stored incore, hence it is a problem
    // for 2 VVVV type matrices.
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
    global_dpd_->buf4_init(&L, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tau <OO|VV>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <AD|CB>");
    // G_ABCD(2) = 1/4\sum_{M,N} T_MN^CB(1) (2T_MN^AD(1) - T_MN^DA(1))
    global_dpd_->contract444(&L, &T, &G, 1, 1, 0.25, 0.0);
    global_dpd_->buf4_sort(&G, PSIF_OCC_DENSITY, prsq, ID("[V,V]"), ID("[V,V]"), "TPDM <AC|BD>");
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <AC|BD>");
    global_dpd_->buf4_sort(&G, PSIF_OCC_DENSITY, prqs, ID("[V,V]"), ID("[V,V]"), "TPDM <VV|VV>");
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
              ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
    // G_ABCD(2) += 1/4\sum_{M,N} T_MN^CD(1) (2T_MN^AB(1) - T_MN^BA(1))
    global_dpd_->contract444(&L, &T, &G, 1, 1, 0.25, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    } // end else if (time4grad == 1)

    /*
    if (time4grad == 1) {
    // Symmetrize the VVVV block
    dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
              ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
    dpd_buf4_sort(&G, PSIF_OCC_DENSITY, prsq, ID("[V,V]"), ID("[V,V]"), "TPDM <AC|DB>");
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <AC|DB>");
    dpd_buf4_sort(&G, PSIF_OCC_DENSITY, prqs, ID("[V,V]"), ID("[V,V]"), "TPDM <AD|CB>");
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
    dpd_buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <AD|CB>");
    dpd_buf4_axpy(&V, &G, 1.0); // 1.0*V + G -> G
    dpd_buf4_close(&V);
    dpd_buf4_scm(&G, 0.5);
    dpd_buf4_close(&G);
    } // enf if (time4grad == 1)
    */

    //Print
    if (print_ > 3) {
      global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
      global_dpd_->buf4_print(&G, "outfile", 1);
      global_dpd_->buf4_close(&G);
    }
 }// end if (reference_ == "RESTRICTED")

 else if (reference_ == "UNRESTRICTED") {
    // Alpha-Alpha spin case
    // G_ABCD(2) = 1/8 \sum_{M,N} T_MN^CD(1) L_AB^MN(1) = 1/8 \sum_{M,N} T_MN^AB(1) T_MN^CD(1)
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
    global_dpd_->buf4_init(&L, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
    global_dpd_->contract444(&L, &T, &G, 1, 1, 0.125, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&G);

    // Beta-Beta spin case
    // G_abcd(2) = 1/8 \sum_{m,n} T_mn^cd(1) L_ab^mn(1) = 1/8 \sum_{m,n} T_mn^ab(1) T_mn^cd(1)
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
    global_dpd_->buf4_init(&L, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "TPDM <vv|vv>");
    global_dpd_->contract444(&L, &T, &G, 1, 1, 0.125, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&G);


    // Alpha-Beta spin case
    // G_AbCd(2) = 1/4 \sum_{M,n} T_Mn^Cd(1) L_Ab^Mn(1) = 1/4 \sum_{M,n} T_Mn^Ab(1) T_Mn^Cd(1)
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
    global_dpd_->buf4_init(&L, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "TPDM <Vv|Vv>");
    global_dpd_->contract444(&L, &T, &G, 1, 1, 0.25, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&G);


    //Print
    if (print_ > 3) {
      global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
      global_dpd_->buf4_print(&G, "outfile", 1);
      global_dpd_->buf4_close(&G);

      global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "TPDM <vv|vv>");
      global_dpd_->buf4_print(&G, "outfile", 1);
      global_dpd_->buf4_close(&G);

      global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "TPDM <Vv|Vv>");
      global_dpd_->buf4_print(&G, "outfile", 1);
      global_dpd_->buf4_close(&G);
    }
 }// end if (reference_ == "UNRESTRICTED")

    psio_->close(PSIF_OCC_DENSITY, 1);
    psio_->close(PSIF_OCC_DPD, 1);
    //outfile->Printf("\n ocepa_tpdm_vvvv done. \n");

} // end of twopdm_vvvv



/*=======================*/
/*  tpdm_ovov()        */
/*=======================*/
void OCCWave::tpdm_ovov()
{
    dpdbuf4 G, T, V;

    psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);

 if (reference_ == "RESTRICTED") {
    // Fully-symmetric
    // G(IA,JB) = -1/4 (V_IAJB +  V_IBJA)
    global_dpd_->buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "V <OV|OV>");
    global_dpd_->buf4_sort(&V, PSIF_OCC_DENSITY , psrq, ID("[O,V]"), ID("[O,V]"), "TPDM <OV|OV>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
    global_dpd_->buf4_axpy(&V, &G, 1.0); // 1.0*V + G -> G
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_scm(&G, -0.25);
    // OMP2.5: G(IA,JB) = -1/8 (V_IAJB +  V_IBJA)
    if (wfn_type_ == "OMP2.5") global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_close(&G);

    /*
    // Partially-symmetric
    // G(IA,JB) = -1/2 V_IAJB
    dpd_buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "V <OV|OV>");
    dpd_buf4_copy(&V, PSIF_OCC_DENSITY , "TPDM <OV|OV>");
    dpd_buf4_close(&V);
    dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
    dpd_buf4_scm(&G, -0.5);
    dpd_buf4_close(&G);
    */
 }// end if (reference_ == "RESTRICTED")

 else if (reference_ == "UNRESTRICTED") {
    // Build G_IAJB
    global_dpd_->buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "V <OV|OV>");
    global_dpd_->buf4_copy(&V, PSIF_OCC_DENSITY, "TPDM <OV|OV>");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
    global_dpd_->buf4_scm(&G, -0.5);
    if (wfn_type_ == "OMP2.5") global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_close(&G);

    // Build G_iajb
    global_dpd_->buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "V <ov|ov>");
    global_dpd_->buf4_copy(&V, PSIF_OCC_DENSITY, "TPDM <ov|ov>");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "TPDM <ov|ov>");
    global_dpd_->buf4_scm(&G, -0.5);
    if (wfn_type_ == "OMP2.5") global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_close(&G);

    // Build G_IaJb
    global_dpd_->buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "V <Ov|Ov>");
    global_dpd_->buf4_copy(&V, PSIF_OCC_DENSITY, "TPDM <Ov|Ov>");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "TPDM <Ov|Ov>");
    global_dpd_->buf4_scm(&G, -0.5);
    if (wfn_type_ == "OMP2.5") global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_close(&G);

 }// end if (reference_ == "UNRESTRICTED")
    psio_->close(PSIF_OCC_DENSITY, 1);

} // end of twopdm_ovov


/*=======================*/
/*  tpdm_vovo()        */
/*=======================*/
void OCCWave::tpdm_vovo()
{
    dpdbuf4 G, T, V;

    psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);

    // G_AIBJ = G_IAJB
    // G_aibj = G_iajb

    // G_AiBj = -1/2 V_iAjB
    global_dpd_->buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "V <oV|oV>");
    global_dpd_->buf4_sort(&V, PSIF_OCC_DENSITY , qpsr, ID("[V,o]"), ID("[V,o]"), "TPDM <Vo|Vo>");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "TPDM <Vo|Vo>");
    global_dpd_->buf4_scm(&G, -0.5);
    if (wfn_type_ == "OMP2.5") global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_close(&G);

    psio_->close(PSIF_OCC_DENSITY, 1);

} // end of twopdm_vovo


/*=======================*/
/*  tpdm_ovvo()        */
/*=======================*/
void OCCWave::tpdm_ovvo()
{
    dpdbuf4 G, T, V;

    psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);

    // G_IABJ = -G_IAJB
    // G_iabj = -G_iajb

    // G_IaBj = 1/2 V_IajB
    global_dpd_->buf4_init(&V, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "V <Ov|oV>");
    global_dpd_->buf4_sort(&V, PSIF_OCC_DENSITY , pqsr, ID("[O,v]"), ID("[V,o]"), "TPDM <Ov|Vo>");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "TPDM <Ov|Vo>");
    global_dpd_->buf4_scm(&G, 0.5);
    if (wfn_type_ == "OMP2.5") global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_close(&G);

    // VoOv block is here!
    // G_AiJb = G_JbAi so I do not need to VoOv block, however for proper contraction in the GFock.cc I need it.
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "TPDM <Ov|Vo>");
    global_dpd_->buf4_sort(&G, PSIF_OCC_DENSITY , rspq, ID("[V,o]"), ID("[O,v]"), "TPDM <Vo|Ov>");
    global_dpd_->buf4_close(&G);


    //Print
    if (print_ > 3) {
      global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "TPDM <Ov|Vo>");
      global_dpd_->buf4_print(&G, "outfile", 1);
      global_dpd_->buf4_close(&G);
    }

    psio_->close(PSIF_OCC_DENSITY, 1);

} // end of twopdm_ovvo

}} // End Namespaces
