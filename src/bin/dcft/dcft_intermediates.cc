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

#include "dcft.h"
#include <libdpd/dpd.h>
#include <libpsio/psio.hpp>
#include <libtrans/integraltransform.h>
#include "defines.h"

namespace psi{ namespace dcft{

/**
 * Builds the intermediate tensors
 */
void
DCFTSolver::build_cumulant_intermediates()
{
    dcft_timer_on("DCFTSolver::build_intermediates()");

    compute_G_intermediate();

    if (exact_tau_) {
        form_density_weighted_fock();
    }

    compute_F_intermediate();

    // Add third-order N-representability terms if needed
    if (options_.get_str("DCFT_FUNCTIONAL") == "ODC-13") {
        compute_V_intermediate();
    }

    dcft_timer_off("DCFTSolver::build_intermediates()");
}

void
DCFTSolver::compute_G_intermediate() {

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdbuf4 I, L, G, T, Taa, Tab, Tbb, Laa, Lab, Lbb;

    /*
     * G_ijab = 1/2 Sum_cd gbar_cdab lambda_ijcd
     */
    if(options_.get_str("AO_BASIS") == "NONE"){
        // G_IJAB += 1/2 Sum_CD gbar_CDAB lambda_IJCD
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V>V]-"), ID("[V>V]-"),
                      ID("[V,V]"), ID("[V,V]"), 1, "MO Ints <VV|VV>");
        global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                      ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
        global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                      ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
        global_dpd_->contract444(&L, &I, &G, 0, 0, 1.0, 0.0);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&L);
        global_dpd_->buf4_close(&G);


        // G_IjAb += Sum_Cd g_CdAb lambda_IjCd
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,v]"), ID("[V,v]"),
                      ID("[V,v]"), ID("[V,v]"), 0, "MO Ints <Vv|Vv>");
        global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
        global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
        global_dpd_->contract444(&L, &I, &G, 0, 0, 1.0, 0.0);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&L);
        global_dpd_->buf4_close(&G);


        // G_ijab += 1/2 Sum_cd gbar_cdab lambda_ijcd
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v>v]-"), ID("[v>v]-"),
                      ID("[v,v]"), ID("[v,v]"), 1, "MO Ints <vv|vv>");
        global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                      ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
        global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                      ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
        global_dpd_->contract444(&L, &I, &G, 0, 0, 1.0, 0.0);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&L);
        global_dpd_->buf4_close(&G);
    }
    else{

        /***********AA***********/
        global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                      ID("[O,O]"), ID("[V,V]"), 0, "tau(temp) <OO|VV>");
        global_dpd_->buf4_copy(&L, PSIF_DCFT_DPD, "G <OO|VV>");
        global_dpd_->buf4_close(&L);

        /***********BB***********/
        global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                      ID("[o,o]"), ID("[v,v]"), 0, "tau(temp) <oo|vv>");
        global_dpd_->buf4_copy(&L, PSIF_DCFT_DPD, "G <oo|vv>");
        global_dpd_->buf4_close(&L);

        /***********AB***********/
        global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "tau(temp) <Oo|Vv>");
        global_dpd_->buf4_copy(&L, PSIF_DCFT_DPD, "G <Oo|Vv>");
        global_dpd_->buf4_close(&L);
    }

    /*
     * G_ijab += 1/2 Sum_kl gbar_ijkl lambda_klab
     */
    // G_IJAB += 1/2 Sum_KL gbar_IJKL lambda_KLAB
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[O>O]-"),
            ID("[O,O]"), ID("[O,O]"), 1, "MO Ints <OO|OO>");
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    global_dpd_->contract444(&I, &L, &G, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&G);

    // G_IjAb += Sum_Kl g_IjKl lambda_KlAb
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "MO Ints <Oo|Oo>");
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
    global_dpd_->contract444(&I, &L, &G, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&G);

    // G_ijab += 1/2 Sum_kl gbar_ijkl lambda_klab
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o>o]-"), ID("[o>o]-"),
                  ID("[o,o]"), ID("[o,o]"), 1, "MO Ints <oo|oo>");
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    global_dpd_->contract444(&I, &L, &G, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&G);


    /*
     * G_ijab -= P(ij)P(ab) Sum_kc gbar_jckb lambda_ikac
     */
    global_dpd_->buf4_init(&Taa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Temp (OV|OV)");
    global_dpd_->buf4_init(&Tab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Temp (Ov|oV)");
    global_dpd_->buf4_init(&Tbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Temp (ov|ov)");
    global_dpd_->buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    global_dpd_->buf4_sort(&Laa, PSIF_DCFT_DPD, prqs, ID("[O,V]"), ID("[O,V]"), "Lambda (OV|OV)");
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Lambda (OV|OV)");
    global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    global_dpd_->buf4_sort(&Lab, PSIF_DCFT_DPD, psqr, ID("[O,v]"), ID("[o,V]"), "Lambda (Ov|oV)");
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Lambda (Ov|oV)");

    global_dpd_->buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    global_dpd_->buf4_sort(&Lbb, PSIF_DCFT_DPD, prqs, ID("[o,v]"),ID("[o,v]"), "Lambda (ov|ov)");
    global_dpd_->buf4_close(&Lbb);
    global_dpd_->buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Lambda (ov|ov)");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "MO Ints <oV|oV>");

    // T_IbjA = -Sum_kC lambda_IbkC g_jAkC
    global_dpd_->contract444(&Lab, &I, &Tab, 0, 0, -1.0, 0.0);

    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "MO Ints <Ov|Ov>") ;

    // T_IbjA -= Sum_Kc g_IbKc lambda_KcjA
    global_dpd_->contract444(&I, &Lab, &Tab, 0, 1, -1.0, 1.0);

    // T_IbjA -> T_IAjb
    global_dpd_->buf4_sort(&Tab, PSIF_DCFT_DPD, psrq, ID("[O,V]"),ID("[o,v]"), "Temp (OV|ov)");
    global_dpd_->buf4_close(&Tab);
    global_dpd_->buf4_init(&Tab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Temp (OV|ov)");

    // Lambda_IbkC -> Lambda_ICkb
    global_dpd_->buf4_sort(&Lab, PSIF_DCFT_DPD, psrq, ID("[O,V]"),ID("[o,v]"), "Lambda (OV|ov)");
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");

    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");

    // T_IAJB = Sum_kc lambda_IAkc g_JBkc
    global_dpd_->contract444(&Lab, &I, &Taa, 0, 0, 1.0, 0.0);

    // T_iajb = Sum_KC lambda_KCia g_KCjb
    global_dpd_->contract444(&Lab, &I, &Tbb, 1, 1, 1.0, 0.0);

    // T_IAjb += Sum_kc g_IAkc lambda_jbkc
    global_dpd_->contract444(&I, &Lbb, &Tab, 0, 0, 1.0, 1.0);

    // T_IAjb += Sum_KC lambda_IAKC g_KCjb
    global_dpd_->contract444(&Laa, &I, &Tab, 0, 1, 1.0, 1.0);


    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>");

    // T_IAJB -= Sum_KC lambda_IAKC g_JBKC
    global_dpd_->contract444(&Laa, &I, &Taa, 0, 0, -1.0, 1.0);

    // T_IAjb -= Sum_KC g_IAKC lambda_KCjb
    global_dpd_->contract444(&I, &Lab, &Tab, 0, 1, -1.0, 1.0);

    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");

    // T_IAJB += Sum_KC lambda_IAKC (JB|KC)
    global_dpd_->contract444(&Laa, &I, &Taa, 0, 0, 1.0, 1.0);

    // T_IAjb += Sum_KC (JB|KC) lambda_KCjb
    global_dpd_->contract444(&I, &Lab, &Tab, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov|ov>");

    // T_iajb -= Sum_kc lambda_iakc g_jbkc
    global_dpd_->contract444(&Lbb, &I, &Tbb, 0, 0, -1.0, 1.0);

    // T_IAjb -= Sum_KC lambda_IAkc g_jbkc
    global_dpd_->contract444(&Lab, &I, &Tab, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov|ov)");

    // T_iajb += Sum_kc lambda_iakc (JB|KC)
    global_dpd_->contract444(&Lbb, &I, &Tbb, 0, 0, 1.0, 1.0);

    // T_IAjb += Sum_KC lambda_IAkc (kc|jb)
    global_dpd_->contract444(&Lab, &I, &Tab, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&Lbb);

    // T_IAJB -> T_IJAB
    global_dpd_->buf4_sort(&Taa, PSIF_DCFT_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "Temp <OO|VV>");
    global_dpd_->buf4_close(&Taa);
    // G_IJAB += T_IJAB
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    dpd_buf4_add(&G, &T, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&Taa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");

    // T_IJAB -> T_JIAB
    global_dpd_->buf4_sort(&Taa, PSIF_DCFT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");



    // G_IJAB -= T_JIAB
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    dpd_buf4_add(&G, &T, -1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&T);

    // T_IJAB -> T_IJBA
    global_dpd_->buf4_sort(&Taa, PSIF_DCFT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    // G_IJAB -= T_IJBA
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    dpd_buf4_add(&G, &T, -1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&T);

    // T_IJAB -> T_JIBA
    global_dpd_->buf4_sort(&Taa, PSIF_DCFT_DPD, qpsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    // G_IJAB += T_JIBA
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    dpd_buf4_add(&G, &T, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_close(&Taa);

    // T_IAjb -> T_IjAb
    global_dpd_->buf4_sort(&Tab, PSIF_DCFT_DPD, prqs, ID("[O,o]"), ID("[V,v]"), "Temp <Oo|Vv>");
    global_dpd_->buf4_close(&Tab);
    // G_IjAb += T_IjAb
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Temp <Oo|Vv>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
    dpd_buf4_add(&G, &T, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&T);

    // T_iajb -> T_ijab
    global_dpd_->buf4_sort(&Tbb, PSIF_DCFT_DPD, prqs, ID("[o,o]"), ID("[v,v]"), "Temp <oo|vv>");
    global_dpd_->buf4_close(&Tbb);
    // G_ijab += T_ijab
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    dpd_buf4_add(&G, &T, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&Tbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");

    // T_ijab -> T_jiab
    global_dpd_->buf4_sort(&Tbb, PSIF_DCFT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    // G_ijab -= T_jiab
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    dpd_buf4_add(&G, &T, -1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&T);

    // T_ijab -> T_ijba
    global_dpd_->buf4_sort(&Tbb, PSIF_DCFT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    // G_ijab -= T_ijba
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    dpd_buf4_add(&G, &T, -1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&T);

    // T_ijab -> T_jiba
    global_dpd_->buf4_sort(&Tbb, PSIF_DCFT_DPD, qpsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    // G_ijab += T_jiba
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    dpd_buf4_add(&G, &T, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_close(&Tbb);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

}

void
DCFTSolver::compute_F_intermediate() {

    dpdfile2 F_OO, F_oo, F_VV, F_vv;
    dpdbuf4 F, T, Laa, Lab, Lbb;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    /*
    * F_ijab += P(ab) F_ca lambda_ijcb - P(ij) F_ki lambda_jkab
    */
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    // Temp_IJAB = lambda_IJCB F_AC
    global_dpd_->file2_init(&F_VV, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");
    global_dpd_->contract244(&F_VV, &Laa, &T, 1, 2, 1, 1.0, 0.0);
    global_dpd_->file2_close(&F_VV);
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&T);
    // F_IJAB = Temp_IJAB
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_copy(&T, PSIF_DCFT_DPD, "F <OO|VV>");
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "F <OO|VV>");
    global_dpd_->buf4_sort(&T, PSIF_DCFT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    global_dpd_->buf4_close(&T);
    // F_IJAB -= Temp_IJBA
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_add(&F, &T, -1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    // Temp_IJAB = -lambda_KJAB F_IK
    global_dpd_->file2_init(&F_OO, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    global_dpd_->contract244(&F_OO, &Laa, &T, 1, 0, 0, -1.0, 0.0);

    global_dpd_->file2_close(&F_OO);
    global_dpd_->buf4_close(&Laa);
    // F_IJAB += Temp_IJAB
    dpd_buf4_add(&F, &T, 1.0);
    global_dpd_->buf4_sort(&T, PSIF_DCFT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    global_dpd_->buf4_close(&T);
    // F_IJAB -= Temp_JIAB
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_add(&F, &T, -1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "F <Oo|Vv>");
    global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    // F_IjAb += lambda_IjCb F_AC
    global_dpd_->file2_init(&F_VV, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");
    global_dpd_->contract244(&F_VV, &Lab, &F, 1, 2, 1, 1.0, 0.0);
    global_dpd_->file2_close(&F_VV);
    // F_IjAb += lambda_IjAc F_bc
    global_dpd_->file2_init(&F_vv, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "F <v|v>");
    global_dpd_->contract424(&Lab, &F_vv, &F, 3, 1, 0, 1.0, 1.0);
    global_dpd_->file2_close(&F_vv);
    // F_IjAb -= lambda_KjAb F_IK
    global_dpd_->file2_init(&F_OO, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    global_dpd_->contract244(&F_OO, &Lab, &F, 1, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&F_OO);
    // F_IjAb -= lambda_IkAb F_jk
    global_dpd_->file2_init(&F_oo, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "F <o|o>");
    global_dpd_->contract424(&Lab, &F_oo, &F, 1, 1, 1, -1.0, 1.0);
    global_dpd_->file2_close(&F_oo);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    // Temp_ijab = lambda_ijcb F_ac
    global_dpd_->file2_init(&F_vv, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "F <v|v>");
    global_dpd_->contract244(&F_vv, &Lbb, &T, 1, 2, 1, 1.0, 0.0);
    global_dpd_->file2_close(&F_vv);
    global_dpd_->buf4_close(&Lbb);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    // F_ijab = Temp_ijab
    global_dpd_->buf4_copy(&T, PSIF_DCFT_DPD, "F <oo|vv>");
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "F <oo|vv>");
    global_dpd_->buf4_sort(&T, PSIF_DCFT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    global_dpd_->buf4_close(&T);
    // F_ijab -= Temp_ijba
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_add(&F, &T, -1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    // Temp_ijab = -lambda_kjab X_ik
    global_dpd_->file2_init(&F_oo, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "F <o|o>");
    global_dpd_->contract244(&F_oo, &Lbb, &T, 1, 0, 0, -1.0, 0.0);
    global_dpd_->file2_close(&F_oo);
    global_dpd_->buf4_close(&Lbb);
    // F_ijab += Temp_ijab
    dpd_buf4_add(&F, &T, 1.0);
    global_dpd_->buf4_sort(&T, PSIF_DCFT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    global_dpd_->buf4_close(&T);
    // F_ijab -= Temp_jiab
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_add(&F, &T, -1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&F);

    psio_->close(PSIF_LIBTRANS_DPD, 1);
}

void
DCFTSolver::form_density_weighted_fock() {

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdfile2 T_OO, T_oo, T_VV, T_vv, F_OO, F_oo, F_VV, F_vv;

    global_dpd_->file2_init(&T_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    global_dpd_->file2_init(&T_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
    global_dpd_->file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    global_dpd_->file2_init(&T_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");
    global_dpd_->file2_mat_init(&T_OO);
    global_dpd_->file2_mat_init(&T_oo);
    global_dpd_->file2_mat_init(&T_VV);
    global_dpd_->file2_mat_init(&T_vv);
    global_dpd_->file2_mat_rd(&T_OO);
    global_dpd_->file2_mat_rd(&T_oo);
    global_dpd_->file2_mat_rd(&T_VV);
    global_dpd_->file2_mat_rd(&T_vv);

    // Copy Tau in MO basis from the DPD library
    SharedMatrix a_tau_mo (new Matrix ("Alpha Tau in the MO basis", nirrep_, nmopi_, nmopi_));
    SharedMatrix b_tau_mo (new Matrix ("Beta Tau in the MO basis", nirrep_, nmopi_, nmopi_));

    for(int h = 0; h < nirrep_; ++h){
        if(nsopi_[h] == 0) continue;

        // Alpha occupied
        for(int p = 0 ; p < naoccpi_[h]; ++p){
            for(int q = 0 ; q <= p; ++q){
                double value = T_OO.matrix[h][p][q];
                a_tau_mo->set(h, p, q, value);
                if (p != q) a_tau_mo->set(h, q, p, value);
            }
        }

        // Beta occupied
        for(int p = 0 ; p < nboccpi_[h]; ++p){
            for(int q = 0 ; q <= p; ++q){
                double value = T_oo.matrix[h][p][q];
                b_tau_mo->set(h, p, q, value);
                if (p != q) b_tau_mo->set(h, q, p, value);
            }
        }

        // Alpha virtual
        for(int p = 0 ; p < navirpi_[h]; ++p){
            for(int q = 0 ; q <= p; ++q){
                double value = T_VV.matrix[h][p][q];
                a_tau_mo->set(h, p + naoccpi_[h], q + naoccpi_[h], value);
                if (p != q) a_tau_mo->set(h, q + naoccpi_[h], p + naoccpi_[h], value);
            }
        }

        // Beta virtual
        for(int p = 0 ; p < nbvirpi_[h]; ++p){
            for(int q = 0 ; q <= p; ++q){
                double value = T_vv.matrix[h][p][q];
                b_tau_mo->set(h, p + nboccpi_[h], q + nboccpi_[h], value);
                if (p != q) b_tau_mo->set(h, q + nboccpi_[h], p + nboccpi_[h], value);
            }
        }
    }

    global_dpd_->file2_close(&T_OO);
    global_dpd_->file2_close(&T_oo);
    global_dpd_->file2_close(&T_VV);
    global_dpd_->file2_close(&T_vv);

    SharedMatrix a_evecs (new Matrix ("Tau Eigenvectors (Alpha)", nirrep_, nmopi_, nmopi_));
    SharedMatrix b_evecs (new Matrix ("Tau Eigenvectors (Beta)", nirrep_, nmopi_, nmopi_));
    SharedVector a_evals (new Vector ("Tau Eigenvalues (Alpha)", nirrep_, nmopi_));
    SharedVector b_evals (new Vector ("Tau Eigenvalues (Beta)", nirrep_, nmopi_));

    // Diagonalize Tau
    a_tau_mo->diagonalize(a_evecs, a_evals);
    a_tau_mo->zero();
    a_tau_mo->set_diagonal(a_evals);
    b_tau_mo->diagonalize(b_evecs, b_evals);
    b_tau_mo->zero();
    b_tau_mo->set_diagonal(b_evals);

    // Transform the Fock matrix to NSO basis
    SharedMatrix nso_Fa (new Matrix ("Alpha Fock in the NSO basis", nirrep_, nmopi_, nmopi_));
    SharedMatrix nso_Fb (new Matrix ("Beta Fock in the NSO basis", nirrep_, nmopi_, nmopi_));

    // Alpha spin
    nso_Fa->transform(moFa_, a_evecs);
    // Beta spin
    nso_Fb->transform(moFb_, b_evecs);

    // Form density-weighted Fock matrix in the NSO basis

    for(int h = 0; h < nirrep_; ++h){
        if(nsopi_[h] == 0) continue;

        // Alpha occupied
        for(int ip = 0 ; ip < naoccpi_[h]; ++ip){
            for(int kp = 0 ; kp < naoccpi_[h]; ++kp){
                double value = nso_Fa->get(h, ip, kp) / (1.0 + a_evals->get(h, ip) + a_evals->get(h, kp));
                nso_Fa->set(h, ip, kp, value);
            }
        }

        // Beta occupied
        for(int ip = 0 ; ip < nboccpi_[h]; ++ip){
            for(int kp = 0 ; kp < nboccpi_[h]; ++kp){
                double value = nso_Fb->get(h, ip, kp) / (1.0 + b_evals->get(h, ip) + b_evals->get(h, kp));
                nso_Fb->set(h, ip, kp, value);
            }
        }

        // Alpha virtual
        for(int ap = naoccpi_[h]; ap < nmopi_[h]; ++ap){
            for(int dp = naoccpi_[h]; dp < nmopi_[h]; ++dp){
                double value = nso_Fa->get(h, dp, ap) / (1.0 - a_evals->get(h, dp) - a_evals->get(h, ap));
                nso_Fa->set(h, dp, ap, value);
            }
        }

        // Beta virtual
        for(int ap = nboccpi_[h]; ap < nmopi_[h]; ++ap){
            for(int dp = nboccpi_[h]; dp < nmopi_[h]; ++dp){
                double value = nso_Fb->get(h, dp, ap) / (1.0 - b_evals->get(h, dp) - b_evals->get(h, ap));
                nso_Fb->set(h, dp, ap, value);
            }
        }
    }

    // Transform density-weighted density matrix back to the original basis

    // Alpha spin
    nso_Fa->back_transform(a_evecs);
    // Beta spin
    nso_Fb->back_transform(b_evecs);

    global_dpd_->file2_init(&F_OO, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    global_dpd_->file2_init(&F_oo, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "F <o|o>");
    global_dpd_->file2_init(&F_VV, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");
    global_dpd_->file2_init(&F_vv, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "F <v|v>");

    global_dpd_->file2_mat_init(&F_OO);
    global_dpd_->file2_mat_init(&F_oo);
    global_dpd_->file2_mat_init(&F_VV);
    global_dpd_->file2_mat_init(&F_vv);

    for(int h = 0; h < nirrep_; ++h){
        if(nsopi_[h] == 0) continue;

        // Alpha occupied
        for(int ip = 0 ; ip < naoccpi_[h]; ++ip){
            for(int kp = 0 ; kp < naoccpi_[h]; ++kp){
                F_OO.matrix[h][ip][kp] = nso_Fa->get(h, ip, kp);
            }
        }

        // Beta occupied
        for(int ip = 0 ; ip < nboccpi_[h]; ++ip){
            for(int kp = 0 ; kp < nboccpi_[h]; ++kp){
                F_oo.matrix[h][ip][kp] = nso_Fb->get(h, ip, kp);
            }
        }

        // Alpha virtual
        for(int ap = 0 ; ap < navirpi_[h]; ++ap){
            for(int dp = 0 ; dp < navirpi_[h]; ++dp){
                F_VV.matrix[h][ap][dp] = nso_Fa->get(h, ap + naoccpi_[h], dp + naoccpi_[h]);
            }
        }

        // Beta virtual
        for(int ap = 0 ; ap < nbvirpi_[h]; ++ap){
            for(int dp = 0 ; dp < nbvirpi_[h]; ++dp){
                F_vv.matrix[h][ap][dp] = nso_Fb->get(h, ap + nboccpi_[h], dp + nboccpi_[h]);
            }
        }

    }

    global_dpd_->file2_mat_wrt(&F_OO);
    global_dpd_->file2_mat_wrt(&F_oo);
    global_dpd_->file2_mat_wrt(&F_VV);
    global_dpd_->file2_mat_wrt(&F_vv);

    global_dpd_->file2_close(&F_OO);
    global_dpd_->file2_close(&F_oo);
    global_dpd_->file2_close(&F_VV);
    global_dpd_->file2_close(&F_vv);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

}

void
DCFTSolver::compute_V_intermediate() {

    dpdbuf4 I, L, V, T, II, J, Taa, Tab, Tbb, Kaa, Kab, Kbb;
    dpdfile2 T_OO, T_oo, T_VV, T_vv, H_OO, H_oo, H_VV, H_vv;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    global_dpd_->file2_init(&T_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "T <O|O>");
    global_dpd_->file2_init(&T_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "T <o|o>");
    global_dpd_->file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "T <V|V>");
    global_dpd_->file2_init(&T_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "T <v|v>");

    /*
     * V_ijab = 1/3 P_(ij) gbar_abik T_kj
     */

    // OOVV

    // V_IJAB = 1/3 gbar_IKAB * T_KJ
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 1, "MO Ints <OO|VV>");
    global_dpd_->contract424(&I, &T_OO, &T, 1, 0, 1, 1.0/3.0, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&T);

    // Temp_IJAB -> Temp_JIAB
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_sort(&T, PSIF_DCFT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    global_dpd_->buf4_scm(&V, 0.0);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // V_IJAB -= 1/3 gbar_JKAB * T_KI
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // OoVv

    // V_IjAb += 1/3 gbar_IkAb * T_kj
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    global_dpd_->contract424(&I, &T_oo, &V, 1, 0, 1, 1.0/3.0, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&V);

    // V_IjAb += 1/3 T_IK * gbar_KjAb
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    global_dpd_->contract244(&T_OO, &I, &V, 1, 0, 0, 1.0/3.0, 1.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&V);

    // oovv

    // V_ijab = 1/3 gbar_ikab * T_kj
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o,o]"), ID("[v,v]"), 1, "MO Ints <oo|vv>");
    global_dpd_->contract424(&I, &T_oo, &T, 1, 0, 1, 1.0/3.0, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&T);

    // Temp_ijab -> Temp_jiab
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_sort(&T, PSIF_DCFT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    global_dpd_->buf4_scm(&V, 0.0);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // V_ijab -= 1/3 gbar_jkab * T_ki
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    /*
     * V_ijab -= 1/3 P_(ab) gbar_ijac T_cb
     */

    // OOVV

    // V_IJAB -= 1/3 gbar_IJAC * T_CB
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 1, "MO Ints <OO|VV>");
    global_dpd_->contract424(&I, &T_VV, &T, 3, 0, 0, -1.0/3.0, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&T);

    // Temp_IJAB -> Temp_IJBA
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_sort(&T, PSIF_DCFT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // V_IJAB += 1/3 gbar_IJBC * T_CA
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // OoVv

    // V_IjAb -= 1/3 gbar_IjAc * T_cb
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    global_dpd_->contract424(&I, &T_vv, &V, 3, 0, 0, -1.0/3.0, 1.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&V);

    // V_IjAb -= 1/3 gbar_IjCb * T_CA
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    global_dpd_->contract244(&T_VV, &I, &V, 1, 2, 1, -1.0/3.0, 1.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&V);

    // oovv

    // V_ijab -= 1/3 gbar_ijac * T_cb
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o,o]"), ID("[v,v]"), 1, "MO Ints <oo|vv>");
    global_dpd_->contract424(&I, &T_vv, &T, 3, 0, 0, -1.0/3.0, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&T);

    // Temp_ijab -> Temp_ijba
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_sort(&T, PSIF_DCFT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // V_ijab += 1/3 gbar_ijbc * T_ca
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    global_dpd_->file2_close(&T_OO);
    global_dpd_->file2_close(&T_oo);
    global_dpd_->file2_close(&T_VV);
    global_dpd_->file2_close(&T_vv);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    /*
     * V_ijab += 1/3 P_(ij) lambda_abik * H_kj
     */

    compute_H_intermediate();

    global_dpd_->file2_init(&H_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "H <O|O>");
    global_dpd_->file2_init(&H_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "H <o|o>");
    global_dpd_->file2_init(&H_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "H <V|V>");
    global_dpd_->file2_init(&H_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "H <v|v>");

    // OOVV

    // V_IJAB += 1/3 lambda_IKAB * H_KJ
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    global_dpd_->contract424(&L, &H_OO, &T, 1, 0, 1, 1.0/3.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);

    // Temp_IJAB -> Temp_JIAB
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_sort(&T, PSIF_DCFT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // V_IJAB -= 1/3 lambda_JKAB * H_KI
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // OoVv

    // V_IjAb += 1/3 lambda_IkAb * H_kj
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    global_dpd_->contract424(&L, &H_oo, &V, 1, 0, 1, 1.0/3.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    // V_IjAb += 1/3 H_IK * lambda_KjAb
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    global_dpd_->contract244(&H_OO, &L, &V, 1, 0, 0, 1.0/3.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    // oovv

    // V_ijab += 1/3 lambda_ikab * H_kj
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    global_dpd_->contract424(&L, &H_oo, &T, 1, 0, 1, 1.0/3.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);

    // Temp_ijab -> Temp_jiab
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_sort(&T, PSIF_DCFT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // V_ijab -= 1/3 lambda_jkab * H_ki
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    /*
     * V_ijab -= 1/3 P_(ab) lambda_ijac H_cb
     */

    // OOVV

    // V_IJAB -= 1/3 lambda_IJAC * H_CB
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    global_dpd_->contract424(&L, &H_VV, &T, 3, 0, 0, -1.0/3.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);

    // Temp_IJAB -> Temp_IJBA
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_sort(&T, PSIF_DCFT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // V_IJAB += 1/3 lambda_IJBC * H_CA
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // OoVv

    // V_IjAb -= 1/3 lambda_IjAc * H_cb
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    global_dpd_->contract424(&L, &H_vv, &V, 3, 0, 0, -1.0/3.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    // V_IjAb -= 1/3 lambda_IjCb * H_CA
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    global_dpd_->contract244(&H_VV, &L, &V, 1, 2, 1, -1.0/3.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    // oovv

    // V_ijab -= 1/3 lambda_ijac * H_cb
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    global_dpd_->contract424(&L, &H_vv, &T, 3, 0, 0, -1.0/3.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);

    // Temp_ijab -> Temp_ijba
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_sort(&T, PSIF_DCFT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // V_ijab += 1/3 lambda_ijbc * H_ca
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    global_dpd_->file2_close(&H_OO);
    global_dpd_->file2_close(&H_oo);
    global_dpd_->file2_close(&H_VV);
    global_dpd_->file2_close(&H_vv);

    /*
     * V_ijab += 1/6 gbar_abkl I_klij
     */

    compute_I_intermediate();

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // V_IJAB += 1/6 * gbar_ABKL * I_KLIJ
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                           ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                           ID("[O,O]"), ID("[V,V]"), 1, "MO Ints <OO|VV>");
    global_dpd_->buf4_init(&II, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[O>O]-"),
                           ID("[O>O]-"), ID("[O>O]-"), 0, "I <OO|OO>");
    global_dpd_->contract444(&II, &I, &V, 0, 1, 1.0/3.0, 1.0);
    global_dpd_->buf4_close(&II);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&V);

    // V_IjAb += 1/3 * gbar_AbKl * I_KlIj
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    global_dpd_->buf4_init(&II, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                           ID("[O,o]"), ID("[O,o]"), 0, "I <Oo|Oo>");
    global_dpd_->contract444(&II, &I, &V, 0, 1, 1.0/3.0, 1.0);
    global_dpd_->buf4_close(&II);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&V);

    // V_ijab += 1/6 * gbar_abkl * I_klij
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                           ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                           ID("[o,o]"), ID("[v,v]"), 1, "MO Ints <oo|vv>");
    global_dpd_->buf4_init(&II, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[o>o]-"),
                           ID("[o>o]-"), ID("[o>o]-"), 0, "I <oo|oo>");
    global_dpd_->contract444(&II, &I, &V, 0, 1, 1.0/3.0, 1.0);
    global_dpd_->buf4_close(&II);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&V);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    /*
     * V_ijab += 1/6 lambda_abkl J_klij
     */

    compute_J_intermediate();

    // V_IJAB += 1/6 * lambda_ABKL * J_KLIJ
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                           ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                           ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    global_dpd_->buf4_init(&J, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[O>O]-"),
                           ID("[O>O]-"), ID("[O>O]-"), 0, "J <OO|OO>");
    global_dpd_->contract444(&J, &L, &V, 0, 1, 1.0/3.0, 1.0);
    global_dpd_->buf4_close(&J);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    // V_IjAb += 1/3 * lambda_AbKl * J_KlIj
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    global_dpd_->buf4_init(&J, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                           ID("[O,o]"), ID("[O,o]"), 0, "J <Oo|Oo>");
    global_dpd_->contract444(&J, &L, &V, 0, 1, 1.0/3.0, 1.0);
    global_dpd_->buf4_close(&J);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    // V_ijab += 1/6 * lambda_abkl * J_klij
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                           ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                           ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    global_dpd_->buf4_init(&J, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[o>o]-"),
                           ID("[o>o]-"), ID("[o>o]-"), 0, "J <oo|oo>");
    global_dpd_->contract444(&J, &L, &V, 0, 1, 1.0/3.0, 1.0);
    global_dpd_->buf4_close(&J);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    /*
     * V_ijab += 2/3 P_(ij) P_(ab) gbar_acik K_kbjc
     */

    compute_K_intermediate();

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // OOVV
    global_dpd_->buf4_init(&Taa, PSIF_DCFT_DPD, 0, ID("[V,O]"), ID("[O,V]"),
                  ID("[V,O]"), ID("[O,V]"), 0, "Temp (VO|OV)");

    global_dpd_->buf4_init(&Kaa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "K (OV|OV)");
    global_dpd_->buf4_init(&Kab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "K (OV|ov)");

    // T_AIJB = 2/3 (AI|KC) K_(KC|JB)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,V]"),
                  ID("[V,O]"), ID("[O,V]"), 0, "MO Ints (VO|OV)");
    global_dpd_->contract444(&I, &Kaa, &Taa, 0, 0, 2.0/3.0, 0.0);
    global_dpd_->buf4_close(&I);

    // T_AIJB -= 2/3 <AI|KC> K_(KC|JB)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,V]"),
                  ID("[V,O]"), ID("[O,V]"), 0, "MO Ints <VO|OV>");
    global_dpd_->contract444(&I, &Kaa, &Taa, 0, 0, -2.0/3.0, 1.0);
    global_dpd_->buf4_close(&I);

    // T_AIJB += 2/3 (AI|kc) K_(kc|JB)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[o,v]"),
                  ID("[V,O]"), ID("[o,v]"), 0, "MO Ints (VO|ov)");
    global_dpd_->contract444(&I, &Kab, &Taa, 0, 0, 2.0/3.0, 1.0);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_close(&Kaa);
    global_dpd_->buf4_close(&Kab);

    // T_AIJB -> T_IJAB
    global_dpd_->buf4_sort(&Taa, PSIF_DCFT_DPD, qrps, ID("[O,O]"), ID("[V,V]"), "Temp <OO|VV>");
    global_dpd_->buf4_close(&Taa);

    // V_IJAB += T_IJAB
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&Taa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");

    // T_IJAB -> T_JIAB
    global_dpd_->buf4_sort(&Taa, PSIF_DCFT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");

    // V_IJAB -= T_JIAB
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // T_IJAB -> T_IJBA
    global_dpd_->buf4_sort(&Taa, PSIF_DCFT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");

    // V_IJAB -= T_IJBA
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // T_IJAB -> T_JIBA
    global_dpd_->buf4_sort(&Taa, PSIF_DCFT_DPD, qpsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");

    // V_IJAB += T_JIBA
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "V <OO|VV>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_print(&V, outfile, 1);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_close(&Taa);

    // OoVv
    global_dpd_->buf4_init(&Kaa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "K (OV|OV)");
    global_dpd_->buf4_init(&Kab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "K (OV|ov)");
    global_dpd_->buf4_init(&Kbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "K (ov|ov)");

    global_dpd_->buf4_init(&Tab, PSIF_DCFT_DPD, 0, ID("[V,O]"), ID("[o,v]"),
                  ID("[V,O]"), ID("[o,v]"), 0, "Temp (VO|ov)");

    // T_AIjb = 2/3 (AI|kc) K_(kc|jb)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[o,v]"),
                  ID("[V,O]"), ID("[o,v]"), 0, "MO Ints (VO|ov)");
    global_dpd_->contract444(&I, &Kbb, &Tab, 0, 1, 2.0/3.0, 0.0);
    global_dpd_->buf4_close(&I);

    // T_AIjb += 2/3 (AI|KC) K_(KC|jb)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,V]"),
                  ID("[V,O]"), ID("[O,V]"), 0, "MO Ints (VO|OV)");
    global_dpd_->contract444(&I, &Kab, &Tab, 0, 1, 2.0/3.0, 1.0);
    global_dpd_->buf4_close(&I);

    // T_AIjb -= 2/3 <AI|KC> K_(KC|jb)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,V]"),
                  ID("[V,O]"), ID("[O,V]"), 0, "MO Ints <VO|OV>");
    global_dpd_->contract444(&I, &Kab, &Tab, 0, 1, -2.0/3.0, 1.0);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_close(&Tab);

    // T_AIjb -> T_IjAb
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[V,O]"), ID("[o,v]"),
                  ID("[V,O]"), ID("[o,v]"), 0, "Temp (VO|ov)");
    global_dpd_->buf4_sort(&T, PSIF_DCFT_DPD, qrps, ID("[O,o]"), ID("[V,v]"), "Temp <Oo|Vv>");
    global_dpd_->buf4_close(&T);

    // V_IjAb += T_IjAb
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Temp <Oo|Vv>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&Tab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[v,o]"),
                  ID("[O,V]"), ID("[v,o]"), 0, "Temp (OV|vo)");

    // T_IAbj = 2/3 (bj|KC) K_(KC|IA)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[O,V]"),
                  ID("[v,o]"), ID("[O,V]"), 0, "MO Ints (vo|OV)");
    global_dpd_->contract444(&Kaa, &I, &Tab, 0, 0, 2.0/3.0, 0.0);
    global_dpd_->buf4_close(&I);

    // T_IAbj += 2/3 (bj|kc) K_(kc|IA)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,v]"),
                  ID("[v,o]"), ID("[o,v]"), 0, "MO Ints (vo|ov)");
    global_dpd_->contract444(&Kab, &I, &Tab, 0, 0, 2.0/3.0, 1.0);
    global_dpd_->buf4_close(&I);

    // T_IAbj -= 2/3 <bj|kc> K_(kc|IA)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,v]"),
                  ID("[v,o]"), ID("[o,v]"), 0, "MO Ints <vo|ov>");
    global_dpd_->contract444(&Kab, &I, &Tab, 0, 0, -2.0/3.0, 1.0);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_close(&Tab);

    // T_IAbj -> T_IjAb
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[v,o]"),
                  ID("[O,V]"), ID("[v,o]"), 0, "Temp (OV|vo)");
    global_dpd_->buf4_sort(&T, PSIF_DCFT_DPD, psqr, ID("[O,o]"), ID("[V,v]"), "Temp <Oo|Vv>");
    global_dpd_->buf4_close(&T);

    // V_IjAb += T_IjAb
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Temp <Oo|Vv>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_close(&Kaa);
    global_dpd_->buf4_close(&Kab);
    global_dpd_->buf4_close(&Kbb);

    // T_IbAj = 2/3 K_(Ib|Kc) <Kc|Aj>
    global_dpd_->buf4_init(&Tab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "Temp (Ov|Vo)");
    global_dpd_->buf4_init(&Kab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "K <Ov|Ov>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "MO Ints <Ov|Vo>");
    global_dpd_->contract444(&Kab, &I, &Tab, 0, 1, 2.0/3.0, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&Kab);
    global_dpd_->buf4_close(&Tab);

    // K_kCjA -> K_CkAj
    global_dpd_->buf4_init(&Kab, PSIF_DCFT_DPD, 0, ID("[o,V]"), ID("[o,V]"),
                           ID("[o,V]"), ID("[o,V]"), 0, "K <oV|oV>");
    global_dpd_->buf4_sort(&Kab, PSIF_DCFT_DPD, qpsr, ID("[V,o]"), ID("[V,o]"), "K <Vo|Vo>");
    global_dpd_->buf4_close(&Kab);

    // T_IbAj += 2/3 <Ib|Ck> K_(Ck|Aj)
    global_dpd_->buf4_init(&Tab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "Temp (Ov|Vo)");
    global_dpd_->buf4_init(&Kab, PSIF_DCFT_DPD, 0, ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "K <Vo|Vo>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "MO Ints <Ov|Vo>");
    global_dpd_->contract444(&I, &Kab, &Tab, 0, 1, 2.0/3.0, 1.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&Kab);
    global_dpd_->buf4_close(&Tab);

    // T_IbAj -> T_IjAb
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "Temp (Ov|Vo)");
    global_dpd_->buf4_sort(&T, PSIF_DCFT_DPD, psrq, ID("[O,o]"), ID("[V,v]"), "Temp <Oo|Vv>");
    global_dpd_->buf4_close(&T);

    // V_IjAb += T_IjAb
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Temp <Oo|Vv>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "V <Oo|Vv>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_print(&V, outfile, 1);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // oovv
    global_dpd_->buf4_init(&Tbb, PSIF_DCFT_DPD, 0, ID("[v,o]"), ID("[o,v]"),
                  ID("[v,o]"), ID("[o,v]"), 0, "Temp (vo|ov)");

    global_dpd_->buf4_init(&Kbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "K (ov|ov)");
    global_dpd_->buf4_init(&Kab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "K (OV|ov)");

    // T_aijb = 2/3 (ai|kc) K_(kc|jb)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,v]"),
                  ID("[v,o]"), ID("[o,v]"), 0, "MO Ints (vo|ov)");
    global_dpd_->contract444(&I, &Kbb, &Tbb, 0, 0, 2.0/3.0, 0.0);
    global_dpd_->buf4_close(&I);

    // T_aijb -= 2/3 <ai|kc> K_(kc|jb)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,v]"),
                  ID("[v,o]"), ID("[o,v]"), 0, "MO Ints <vo|ov>");
    global_dpd_->contract444(&I, &Kbb, &Tbb, 0, 0, -2.0/3.0, 1.0);
    global_dpd_->buf4_close(&I);

    // T_aijb += 2/3 (ai|KC) K_(KC|jb)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[O,V]"),
                  ID("[v,o]"), ID("[O,V]"), 0, "MO Ints (vo|OV)");
    global_dpd_->contract444(&I, &Kab, &Tbb, 0, 1, 2.0/3.0, 1.0);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_close(&Kab);
    global_dpd_->buf4_close(&Kbb);

    // T_aijb -> T_ijab
    global_dpd_->buf4_sort(&Tbb, PSIF_DCFT_DPD, qrps, ID("[o,o]"), ID("[v,v]"), "Temp <oo|vv>");
    global_dpd_->buf4_close(&Tbb);

    // V_ijab += T_ijab
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&Tbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");

    // T_ijab -> T_jiab
    global_dpd_->buf4_sort(&Tbb, PSIF_DCFT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");

    // V_ijab -= T_jiab
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // T_ijab -> T_ijba
    global_dpd_->buf4_sort(&Tbb, PSIF_DCFT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");

    // V_ijab -= T_ijba
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    // T_ijab -> T_jiba
    global_dpd_->buf4_sort(&Tbb, PSIF_DCFT_DPD, qpsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");

    // V_ijab += T_jiba
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    global_dpd_->buf4_init(&V, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "V <oo|vv>");
    dpd_buf4_add(&V, &T, 1.0);
    global_dpd_->buf4_print(&V, outfile, 1);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_close(&Tbb);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    /*
     * V_ijab += 1/3 P_(ij) P_(ab) lambda_acik L_kbjc
     */

    compute_L_intermediate();


}

void
DCFTSolver::compute_H_intermediate() {

    dpdbuf4 L, I;
    dpdfile2 H_OO, H_oo, H_VV, H_vv, HH_OO, HH_oo, HH_VV, HH_vv;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    global_dpd_->file2_init(&H_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "H <O|O>");
    global_dpd_->file2_init(&H_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "H <o|o>");
    global_dpd_->file2_init(&H_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "H <V|V>");
    global_dpd_->file2_init(&H_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "H <v|v>");

    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 1, "MO Ints <OO|VV>");
    /*
     * H_IJ = -1/2 Lambda_IKAB gbar_JKAB
     */
    global_dpd_->contract442(&L, &I, &H_OO, 0, 0, -1.0, 0.0);
    /*
     * H_AB = +1/2 Lambda_IJAC gbar_IJBC
     */
    global_dpd_->contract442(&L, &I, &H_VV, 2, 2, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o,o]"), ID("[v,v]"), 1, "MO Ints <oo|vv>");
    /*
     * H_ij = -1/2 Lambda_ikab gbar_jkab
     */
    global_dpd_->contract442(&L, &I, &H_oo, 0, 0, -1.0, 0.0);
    /*
     * H_ab = +1/2 Lambda_ijac gbar_ijbc
     */
    global_dpd_->contract442(&L, &I, &H_vv, 2, 2, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    /*
     * H_IJ -= Lambda_IkAb gbar_JkAb - Lambda_IkaB gbar_JkaB
     */
    global_dpd_->contract442(&L, &I, &H_OO, 0, 0, -2.0, 1.0);
    /*
     * H_ij -= Lambda_KiAb gbar_KjAb - Lambda_KiaB gbar_KjaB
     */
    global_dpd_->contract442(&L, &I, &H_oo, 1, 1, -2.0, 1.0);
    /*
     * H_AB += Lambda_IjAc gbar_IjBc + Lambda_iJAc gbar_iJBc
     */
    global_dpd_->contract442(&L, &I, &H_VV, 2, 2, 2.0, 1.0);
    /*
     * H_ab += Lambda_IjCa gbar_IjCb + Lambda_iJCa gbar_iJCb
     */
    global_dpd_->contract442(&L, &I, &H_vv, 3, 3, 2.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&I);

    global_dpd_->file2_close(&H_OO);
    global_dpd_->file2_close(&H_oo);
    global_dpd_->file2_close(&H_VV);
    global_dpd_->file2_close(&H_vv);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    // Symmetrize H
    global_dpd_->file2_init(&H_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "H <O|O>");
    global_dpd_->file2_init(&H_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "H <o|o>");
    global_dpd_->file2_init(&H_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "H <V|V>");
    global_dpd_->file2_init(&H_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "H <v|v>");
    global_dpd_->file2_init(&HH_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "H <O|O>");
    global_dpd_->file2_init(&HH_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "H <o|o>");
    global_dpd_->file2_init(&HH_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "H <V|V>");
    global_dpd_->file2_init(&HH_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "H <v|v>");

    global_dpd_->file2_axpy(&HH_OO, &H_OO, 1.0, 1);
    global_dpd_->file2_axpy(&HH_oo, &H_oo, 1.0, 1);
    global_dpd_->file2_axpy(&HH_VV, &H_VV, 1.0, 1);
    global_dpd_->file2_axpy(&HH_vv, &H_vv, 1.0, 1);

    global_dpd_->file2_close(&HH_OO);
    global_dpd_->file2_close(&HH_oo);
    global_dpd_->file2_close(&HH_VV);
    global_dpd_->file2_close(&HH_vv);
    global_dpd_->file2_close(&H_OO);
    global_dpd_->file2_close(&H_oo);
    global_dpd_->file2_close(&H_VV);
    global_dpd_->file2_close(&H_vv);

}

void
DCFTSolver::compute_I_intermediate() {

    dpdbuf4 LLaa, LLab, LLbb, Laa, Lab, Lbb, Iaa, Iab, Ibb;

    // I_ijkl = Lambda_ijab * Lambda_klab
    global_dpd_->buf4_init(&Iaa, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[O>O]-"),
                           ID("[O>O]-"), ID("[O>O]-"), 0, "I <OO|OO>");
    global_dpd_->buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                           ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    global_dpd_->buf4_init(&LLaa, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                           ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    global_dpd_->contract444(&Laa, &LLaa, &Iaa, 0, 0, 2.0, 0.0);
    global_dpd_->buf4_close(&LLaa);
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&Iaa);

    global_dpd_->buf4_init(&Iab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                           ID("[O,o]"), ID("[O,o]"), 0, "I <Oo|Oo>");
    global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    global_dpd_->buf4_init(&LLab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    global_dpd_->contract444(&Lab, &LLab, &Iab, 0, 0, 2.0, 0.0);
    global_dpd_->buf4_close(&LLab);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&Iab);

    global_dpd_->buf4_init(&Ibb, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[o>o]-"),
                           ID("[o>o]-"), ID("[o>o]-"), 0, "I <oo|oo>");
    global_dpd_->buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                           ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    global_dpd_->buf4_init(&LLbb, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                           ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    global_dpd_->contract444(&Lbb, &LLbb, &Ibb, 0, 0, 2.0, 0.0);
    global_dpd_->buf4_close(&LLbb);
    global_dpd_->buf4_close(&Lbb);
    global_dpd_->buf4_close(&Ibb);

}

void
DCFTSolver::compute_J_intermediate() {

    dpdbuf4 Iaa, Iab, Ibb, Laa, Lab, Lbb, Jaa, Jab, Jbb;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // J_ijkl = Lambda_ijab * gbar_klab + gbar_ijab * lambda_klab
    global_dpd_->buf4_init(&Jaa, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[O>O]-"),
                           ID("[O>O]-"), ID("[O>O]-"), 0, "J <OO|OO>");
    global_dpd_->buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                           ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    global_dpd_->buf4_init(&Iaa, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                           ID("[O,O]"), ID("[V,V]"), 1, "MO Ints <OO|VV>");
    global_dpd_->contract444(&Laa, &Iaa, &Jaa, 0, 0, 4.0, 0.0);
    global_dpd_->buf4_close(&Iaa);
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&Jaa);

    global_dpd_->buf4_init(&Jab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                           ID("[O,o]"), ID("[O,o]"), 0, "J <Oo|Oo>");
    global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    global_dpd_->buf4_init(&Iab, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    global_dpd_->contract444(&Lab, &Iab, &Jab, 0, 0, 4.0, 0.0);
    global_dpd_->buf4_close(&Iab);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&Jab);

    global_dpd_->buf4_init(&Jbb, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[o>o]-"),
                           ID("[o>o]-"), ID("[o>o]-"), 0, "J <oo|oo>");
    global_dpd_->buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                           ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    global_dpd_->buf4_init(&Ibb, PSIF_LIBTRANS_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                           ID("[o,o]"), ID("[v,v]"), 1, "MO Ints <oo|vv>");
    global_dpd_->contract444(&Lbb, &Ibb, &Jbb, 0, 0, 4.0, 0.0);
    global_dpd_->buf4_close(&Ibb);
    global_dpd_->buf4_close(&Lbb);
    global_dpd_->buf4_close(&Jbb);

    // Bra-ket symmetrize J
    global_dpd_->buf4_init(&Jaa, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[O>O]-"),
                           ID("[O>O]-"), ID("[O>O]-"), 0, "J <OO|OO>");
    global_dpd_->buf4_symm(&Jaa);
    global_dpd_->buf4_close(&Jaa);

    global_dpd_->buf4_init(&Jab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                           ID("[O,o]"), ID("[O,o]"), 0, "J <Oo|Oo>");
    global_dpd_->buf4_symm(&Jab);
    global_dpd_->buf4_close(&Jab);

    global_dpd_->buf4_init(&Jbb, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[o>o]-"),
                           ID("[o>o]-"), ID("[o>o]-"), 0, "J <oo|oo>");
    global_dpd_->buf4_symm(&Jbb);
    global_dpd_->buf4_close(&Jbb);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

}

void
DCFTSolver::compute_K_intermediate() {

    dpdbuf4 LLaa, LLab, LLbb, Laa, Lab, Lbb, Kaa, Kab, Kba, Kbb;

    // There are five unique spin cases: K<IAJB>, K<iajb>, K<IaJb>, K<iAjB>, K<IajB>

    // K<IAJB> spin case

    global_dpd_->buf4_init(&Kaa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                           ID("[O,V]"), ID("[O,V]"), 0, "K (OV|OV)");
    global_dpd_->buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                           ID("[O,V]"), ID("[O,V]"), 0, "Lambda (OV|OV)");
    global_dpd_->buf4_init(&LLaa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                           ID("[O,V]"), ID("[O,V]"), 0, "Lambda (OV|OV)");
    global_dpd_->contract444(&Laa, &LLaa, &Kaa, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&LLaa);
    global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                           ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");
    global_dpd_->buf4_init(&LLab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                           ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");
    global_dpd_->contract444(&Lab, &LLab, &Kaa, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&LLab);
    global_dpd_->buf4_close(&Kaa);

    // Resort K(OV|OV) to the K<OV|OV>
    global_dpd_->buf4_init(&Kaa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                           ID("[O,V]"), ID("[O,V]"), 0, "K (OV|OV)");
    global_dpd_->buf4_sort(&Kaa, PSIF_DCFT_DPD, psrq, ID("[O,V]"),ID("[O,V]"), "K <OV|OV>");
    global_dpd_->buf4_close(&Kaa);

    // K<IaJb> and K<iAjB> spin cases:

    // Although we denote K <Ov|Ov> and K <Vo|Vo> as in physist's notation, it is actually in chemist's notation.
    // However, it is convenient to store that intermediate in this form to avoid unnecessary tensor sorts for the density computation
    // E.g. K_(Kb|Ic) <- Lambda_(Kb|mE) Lambda_(Ic|mE) = Lambda_<Km|Eb> Lambda_<Im|Ec>
    global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                           ID("[O,v]"), ID("[o,V]"), 0, "Lambda (Ov|oV)");
    global_dpd_->buf4_init(&LLab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                           ID("[O,v]"), ID("[o,V]"), 0, "Lambda (Ov|oV)");
    global_dpd_->buf4_init(&Kab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                           ID("[O,v]"), ID("[O,v]"), 0, "K <Ov|Ov>");
    global_dpd_->contract444(&Lab, &LLab, &Kab, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&Kab);
    global_dpd_->buf4_init(&Kba, PSIF_DCFT_DPD, 0, ID("[o,V]"), ID("[o,V]"),
                           ID("[o,V]"), ID("[o,V]"), 0, "K <oV|oV>");
    global_dpd_->contract444(&Lab, &LLab, &Kba, 1, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&Kba);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&LLab);

    // K<IajB> spin case:

    global_dpd_->buf4_init(&Kab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                           ID("[O,V]"), ID("[o,v]"), 0, "K (OV|ov)");
    global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                           ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");
    global_dpd_->buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                           ID("[O,V]"), ID("[O,V]"), 0, "Lambda (OV|OV)");
    global_dpd_->contract444(&Laa, &Lab, &Kab, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                           ID("[o,v]"), ID("[o,v]"), 0, "Lambda (ov|ov)");
    global_dpd_->contract444(&Lab, &Lbb, &Kab, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&Lbb);
    global_dpd_->buf4_close(&Kab);
    global_dpd_->buf4_init(&Kab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                           ID("[O,V]"), ID("[o,v]"), 0, "K (OV|ov)");
    global_dpd_->buf4_sort(&Kab, PSIF_DCFT_DPD, psrq, ID("[O,v]"), ID("[o,V]"), "K <Ov|oV>");
    // Resort to get the K_oVOv. Used for the MO Lagrangian
    global_dpd_->buf4_sort(&Kab, PSIF_DCFT_DPD, rqps, ID("[o,V]"), ID("[O,v]"), "K <oV|Ov>");
    global_dpd_->buf4_close(&Kab);
    global_dpd_->buf4_close(&Lab);

    // K<iajb> spin case:

    global_dpd_->buf4_init(&Kbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                           ID("[o,v]"), ID("[o,v]"), 0, "K (ov|ov)");
    global_dpd_->buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                           ID("[o,v]"), ID("[o,v]"), 0, "Lambda (ov|ov)");
    global_dpd_->buf4_init(&LLbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                           ID("[o,v]"), ID("[o,v]"), 0, "Lambda (ov|ov)");
    global_dpd_->contract444(&Lbb, &LLbb, &Kbb, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&Lbb);
    global_dpd_->buf4_close(&LLbb);
    global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                           ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");
    global_dpd_->buf4_init(&LLab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                           ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");
    global_dpd_->contract444(&Lab, &LLab, &Kbb, 1, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&LLab);
    global_dpd_->buf4_close(&Kbb);

    // Resort K(ov|ov) to the K<ov|ov>
    global_dpd_->buf4_init(&Kbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                           ID("[o,v]"), ID("[o,v]"), 0, "K (ov|ov)");
    global_dpd_->buf4_sort(&Kbb, PSIF_DCFT_DPD, psrq, ID("[o,v]"),ID("[o,v]"), "K <ov|ov>");
    global_dpd_->buf4_close(&Kbb);

}

void
DCFTSolver::compute_L_intermediate() {

    dpdbuf4 I, Iaa, Iab, Ibb, Lambda_aa, Lambda_ab, Lambda_bb, Laa, Lab, Lba, Lbb;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // There are six unique spin cases: L<IAJB>, L<iajb>, L<IaJb>, L<iAjB>, L<IajB>, L<iAJb>

    // L<IAJB> spin case
    global_dpd_->buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[V,O]"),
                           ID("[O,V]"), ID("[V,O]"), 0, "L (OV|VO)");

    global_dpd_->buf4_init(&Lambda_aa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                           ID("[O,V]"), ID("[O,V]"), 0, "Lambda (OV|OV)");

    // L_(IB|AJ) = Lambda_(IB|KC) * (AJ|KC)
    global_dpd_->buf4_init(&Iaa, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,V]"),
                           ID("[V,O]"), ID("[O,V]"), 0, "MO Ints (VO|OV)");
    global_dpd_->contract444(&Lambda_aa, &Iaa, &Laa, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&Iaa);

    // L_(IB|AJ) -= Lambda_(IB|KC) * <AJ|KC>
    global_dpd_->buf4_init(&Iaa, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,V]"),
                           ID("[V,O]"), ID("[O,V]"), 0, "MO Ints <VO|OV>");
    global_dpd_->contract444(&Lambda_aa, &Iaa, &Laa, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&Iaa);

    global_dpd_->buf4_close(&Lambda_aa);

    // L_(IB|AJ) += Lambda_(IB|kc) * (AJ|kc)
    global_dpd_->buf4_init(&Lambda_ab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                           ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");
    global_dpd_->buf4_init(&Iab, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[o,v]"),
                           ID("[V,O]"), ID("[o,v]"), 0, "MO Ints (VO|ov)");
    global_dpd_->contract444(&Lambda_ab, &Iab, &Laa, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&Iab);
    global_dpd_->buf4_close(&Lambda_ab);

    global_dpd_->buf4_print(&Laa, outfile, 1);
    global_dpd_->buf4_close(&Laa);

    // L<IaJb> and L<iAjB> spin cases:

    // Although we denote L <Ov|Ov> and L <Vo|Vo> as in physist's notation, it is actually in chemist's notation.
    // However, it is convenient to store that intermediate in this form to avoid unnecessary tensor sorts
    global_dpd_->buf4_init(&Lambda_ab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                           ID("[O,v]"), ID("[o,V]"), 0, "Lambda (Ov|oV)");
    global_dpd_->buf4_sort(&Lambda_ab, PSIF_DCFT_DPD, pqsr, ID("[O,v]"), ID("[V,o]"), "Lambda (Ov|Vo)");
    global_dpd_->buf4_close(&Lambda_ab);

    global_dpd_->buf4_init(&Lambda_ab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[V,o]"),
                           ID("[O,v]"), ID("[V,o]"), 0, "Lambda (Ov|Vo)");
    global_dpd_->buf4_init(&Iab, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,o]"),
                           ID("[O,v]"), ID("[V,o]"), 0, "MO Ints <Ov|Vo>");

    // L_<Ib|Ja> = Lambda_(Ib|Ck) * <Ja|Ck>
    global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                           ID("[O,v]"), ID("[O,v]"), 0, "L <Ov|Ov>");
    global_dpd_->contract444(&Lambda_ab, &Iab, &Lab, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_print(&Lab, outfile, 1);
    global_dpd_->buf4_close(&Lab);

    // L_<Bi|Aj> = Lambda_(Kc|Bi) * <Kc|Aj>
    global_dpd_->buf4_init(&Lba, PSIF_DCFT_DPD, 0, ID("[V,o]"), ID("[V,o]"),
                           ID("[V,o]"), ID("[V,o]"), 0, "L <Vo|Vo>");
    global_dpd_->contract444(&Lambda_ab, &Iab, &Lba, 1, 1, 1.0, 0.0);
    global_dpd_->buf4_print(&Lba, outfile, 1);
    global_dpd_->buf4_close(&Lba);

    global_dpd_->buf4_close(&Iab);
    global_dpd_->buf4_close(&Lambda_ab);

    // L<IajB> spin case:
    global_dpd_->buf4_init(&Lambda_aa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                           ID("[O,V]"), ID("[O,V]"), 0, "Lambda (OV|OV)");
    global_dpd_->buf4_init(&Lambda_ab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                           ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");
    global_dpd_->buf4_init(&Lambda_bb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                           ID("[o,v]"), ID("[o,v]"), 0, "Lambda (ov|ov)");

    global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[v,o]"),
                           ID("[O,V]"), ID("[v,o]"), 0, "L (OV|vo)");

    // L_(IB|aj) = Lambda_(IB|kc) * (aj|kc)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,v]"),
                           ID("[v,o]"), ID("[o,v]"), 0, "MO Ints (vo|ov)");
    global_dpd_->contract444(&Lambda_ab, &I, &Lab, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&I);

    // L_(IB|aj) -= Lambda_(IB|kc) * <aj|kc>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,v]"),
                           ID("[v,o]"), ID("[o,v]"), 0, "MO Ints <vo|ov>");
    global_dpd_->contract444(&Lambda_ab, &I, &Lab, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&I);

    // L_(IB|aj) += Lambda_(IB|KC) * (aj|KC)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[O,V]"),
                           ID("[v,o]"), ID("[O,V]"), 0, "MO Ints (vo|OV)");
    global_dpd_->contract444(&Lambda_aa, &I, &Lab, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_print(&Lab, outfile, 1);
    global_dpd_->buf4_close(&Lab);

    // L<iAJb> spin case:
    global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[V,O]"), ID("[o,v]"),
                           ID("[V,O]"), ID("[o,v]"), 0, "L (VO|ov)");

    // L_(AJ|ib) += (AJ|KC) * Lambda_(KC|ib)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,V]"),
                           ID("[V,O]"), ID("[O,V]"), 0, "MO Ints (VO|OV)");
    global_dpd_->contract444(&I, &Lambda_ab, &Lab, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&I);

    // L_(AJ|ib) -= <AJ|KC> * Lambda_(KC|ib)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,V]"),
                           ID("[V,O]"), ID("[O,V]"), 0, "MO Ints <VO|OV>");
    global_dpd_->contract444(&I, &Lambda_ab, &Lab, 0, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&I);

    // L_(AJ|ib) += (AJ|kc) * Lambda_(ib|kc)
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[o,v]"),
                           ID("[V,O]"), ID("[o,v]"), 0, "MO Ints (VO|ov)");
    global_dpd_->contract444(&I, &Lambda_bb, &Lab, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_print(&Lab, outfile, 1);
    global_dpd_->buf4_close(&Lab);

    global_dpd_->buf4_close(&Lambda_bb);
    global_dpd_->buf4_close(&Lambda_ab);
    global_dpd_->buf4_close(&Lambda_aa);

    // L<iajb> spin case
    global_dpd_->buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[v,o]"),
                           ID("[o,v]"), ID("[v,o]"), 0, "L (ov|vo)");

    global_dpd_->buf4_init(&Lambda_bb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                           ID("[o,v]"), ID("[o,v]"), 0, "Lambda (ov|ov)");

    // L_(ib|aj) = Lambda_(ib|kc) * (aj|kc)
    global_dpd_->buf4_init(&Ibb, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,v]"),
                           ID("[v,o]"), ID("[o,v]"), 0, "MO Ints (vo|ov)");
    global_dpd_->contract444(&Lambda_bb, &Ibb, &Lbb, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&Ibb);

    // L_(ib|aj) -= Lambda_(ib|kc) * <aj|kc>
    global_dpd_->buf4_init(&Ibb, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,v]"),
                           ID("[v,o]"), ID("[o,v]"), 0, "MO Ints <vo|ov>");
    global_dpd_->contract444(&Lambda_bb, &Ibb, &Lbb, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&Ibb);

    global_dpd_->buf4_close(&Lambda_bb);

    // L_(ib|aj) += Lambda_(ib|KC) * (aj|KC)
    global_dpd_->buf4_init(&Lambda_ab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                           ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");
    global_dpd_->buf4_init(&Iab, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[O,V]"),
                           ID("[v,o]"), ID("[O,V]"), 0, "MO Ints (vo|OV)");
    global_dpd_->contract444(&Lambda_ab, &Iab, &Lbb, 1, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&Iab);
    global_dpd_->buf4_close(&Lambda_ab);

    global_dpd_->buf4_print(&Lbb, outfile, 1);
    global_dpd_->buf4_close(&Lbb);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

}

}} // Namespaces



