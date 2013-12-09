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

}} // Namespaces



