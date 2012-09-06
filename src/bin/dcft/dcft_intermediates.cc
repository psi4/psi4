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
DCFTSolver::build_intermediates()
{
    dcft_timer_on("DCFTSolver::build_intermediates()");

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdbuf4 I, L, G, T, Taa, Tab, Tbb, Laa, Lab, Lbb;

    /*
     * G_ijab = <ij||ab>
     */
    // G_IJAB = <IJ||AB>
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O,O]"), ID("[V,V]"), 1, "MO Ints <OO|VV>");
    dpd_buf4_copy(&I, PSIF_DCFT_DPD, "G <OO|VV>");
    dpd_buf4_close(&I);

    // G_IjAb = <Ij|Ab>
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    dpd_buf4_copy(&I, PSIF_DCFT_DPD, "G <Oo|Vv>");
    dpd_buf4_close(&I);

    // G_ijab = <ij||ab>
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o,o]"), ID("[v,v]"), 1, "MO Ints <oo|vv>");
    dpd_buf4_copy(&I, PSIF_DCFT_DPD, "G <oo|vv>");
    dpd_buf4_close(&I);


    /*
     * G_ijab += 1/2 Sum_cd gbar_cdab lambda_ijcd
     */
    if(options_.get_str("AO_BASIS") == "NONE"){
        // G_IJAB += 1/2 Sum_CD gbar_CDAB lambda_IJCD
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V>V]-"), ID("[V>V]-"),
                      ID("[V,V]"), ID("[V,V]"), 1, "MO Ints <VV|VV>");
        dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                      ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
        dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                      ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
        dpd_contract444(&L, &I, &G, 0, 0, 1.0, 1.0);
        dpd_buf4_close(&I);
        dpd_buf4_close(&L);
        dpd_buf4_close(&G);


        // G_IjAb += Sum_Cd g_CdAb lambda_IjCd
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,v]"), ID("[V,v]"),
                      ID("[V,v]"), ID("[V,v]"), 0, "MO Ints <Vv|Vv>");
        dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
        dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
        dpd_contract444(&L, &I, &G, 0, 0, 1.0, 1.0);
        dpd_buf4_close(&I);
        dpd_buf4_close(&L);
        dpd_buf4_close(&G);


        // G_ijab += 1/2 Sum_cd gbar_cdab lambda_ijcd
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v>v]-"), ID("[v>v]-"),
                      ID("[v,v]"), ID("[v,v]"), 1, "MO Ints <vv|vv>");
        dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                      ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
        dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                      ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
        dpd_contract444(&L, &I, &G, 0, 0, 1.0, 1.0);
        dpd_buf4_close(&I);
        dpd_buf4_close(&L);
        dpd_buf4_close(&G);
    }
    else{

        /***********AA***********/
        dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
        dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "tau(temp) <OO|VV>");
        dpd_buf4_axpy(&L, &G, 1.0);
        dpd_buf4_close(&L);
        dpd_buf4_close(&G);

        /***********BB***********/
        dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
        dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o,o]"), ID("[v,v]"), 0, "tau(temp) <oo|vv>");
        dpd_buf4_axpy(&L, &G, 1.0);
        dpd_buf4_close(&L);
        dpd_buf4_close(&G);

        /***********AB***********/
        dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
        dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "tau(temp) <Oo|Vv>");
        dpd_buf4_axpy(&L, &G, 1.0);
        dpd_buf4_close(&L);
        dpd_buf4_close(&G);
    }

    /*
     * G_ijab += 1/2 Sum_kl gbar_ijkl lambda_klab
     */
    // G_IJAB += 1/2 Sum_KL gbar_IJKL lambda_KLAB
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[O>O]-"),
            ID("[O,O]"), ID("[O,O]"), 1, "MO Ints <OO|OO>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    dpd_contract444(&I, &L, &G, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);
    dpd_buf4_close(&G);

    // G_IjAb += Sum_Kl g_IjKl lambda_KlAb
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "MO Ints <Oo|Oo>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
    dpd_contract444(&I, &L, &G, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);
    dpd_buf4_close(&G);

    // G_ijab += 1/2 Sum_kl gbar_ijkl lambda_klab
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o>o]-"), ID("[o>o]-"),
                  ID("[o,o]"), ID("[o,o]"), 1, "MO Ints <oo|oo>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    dpd_contract444(&I, &L, &G, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);
    dpd_buf4_close(&G);


    /*
     * G_ijab -= P(ij)P(ab) Sum_kc gbar_jckb lambda_ikac
     */
    dpd_buf4_init(&Taa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Temp (OV|OV)");
    dpd_buf4_init(&Tab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Temp (Ov|oV)");
    dpd_buf4_init(&Tbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Temp (ov|ov)");
    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    dpd_buf4_sort(&Laa, PSIF_DCFT_DPD, prqs, ID("[O,V]"), ID("[O,V]"), "Lambda (OV|OV)");
    dpd_buf4_close(&Laa);
    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Lambda (OV|OV)");
    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_sort(&Lab, PSIF_DCFT_DPD, psqr, ID("[O,v]"), ID("[o,V]"), "Lambda (Ov|oV)");
    dpd_buf4_close(&Lab);
    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Lambda (Ov|oV)");

    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    dpd_buf4_sort(&Lbb, PSIF_DCFT_DPD, prqs, ID("[o,v]"),ID("[o,v]"), "Lambda (ov|ov)");
    dpd_buf4_close(&Lbb);
    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Lambda (ov|ov)");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "MO Ints <oV|oV>");

    // T_IbjA = -Sum_kC lambda_IbkC g_jAkC
    dpd_contract444(&Lab, &I, &Tab, 0, 0, -1.0, 0.0);

    dpd_buf4_close(&I);
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "MO Ints <Ov|Ov>") ;

    // T_IbjA -= Sum_Kc g_IbKc lambda_KcjA
    dpd_contract444(&I, &Lab, &Tab, 0, 1, -1.0, 1.0);

    // T_IbjA -> T_IAjb
    dpd_buf4_sort(&Tab, PSIF_DCFT_DPD, psrq, ID("[O,V]"),ID("[o,v]"), "Temp (OV|ov)");
    dpd_buf4_close(&Tab);
    dpd_buf4_init(&Tab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Temp (OV|ov)");

    // Lambda_IbkC -> Lambda_ICkb
    dpd_buf4_sort(&Lab, PSIF_DCFT_DPD, psrq, ID("[O,V]"),ID("[o,v]"), "Lambda (OV|ov)");
    dpd_buf4_close(&Lab);
    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");

    dpd_buf4_close(&I);
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");

    // T_IAJB = Sum_kc lambda_IAkc g_JBkc
    dpd_contract444(&Lab, &I, &Taa, 0, 0, 1.0, 0.0);

    // T_iajb = Sum_KC lambda_KCia g_KCjb
    dpd_contract444(&Lab, &I, &Tbb, 1, 1, 1.0, 0.0);

    // T_IAjb += Sum_kc g_IAkc lambda_jbkc
    dpd_contract444(&I, &Lbb, &Tab, 0, 0, 1.0, 1.0);

    // T_IAjb += Sum_KC lambda_IAKC g_KCjb
    dpd_contract444(&Laa, &I, &Tab, 0, 1, 1.0, 1.0);


    dpd_buf4_close(&I);
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>");

    // T_IAJB -= Sum_KC lambda_IAKC g_JBKC
    dpd_contract444(&Laa, &I, &Taa, 0, 0, -1.0, 1.0);

    // T_IAjb -= Sum_KC g_IAKC lambda_KCjb
    dpd_contract444(&I, &Lab, &Tab, 0, 1, -1.0, 1.0);

    dpd_buf4_close(&I);
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");

    // T_IAJB += Sum_KC lambda_IAKC (JB|KC)
    dpd_contract444(&Laa, &I, &Taa, 0, 0, 1.0, 1.0);

    // T_IAjb += Sum_KC (JB|KC) lambda_KCjb
    dpd_contract444(&I, &Lab, &Tab, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov|ov>");

    // T_iajb -= Sum_kc lambda_iakc g_jbkc
    dpd_contract444(&Lbb, &I, &Tbb, 0, 0, -1.0, 1.0);

    // T_IAjb -= Sum_KC lambda_IAkc g_jbkc
    dpd_contract444(&Lab, &I, &Tab, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov|ov)");

    // T_iajb += Sum_kc lambda_iakc (JB|KC)
    dpd_contract444(&Lbb, &I, &Tbb, 0, 0, 1.0, 1.0);

    // T_IAjb += Sum_KC lambda_IAkc (kc|jb)
    dpd_contract444(&Lab, &I, &Tab, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);

    dpd_buf4_close(&Laa);
    dpd_buf4_close(&Lab);
    dpd_buf4_close(&Lbb);

    // T_IAJB -> T_IJAB
    dpd_buf4_sort(&Taa, PSIF_DCFT_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "Temp <OO|VV>");
    dpd_buf4_close(&Taa);
    // G_IJAB += T_IJAB
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    dpd_buf4_add(&G, &T, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    dpd_buf4_init(&Taa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");

    // T_IJAB -> T_JIAB
    dpd_buf4_sort(&Taa, PSIF_DCFT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");



    // G_IJAB -= T_JIAB
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    dpd_buf4_add(&G, &T, -1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    // T_IJAB -> T_IJBA
    dpd_buf4_sort(&Taa, PSIF_DCFT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    // G_IJAB -= T_IJBA
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    dpd_buf4_add(&G, &T, -1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    // T_IJAB -> T_JIBA
    dpd_buf4_sort(&Taa, PSIF_DCFT_DPD, qpsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    // G_IJAB += T_JIBA
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    dpd_buf4_add(&G, &T, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    dpd_buf4_close(&Taa);

    // T_IAjb -> T_IjAb
    dpd_buf4_sort(&Tab, PSIF_DCFT_DPD, prqs, ID("[O,o]"), ID("[V,v]"), "Temp <Oo|Vv>");
    dpd_buf4_close(&Tab);
    // G_IjAb += T_IjAb
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Temp <Oo|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
    dpd_buf4_add(&G, &T, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    // T_iajb -> T_ijab
    dpd_buf4_sort(&Tbb, PSIF_DCFT_DPD, prqs, ID("[o,o]"), ID("[v,v]"), "Temp <oo|vv>");
    dpd_buf4_close(&Tbb);
    // G_ijab += T_ijab
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    dpd_buf4_add(&G, &T, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    dpd_buf4_init(&Tbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");

    // T_ijab -> T_jiab
    dpd_buf4_sort(&Tbb, PSIF_DCFT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    // G_ijab -= T_jiab
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    dpd_buf4_add(&G, &T, -1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    // T_ijab -> T_ijba
    dpd_buf4_sort(&Tbb, PSIF_DCFT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    // G_ijab -= T_ijba
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    dpd_buf4_add(&G, &T, -1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    // T_ijab -> T_jiba
    dpd_buf4_sort(&Tbb, PSIF_DCFT_DPD, qpsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    // G_ijab += T_jiba
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    dpd_buf4_add(&G, &T, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    dpd_buf4_close(&Tbb);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    if (options_.get_str("TAU") == "APPROXIMATE") {
        compute_F_intermediate();
    }
    else {
        compute_refined_F_intermediate();
    }


    dcft_timer_off("DCFTSolver::build_intermediates()");
}

void
DCFTSolver::compute_F_intermediate() {

    dpdfile2 F_OO, F_oo, F_VV, F_vv;
    dpdbuf4 F, T, Laa, Lab, Lbb;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    /*
    * F_ijab += P(ab) F_ca lambda_ijcb - P(ij) F_ki lambda_jkab
    */
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    // Temp_IJAB = lambda_IJCB F_AC
    dpd_file2_init(&F_VV, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");
    dpd_contract244(&F_VV, &Laa, &T, 1, 2, 1, 1.0, 0.0);
    dpd_file2_close(&F_VV);
    dpd_buf4_close(&Laa);
    dpd_buf4_close(&T);
    // F_IJAB = Temp_IJAB
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_copy(&T, PSIF_DCFT_DPD, "F <OO|VV>");
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "F <OO|VV>");
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    dpd_buf4_close(&T);
    // F_IJAB -= Temp_IJBA
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_add(&F, &T, -1.0);
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    // Temp_IJAB = -lambda_KJAB F_IK
    dpd_file2_init(&F_OO, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    dpd_contract244(&F_OO, &Laa, &T, 1, 0, 0, -1.0, 0.0);

    dpd_file2_close(&F_OO);
    dpd_buf4_close(&Laa);
    // F_IJAB += Temp_IJAB
    dpd_buf4_add(&F, &T, 1.0);
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    dpd_buf4_close(&T);
    // F_IJAB -= Temp_JIAB
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_add(&F, &T, -1.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "F <Oo|Vv>");
    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    // F_IjAb += lambda_IjCb F_AC
    dpd_file2_init(&F_VV, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");
    dpd_contract244(&F_VV, &Lab, &F, 1, 2, 1, 1.0, 0.0);
    dpd_file2_close(&F_VV);
    // F_IjAb += lambda_IjAc F_bc
    dpd_file2_init(&F_vv, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "F <v|v>");
    dpd_contract424(&Lab, &F_vv, &F, 3, 1, 0, 1.0, 1.0);
    dpd_file2_close(&F_vv);
    // F_IjAb -= lambda_KjAb F_IK
    dpd_file2_init(&F_OO, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    dpd_contract244(&F_OO, &Lab, &F, 1, 0, 0, -1.0, 1.0);
    dpd_file2_close(&F_OO);
    // F_IjAb -= lambda_IkAb F_jk
    dpd_file2_init(&F_oo, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "F <o|o>");
    dpd_contract424(&Lab, &F_oo, &F, 1, 1, 1, -1.0, 1.0);
    dpd_file2_close(&F_oo);
    dpd_buf4_close(&Lab);
    dpd_buf4_close(&F);

    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    // Temp_ijab = lambda_ijcb F_ac
    dpd_file2_init(&F_vv, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "F <v|v>");
    dpd_contract244(&F_vv, &Lbb, &T, 1, 2, 1, 1.0, 0.0);
    dpd_file2_close(&F_vv);
    dpd_buf4_close(&Lbb);
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    // F_ijab = Temp_ijab
    dpd_buf4_copy(&T, PSIF_DCFT_DPD, "F <oo|vv>");
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    dpd_buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "F <oo|vv>");
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    dpd_buf4_close(&T);
    // F_ijab -= Temp_ijba
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_add(&F, &T, -1.0);
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    // Temp_ijab = -lambda_kjab X_ik
    dpd_file2_init(&F_oo, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "F <o|o>");
    dpd_contract244(&F_oo, &Lbb, &T, 1, 0, 0, -1.0, 0.0);
    dpd_file2_close(&F_oo);
    dpd_buf4_close(&Lbb);
    // F_ijab += Temp_ijab
    dpd_buf4_add(&F, &T, 1.0);
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    dpd_buf4_close(&T);
    // F_ijab -= Temp_jiab
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_add(&F, &T, -1.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&F);

    psio_->close(PSIF_LIBTRANS_DPD, 1);
}

void
DCFTSolver::compute_refined_F_intermediate() {

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdfile2 T_OO, T_oo, T_VV, T_vv, tF_OO, tF_oo, tF_VV, tF_vv, U_OO, U_oo, U_VV, U_vv;
    dpdbuf4 F, T, Laa, Lab, Lba, Lbb;


    dpd_file2_init(&T_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    dpd_file2_init(&T_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
    dpd_file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    dpd_file2_init(&T_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");
    dpd_file2_mat_init(&T_OO);
    dpd_file2_mat_init(&T_oo);
    dpd_file2_mat_init(&T_VV);
    dpd_file2_mat_init(&T_vv);
    dpd_file2_mat_rd(&T_OO);
    dpd_file2_mat_rd(&T_oo);
    dpd_file2_mat_rd(&T_VV);
    dpd_file2_mat_rd(&T_vv);

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

    dpd_file2_close(&T_OO);
    dpd_file2_close(&T_oo);
    dpd_file2_close(&T_VV);
    dpd_file2_close(&T_vv);

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
    SharedMatrix temp (new Matrix ("Temp matrix", nirrep_, nmopi_, nmopi_));

    // Alpha spin
    temp->gemm(true, false, 1.0, a_evecs, moFa_, 0.0);
    nso_Fa->gemm(false, false, 1.0, temp, a_evecs, 0.0);
    // Beta spin
    temp->gemm(true, false, 1.0, b_evecs, moFb_, 0.0);
    nso_Fb->gemm(false, false, 1.0, temp, b_evecs, 0.0);

    // Form half-transformed modified Fock matrix and copy Tau eigenvectors to the DPD file
    dpd_file2_init(&tF_OO, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "tF <O|O'>");
    dpd_file2_init(&tF_oo, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "tF <o|o'>");
    dpd_file2_init(&tF_VV, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "tF <V|V'>");
    dpd_file2_init(&tF_vv, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "tF <v|v'>");
    dpd_file2_init(&U_OO, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "U <O|O'>");
    dpd_file2_init(&U_oo, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "U <o|o'>");
    dpd_file2_init(&U_VV, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "U <V|V'>");
    dpd_file2_init(&U_vv, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "U <v|v'>");
    dpd_file2_mat_init(&tF_OO);
    dpd_file2_mat_init(&tF_oo);
    dpd_file2_mat_init(&tF_VV);
    dpd_file2_mat_init(&tF_vv);
    dpd_file2_mat_init(&U_OO);
    dpd_file2_mat_init(&U_oo);
    dpd_file2_mat_init(&U_VV);
    dpd_file2_mat_init(&U_vv);

    for(int h = 0; h < nirrep_; ++h){
        if(nsopi_[h] == 0) continue;

        // Alpha occupied
        for(int i = 0 ; i < naoccpi_[h]; ++i){
            for(int kp = 0 ; kp < naoccpi_[h]; ++kp){
                double value = 0.0;
                for(int ip = 0 ; ip < naoccpi_[h]; ++ip){
                    value += nso_Fa->get(h, ip, kp) * a_evecs->get(h, i, ip) / (1.0 + a_evals->get(h, ip) + a_evals->get(h, kp));
                }
                tF_OO.matrix[h][i][kp] = value;
                U_OO.matrix[h][i][kp] = a_evecs->get(h, i, kp);
            }
        }

        // Beta occupied
        for(int i = 0 ; i < nboccpi_[h]; ++i){
            for(int kp = 0 ; kp < nboccpi_[h]; ++kp){
                double value = 0.0;
                for(int ip = 0 ; ip < nboccpi_[h]; ++ip){
                    value += nso_Fb->get(h, ip, kp) * b_evecs->get(h, i, ip) / (1.0 + b_evals->get(h, ip) + b_evals->get(h, kp));
                }
                tF_oo.matrix[h][i][kp] = value;
                U_oo.matrix[h][i][kp] = b_evecs->get(h, i, kp);
            }
        }

        // Alpha virtual
        for(int a = 0 ; a < navirpi_[h]; ++a){
            for(int dp = 0 ; dp < navirpi_[h]; ++dp){
                double value = 0.0;
                for(int ap = 0 ; ap < navirpi_[h]; ++ap){
                    value += nso_Fa->get(h, dp + naoccpi_[h], ap + naoccpi_[h]) * a_evecs->get(h, a + naoccpi_[h], ap + naoccpi_[h])
                            / (1.0 - a_evals->get(h, dp + naoccpi_[h]) - a_evals->get(h, ap + naoccpi_[h]));
                }
                tF_VV.matrix[h][a][dp] = value;
                U_VV.matrix[h][a][dp] = a_evecs->get(h, a + naoccpi_[h], dp + naoccpi_[h]);
            }
        }

        // Beta virtual
        for(int a = 0 ; a < nbvirpi_[h]; ++a){
            for(int dp = 0 ; dp < nbvirpi_[h]; ++dp){
                double value = 0.0;
                for(int ap = 0 ; ap < nbvirpi_[h]; ++ap){
                    value += nso_Fb->get(h, dp + nboccpi_[h], ap + nboccpi_[h]) * b_evecs->get(h, a + nboccpi_[h], ap + nboccpi_[h])
                            / (1.0 - b_evals->get(h, dp + nboccpi_[h]) - b_evals->get(h, ap + nboccpi_[h]));
                }
                tF_vv.matrix[h][a][dp] = value;
                U_vv.matrix[h][a][dp] = b_evecs->get(h, a + nboccpi_[h], dp + nboccpi_[h]);
            }
        }

    }

    dpd_file2_mat_wrt(&U_OO);
    dpd_file2_mat_wrt(&U_oo);
    dpd_file2_mat_wrt(&U_VV);
    dpd_file2_mat_wrt(&U_vv);
    dpd_file2_mat_wrt(&tF_OO);
    dpd_file2_mat_wrt(&tF_oo);
    dpd_file2_mat_wrt(&tF_VV);
    dpd_file2_mat_wrt(&tF_vv);

    dpd_file2_close(&tF_OO);
    dpd_file2_close(&tF_oo);
    dpd_file2_close(&tF_VV);
    dpd_file2_close(&tF_vv);
    dpd_file2_close(&U_OO);
    dpd_file2_close(&U_oo);
    dpd_file2_close(&U_VV);
    dpd_file2_close(&U_vv);

    // Perform one-index transformation of the density cumulant to NSO basis
    one_index_transform();

    // Compute the F intermediate using the half-transformed modified Fock matrix and one-index transformed density cumulant

    /*
    * F_ijab += P(ab) tF_ca lambda_ijcb - P(ij) tF_ki lambda_jkab
    */
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V,V]"), 0, "Lambda <OO|V'V>");
    // Temp_IJAB = lambda_IJCB tF_AC
    dpd_file2_init(&tF_VV, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "tF <V|V'>");
    dpd_contract244(&tF_VV, &Laa, &T, 1, 2, 1, 1.0, 0.0);
    dpd_file2_close(&tF_VV);
    dpd_buf4_close(&Laa);
    dpd_buf4_close(&T);
    // F_IJAB = Temp_IJAB
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_copy(&T, PSIF_DCFT_DPD, "F <OO|VV>");
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "F <OO|VV>");
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    dpd_buf4_close(&T);
    // F_IJAB -= Temp_IJBA
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_add(&F, &T, -1.0);
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V>V]-"), 0, "Lambda <O'O|VV>");
    // Temp_IJAB = -lambda_KJAB tF_IK
    dpd_file2_init(&tF_OO, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "tF <O|O'>");
    dpd_contract244(&tF_OO, &Laa, &T, 1, 0, 0, -1.0, 0.0);

    dpd_file2_close(&tF_OO);
    dpd_buf4_close(&Laa);
    // F_IJAB += Temp_IJAB
    dpd_buf4_add(&F, &T, 1.0);
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    dpd_buf4_close(&T);
    // F_IJAB -= Temp_JIAB
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_add(&F, &T, -1.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "F <Oo|Vv>");
    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|V'v>");
    dpd_buf4_init(&Lba, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv'>");
    // F_IjAb += lambda_IjCb tF_AC
    dpd_file2_init(&tF_VV, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "tF <V|V'>");
    dpd_contract244(&tF_VV, &Lab, &F, 1, 2, 1, 1.0, 0.0);
    dpd_file2_close(&tF_VV);
    // F_IjAb += lambda_IjAc tF_bc
    dpd_file2_init(&tF_vv, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "tF <v|v'>");
    dpd_contract424(&Lba, &tF_vv, &F, 3, 1, 0, 1.0, 1.0);
    dpd_file2_close(&tF_vv);
    dpd_buf4_close(&Lba);
    dpd_buf4_close(&Lab);
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "F <Oo|Vv>");
    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <O'o|Vv>");
    dpd_buf4_init(&Lba, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo'|Vv>");
    // F_IjAb -= lambda_KjAb tF_IK
    dpd_file2_init(&tF_OO, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "tF <O|O'>");
    dpd_contract244(&tF_OO, &Lab, &F, 1, 0, 0, -1.0, 1.0);
    dpd_file2_close(&tF_OO);
    // F_IjAb -= lambda_IkAb tF_jk
    dpd_file2_init(&tF_oo, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "tF <o|o'>");
    dpd_contract424(&Lba, &tF_oo, &F, 1, 1, 1, -1.0, 1.0);
    dpd_file2_close(&tF_oo);
    dpd_buf4_close(&Lba);
    dpd_buf4_close(&Lab);
    dpd_buf4_close(&F);

    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v,v]"), 0, "Lambda <oo|v'v>");
    // Temp_ijab = lambda_ijcb tF_ac
    dpd_file2_init(&tF_vv, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "tF <v|v'>");
    dpd_contract244(&tF_vv, &Lbb, &T, 1, 2, 1, 1.0, 0.0);
    dpd_file2_close(&tF_vv);
    dpd_buf4_close(&Lbb);
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    // F_ijab = Temp_ijab
    dpd_buf4_copy(&T, PSIF_DCFT_DPD, "F <oo|vv>");
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    dpd_buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "F <oo|vv>");
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    dpd_buf4_close(&T);
    // F_ijab -= Temp_ijba
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_add(&F, &T, -1.0);
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v>v]-"), 0, "Lambda <o'o|vv>");
    // Temp_ijab = -lambda_kjab tF_ik
    dpd_file2_init(&tF_oo, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "tF <o|o'>");
    dpd_contract244(&tF_oo, &Lbb, &T, 1, 0, 0, -1.0, 0.0);
    dpd_file2_close(&tF_oo);
    dpd_buf4_close(&Lbb);
    // F_ijab += Temp_ijab
    dpd_buf4_add(&F, &T, 1.0);
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    dpd_buf4_close(&T);
    // F_ijab -= Temp_jiab
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_add(&F, &T, -1.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&F);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

}

void
DCFTSolver::one_index_transform() {

    dpdfile2 U;
    dpdbuf4 L, Lt;

    // Lambda_I'JAB = U_KI' Lambda_KJAB
    dpd_buf4_init(&Lt, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V>V]-"), 0, "Lambda <O'O|VV>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    dpd_file2_init(&U, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "U <O|O'>");
    dpd_contract244(&U, &L, &Lt, 0, 0, 0, 1.0, 0.0);
    dpd_file2_close(&U);
    dpd_buf4_close(&L);
    dpd_buf4_close(&Lt);

    // Lambda_i'jab = Lambda_kjab U_ki'
    dpd_buf4_init(&Lt, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v>v]-"), 0, "Lambda <o'o|vv>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    dpd_file2_init(&U, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "U <o|o'>");
    dpd_contract244(&U, &L, &Lt, 0, 0, 0, 1.0, 0.0);
    dpd_file2_close(&U);
    dpd_buf4_close(&L);
    dpd_buf4_close(&Lt);

    // Lambda_I'jAb = Lambda_KjAb U_KI'
    dpd_buf4_init(&Lt, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <O'o|Vv>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_file2_init(&U, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "U <O|O'>");
    dpd_contract244(&U, &L, &Lt, 0, 0, 0, 1.0, 0.0);
    dpd_file2_close(&U);
    dpd_buf4_close(&L);
    dpd_buf4_close(&Lt);

    // Lambda_Ij'Ab = Lambda_IkAb U_kj'
    dpd_buf4_init(&Lt, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo'|Vv>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_file2_init(&U, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "U <o|o'>");
    dpd_contract424(&L, &U, &Lt, 1, 0, 1, 1.0, 0.0);
    dpd_file2_close(&U);
    dpd_buf4_close(&L);
    dpd_buf4_close(&Lt);

    // Lambda_IJA'B = Lambda_IJCB U_CA'
    dpd_buf4_init(&Lt, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V,V]"), 0, "Lambda <OO|V'V>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    dpd_file2_init(&U, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "U <V|V'>");
    dpd_contract244(&U, &L, &Lt, 0, 2, 1, 1.0, 0.0);
    dpd_file2_close(&U);
    dpd_buf4_close(&L);
    dpd_buf4_close(&Lt);

    // Lambda_ija'b = Lambda_ijcb U_ca'
    dpd_buf4_init(&Lt, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v,v]"), 0, "Lambda <oo|v'v>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    dpd_file2_init(&U, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "U <v|v'>");
    dpd_contract244(&U, &L, &Lt, 0, 2, 1, 1.0, 0.0);
    dpd_file2_close(&U);
    dpd_buf4_close(&L);
    dpd_buf4_close(&Lt);

    // Lambda_IjA'b = Lambda_IjCb U_CA'
    dpd_buf4_init(&Lt, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|V'v>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_file2_init(&U, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "U <V|V'>");
    dpd_contract244(&U, &L, &Lt, 0, 2, 1, 1.0, 0.0);
    dpd_file2_close(&U);
    dpd_buf4_close(&L);
    dpd_buf4_close(&Lt);

    // Lambda_IjAb' = Lambda_IjAc U_cb'
    dpd_buf4_init(&Lt, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv'>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_file2_init(&U, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "U <V|V'>");
    dpd_contract424(&L, &U, &Lt, 3, 0, 0, 1.0, 0.0);
    dpd_file2_close(&U);
    dpd_buf4_close(&L);
    dpd_buf4_close(&Lt);

}

}} // Namespaces



