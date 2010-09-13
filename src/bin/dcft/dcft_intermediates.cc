#include "dcft.h"
#include "libdpd/dpd.h"
#include "defines.h"

namespace psi{ namespace dcft{

/**
 * Builds the intermediate tensors
 */
void
DCFTSolver::build_intermediates()
{
    _psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdfile2 X_OO, X_oo, X_VV, X_vv,
             T_OO, T_oo, T_VV, T_vv,
             F_OO, F_oo, F_VV, F_vv, tmp;
    dpdbuf4 I, L, G, T, A, D,
            Taa, Tab, Tbb,
            Laa, Lab, Lbb;
    if(!_options.get_bool("IGNORE_TAU")){
        /*
         * Form the X intermediates from the MO integrals
         */
        dpd_file2_init(&T_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
        dpd_file2_init(&T_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
        dpd_file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
        dpd_file2_init(&T_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");

        dpd_file2_init(&X_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "X <V|V>");

        if(_options.get_str("AO_BASIS") == "NONE"){
            // X_AB = (AB|CD) Tau_CD
            dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                    ID("[V>=V]+"), ID("[V>=V]+"), 0, "MO Ints (VV|VV)");
            dpd_contract422(&I, &T_VV, &X_VV, 0, 0, 1.0, 0.0);
            dpd_buf4_close(&I);
            // X_AB -= <AB|CD> Tau_CD
            dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                    ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV|VV>");
            dpd_contract422(&I, &T_VV, &X_VV, 0, 0, -1.0, 1.0);
            dpd_buf4_close(&I);
            // X_AB += (AB|cd) Tau_cd
            dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[v,v]"),
                    ID("[V>=V]+"), ID("[v>=v]+"), 0, "MO Ints (VV|vv)");
            dpd_contract422(&I, &T_vv, &X_VV, 0, 0, 1.0, 1.0);
            dpd_buf4_close(&I);
        }
        else{ // if(_options.get_str("AO_BASIS") == "DISK"){
            dpd_file2_init(&tmp, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "s(add)A <V|V>");
            dpd_file2_scm(&X_VV, 0.0);
            dpd_file2_axpy(&tmp, &X_VV, 1.0, 0);
            dpd_file2_close(&tmp);
        }

        // X_AB = +(AB|IJ) T_IJ
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,O]"),
                      ID("[V>=V]+"), ID("[O>=O]+"), 0, "MO Ints (VV|OO)");
        dpd_contract422(&I, &T_OO, &X_VV, 0, 0, 1.0, 1.0);
        dpd_buf4_close(&I);
        // X_AB -= <AB|IJ> T_IJ
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,O]"),
                      ID("[V,V]"), ID("[O,O]"), 0, "MO Ints <VV|OO>");
        dpd_contract422(&I, &T_OO, &X_VV, 0, 0, -1.0, 1.0);
        dpd_buf4_close(&I);
        // X_AB += (AB|ij) T_ij
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[o,o]"),
                      ID("[V>=V]+"), ID("[o>=o]+"), 0, "MO Ints (VV|oo)");
        dpd_contract422(&I, &T_oo, &X_VV, 0, 0, 1.0, 1.0);
        dpd_buf4_close(&I);
        dpd_file2_close(&X_VV);


        dpd_file2_init(&X_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "X <v|v>");

        if(_options.get_str("AO_BASIS") == "NONE"){
            // X_ab = +(ab|cd) Tau_cd
            dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                    ID("[v>=v]+"), ID("[v>=v]+"), 0, "MO Ints (vv|vv)");
            dpd_contract422(&I, &T_vv, &X_vv, 0, 0, 1.0, 0.0);
            dpd_buf4_close(&I);
            // X_ab -= <ab|cd> Tau_cd
            dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                    ID("[v,v]"), ID("[v,v]"), 0, "MO Ints <vv|vv>");
            dpd_contract422(&I, &T_vv, &X_vv, 0, 0, -1.0, 1.0);
            dpd_buf4_close(&I);
            // X_ab += (ab|CD) Tau_CD
            dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[V,V]"),
                    ID("[v,v]"), ID("[V,V]"), 0, "MO Ints (vv|VV)");
            dpd_contract422(&I, &T_VV, &X_vv, 0, 0, 1.0, 1.0);
            dpd_buf4_close(&I);
        }
        else{ // if(_options.get_str("AO_BASIS") == "DISK"){
            dpd_file2_init(&tmp, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "s(add)B <v|v>");
            dpd_file2_scm(&X_vv, 0.0);
            dpd_file2_axpy(&tmp, &X_vv, 1.0, 0);
            dpd_file2_close(&tmp);
        }

        // X_ab = +(ab|ij) T_ij
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[o,o]"),
                      ID("[v>=v]+"), ID("[o>=o]+"), 0, "MO Ints (vv|oo)");
        dpd_contract422(&I, &T_oo, &X_vv, 0, 0, 1.0, 1.0);
        dpd_buf4_close(&I);
        // X_ab -= <ab|ij> T_ij
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[o,o]"),
                      ID("[v,v]"), ID("[o,o]"), 0, "MO Ints <vv|oo>");
        dpd_contract422(&I, &T_oo, &X_vv, 0, 0, -1.0, 1.0);
        dpd_buf4_close(&I);
        // X_ab += (ab|IJ) T_IJ
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[O,O]"),
                      ID("[v,v]"), ID("[O,O]"), 0, "MO Ints (vv|OO)");
        dpd_contract422(&I, &T_OO, &X_vv, 0, 0, 1.0, 1.0);
        dpd_buf4_close(&I);
        dpd_file2_close(&X_vv);

        dpd_file2_init(&X_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
        // X_IJ = +(IJ|AB) T_AB
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O>=O]+"), ID("[V>=V]+"), 0, "MO Ints (OO|VV)");
        dpd_contract422(&I, &T_VV, &X_OO, 0, 0, 1.0, 0.0);
        dpd_buf4_close(&I);

        // X_IJ -= <IJ|AB> T_AB
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
        dpd_contract422(&I, &T_VV, &X_OO, 0, 0, -1.0, 1.0);
        dpd_buf4_close(&I);
        // X_IJ += (IJ|ab) T_ab
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[v,v]"),
                      ID("[O>=O]+"), ID("[v>=v]+"), 0, "MO Ints (OO|vv)");
        dpd_contract422(&I, &T_vv, &X_OO, 0, 0, 1.0, 1.0);
        dpd_buf4_close(&I);
        // X_IJ = +(IJ|KL) T_KL
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                      ID("[O>=O]+"), ID("[O>=O]+"), 0, "MO Ints (OO|OO)");
        dpd_contract422(&I, &T_OO, &X_OO, 0, 0, 1.0, 1.0);
        dpd_buf4_close(&I);
        // X_IJ -= <IJ|KL> T_KL
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                      ID("[O,O]"), ID("[O,O]"), 0, "MO Ints <OO|OO>");
        dpd_contract422(&I, &T_OO, &X_OO, 0, 0, -1.0, 1.0);
        dpd_buf4_close(&I);
        // X_IJ += (IJ|kl) T_kl
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[o,o]"),
                      ID("[O>=O]+"), ID("[o>=o]+"), 0, "MO Ints (OO|oo)");
        dpd_contract422(&I, &T_oo, &X_OO, 0, 0, 1.0, 1.0);
        dpd_buf4_close(&I);
        dpd_file2_close(&X_OO);


        dpd_file2_init(&X_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "X <o|o>");
        // X_ij = +(ij|ab) T_ab
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o>=o]+"), ID("[v>=v]+"), 0, "MO Ints (oo|vv)");
        dpd_contract422(&I, &T_vv, &X_oo, 0, 0, 1.0, 0.0);
        dpd_buf4_close(&I);
        // X_ij -= <ij|ab> T_ab
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo|vv>");
        dpd_contract422(&I, &T_vv, &X_oo, 0, 0, -1.0, 1.0);
        dpd_buf4_close(&I);
        // X_ij += (ij|AB) T_AB
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[V,V]"),
                      ID("[o,o]"), ID("[V,V]"), 0, "MO Ints (oo|VV)");
        dpd_contract422(&I, &T_VV, &X_oo, 0, 0, 1.0, 1.0);
        dpd_buf4_close(&I);
        // X_ij = +(ij|kl) T_kl
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                      ID("[o>=o]+"), ID("[o>=o]+"), 0, "MO Ints (oo|oo)");
        dpd_contract422(&I, &T_oo, &X_oo, 0, 0, 1.0, 1.0);
        dpd_buf4_close(&I);
        // X_ij -= <ij|kl> T_kl
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                      ID("[o,o]"), ID("[o,o]"), 0, "MO Ints <oo|oo>");
        dpd_contract422(&I, &T_oo, &X_oo, 0, 0, -1.0, 1.0);
        dpd_buf4_close(&I);
        // X_IJ += (ij|KL) T_KL
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[O,O]"),
                      ID("[o,o]"), ID("[O,O]"), 0, "MO Ints (oo|OO)");
        dpd_contract422(&I, &T_OO, &X_oo, 0, 0, 1.0, 1.0);
        dpd_buf4_close(&I);
        dpd_file2_close(&X_oo);

        dpd_file2_close(&T_OO);
        dpd_file2_close(&T_oo);
        dpd_file2_close(&T_VV);
        dpd_file2_close(&T_vv);
        // End of the X intermediates


        /*
         * The T intermediates
         * T_ijab = P(ab) X_ca lambda_ijcb + P(ij) X_kj lambda_ikab
         */
        dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
        dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
        // Temp_IJAB = lambda_IJCB X_AC
        dpd_file2_init(&X_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "X <V|V>");
        dpd_contract244(&X_VV, &Laa, &T, 1, 2, 1, 1.0, 0.0);
        dpd_file2_close(&X_VV);
        dpd_buf4_close(&Laa);
        // T_IJAB = Temp_IJAB
        dpd_buf4_copy(&T, PSIF_DCFT_DPD, "T <OO|VV>");
        dpd_buf4_sort(&T, PSIF_DCFT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
        dpd_buf4_close(&T);
        // T_IJAB -= Temp_IJBA
        dpd_buf4_init(&Taa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "T <OO|VV>");
        dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
        dpd_buf4_add(&Taa, &T, -1.0);
        dpd_buf4_close(&T);
        dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
        dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
        // Temp_IJAB = -lambda_KJAB X_IK
        dpd_file2_init(&X_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
        dpd_contract244(&X_OO, &Laa, &T, 1, 0, 0, -1.0, 0.0);
        dpd_file2_close(&X_OO);
        dpd_buf4_close(&Laa);
        // T_IJAB += TEMP_IJAB
        dpd_buf4_add(&Taa, &T, 1.0);
        dpd_buf4_sort(&T, PSIF_DCFT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
        dpd_buf4_close(&T);
        // T_IJAB -= Temp_JIAB
        dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
        dpd_buf4_add(&Taa, &T, -1.0);
        dpd_buf4_close(&T);
        dpd_buf4_close(&Taa);

        // The Tab intermediate
        dpd_buf4_init(&Tab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "T <Oo|Vv>");
        dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
        // T_IjAb = lambda_IjCb X_AC
        dpd_file2_init(&X_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "X <V|V>");
        dpd_contract244(&X_VV, &Lab, &Tab, 1, 2, 1, 1.0, 0.0);
        dpd_file2_close(&X_VV);
        // T_IjAb += lambda_IjAc X_bc
        dpd_file2_init(&X_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "X <v|v>");
        dpd_contract424(&Lab, &X_vv, &Tab, 3, 1, 0, 1.0, 1.0);
        dpd_file2_close(&X_vv);
        // T_IjAb -= lambda_KjAb X_IK
        dpd_file2_init(&X_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
        dpd_contract244(&X_OO, &Lab, &Tab, 1, 0, 0, -1.0, 1.0);
        dpd_file2_close(&X_OO);
        // T_IjAb -= lambda_IkAb X_jk
        dpd_file2_init(&X_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "X <o|o>");
        dpd_contract424(&Lab, &X_oo, &Tab, 1, 1, 1, -1.0, 1.0);
        dpd_file2_close(&X_oo);
        dpd_buf4_close(&Lab);
        dpd_buf4_close(&Tab);

        // The Tbb intermediate
        dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
        dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
        // Temp_ijab = lambda_ijcb X_ac
        dpd_file2_init(&X_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "X <v|v>");
        dpd_contract244(&X_vv, &Lbb, &T, 1, 2, 1, 1.0, 0.0);
        dpd_file2_close(&X_vv);
        dpd_buf4_close(&Lbb);
        // T_ijab = Temp_ijab
        dpd_buf4_copy(&T, PSIF_DCFT_DPD, "T <oo|vv>");
        dpd_buf4_sort(&T, PSIF_DCFT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
        dpd_buf4_close(&T);

        // T_ijab -= Temp_ijba
        dpd_buf4_init(&Tbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o,o]"), ID("[v,v]"), 0, "T <oo|vv>");
        dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
        dpd_buf4_add(&Tbb, &T, -1.0);
        dpd_buf4_close(&T);
        dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
        dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
        // Temp_ijab = -lambda_kjab X_ik
        dpd_file2_init(&X_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "X <o|o>");
        dpd_contract244(&X_oo, &Lbb, &T, 1, 0, 0, -1.0, 0.0);
        dpd_file2_close(&X_oo);
        dpd_buf4_close(&Lbb);
        // T_ijab += TEMP_ijab
        dpd_buf4_add(&Tbb, &T, 1.0);
        dpd_buf4_sort(&T, PSIF_DCFT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
        dpd_buf4_close(&T);
        // T_ijab -= Temp_jiab
        dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
        dpd_buf4_add(&Tbb, &T, -1.0);
        dpd_buf4_close(&T);
        dpd_buf4_close(&Tbb);
        //End of the T intermediates

    }// End "if ignore tau"





    /*
     * G_ijab = <ij||ab>
     */
    // G_IJAB = <IJ||AB>
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 1, "MO Ints <OO|VV>");
    dpd_buf4_copy(&I, PSIF_DCFT_DPD, "G <OO|VV>");
    dpd_buf4_close(&I);

    // G_IjAb = <Ij|Ab>
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    dpd_buf4_copy(&I, PSIF_DCFT_DPD, "G <Oo|Vv>");
    dpd_buf4_close(&I);

    // G_ijab = <ij||ab>
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 1, "MO Ints <oo|vv>");
    dpd_buf4_copy(&I, PSIF_DCFT_DPD, "G <oo|vv>");
    dpd_buf4_close(&I);


    /*
     * G_ijab += 1/2 Sum_cd gbar_cdab lambda_ijcd
     */
    if(_options.get_str("AO_BASIS") == "NONE"){
        // G_IJAB += 1/2 Sum_CD gbar_CDAB lambda_IJCD
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                      ID("[V,V]"), ID("[V,V]"), 1, "MO Ints <VV|VV>");
        dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
        dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "G <OO|VV>");
        dpd_contract444(&L, &I, &G, 0, 0, 0.5, 1.0);
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
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                      ID("[v,v]"), ID("[v,v]"), 1, "MO Ints <vv|vv>");
        dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
        dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o,o]"), ID("[v,v]"), 0, "G <oo|vv>");
        dpd_contract444(&L, &I, &G, 0, 0, 0.5, 1.0);
        dpd_buf4_close(&I);
        dpd_buf4_close(&L);
        dpd_buf4_close(&G);


    }
    else{

        /***********AA***********/
        dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "G <OO|VV>");
        dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "tau(temp) <OO|VV>");
        dpd_buf4_axpy(&L, &G, 1.0);
        dpd_buf4_close(&L);
        dpd_buf4_close(&G);

        /***********BB***********/
        dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o,o]"), ID("[v,v]"), 0, "G <oo|vv>");
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
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
            ID("[O,O]"), ID("[O,O]"), 1, "MO Ints <OO|OO>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "G <OO|VV>");
    dpd_contract444(&I, &L, &G, 0, 1, 0.5, 1.0);
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
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 1, "MO Ints <oo|oo>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "G <oo|vv>");
    dpd_contract444(&I, &L, &G, 0, 1, 0.5, 1.0);
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
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
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
                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
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
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV> - (OV|OV)");

    // T_IAJB -= Sum_KC lambda_IAKC gbar_JBKC
    dpd_contract444(&Laa, &I, &Taa, 0, 0, -1.0, 1.0);

    // T_IAjb -= Sum_KC gbar_IAKC lambda_KCjb
    dpd_contract444(&I, &Lab, &Tab, 0, 1, -1.0, 1.0);

    dpd_buf4_close(&I);
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov|ov> - (ov|ov)");

    // T_iajb -= Sum_kc lambda_iakc gbar_jbkc
    dpd_contract444(&Lbb, &I, &Tbb, 0, 0, -1.0, 1.0);

    // T_IAjb -= Sum_KC lambda_IAkc gbar_jbkc
    dpd_contract444(&Lab, &I, &Tab, 0, 0, -1.0, 1.0);

    dpd_buf4_close(&I);
    dpd_buf4_close(&Laa);
    dpd_buf4_close(&Lab);
    dpd_buf4_close(&Lbb);

    // T_IAJB -> T_IJAB
    dpd_buf4_sort(&Taa, PSIF_DCFT_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "Temp <OO|VV>");
    dpd_buf4_close(&Taa);
    // G_IJAB += T_IJAB
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "G <OO|VV>");
    dpd_buf4_add(&G, &T, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    dpd_buf4_init(&Taa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");

    // T_IJAB -> T_JIAB
    dpd_buf4_sort(&Taa, PSIF_DCFT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    // G_IJAB -= T_JIAB
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "G <OO|VV>");
    dpd_buf4_add(&G, &T, -1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    // T_IJAB -> T_IJBA
    dpd_buf4_sort(&Taa, PSIF_DCFT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    // G_IJAB -= T_IJBA
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "G <OO|VV>");
    dpd_buf4_add(&G, &T, -1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    // T_IJAB -> T_JIBA
    dpd_buf4_sort(&Taa, PSIF_DCFT_DPD, qpsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    // G_IJAB += T_JIBA
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "G <OO|VV>");
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
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "G <oo|vv>");
    dpd_buf4_add(&G, &T, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    dpd_buf4_init(&Tbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");

    // T_ijab -> T_jiab
    dpd_buf4_sort(&Tbb, PSIF_DCFT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    // G_ijab -= T_jiab
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "G <oo|vv>");
    dpd_buf4_add(&G, &T, -1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    // T_ijab -> T_ijba
    dpd_buf4_sort(&Tbb, PSIF_DCFT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    // G_ijab -= T_ijba
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "G <oo|vv>");
    dpd_buf4_add(&G, &T, -1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    // T_ijab -> T_jiba
    dpd_buf4_sort(&Tbb, PSIF_DCFT_DPD, qpsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    // G_ijab += T_jiba
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "G <oo|vv>");
    dpd_buf4_add(&G, &T, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    dpd_buf4_close(&Tbb);

    /*
     * G_ijab += P(ab) F0_ca lambda_ijcb - P(ij) F_ki lambda_jkab
     */
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
    // Temp_IJAB = lambda_IJCB F0_AC
    dpd_file2_init(&F_VV, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F0 <V|V>");
    dpd_contract244(&F_VV, &Laa, &T, 1, 2, 1, 1.0, 0.0);
    dpd_file2_close(&F_VV);
    dpd_buf4_close(&Laa);
    // G_IJAB += Temp_IJAB
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "G <OO|VV>");
    dpd_buf4_add(&G, &T, 1.0);
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    dpd_buf4_close(&T);
    // G_IJAB -= Temp_IJBA
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_add(&G, &T, -1.0);
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
              ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
    // Temp_IJAB = -lambda_KJAB X_IK
    dpd_file2_init(&F_OO, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F0 <O|O>");
    dpd_contract244(&F_OO, &Laa, &T, 1, 0, 0, -1.0, 0.0);
    dpd_file2_close(&F_OO);
    dpd_buf4_close(&Laa);
    // G_IJAB += Temp_IJAB
    dpd_buf4_add(&G, &T, 1.0);
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    dpd_buf4_close(&T);
    // G_IJAB -= Temp_JIAB
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_add(&G, &T, -1.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    // G_IjAb += lambda_IjCb F0_AC
    dpd_file2_init(&F_VV, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F0 <V|V>");
    dpd_contract244(&F_VV, &Lab, &G, 1, 2, 1, 1.0, 1.0);
    dpd_file2_close(&F_VV);
    // G_IjAb += lambda_IjAc F0_bc
    dpd_file2_init(&F_vv, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "F0 <v|v>");
    dpd_contract424(&Lab, &F_vv, &G, 3, 1, 0, 1.0, 1.0);
    dpd_file2_close(&F_vv);
    // G_IjAb -= lambda_KjAb F0_IK
    dpd_file2_init(&F_OO, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F0 <O|O>");
    dpd_contract244(&F_OO, &Lab, &G, 1, 0, 0, -1.0, 1.0);
    dpd_file2_close(&F_OO);
    // G_IjAb -= lambda_IkAb F0_jk
    dpd_file2_init(&F_oo, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "F0 <o|o>");
    dpd_contract424(&Lab, &F_oo, &G, 1, 1, 1, -1.0, 1.0);
    dpd_file2_close(&F_oo);
    dpd_buf4_close(&Lab);
    dpd_buf4_close(&G);

    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
    // Temp_ijab = lambda_ijcb F0_ac
    dpd_file2_init(&F_vv, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "F0 <v|v>");
    dpd_contract244(&F_vv, &Lbb, &T, 1, 2, 1, 1.0, 0.0);
    dpd_file2_close(&F_vv);
    dpd_buf4_close(&Lbb);
    // G_ijab = Temp_ijab
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "G <oo|vv>");
    dpd_buf4_add(&G, &T, 1.0);
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    dpd_buf4_close(&T);
    // G_ijab -= Temp_ijba
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_add(&G, &T, -1.0);
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
              ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
    // Temp_ijab = -lambda_kjab X_ik
    dpd_file2_init(&F_oo, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "F0 <o|o>");
    dpd_contract244(&F_oo, &Lbb, &T, 1, 0, 0, -1.0, 0.0);
    dpd_file2_close(&F_oo);
    dpd_buf4_close(&Lbb);
    // G_ijab += Temp_ijab
    dpd_buf4_add(&G, &T, 1.0);
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    dpd_buf4_close(&T);
    // G_ijab -= Temp_jiab
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_add(&G, &T, -1.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&G);
    


    /*
     * A_ijab = lambda_ijab / D_ijab
     */
    // A_IJAB = lambda_IJAB / D_IJAB
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "D <OO|VV>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
    dpd_buf4_init(&A, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "A <OO|VV>");
    for(int h = 0; h < _nIrreps; ++h){
        dpd_buf4_mat_irrep_init(&D, h);
        dpd_buf4_mat_irrep_init(&L, h);
        dpd_buf4_mat_irrep_init(&A, h);
        dpd_buf4_mat_irrep_rd(&L, h);
        dpd_buf4_mat_irrep_rd(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            for(int col = 0; col < D.params->coltot[h]; ++col){
                A.matrix[h][row][col] = L.matrix[h][row][col] / D.matrix[h][row][col];
            }
        }
        dpd_buf4_mat_irrep_wrt(&A, h);
        dpd_buf4_mat_irrep_close(&D, h);
        dpd_buf4_mat_irrep_close(&L, h);
        dpd_buf4_mat_irrep_close(&A, h);
    }
    dpd_buf4_close(&D);
    dpd_buf4_close(&L);
    dpd_buf4_close(&A);

    // A_IjAb = lambda_IjAb / D_IjAb
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "D <Oo|Vv>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_init(&A, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "A <Oo|Vv>");
    for(int h = 0; h < _nIrreps; ++h){
        dpd_buf4_mat_irrep_init(&D, h);
        dpd_buf4_mat_irrep_init(&L, h);
        dpd_buf4_mat_irrep_init(&A, h);
        dpd_buf4_mat_irrep_rd(&L, h);
        dpd_buf4_mat_irrep_rd(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            for(int col = 0; col < D.params->coltot[h]; ++col){
                A.matrix[h][row][col] = L.matrix[h][row][col] / D.matrix[h][row][col];
            }
        }
        dpd_buf4_mat_irrep_wrt(&A, h);
        dpd_buf4_mat_irrep_close(&D, h);
        dpd_buf4_mat_irrep_close(&L, h);
        dpd_buf4_mat_irrep_close(&A, h);
    }
    dpd_buf4_close(&D);
    dpd_buf4_close(&L);
    dpd_buf4_close(&A);

    // A_ijab = lambda_ijab / D_ijab
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "D <oo|vv>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
    dpd_buf4_init(&A, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "A <oo|vv>");
    for(int h = 0; h < _nIrreps; ++h){
        dpd_buf4_mat_irrep_init(&D, h);
        dpd_buf4_mat_irrep_init(&L, h);
        dpd_buf4_mat_irrep_init(&A, h);
        dpd_buf4_mat_irrep_rd(&L, h);
        dpd_buf4_mat_irrep_rd(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            for(int col = 0; col < D.params->coltot[h]; ++col){
                A.matrix[h][row][col] = L.matrix[h][row][col] / D.matrix[h][row][col];
            }
        }
        dpd_buf4_mat_irrep_wrt(&A, h);
        dpd_buf4_mat_irrep_close(&D, h);
        dpd_buf4_mat_irrep_close(&L, h);
        dpd_buf4_mat_irrep_close(&A, h);
    }
    dpd_buf4_close(&D);
    dpd_buf4_close(&L);
    dpd_buf4_close(&A);

    _psio->close(PSIF_LIBTRANS_DPD, 1);
}

}} // Namespaces



