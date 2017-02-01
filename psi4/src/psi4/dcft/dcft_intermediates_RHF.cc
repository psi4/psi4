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

#include "dcft.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libtrans/integraltransform.h"
#include "defines.h"

namespace psi{ namespace dcft{

/**
 * Builds the intermediate tensors
 */
void
DCFTSolver::build_cumulant_intermediates_RHF()
{
    dcft_timer_on("DCFTSolver::build_intermediates()");

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    dpdbuf4 I, L, G, T;

    /*
     * G_IjAb = <Ij||Ab>
     */
    dcft_timer_on("DCFTSolver::copy <Ij|Ab>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>"); // MO Ints <Oo|Vv>
    global_dpd_->buf4_copy(&I, PSIF_DCFT_DPD, "G <OO|VV>"); // G <Oo|Vv>
    global_dpd_->buf4_close(&I);
    dcft_timer_off("DCFTSolver::copy <Ij|Ab>");

    /*
     * G_IjAb += Sum_Cd gbar_CdAb lambda_IjCd
     */
    dcft_timer_on("DCFTSolver::g_AbCd lambda_IjCd");
    if (options_.get_str("AO_BASIS") == "NONE" && options_.get_str("DCFT_TYPE") == "CONV"){
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                                ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV|VV>"); // MO Ints <Vv|Vv>
        global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                ID("[O,O]"), ID("[V,V]"), 0, "Lambda SF <OO|VV>"); // Lambda <Oo|Vv>
        global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                ID("[O,O]"), ID("[V,V]"), 0, "G <OO|VV>"); // G <Oo|Vv>
        global_dpd_->contract444(&L, &I, &G, 0, 0, 1.0, 1.0);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&L);
        global_dpd_->buf4_close(&G);
    }
    else{
        global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                               ID("[O,O]"), ID("[V,V]"), 0, "G <OO|VV>"); // G <Oo|Vv>
        global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                               ID("[O,O]"), ID("[V,V]"), 0, "tau(temp) SF <OO|VV>"); // tau(temp) <Oo|Vv>
        global_dpd_->buf4_axpy(&L, &G, 1.0);
        global_dpd_->buf4_close(&L);
        global_dpd_->buf4_close(&G);
    }
    dcft_timer_off("DCFTSolver::g_AbCd lambda_IjCd");

    /*
     * G_IjAb += Sum_Kl gbar_IjKl lambda_KlAb
     */
    dcft_timer_on("DCFTSolver::g_IjKl lambda_KlAb");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "MO Ints <OO|OO>"); // MO Ints <Oo|Oo>
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "Lambda SF <OO|VV>"); // Lambda <Oo|Vv>
    global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "G <OO|VV>"); // G <Oo|Vv>
    global_dpd_->contract444(&I, &L, &G, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&G);
    dcft_timer_off("DCFTSolver::g_IjKl lambda_KlAb");

    /*
     * G_ijab -= P(ij)P(ab) Sum_kc gbar_jckb lambda_ikac
     */
    dcft_timer_on("DCFTSolver::g_JcKb lambda_IkAc (4 times)")
    dpdbuf4 Laa, Lbb, Lab, Tab;

    global_dpd_->buf4_init(&Tab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Temp SF (OV|OV):(Ov|oV)");

    global_dpd_->buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
    global_dpd_->buf4_sort(&Laa, PSIF_DCFT_DPD, prqs, ID("[O,V]"), ID("[O,V]"), "Lambda (OV|OV)");
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Lambda (OV|OV)");

    global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "Lambda SF <OO|VV>");
    global_dpd_->buf4_sort(&Lab, PSIF_DCFT_DPD, psqr, ID("[O,V]"), ID("[O,V]"), "Lambda SF (OV|OV):(Ov|oV)");
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                           ID("[O,V]"), ID("[O,V]"), 0, "Lambda SF (OV|OV):(Ov|oV)");

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                           ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>");
    // In RHF, g_jAkC = g_JAKC
    // T_IbjA = -Sum_kC lambda_IbkC g_jAkC
    global_dpd_->contract444(&Lab, &I, &Tab, 0, 0, -1.0, 0.0);

    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>") ; // MO Ints <Ov|Ov>

    // T_IbjA -= Sum_Kc g_IbKc lambda_KcjA
    global_dpd_->contract444(&I, &Lab, &Tab, 0, 1, -1.0, 1.0);

    // T_IbjA -> T_IAjb
    global_dpd_->buf4_sort(&Tab, PSIF_DCFT_DPD, psrq, ID("[O,V]"),ID("[O,V]"), "Temp SF (OV|OV):(OV|ov)");
    global_dpd_->buf4_close(&Tab);
    global_dpd_->buf4_init(&Tab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Temp SF (OV|OV):(OV|ov)");

    // Lambda_IbkC -> Lambda_ICkb
    global_dpd_->buf4_sort(&Lab, PSIF_DCFT_DPD, psrq, ID("[O,V]"),ID("[O,V]"), "Lambda SF (OV|OV):(OV|ov)");

    global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                           ID("[O,V]"), ID("[O,V]"), 0, "Lambda SF (OV|OV):(OV|ov)");

    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");

    // In RHF, g_IAkc = g_IAKC

    // T_IAjb += Sum_kc g_IAkc lambda_jbkc
    global_dpd_->contract444(&I, &Laa, &Tab, 0, 0, 1.0, 1.0);

    // T_IAjb += Sum_KC lambda_IAKC g_KCjb
    global_dpd_->contract444(&Laa, &I, &Tab, 0, 1, 1.0, 1.0);

    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>");

    // T_IAjb -= Sum_KC g_IAKC lambda_KCjb
    global_dpd_->contract444(&I, &Lab, &Tab, 0, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");

    // T_IAjb += Sum_KC (JB|KC) lambda_KCjb
    global_dpd_->contract444(&I, &Lab, &Tab, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>"); // MO Ints <ov|ov>

    // T_IAjb -= Sum_KC lambda_IAkc g_jbkc
    global_dpd_->contract444(&Lab, &I, &Tab, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");

    // In RHF, (kc|jb) = (KC|JB)

    // T_IAjb += Sum_KC lambda_IAkc (kc|jb)
    global_dpd_->contract444(&Lab, &I, &Tab, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&Lab);

    // T_IAjb -> T_IjAb
    global_dpd_->buf4_sort(&Tab, PSIF_DCFT_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "Temp <OO|VV>"); // Temp <Oo|Vv>

    global_dpd_->buf4_close(&Tab);
    // G_IjAb += T_IjAb
    global_dpd_->buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>"); // Temp <Oo|Vv>
    global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "G <OO|VV>"); // G <Oo|Vv>
    dpd_buf4_add(&G, &T, 1.0);

    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&T);
    dcft_timer_off("DCFTSolver::g_JcKb lambda_IkAc (4 times)")

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    if (exact_tau_) {
        form_density_weighted_fock_RHF();
    }
    compute_F_intermediate_RHF();

    dcft_timer_off("DCFTSolver::build_intermediates()");

}

void
DCFTSolver::compute_F_intermediate_RHF() {
    dpdfile2 F_OO, F_oo, F_VV, F_vv;
    dpdbuf4 F, Lab;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    /*
    * F_ijab += P(ab) F_ca lambda_ijcb - P(ij) F_ki lambda_jkab
    */
    global_dpd_->buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "F <OO|VV>"); // F <Oo|Vv>
    global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda SF <OO|VV>"); // Lambda <Oo|Vv>

    // F_IjAb += lambda_IjCb F_AC
    global_dpd_->file2_init(&F_VV, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");
    global_dpd_->contract244(&F_VV, &Lab, &F, 1, 2, 1, 1.0, 0.0);
    global_dpd_->file2_close(&F_VV);
    // F_IjAb += lambda_IjAc F_bc
    global_dpd_->file2_init(&F_vv, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>"); // F <v|v>
    global_dpd_->contract424(&Lab, &F_vv, &F, 3, 1, 0, 1.0, 1.0);
    global_dpd_->file2_close(&F_vv);
    // F_IjAb -= lambda_KjAb F_IK
    global_dpd_->file2_init(&F_OO, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    global_dpd_->contract244(&F_OO, &Lab, &F, 1, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&F_OO);
    // F_IjAb -= lambda_IkAb F_jk
    global_dpd_->file2_init(&F_oo, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>"); // F <o|o>
    global_dpd_->contract424(&Lab, &F_oo, &F, 1, 1, 1, -1.0, 1.0);
    global_dpd_->file2_close(&F_oo);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&F);

    psio_->close(PSIF_LIBTRANS_DPD, 1);
}

void
DCFTSolver::form_density_weighted_fock_RHF() {
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdfile2 T_OO, T_VV, F_OO, F_VV;

    global_dpd_->file2_init(&T_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    global_dpd_->file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    global_dpd_->file2_mat_init(&T_OO);
    global_dpd_->file2_mat_init(&T_VV);
    global_dpd_->file2_mat_rd(&T_OO);
    global_dpd_->file2_mat_rd(&T_VV);

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

        // Alpha virtual
        for(int p = 0 ; p < navirpi_[h]; ++p){
            for(int q = 0 ; q <= p; ++q){
                double value = T_VV.matrix[h][p][q];
                a_tau_mo->set(h, p + naoccpi_[h], q + naoccpi_[h], value);
                if (p != q) a_tau_mo->set(h, q + naoccpi_[h], p + naoccpi_[h], value);
            }
        }
    }

    b_tau_mo->copy(a_tau_mo);

    global_dpd_->file2_close(&T_OO);
    global_dpd_->file2_close(&T_VV);

    SharedMatrix a_evecs (new Matrix ("Tau Eigenvectors (Alpha)", nirrep_, nmopi_, nmopi_));
    SharedVector a_evals (new Vector ("Tau Eigenvalues (Alpha)", nirrep_, nmopi_));

    // Diagonalize Tau
    a_tau_mo->diagonalize(a_evecs, a_evals);
    a_tau_mo->zero();
    a_tau_mo->set_diagonal(a_evals);

    // Transform the Fock matrix to NSO basis
    SharedMatrix nso_Fa (new Matrix ("Alpha Fock in the NSO basis", nirrep_, nmopi_, nmopi_));

    // Alpha spin
    nso_Fa->transform(moFa_, a_evecs);

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

        // Alpha virtual
        for(int ap = naoccpi_[h]; ap < nmopi_[h]; ++ap){
            for(int dp = naoccpi_[h]; dp < nmopi_[h]; ++dp){
                double value = nso_Fa->get(h, dp, ap) / (1.0 - a_evals->get(h, dp) - a_evals->get(h, ap));
                nso_Fa->set(h, dp, ap, value);
            }
        }
    }

    // Transform density-weighted density matrix back to the original basis

    // Alpha spin
    nso_Fa->back_transform(a_evecs);

    global_dpd_->file2_init(&F_OO, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    global_dpd_->file2_init(&F_VV, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");

    global_dpd_->file2_mat_init(&F_OO);
    global_dpd_->file2_mat_init(&F_VV);

    for(int h = 0; h < nirrep_; ++h){
        if(nsopi_[h] == 0) continue;

        // Alpha occupied
        for(int ip = 0 ; ip < naoccpi_[h]; ++ip){
            for(int kp = 0 ; kp < naoccpi_[h]; ++kp){
                F_OO.matrix[h][ip][kp] = nso_Fa->get(h, ip, kp);
            }
        }

        // Alpha virtual
        for(int ap = 0 ; ap < navirpi_[h]; ++ap){
            for(int dp = 0 ; dp < navirpi_[h]; ++dp){
                F_VV.matrix[h][ap][dp] = nso_Fa->get(h, ap + naoccpi_[h], dp + naoccpi_[h]);
                if (ap == dp) F_VV.matrix[h][ap][dp] += energy_level_shift_;
            }
        }
    }

    global_dpd_->file2_mat_wrt(&F_OO);
    global_dpd_->file2_mat_wrt(&F_VV);
    global_dpd_->file2_close(&F_OO);
    global_dpd_->file2_close(&F_VV);

    psio_->close(PSIF_LIBTRANS_DPD, 1);
}

}}// Namespace
