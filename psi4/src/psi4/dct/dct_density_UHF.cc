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

#include "dct.h"
#include "psi4/psifiles.h"

#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/molecule.h"

namespace psi {
namespace dct {

void DCTSolver::compute_unrelaxed_density_OOOO(bool cumulant_only) {
    dpdbuf4 Iaa, Iab, Ibb, Gaa, Gab, Gbb;

    // Compute the N^6 terms for Gamma OOOO

    if (options_.get_str("DCT_FUNCTIONAL") != "ODC-13") {
        compute_I_intermediate();
    }

    const std::string density_variable = cumulant_only ? "Lambda " : "Gamma ";
    auto varname = [&density_variable](const std::string& x) { return (density_variable + x); };

    psio_->open(PSIF_DCT_DENSITY, PSIO_OPEN_OLD);

    // Gamma_ijkl = 1/8 * I_ijkl
    global_dpd_->buf4_init(&Iaa, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[O>O]-"), ID("[O>O]-"), ID("[O>O]-"), 0,
                           "I <OO|OO>");
    global_dpd_->buf4_copy(&Iaa, PSIF_DCT_DENSITY, varname("<OO|OO>"));
    global_dpd_->buf4_close(&Iaa);

    global_dpd_->buf4_init(&Iab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), 0, "I <Oo|Oo>");
    global_dpd_->buf4_copy(&Iab, PSIF_DCT_DENSITY, varname("<Oo|Oo>"));
    global_dpd_->buf4_close(&Iab);

    global_dpd_->buf4_init(&Ibb, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[o>o]-"), ID("[o>o]-"), ID("[o>o]-"), 0,
                           "I <oo|oo>");
    global_dpd_->buf4_copy(&Ibb, PSIF_DCT_DENSITY, varname("<oo|oo>"));
    global_dpd_->buf4_close(&Ibb);

    global_dpd_->buf4_init(&Gaa, PSIF_DCT_DENSITY, 0, ID("[O>O]-"), ID("[O>O]-"), ID("[O>O]-"), ID("[O>O]-"), 0,
                           varname("<OO|OO>"));
    global_dpd_->buf4_scm(&Gaa, 1.0 / 8.0);
    global_dpd_->buf4_close(&Gaa);

    global_dpd_->buf4_init(&Gab, PSIF_DCT_DENSITY, 0, ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), 0,
                           varname("<Oo|Oo>"));
    global_dpd_->buf4_scm(&Gab, 1.0 / 8.0);
    global_dpd_->buf4_close(&Gab);

    global_dpd_->buf4_init(&Gbb, PSIF_DCT_DENSITY, 0, ID("[o>o]-"), ID("[o>o]-"), ID("[o>o]-"), ID("[o>o]-"), 0,
                           varname("<oo|oo>"));
    global_dpd_->buf4_scm(&Gbb, 1.0 / 8.0);
    global_dpd_->buf4_close(&Gbb);

    if (!cumulant_only) {
        compute_unrelaxed_separable_density_OOOO();
    }

    psio_->close(PSIF_DCT_DENSITY, 1);
}

void DCTSolver::compute_unrelaxed_separable_density_OOOO() {
    dpdbuf4 Gaa, Gab, Gbb;

    global_dpd_->buf4_init(&Gaa, PSIF_DCT_DENSITY, 0, ID("[O,O]"), ID("[O,O]"), ID("[O>O]-"), ID("[O>O]-"), 0,
                           "Gamma <OO|OO>");
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&Gaa, h);
        global_dpd_->buf4_mat_irrep_rd(&Gaa, h);

#pragma omp parallel for
        for (size_t ij = 0; ij < Gaa.params->rowtot[h]; ++ij) {
            size_t i = Gaa.params->roworb[h][ij][0];
            int Gi = Gaa.params->psym[i];
            i -= Gaa.params->poff[Gi];
            size_t j = Gaa.params->roworb[h][ij][1];
            int Gj = Gaa.params->qsym[j];
            j -= Gaa.params->qoff[Gj];
            for (size_t kl = 0; kl < Gaa.params->coltot[h]; ++kl) {
                double tpdm = 0.0;
                size_t k = Gaa.params->colorb[h][kl][0];
                int Gk = Gaa.params->rsym[k];
                k -= Gaa.params->roff[Gk];
                size_t l = Gaa.params->colorb[h][kl][1];
                int Gl = Gaa.params->ssym[l];
                l -= Gaa.params->soff[Gl];

                if (Gi == Gk && Gj == Gl) tpdm += 0.25 * kappa_mo_a_->get(Gi, i, k) * kappa_mo_a_->get(Gj, j, l);
                if (Gi == Gl && Gj == Gk) tpdm -= 0.25 * kappa_mo_a_->get(Gi, i, l) * kappa_mo_a_->get(Gj, j, k);

                if (Gi == Gk && Gj == Gl) tpdm += 0.25 * kappa_mo_a_->get(Gi, i, k) * aocc_tau_.get(Gj, j, l);
                if (Gi == Gl && Gj == Gk) tpdm -= 0.25 * kappa_mo_a_->get(Gi, i, l) * aocc_tau_.get(Gj, j, k);
                if (Gj == Gk && Gi == Gl) tpdm -= 0.25 * kappa_mo_a_->get(Gj, j, k) * aocc_tau_.get(Gi, i, l);
                if (Gj == Gl && Gi == Gk) tpdm += 0.25 * kappa_mo_a_->get(Gj, j, l) * aocc_tau_.get(Gi, i, k);

                if (Gi == Gk && Gj == Gl) tpdm += 0.25 * aocc_tau_.get(Gi, i, k) * aocc_tau_.get(Gj, j, l);
                if (Gi == Gl && Gj == Gk) tpdm -= 0.25 * aocc_tau_.get(Gi, i, l) * aocc_tau_.get(Gj, j, k);

                Gaa.matrix[h][ij][kl] += tpdm;
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&Gaa, h);
        global_dpd_->buf4_mat_irrep_close(&Gaa, h);
    }

    global_dpd_->buf4_close(&Gaa);

    global_dpd_->buf4_init(&Gab, PSIF_DCT_DENSITY, 0, ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), 0,
                           "Gamma <Oo|Oo>");
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&Gab, h);
        global_dpd_->buf4_mat_irrep_rd(&Gab, h);

#pragma omp parallel for
        for (size_t ij = 0; ij < Gab.params->rowtot[h]; ++ij) {
            size_t i = Gab.params->roworb[h][ij][0];
            int Gi = Gab.params->psym[i];
            i -= Gab.params->poff[Gi];
            size_t j = Gab.params->roworb[h][ij][1];
            int Gj = Gab.params->qsym[j];
            j -= Gab.params->qoff[Gj];
            for (size_t kl = 0; kl < Gab.params->coltot[h]; ++kl) {
                double tpdm = 0.0;
                size_t k = Gab.params->colorb[h][kl][0];
                int Gk = Gab.params->rsym[k];
                k -= Gab.params->roff[Gk];
                size_t l = Gab.params->colorb[h][kl][1];
                int Gl = Gab.params->ssym[l];
                l -= Gab.params->soff[Gl];
                if (Gi == Gk && Gj == Gl) tpdm += 0.25 * kappa_mo_a_->get(Gi, i, k) * kappa_mo_b_->get(Gj, j, l);

                if (Gi == Gk && Gj == Gl) tpdm += 0.25 * kappa_mo_a_->get(Gi, i, k) * bocc_tau_.get(Gj, j, l);
                if (Gj == Gl && Gi == Gk) tpdm += 0.25 * kappa_mo_b_->get(Gj, j, l) * aocc_tau_.get(Gi, i, k);

                if (Gi == Gk && Gj == Gl) tpdm += 0.25 * aocc_tau_.get(Gi, i, k) * bocc_tau_.get(Gj, j, l);

                Gab.matrix[h][ij][kl] += tpdm;
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&Gab, h);
        global_dpd_->buf4_mat_irrep_close(&Gab, h);
    }

    global_dpd_->buf4_close(&Gab);

    global_dpd_->buf4_init(&Gbb, PSIF_DCT_DENSITY, 0, ID("[o,o]"), ID("[o,o]"), ID("[o>o]-"), ID("[o>o]-"), 0,
                           "Gamma <oo|oo>");
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&Gbb, h);
        global_dpd_->buf4_mat_irrep_rd(&Gbb, h);

#pragma omp parallel for
        for (size_t ij = 0; ij < Gbb.params->rowtot[h]; ++ij) {
            size_t i = Gbb.params->roworb[h][ij][0];
            int Gi = Gbb.params->psym[i];
            i -= Gbb.params->poff[Gi];
            size_t j = Gbb.params->roworb[h][ij][1];
            int Gj = Gbb.params->qsym[j];
            j -= Gbb.params->qoff[Gj];
            for (size_t kl = 0; kl < Gbb.params->coltot[h]; ++kl) {
                double tpdm = 0.0;
                size_t k = Gbb.params->colorb[h][kl][0];
                int Gk = Gbb.params->rsym[k];
                k -= Gbb.params->roff[Gk];
                size_t l = Gbb.params->colorb[h][kl][1];
                int Gl = Gbb.params->ssym[l];
                l -= Gbb.params->soff[Gl];
                if (Gi == Gk && Gj == Gl) tpdm += 0.25 * kappa_mo_b_->get(Gi, i, k) * kappa_mo_b_->get(Gj, j, l);
                if (Gi == Gl && Gj == Gk) tpdm -= 0.25 * kappa_mo_b_->get(Gi, i, l) * kappa_mo_b_->get(Gj, j, k);

                if (Gi == Gk && Gj == Gl) tpdm += 0.25 * kappa_mo_b_->get(Gi, i, k) * bocc_tau_.get(Gj, j, l);
                if (Gi == Gl && Gj == Gk) tpdm -= 0.25 * kappa_mo_b_->get(Gi, i, l) * bocc_tau_.get(Gj, j, k);
                if (Gj == Gk && Gi == Gl) tpdm -= 0.25 * kappa_mo_b_->get(Gj, j, k) * bocc_tau_.get(Gi, i, l);
                if (Gj == Gl && Gi == Gk) tpdm += 0.25 * kappa_mo_b_->get(Gj, j, l) * bocc_tau_.get(Gi, i, k);

                if (Gi == Gk && Gj == Gl) tpdm += 0.25 * bocc_tau_.get(Gi, i, k) * bocc_tau_.get(Gj, j, l);
                if (Gi == Gl && Gj == Gk) tpdm -= 0.25 * bocc_tau_.get(Gi, i, l) * bocc_tau_.get(Gj, j, k);

                Gbb.matrix[h][ij][kl] += tpdm;
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&Gbb, h);
        global_dpd_->buf4_mat_irrep_close(&Gbb, h);
    }

    global_dpd_->buf4_close(&Gbb);
}

void DCTSolver::compute_unrelaxed_density_OOVV(bool cumulant_only) {
    psio_->open(PSIF_DCT_DENSITY, PSIO_OPEN_OLD);

    dpdbuf4 Laa, Lab, Lbb, Gaa, Gab, Gbb;
    dpdbuf4 L, G, T, II, Taa, Tab, Tbb, Kaa, Kab, Kbb;
    dpdfile2 T_OO, T_oo, T_VV, T_vv;

    const std::string density_variable = cumulant_only ? "Lambda " : "Gamma ";
    auto varname = [&density_variable](const std::string& x) { return (density_variable + x); };

    /*
     * The OOVV and VVOO blocks
     */

    // First-order density contribution

    // OOVV
    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->buf4_copy(&Laa, PSIF_DCT_DENSITY, varname("<OO|VV>"));
    global_dpd_->buf4_close(&Laa);

    global_dpd_->buf4_init(&Gaa, PSIF_DCT_DENSITY, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           varname("<OO|VV>"));
    global_dpd_->buf4_scm(&Gaa, 0.5);
    global_dpd_->buf4_close(&Gaa);

    // OoVv
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->buf4_copy(&Lab, PSIF_DCT_DENSITY, varname("<Oo|Vv>"));
    global_dpd_->buf4_close(&Lab);

    global_dpd_->buf4_init(&Gab, PSIF_DCT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           varname("<Oo|Vv>"));
    global_dpd_->buf4_scm(&Gab, 0.5);
    global_dpd_->buf4_close(&Gab);

    // oovv
    global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->buf4_copy(&Lbb, PSIF_DCT_DENSITY, varname("<oo|vv>"));
    global_dpd_->buf4_close(&Lbb);

    global_dpd_->buf4_init(&Gbb, PSIF_DCT_DENSITY, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           varname("<oo|vv>"));
    global_dpd_->buf4_scm(&Gbb, 0.5);
    global_dpd_->buf4_close(&Gbb);

    // Add third-order terms for the oovv density
    if (options_.get_str("DCT_FUNCTIONAL") == "ODC-13") {
        global_dpd_->file2_init(&T_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "T <O|O>");
        global_dpd_->file2_init(&T_oo, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "T <o|o>");
        global_dpd_->file2_init(&T_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "T <V|V>");
        global_dpd_->file2_init(&T_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "T <v|v>");

        /*
         * Gamma_ijab = 1/6 P_(ij) lambda_abik T_kj
         */

        // OOVV

        // G_IJAB = 1/6 lambda_IKAB * T_KJ
        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "Temp <OO|VV>");
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                               "Amplitude <OO|VV>");
        global_dpd_->contract424(&L, &T_OO, &T, 1, 0, 1, 1.0 / 6.0, 0.0);
        global_dpd_->buf4_close(&L);
        global_dpd_->buf4_close(&T);

        // Temp_IJAB -> Temp_JIAB
        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "Temp <OO|VV>");
        global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
        global_dpd_->buf4_close(&T);

        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "Temp <OO|VV>");
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                               varname("<OO|VV>"));
        dpd_buf4_add(&G, &T, 1.0);
        global_dpd_->buf4_close(&G);
        global_dpd_->buf4_close(&T);

        // G_IJAB -= 1/6 lambda_JKAB * T_KI
        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "P(Temp) <OO|VV>");
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                               varname("<OO|VV>"));
        dpd_buf4_add(&G, &T, -1.0);
        global_dpd_->buf4_close(&G);
        global_dpd_->buf4_close(&T);

        // OoVv

        // G_IjAb += 1/6 lambda_IkAb * T_kj
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               varname("<Oo|Vv>"));
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "Amplitude <Oo|Vv>");
        global_dpd_->contract424(&L, &T_oo, &G, 1, 0, 1, 1.0 / 6.0, 1.0);
        global_dpd_->buf4_close(&L);
        global_dpd_->buf4_close(&G);

        // G_IjAb += 1/6 T_IK * lambda_KjAb
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               varname("<Oo|Vv>"));
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "Amplitude <Oo|Vv>");
        global_dpd_->contract244(&T_OO, &L, &G, 1, 0, 0, 1.0 / 6.0, 1.0);
        global_dpd_->buf4_close(&L);
        global_dpd_->buf4_close(&G);

        // oovv

        // G_ijab = 1/6 lambda_ikab * T_kj
        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "Temp <oo|vv>");
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                               "Amplitude <oo|vv>");
        global_dpd_->contract424(&L, &T_oo, &T, 1, 0, 1, 1.0 / 6.0, 0.0);
        global_dpd_->buf4_close(&L);
        global_dpd_->buf4_close(&T);

        // Temp_ijab -> Temp_jiab
        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "Temp <oo|vv>");
        global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
        global_dpd_->buf4_close(&T);

        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "Temp <oo|vv>");
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                               varname("<oo|vv>"));
        dpd_buf4_add(&G, &T, 1.0);
        global_dpd_->buf4_close(&G);
        global_dpd_->buf4_close(&T);

        // G_ijab -= 1/6 lambda_jkab * T_ki
        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "P(Temp) <oo|vv>");
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                               varname("<oo|vv>"));
        dpd_buf4_add(&G, &T, -1.0);
        global_dpd_->buf4_close(&G);
        global_dpd_->buf4_close(&T);

        /*
         * Gamma_ijab -= 1/6 P_(ab) lambda_ijac T_cb
         */

        // OOVV

        // G_IJAB -= 1/6 lambda_IJAC * T_CB
        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "Temp <OO|VV>");
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                               "Amplitude <OO|VV>");
        global_dpd_->contract424(&L, &T_VV, &T, 3, 0, 0, -1.0 / 6.0, 0.0);
        global_dpd_->buf4_close(&L);
        global_dpd_->buf4_close(&T);

        // Temp_IJAB -> Temp_IJBA
        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "Temp <OO|VV>");
        global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
        global_dpd_->buf4_close(&T);

        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "Temp <OO|VV>");
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                               varname("<OO|VV>"));
        dpd_buf4_add(&G, &T, 1.0);
        global_dpd_->buf4_close(&G);
        global_dpd_->buf4_close(&T);

        // G_IJAB += 1/6 lambda_IJBC * T_CA
        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "P(Temp) <OO|VV>");
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                               varname("<OO|VV>"));
        dpd_buf4_add(&G, &T, -1.0);
        global_dpd_->buf4_close(&G);
        global_dpd_->buf4_close(&T);

        // OoVv

        // G_IjAb -= 1/6 lambda_IjAc * T_cb
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               varname("<Oo|Vv>"));
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "Amplitude <Oo|Vv>");
        global_dpd_->contract424(&L, &T_vv, &G, 3, 0, 0, -1.0 / 6.0, 1.0);
        global_dpd_->buf4_close(&L);
        global_dpd_->buf4_close(&G);

        // G_IjAb -= 1/6 lambda_IjCb * T_CA
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               varname("<Oo|Vv>"));
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "Amplitude <Oo|Vv>");
        global_dpd_->contract244(&T_VV, &L, &G, 1, 2, 1, -1.0 / 6.0, 1.0);
        global_dpd_->buf4_close(&L);
        global_dpd_->buf4_close(&G);

        // oovv

        // G_ijab -= 1/6 lambda_ijac * T_cb
        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "Temp <oo|vv>");
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                               "Amplitude <oo|vv>");
        global_dpd_->contract424(&L, &T_vv, &T, 3, 0, 0, -1.0 / 6.0, 0.0);
        global_dpd_->buf4_close(&L);
        global_dpd_->buf4_close(&T);

        // Temp_ijab -> Temp_ijba
        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "Temp <oo|vv>");
        global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
        global_dpd_->buf4_close(&T);

        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "Temp <oo|vv>");
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                               varname("<oo|vv>"));
        dpd_buf4_add(&G, &T, 1.0);
        global_dpd_->buf4_close(&G);
        global_dpd_->buf4_close(&T);

        // G_ijab += 1/6 lambda_ijbc * T_ca
        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "P(Temp) <oo|vv>");
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                               varname("<oo|vv>"));
        dpd_buf4_add(&G, &T, -1.0);
        global_dpd_->buf4_close(&G);
        global_dpd_->buf4_close(&T);

        global_dpd_->file2_close(&T_OO);
        global_dpd_->file2_close(&T_oo);
        global_dpd_->file2_close(&T_VV);
        global_dpd_->file2_close(&T_vv);

        /*
         * Gamma_ijab += 1/12 lambda_abkl I_klij
         */

        // Gamma_IJAB += 1/12 * lambda_ABKL * I_KLIJ
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                               varname("<OO|VV>"));
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                               "Amplitude <OO|VV>");
        global_dpd_->buf4_init(&II, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[O>O]-"), ID("[O>O]-"), ID("[O>O]-"), 0,
                               "I <OO|OO>");
        global_dpd_->contract444(&II, &L, &G, 0, 1, 1.0 / 6.0, 1.0);
        global_dpd_->buf4_close(&II);
        global_dpd_->buf4_close(&L);
        global_dpd_->buf4_close(&G);

        // Gamma_IjAb += 1/6 * lambda_AbKl * I_KlIj
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               varname("<Oo|Vv>"));
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "Amplitude <Oo|Vv>");
        global_dpd_->buf4_init(&II, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), 0,
                               "I <Oo|Oo>");
        global_dpd_->contract444(&II, &L, &G, 0, 1, 1.0 / 6.0, 1.0);
        global_dpd_->buf4_close(&II);
        global_dpd_->buf4_close(&L);
        global_dpd_->buf4_close(&G);

        // Gamma_ijab += 1/12 * lambda_abkl * I_klij
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                               varname("<oo|vv>"));
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                               "Amplitude <oo|vv>");
        global_dpd_->buf4_init(&II, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[o>o]-"), ID("[o>o]-"), ID("[o>o]-"), 0,
                               "I <oo|oo>");
        global_dpd_->contract444(&II, &L, &G, 0, 1, 1.0 / 6.0, 1.0);
        global_dpd_->buf4_close(&II);
        global_dpd_->buf4_close(&L);
        global_dpd_->buf4_close(&G);

        /*
         * Gamma_ijab += 1/6 P_(ij) P_(ab) lambda_acik K_kbjc
         */

        // OOVV
        global_dpd_->buf4_init(&Taa, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "Temp (OV|OV)");

        global_dpd_->buf4_init(&Kaa, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "K (OV|OV)");
        global_dpd_->buf4_init(&Kab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                               "K (OV|ov)");

        // T_IAJB = 1/6 Amplitude_(IA|KC) K_(KC|JB)
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "Amplitude (OV|OV)");
        global_dpd_->contract444(&L, &Kaa, &Taa, 0, 0, 1.0 / 6.0, 0.0);
        global_dpd_->buf4_close(&L);

        // T_IAJB += 1/6 Amplitude_(IA|kc) K_(JB|kc)
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                               "Amplitude (OV|ov)");
        global_dpd_->contract444(&L, &Kab, &Taa, 0, 0, 1.0 / 6.0, 1.0);
        global_dpd_->buf4_close(&L);

        global_dpd_->buf4_close(&Kaa);
        global_dpd_->buf4_close(&Kab);

        // T_IAJB -> T_IJAB
        global_dpd_->buf4_sort(&Taa, PSIF_DCT_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "Temp <OO|VV>");
        global_dpd_->buf4_close(&Taa);

        // Gamma_IJAB += T_IJAB
        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 0,
                               "Temp <OO|VV>");
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                               varname("<OO|VV>"));
        dpd_buf4_add(&G, &T, 1.0);
        global_dpd_->buf4_close(&G);
        global_dpd_->buf4_close(&T);

        global_dpd_->buf4_init(&Taa, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "Temp <OO|VV>");

        // T_IJAB -> T_JIAB
        global_dpd_->buf4_sort(&Taa, PSIF_DCT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");

        // Gamma_IJAB -= T_JIAB
        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 0,
                               "P(Temp) <OO|VV>");
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                               varname("<OO|VV>"));
        dpd_buf4_add(&G, &T, -1.0);
        global_dpd_->buf4_close(&G);
        global_dpd_->buf4_close(&T);

        // T_IJAB -> T_IJBA
        global_dpd_->buf4_sort(&Taa, PSIF_DCT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");

        // Gamma_IJAB -= T_IJBA
        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 0,
                               "P(Temp) <OO|VV>");
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                               varname("<OO|VV>"));
        dpd_buf4_add(&G, &T, -1.0);
        global_dpd_->buf4_close(&G);
        global_dpd_->buf4_close(&T);

        // T_IJAB -> T_JIBA
        global_dpd_->buf4_sort(&Taa, PSIF_DCT_DPD, qpsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");

        // Gamma_IJAB += T_JIBA
        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O,O]"), ID("[V,V]"), 0,
                               "P(Temp) <OO|VV>");
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                               varname("<OO|VV>"));
        dpd_buf4_add(&G, &T, 1.0);
        global_dpd_->buf4_close(&G);
        global_dpd_->buf4_close(&T);

        global_dpd_->buf4_close(&Taa);

        // OoVv
        global_dpd_->buf4_init(&Kaa, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "K (OV|OV)");
        global_dpd_->buf4_init(&Kab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                               "K (OV|ov)");
        global_dpd_->buf4_init(&Kbb, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                               "K (ov|ov)");

        global_dpd_->buf4_init(&Tab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                               "Temp (OV|ov)");

        // T_IAjb = 1/6 Amplitude_(IA|kc) K_(kc|jb)
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                               "Amplitude (OV|ov)");
        global_dpd_->contract444(&L, &Kbb, &Tab, 0, 1, 1.0 / 6.0, 0.0);
        global_dpd_->buf4_close(&L);

        // T_IAjb += 1/6 Amplitude_(IA|KC) K_(KC|jb)
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "Amplitude (OV|OV)");
        global_dpd_->contract444(&L, &Kab, &Tab, 0, 1, 1.0 / 6.0, 1.0);
        global_dpd_->buf4_close(&L);

        // T_IAjb += 1/6 K_(IA|KC) Amplitude_(KC|jb)
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                               "Amplitude (OV|ov)");
        global_dpd_->contract444(&Kaa, &L, &Tab, 0, 1, 1.0 / 6.0, 1.0);
        global_dpd_->buf4_close(&L);

        // T_IAjb += 1/6 K_(IA|kc) Amplitude_(jb|kc)
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                               "Amplitude (ov|ov)");
        global_dpd_->contract444(&Kab, &L, &Tab, 0, 0, 1.0 / 6.0, 1.0);
        global_dpd_->buf4_close(&L);

        global_dpd_->buf4_close(&Tab);

        // T_IAjb -> T_IjAb
        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                               "Temp (OV|ov)");
        global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, prqs, ID("[O,o]"), ID("[V,v]"), "Temp <Oo|Vv>");
        global_dpd_->buf4_close(&T);

        // Gamma_IjAb += T_IjAb
        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "Temp <Oo|Vv>");
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               varname("<Oo|Vv>"));
        dpd_buf4_add(&G, &T, 1.0);
        global_dpd_->buf4_close(&G);
        global_dpd_->buf4_close(&T);

        global_dpd_->buf4_close(&Kaa);
        global_dpd_->buf4_close(&Kab);
        global_dpd_->buf4_close(&Kbb);

        // T_IbAj = 1/6 K_(Ib|Kc) Amplitude_(Kc|Aj)
        global_dpd_->buf4_init(&Tab, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[V,o]"), ID("[O,v]"), ID("[V,o]"), 0,
                               "Temp (Ov|Vo)");
        global_dpd_->buf4_init(&Kab, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), 0,
                               "K <Ov|Ov>");
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[V,o]"), ID("[O,v]"), ID("[V,o]"), 0,
                               "Amplitude (Ov|Vo)");
        global_dpd_->contract444(&Kab, &L, &Tab, 0, 1, 1.0 / 6.0, 0.0);
        global_dpd_->buf4_close(&L);
        global_dpd_->buf4_close(&Kab);
        global_dpd_->buf4_close(&Tab);

        // T_IbAj += 1/6 Amplitude_(Ib|Ck) K_(Ck|Aj)
        global_dpd_->buf4_init(&Tab, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[V,o]"), ID("[O,v]"), ID("[V,o]"), 0,
                               "Temp (Ov|Vo)");
        global_dpd_->buf4_init(&Kab, PSIF_DCT_DPD, 0, ID("[V,o]"), ID("[V,o]"), ID("[V,o]"), ID("[V,o]"), 0,
                               "K <Vo|Vo>");
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[V,o]"), ID("[O,v]"), ID("[V,o]"), 0,
                               "Amplitude (Ov|Vo)");
        global_dpd_->contract444(&L, &Kab, &Tab, 0, 1, 1.0 / 6.0, 1.0);
        global_dpd_->buf4_close(&L);
        global_dpd_->buf4_close(&Kab);
        global_dpd_->buf4_close(&Tab);

        // T_IbAj -> T_IjAb
        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[V,o]"), ID("[O,v]"), ID("[V,o]"), 0,
                               "Temp (Ov|Vo)");
        global_dpd_->buf4_sort(&T, PSIF_DCT_DPD, psrq, ID("[O,o]"), ID("[V,v]"), "Temp <Oo|Vv>");
        global_dpd_->buf4_close(&T);

        // Gamma_IjAb += T_IjAb
        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "Temp <Oo|Vv>");
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               varname("<Oo|Vv>"));
        dpd_buf4_add(&G, &T, 1.0);
        global_dpd_->buf4_close(&G);
        global_dpd_->buf4_close(&T);

        // oovv
        global_dpd_->buf4_init(&Tbb, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                               "Temp (ov|ov)");

        global_dpd_->buf4_init(&Kbb, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                               "K (ov|ov)");
        global_dpd_->buf4_init(&Kab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                               "K (OV|ov)");

        // T_iajb = 1/6 Amplitude_(ia|kc) K_(kc|jb)
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                               "Amplitude (ov|ov)");
        global_dpd_->contract444(&L, &Kbb, &Tbb, 0, 0, 1.0 / 6.0, 0.0);
        global_dpd_->buf4_close(&L);

        // T_iajb += 1/6 Amplitude_(KC|ia) K_(KC|jb)
        global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                               "Amplitude (OV|ov)");
        global_dpd_->contract444(&L, &Kab, &Tbb, 1, 1, 1.0 / 6.0, 1.0);
        global_dpd_->buf4_close(&L);

        global_dpd_->buf4_close(&Kab);
        global_dpd_->buf4_close(&Kbb);

        // T_iajb -> T_ijab
        global_dpd_->buf4_sort(&Tbb, PSIF_DCT_DPD, prqs, ID("[o,o]"), ID("[v,v]"), "Temp <oo|vv>");
        global_dpd_->buf4_close(&Tbb);

        // Gamma_ijab += T_ijab
        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 0,
                               "Temp <oo|vv>");
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                               varname("<oo|vv>"));
        dpd_buf4_add(&G, &T, 1.0);
        global_dpd_->buf4_close(&G);
        global_dpd_->buf4_close(&T);

        global_dpd_->buf4_init(&Tbb, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "Temp <oo|vv>");

        // T_ijab -> T_jiab
        global_dpd_->buf4_sort(&Tbb, PSIF_DCT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");

        // Gamma_ijab -= T_jiab
        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 0,
                               "P(Temp) <oo|vv>");
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                               varname("<oo|vv>"));
        dpd_buf4_add(&G, &T, -1.0);
        global_dpd_->buf4_close(&G);
        global_dpd_->buf4_close(&T);

        // T_ijab -> T_ijba
        global_dpd_->buf4_sort(&Tbb, PSIF_DCT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");

        // Gamma_ijab -= T_ijba
        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 0,
                               "P(Temp) <oo|vv>");
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                               varname("<oo|vv>"));
        dpd_buf4_add(&G, &T, -1.0);
        global_dpd_->buf4_close(&G);
        global_dpd_->buf4_close(&T);

        // T_ijab -> T_jiba
        global_dpd_->buf4_sort(&Tbb, PSIF_DCT_DPD, qpsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");

        // Gamma_ijab += T_jiba
        global_dpd_->buf4_init(&T, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o,o]"), ID("[v,v]"), 0,
                               "P(Temp) <oo|vv>");
        global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                               varname("<oo|vv>"));
        dpd_buf4_add(&G, &T, 1.0);
        global_dpd_->buf4_close(&G);
        global_dpd_->buf4_close(&T);

        global_dpd_->buf4_close(&Tbb);
    }

    // Sort OOVV density to VVOO
    global_dpd_->buf4_init(&Gaa, PSIF_DCT_DENSITY, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           varname("<OO|VV>"));
    global_dpd_->buf4_sort(&Gaa, PSIF_DCT_DENSITY, rspq, ID("[V>V]-"), ID("[O>O]-"), varname("<VV|OO>"));
    global_dpd_->buf4_close(&Gaa);

    // Sort OoVv density to VvOo
    global_dpd_->buf4_init(&Gab, PSIF_DCT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           varname("<Oo|Vv>"));
    global_dpd_->buf4_sort(&Gab, PSIF_DCT_DENSITY, rspq, ID("[V,v]"), ID("[O,o]"), varname("<Vv|Oo>"));
    global_dpd_->buf4_close(&Gab);

    // Sort oovv density to vvoo
    global_dpd_->buf4_init(&Gbb, PSIF_DCT_DENSITY, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           varname("<oo|vv>"));
    global_dpd_->buf4_sort(&Gbb, PSIF_DCT_DENSITY, rspq, ID("[v>v]-"), ID("[o>o]-"), varname("<vv|oo>"));
    global_dpd_->buf4_close(&Gbb);

    psio_->close(PSIF_DCT_DENSITY, 1);
}

void DCTSolver::compute_unrelaxed_density_OVOV(bool cumulant_only) {
    /*
     * The OVOV block
     */

    dpdbuf4 Kaa, Kab, Kba, Kbb, Gaa, Gab, Gba, Gbb;

    if (options_.get_str("DCT_FUNCTIONAL") != "ODC-13") {
        compute_K_intermediate();
    }

    psio_->open(PSIF_DCT_DENSITY, PSIO_OPEN_OLD);

    const std::string density_variable = cumulant_only ? "Lambda " : "Gamma ";
    auto varname = [&density_variable](const std::string& x) { return (density_variable + x); };

    // There are five unique spin cases: Г<IAJB>, Г<iajb>, Г<IaJb>, Г<iAjB>, Г<IajB>

    // Г<IAJB> spin case

    // Gamma_IAJB = -1.0 * K_IAJB
    global_dpd_->buf4_init(&Kaa, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "K <OV|OV>");
    global_dpd_->buf4_copy(&Kaa, PSIF_DCT_DENSITY, varname("<OV|OV>"));
    global_dpd_->buf4_close(&Kaa);

    global_dpd_->buf4_init(&Gaa, PSIF_DCT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           varname("<OV|OV>"));
    global_dpd_->buf4_scm(&Gaa, -1.0);
    global_dpd_->buf4_close(&Gaa);

    // Г<IaJb> and Г<iAjB> spin cases:
    // Gamma_IaJb = -1.0 * K_IaJb
    // Gamma_iAjB = -1.0 * K_iAjB
    global_dpd_->buf4_init(&Kab, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), 0, "K <Ov|Ov>");
    global_dpd_->buf4_copy(&Kab, PSIF_DCT_DENSITY, varname("<Ov|Ov>"));
    global_dpd_->buf4_close(&Kab);

    global_dpd_->buf4_init(&Kba, PSIF_DCT_DPD, 0, ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), 0, "K <oV|oV>");
    global_dpd_->buf4_copy(&Kba, PSIF_DCT_DENSITY, varname("<oV|oV>"));
    global_dpd_->buf4_close(&Kba);

    global_dpd_->buf4_init(&Gab, PSIF_DCT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), 0,
                           varname("<Ov|Ov>"));
    global_dpd_->buf4_scm(&Gab, -1.0);
    global_dpd_->buf4_close(&Gab);

    global_dpd_->buf4_init(&Gba, PSIF_DCT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), 0,
                           varname("<oV|oV>"));
    global_dpd_->buf4_scm(&Gba, -1.0);
    global_dpd_->buf4_close(&Gba);

    // Г<IajB> spin case:
    // Gamma_IajB = -1.0 * K_IajB
    global_dpd_->buf4_init(&Kab, PSIF_DCT_DPD, 0, ID("[O,v]"), ID("[o,V]"), ID("[O,v]"), ID("[o,V]"), 0, "K <Ov|oV>");
    global_dpd_->buf4_copy(&Kab, PSIF_DCT_DENSITY, varname("<Ov|oV>"));
    global_dpd_->buf4_close(&Kab);

    global_dpd_->buf4_init(&Kab, PSIF_DCT_DPD, 0, ID("[o,V]"), ID("[O,v]"), ID("[o,V]"), ID("[O,v]"), 0, "K <oV|Ov>");
    global_dpd_->buf4_copy(&Kab, PSIF_DCT_DENSITY, varname("<oV|Ov>"));
    global_dpd_->buf4_close(&Kab);

    global_dpd_->buf4_init(&Gab, PSIF_DCT_DENSITY, 0, ID("[O,v]"), ID("[o,V]"), ID("[O,v]"), ID("[o,V]"), 0,
                           varname("<Ov|oV>"));
    global_dpd_->buf4_scm(&Gab, -1.0);
    global_dpd_->buf4_close(&Gab);

    global_dpd_->buf4_init(&Gab, PSIF_DCT_DENSITY, 0, ID("[o,V]"), ID("[O,v]"), ID("[o,V]"), ID("[O,v]"), 0,
                           varname("<oV|Ov>"));
    global_dpd_->buf4_scm(&Gab, -1.0);
    global_dpd_->buf4_close(&Gab);

    // Г<iajb> spin case:
    // Gamma_iajb = -1.0 * K_iajb
    global_dpd_->buf4_init(&Kbb, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0, "K <ov|ov>");
    global_dpd_->buf4_copy(&Kbb, PSIF_DCT_DENSITY, varname("<ov|ov>"));
    global_dpd_->buf4_close(&Kbb);

    global_dpd_->buf4_init(&Gbb, PSIF_DCT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                           varname("<ov|ov>"));
    global_dpd_->buf4_scm(&Gbb, -1.0);
    global_dpd_->buf4_close(&Gbb);

    if (!cumulant_only) {
        compute_unrelaxed_separable_density_OVOV();
    };

    psio_->close(PSIF_DCT_DENSITY, 1);
}

void DCTSolver::compute_unrelaxed_separable_density_OVOV() {
    dpdbuf4 Gaa, Gab, Gba, Gbb;

    global_dpd_->buf4_init(&Gaa, PSIF_DCT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Gamma <OV|OV>");
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&Gaa, h);
        global_dpd_->buf4_mat_irrep_rd(&Gaa, h);

#pragma omp parallel for
        for (size_t ia = 0; ia < Gaa.params->rowtot[h]; ++ia) {
            size_t i = Gaa.params->roworb[h][ia][0];
            int Gi = Gaa.params->psym[i];
            i -= Gaa.params->poff[Gi];
            size_t a = Gaa.params->roworb[h][ia][1];
            int Ga = Gaa.params->qsym[a];
            a -= Gaa.params->qoff[Ga];
            for (size_t jb = 0; jb < Gaa.params->coltot[h]; ++jb) {
                size_t j = Gaa.params->colorb[h][jb][0];
                int Gj = Gaa.params->rsym[j];
                j -= Gaa.params->roff[Gj];
                size_t b = Gaa.params->colorb[h][jb][1];
                int Gb = Gaa.params->ssym[b];
                b -= Gaa.params->soff[Gb];
                if (Gi == Gj && Ga == Gb) {
                    Gaa.matrix[h][ia][jb] +=
                        (kappa_mo_a_->get(Gi, i, j) + aocc_tau_.get(Gi, i, j)) * avir_tau_.get(Ga, a, b);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&Gaa, h);
        global_dpd_->buf4_mat_irrep_close(&Gaa, h);
    }

    global_dpd_->buf4_close(&Gaa);

    global_dpd_->buf4_init(&Gab, PSIF_DCT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), 0,
                           "Gamma <Ov|Ov>");
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&Gab, h);
        global_dpd_->buf4_mat_irrep_rd(&Gab, h);

#pragma omp parallel for
        for (size_t ia = 0; ia < Gab.params->rowtot[h]; ++ia) {
            size_t i = Gab.params->roworb[h][ia][0];
            int Gi = Gab.params->psym[i];
            i -= Gab.params->poff[Gi];
            size_t a = Gab.params->roworb[h][ia][1];
            int Ga = Gab.params->qsym[a];
            a -= Gab.params->qoff[Ga];
            for (size_t jb = 0; jb < Gab.params->coltot[h]; ++jb) {
                size_t j = Gab.params->colorb[h][jb][0];
                int Gj = Gab.params->rsym[j];
                j -= Gab.params->roff[Gj];
                size_t b = Gab.params->colorb[h][jb][1];
                int Gb = Gab.params->ssym[b];
                b -= Gab.params->soff[Gb];
                if (Gi == Gj && Ga == Gb) {
                    Gab.matrix[h][ia][jb] +=
                        (kappa_mo_a_->get(Gi, i, j) + aocc_tau_.get(Gi, i, j)) * bvir_tau_.get(Ga, a, b);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&Gab, h);
        global_dpd_->buf4_mat_irrep_close(&Gab, h);
    }

    global_dpd_->buf4_close(&Gab);

    global_dpd_->buf4_init(&Gba, PSIF_DCT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), 0,
                           "Gamma <oV|oV>");
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&Gba, h);
        global_dpd_->buf4_mat_irrep_rd(&Gba, h);

#pragma omp parallel for
        for (size_t ia = 0; ia < Gba.params->rowtot[h]; ++ia) {
            size_t i = Gba.params->roworb[h][ia][0];
            int Gi = Gba.params->psym[i];
            i -= Gba.params->poff[Gi];
            size_t a = Gba.params->roworb[h][ia][1];
            int Ga = Gba.params->qsym[a];
            a -= Gba.params->qoff[Ga];
            for (size_t jb = 0; jb < Gba.params->coltot[h]; ++jb) {
                size_t j = Gba.params->colorb[h][jb][0];
                int Gj = Gba.params->rsym[j];
                j -= Gba.params->roff[Gj];
                size_t b = Gba.params->colorb[h][jb][1];
                int Gb = Gba.params->ssym[b];
                b -= Gba.params->soff[Gb];
                if (Gi == Gj && Ga == Gb) {
                    Gba.matrix[h][ia][jb] +=
                        (kappa_mo_b_->get(Gi, i, j) + bocc_tau_.get(Gi, i, j)) * avir_tau_.get(Ga, a, b);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&Gba, h);
        global_dpd_->buf4_mat_irrep_close(&Gba, h);
    }

    global_dpd_->buf4_close(&Gba);

    global_dpd_->buf4_init(&Gbb, PSIF_DCT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                           "Gamma <ov|ov>");

    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&Gbb, h);
        global_dpd_->buf4_mat_irrep_rd(&Gbb, h);

#pragma omp parallel for
        for (size_t ia = 0; ia < Gbb.params->rowtot[h]; ++ia) {
            size_t i = Gbb.params->roworb[h][ia][0];
            int Gi = Gbb.params->psym[i];
            i -= Gbb.params->poff[Gi];
            size_t a = Gbb.params->roworb[h][ia][1];
            int Ga = Gbb.params->qsym[a];
            a -= Gbb.params->qoff[Ga];
            for (size_t jb = 0; jb < Gbb.params->coltot[h]; ++jb) {
                size_t j = Gbb.params->colorb[h][jb][0];
                int Gj = Gbb.params->rsym[j];
                j -= Gbb.params->roff[Gj];
                size_t b = Gbb.params->colorb[h][jb][1];
                int Gb = Gbb.params->ssym[b];
                b -= Gbb.params->soff[Gb];
                if (Gi == Gj && Ga == Gb) {
                    Gbb.matrix[h][ia][jb] +=
                        (kappa_mo_b_->get(Gi, i, j) + bocc_tau_.get(Gi, i, j)) * bvir_tau_.get(Ga, a, b);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&Gbb, h);
        global_dpd_->buf4_mat_irrep_close(&Gbb, h);
    }

    global_dpd_->buf4_close(&Gbb);
}

void DCTSolver::compute_unrelaxed_density_VVVV(bool cumulant_only) {
    psio_->open(PSIF_DCT_DENSITY, PSIO_OPEN_OLD);

    dpdbuf4 LLaa, LLab, LLbb, Laa, Lab, Lbb, Gaa, Gab, Gbb;

    const std::string density_variable = cumulant_only ? "Lambda " : "Gamma ";
    auto varname = [&density_variable](const std::string& x) { return (density_variable + x); };

    /*
     * The VVVV block
     */

    // Gamma_abcd = 1/16 (Amplitude_ijab * Amplitude_ijcd + Amplitude_ijab * Amplitude_ijcd)
    global_dpd_->buf4_init(&Gaa, PSIF_DCT_DENSITY, 0, ID("[V>V]-"), ID("[V>V]-"), ID("[V>V]-"), ID("[V>V]-"), 0,
                           varname("<VV|VV>"));
    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->buf4_init(&LLaa, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->contract444(&Laa, &LLaa, &Gaa, 1, 1, 0.25, 0.0);
    global_dpd_->buf4_close(&LLaa);
    global_dpd_->buf4_close(&Gaa);
    global_dpd_->buf4_close(&Laa);

    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->buf4_init(&LLab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->buf4_init(&Gab, PSIF_DCT_DENSITY, 0, ID("[V,v]"), ID("[V,v]"), ID("[V,v]"), ID("[V,v]"), 0,
                           varname("<Vv|Vv>"));
    global_dpd_->contract444(&Lab, &LLab, &Gab, 1, 1, 0.25, 0.0);
    global_dpd_->buf4_close(&Gab);
    global_dpd_->buf4_close(&LLab);
    global_dpd_->buf4_close(&Lab);

    global_dpd_->buf4_init(&Gbb, PSIF_DCT_DENSITY, 0, ID("[v>v]-"), ID("[v>v]-"), ID("[v>v]-"), ID("[v>v]-"), 0,
                           varname("<vv|vv>"));
    global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->buf4_init(&LLbb, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->contract444(&Lbb, &LLbb, &Gbb, 1, 1, 0.25, 0.0);
    global_dpd_->buf4_close(&LLbb);
    global_dpd_->buf4_close(&Gbb);
    global_dpd_->buf4_close(&Lbb);

    if (!cumulant_only) {
        compute_unrelaxed_separable_density_VVVV();
    }

    psio_->close(PSIF_DCT_DENSITY, 1);
}

void DCTSolver::compute_unrelaxed_separable_density_VVVV() {
    dpdbuf4 Gaa, Gab, Gbb;

    global_dpd_->buf4_init(&Gaa, PSIF_DCT_DENSITY, 0, ID("[V,V]"), ID("[V,V]"), ID("[V>V]-"), ID("[V>V]-"), 0,
                           "Gamma <VV|VV>");
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&Gaa, h);
        global_dpd_->buf4_mat_irrep_rd(&Gaa, h);

#pragma omp parallel for
        for (size_t ab = 0; ab < Gaa.params->rowtot[h]; ++ab) {
            size_t a = Gaa.params->roworb[h][ab][0];
            int Ga = Gaa.params->psym[a];
            a -= Gaa.params->poff[Ga];
            size_t b = Gaa.params->roworb[h][ab][1];
            int Gb = Gaa.params->qsym[b];
            b -= Gaa.params->qoff[Gb];
            for (size_t cd = 0; cd < Gaa.params->coltot[h]; ++cd) {
                double tpdm = 0.0;
                size_t c = Gaa.params->colorb[h][cd][0];
                int Gc = Gaa.params->rsym[c];
                c -= Gaa.params->roff[Gc];
                size_t d = Gaa.params->colorb[h][cd][1];
                int Gd = Gaa.params->ssym[d];
                d -= Gaa.params->soff[Gd];
                if (Ga == Gc && Gb == Gd) tpdm += 0.25 * avir_tau_.get(Ga, a, c) * avir_tau_.get(Gb, b, d);
                if (Ga == Gd && Gb == Gc) tpdm -= 0.25 * avir_tau_.get(Ga, a, d) * avir_tau_.get(Gb, b, c);

                Gaa.matrix[h][ab][cd] += tpdm;
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&Gaa, h);
        global_dpd_->buf4_mat_irrep_close(&Gaa, h);
    }

    global_dpd_->buf4_close(&Gaa);

    global_dpd_->buf4_init(&Gab, PSIF_DCT_DENSITY, 0, ID("[V,v]"), ID("[V,v]"), ID("[V,v]"), ID("[V,v]"), 0,
                           "Gamma <Vv|Vv>");
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&Gab, h);
        global_dpd_->buf4_mat_irrep_rd(&Gab, h);

#pragma omp parallel for
        for (size_t ab = 0; ab < Gab.params->rowtot[h]; ++ab) {
            size_t a = Gab.params->roworb[h][ab][0];
            int Ga = Gab.params->psym[a];
            a -= Gab.params->poff[Ga];
            size_t b = Gab.params->roworb[h][ab][1];
            int Gb = Gab.params->qsym[b];
            b -= Gab.params->qoff[Gb];
            for (size_t cd = 0; cd < Gab.params->coltot[h]; ++cd) {
                double tpdm = 0.0;
                size_t c = Gab.params->colorb[h][cd][0];
                int Gc = Gab.params->rsym[c];
                c -= Gab.params->roff[Gc];
                size_t d = Gab.params->colorb[h][cd][1];
                int Gd = Gab.params->ssym[d];
                d -= Gab.params->soff[Gd];
                if (Ga == Gc && Gb == Gd) tpdm += 0.25 * avir_tau_.get(Ga, a, c) * bvir_tau_.get(Gb, b, d);
                Gab.matrix[h][ab][cd] += tpdm;
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&Gab, h);
        global_dpd_->buf4_mat_irrep_close(&Gab, h);
    }

    global_dpd_->buf4_close(&Gab);

    global_dpd_->buf4_init(&Gbb, PSIF_DCT_DENSITY, 0, ID("[v,v]"), ID("[v,v]"), ID("[v>v]-"), ID("[v>v]-"), 0,
                           "Gamma <vv|vv>");
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&Gbb, h);
        global_dpd_->buf4_mat_irrep_rd(&Gbb, h);

#pragma omp parallel for
        for (size_t ab = 0; ab < Gbb.params->rowtot[h]; ++ab) {
            size_t a = Gbb.params->roworb[h][ab][0];
            int Ga = Gbb.params->psym[a];
            a -= Gbb.params->poff[Ga];
            size_t b = Gbb.params->roworb[h][ab][1];
            int Gb = Gbb.params->qsym[b];
            b -= Gbb.params->qoff[Gb];
            for (size_t cd = 0; cd < Gbb.params->coltot[h]; ++cd) {
                double tpdm = 0.0;
                size_t c = Gbb.params->colorb[h][cd][0];
                int Gc = Gbb.params->rsym[c];
                c -= Gbb.params->roff[Gc];
                size_t d = Gbb.params->colorb[h][cd][1];
                int Gd = Gbb.params->ssym[d];
                d -= Gbb.params->soff[Gd];
                if (Ga == Gc && Gb == Gd) tpdm += 0.25 * bvir_tau_.get(Ga, a, c) * bvir_tau_.get(Gb, b, d);
                if (Ga == Gd && Gb == Gc) tpdm -= 0.25 * bvir_tau_.get(Ga, a, d) * bvir_tau_.get(Gb, b, c);
                Gbb.matrix[h][ab][cd] += tpdm;
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&Gbb, h);
        global_dpd_->buf4_mat_irrep_close(&Gbb, h);
    }

    global_dpd_->buf4_close(&Gbb);
}

Matrix DCTSolver::construct_oo_density(const Matrix& occtau, const Matrix& virtau, const Matrix& kappa,
                                       const Matrix& C) {
    auto opdm = Matrix("MO basis OPDM", nmopi_, nmopi_);

    auto occdim = occtau.rowspi();
    auto virdim = virtau.rowspi();

    for (int h = 0; h < nirrep_; ++h) {
        for (int i = 0; i < occdim[h]; ++i) {
            for (int j = 0; j <= i; ++j) {
                opdm.set(h, i, j, occtau.get(h, i, j) + kappa.get(h, i, j));
                if (i != j) opdm.set(h, j, i, occtau.get(h, i, j) + kappa.get(h, i, j));
            }
        }
        for (int a = 0; a < virdim[h]; ++a) {
            for (int b = 0; b <= a; ++b) {
                opdm.set(h, a + occdim[h], b + occdim[h], virtau.get(h, a, b));
                if (a != b) opdm.set(h, b + occdim[h], a + occdim[h], virtau.get(h, a, b));
            }
        }
    }

    return linalg::triplet(C, opdm, C, false, false, true);
}

void DCTSolver::compute_TPDM_trace(bool cumulant_only) {
    dpdbuf4 G;

    // If we only have the cumulant saved, we'll need to figure out the rest of the trace from the 1RDM.
    const std::string density_variable = cumulant_only ? "Lambda " : "Gamma ";
    auto varname = [&density_variable](const std::string& x) { return (density_variable + x); };

    psio_->open(PSIF_DCT_DENSITY, PSIO_OPEN_OLD);

    double tpdm_trace = 0.0;

    // OOOO density
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O>O]-"), ID("[O>O]-"), ID("[O>O]-"), ID("[O>O]-"), 0,
                           varname("<OO|OO>"));
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);

        for (size_t ij = 0; ij < G.params->rowtot[h]; ++ij) {
            tpdm_trace += 8.0 * G.matrix[h][ij][ij];
        }

        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), 0,
                           varname("<Oo|Oo>"));
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);

        for (size_t ij = 0; ij < G.params->rowtot[h]; ++ij) {
            tpdm_trace += 8.0 * G.matrix[h][ij][ij];
        }

        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[o>o]-"), ID("[o>o]-"), ID("[o>o]-"), ID("[o>o]-"), 0,
                           varname("<oo|oo>"));
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);

        for (size_t ij = 0; ij < G.params->rowtot[h]; ++ij) {
            tpdm_trace += 8.0 * G.matrix[h][ij][ij];
        }

        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    // VVVV density
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[V>V]-"), ID("[V>V]-"), ID("[V>V]-"), ID("[V>V]-"), 0,
                           varname("<VV|VV>"));
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);

        for (size_t ab = 0; ab < G.params->rowtot[h]; ++ab) {
            tpdm_trace += 8.0 * G.matrix[h][ab][ab];
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }

    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[V,v]"), ID("[V,v]"), ID("[V,v]"), ID("[V,v]"), 0,
                           varname("<Vv|Vv>"));
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);

        for (size_t ab = 0; ab < G.params->rowtot[h]; ++ab) {
            tpdm_trace += 8.0 * G.matrix[h][ab][ab];
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }

    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[v>v]-"), ID("[v>v]-"), ID("[v>v]-"), ID("[v>v]-"), 0,
                           varname("<vv|vv>"));
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);

        for (size_t ab = 0; ab < G.params->rowtot[h]; ++ab) {
            tpdm_trace += 8.0 * G.matrix[h][ab][ab];
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }

    global_dpd_->buf4_close(&G);

    // OVOV density
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           varname("<OV|OV>"));
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);

        for (size_t ia = 0; ia < G.params->rowtot[h]; ++ia) {
            tpdm_trace += 2.0 * G.matrix[h][ia][ia];
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }

    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), 0,
                           varname("<Ov|Ov>"));
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);

        for (size_t ia = 0; ia < G.params->rowtot[h]; ++ia) {
            tpdm_trace += 2.0 * G.matrix[h][ia][ia];
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }

    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), 0,
                           varname("<oV|oV>"));
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);

        for (size_t ia = 0; ia < G.params->rowtot[h]; ++ia) {
            tpdm_trace += 2.0 * G.matrix[h][ia][ia];
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }

    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                           varname("<ov|ov>"));

    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);

        for (size_t ia = 0; ia < G.params->rowtot[h]; ++ia) {
            tpdm_trace += 2.0 * G.matrix[h][ia][ia];
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }

    global_dpd_->buf4_close(&G);

    // Compute OPDM trace
    double opdm_trace = (kappa_mo_a_->trace() + kappa_mo_b_->trace());
    opdm_trace += (aocc_tau_.trace() + bocc_tau_.trace() + avir_tau_.trace() + bvir_tau_.trace());

    if (cumulant_only) {
        // Add the partial trace of the antisymmetrized product of 1RDMs.
        tpdm_trace += opdm_trace * opdm_trace;
        auto alpha = linalg::doublet(mo_gammaA_, mo_gammaA_, false, false);
        tpdm_trace -= alpha.trace();
        auto beta = linalg::doublet(mo_gammaB_, mo_gammaB_, false, false);
        tpdm_trace -= beta.trace();
    }

    // Compute deviations from N-representability
    auto N = (double)(nalpha_ + nbeta_);
    double opdm_dev = N - opdm_trace;
    double tpdm_dev = N * (N - 1.0) - tpdm_trace;

    outfile->Printf("\t OPDM trace: \t%10.6f\t\tDeviation: \t%8.4e\n", opdm_trace, opdm_dev);
    outfile->Printf("\t TPDM trace: \t%10.6f\t\tDeviation: \t%8.4e\n", tpdm_trace, tpdm_dev);

    psio_->close(PSIF_DCT_DENSITY, 1);
}

}  // namespace dct
}  // namespace psi
