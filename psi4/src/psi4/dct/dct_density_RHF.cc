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

#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"
#include "dct.h"
#include "psi4/psifiles.h"

namespace psi {
namespace dct {

void DCTSolver::compute_unrelaxed_density_OOOO_RHF(bool cumulant_only) {
    psio_->open(PSIF_DCT_DENSITY, PSIO_OPEN_OLD);

    dpdbuf4 LLab, Lab, Gab;

    const std::string density_variable = cumulant_only ? "Lambda " : "Gamma ";
    auto varname = [&density_variable](const std::string& x) { return (density_variable + x); };

    // Compute the N^6 terms for Gamma OOOO

    // The Alpha - Beta case
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Amplitude SF <OO|VV>");  // Amplitude <Oo|Vv>
    global_dpd_->buf4_init(&LLab, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Amplitude SF <OO|VV>");  // Amplitude <Oo|Vv>
    global_dpd_->buf4_init(&Gab, PSIF_DCT_DENSITY, 0, ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), 0,
                           varname("SF <OO|OO>"));  // Gamma <Oo|Oo>
    global_dpd_->contract444(&Lab, &LLab, &Gab, 0, 0, 0.25, 0.0);
    global_dpd_->buf4_close(&Gab);
    global_dpd_->buf4_close(&LLab);
    global_dpd_->buf4_close(&Lab);

    if (!cumulant_only) {
        compute_unrelaxed_separable_density_OOOO_RHF();
    };

    // Form Alpha-Alpha Gamma_OOOO case for later use
    // Gamma_IJKL = Gamma_IjKl - Gamma_JiKl
    global_dpd_->buf4_init(&Gab, PSIF_DCT_DENSITY, 0, ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), 1,
                           varname("SF <OO|OO>"));
    global_dpd_->buf4_copy(&Gab, PSIF_DCT_DENSITY, varname("<OO|OO>"));
    global_dpd_->buf4_close(&Gab);

    psio_->close(PSIF_DCT_DENSITY, 1);
}

void DCTSolver::compute_unrelaxed_separable_density_OOOO_RHF() {
    dpdbuf4 Gab;

    global_dpd_->buf4_init(&Gab, PSIF_DCT_DENSITY, 0, ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), ID("[O,O]"), 0,
                           "Gamma SF <OO|OO>");  // Gamma <Oo|Oo>
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&Gab, h);
        global_dpd_->buf4_mat_irrep_rd(&Gab, h);

#pragma omp parallel for
        for (long int ij = 0; ij < Gab.params->rowtot[h]; ++ij) {
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
                if (Gi == Gk && Gj == Gl) tpdm += 0.25 * kappa_mo_a_->get(Gi, i, k) * kappa_mo_a_->get(Gj, j, l);

                if (Gi == Gk && Gj == Gl) tpdm += 0.25 * kappa_mo_a_->get(Gi, i, k) * aocc_tau_.get(Gj, j, l);
                if (Gj == Gl && Gi == Gk) tpdm += 0.25 * kappa_mo_a_->get(Gj, j, l) * aocc_tau_.get(Gi, i, k);

                if (Gi == Gk && Gj == Gl) tpdm += 0.25 * aocc_tau_.get(Gi, i, k) * aocc_tau_.get(Gj, j, l);

                Gab.matrix[h][ij][kl] += tpdm;
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&Gab, h);
        global_dpd_->buf4_mat_irrep_close(&Gab, h);
    }
    global_dpd_->buf4_close(&Gab);
}

void DCTSolver::compute_unrelaxed_density_OOVV_RHF(bool cumulant_only) {
    psio_->open(PSIF_DCT_DENSITY, PSIO_OPEN_OLD);

    dpdbuf4 Lab, Gab;

    const std::string density_variable = cumulant_only ? "Lambda " : "Gamma ";
    auto varname = [&density_variable](const std::string& x) { return (density_variable + x); };

    /*
     * The OOVV and VVOO blocks
     */

    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Amplitude SF <OO|VV>");                         // Amplitude <Oo|Vv>
    global_dpd_->buf4_copy(&Lab, PSIF_DCT_DENSITY, varname("SF <OO|VV>"));  // Gamma <Oo|Vv>
    global_dpd_->buf4_sort(&Lab, PSIF_DCT_DENSITY, rspq, ID("[V,V]"), ID("[O,O]"),
                           varname("SF <VV|OO>"));  // Gamma <Vv|Oo>
    global_dpd_->buf4_close(&Lab);

    global_dpd_->buf4_init(&Gab, PSIF_DCT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           varname("SF <OO|VV>"));  // Gamma <Oo|Vv>
    global_dpd_->buf4_scm(&Gab, 0.5);
    global_dpd_->buf4_close(&Gab);
    global_dpd_->buf4_init(&Gab, PSIF_DCT_DENSITY, 0, ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), 0,
                           varname("SF <VV|OO>"));  // Gamma <Vv|Oo>
    global_dpd_->buf4_scm(&Gab, 0.5);
    global_dpd_->buf4_close(&Gab);

    // Form alpha-alpha case
    // Gamma_IJAB = Gamma_IjAb - Gamma_JiAb
    dpdbuf4 Gaa, T;

    global_dpd_->buf4_init(&Gaa, PSIF_DCT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 1,
                           varname("SF <OO|VV>"));
    global_dpd_->buf4_copy(&Gaa, PSIF_DCT_DENSITY, varname("<OO|VV>"));
    global_dpd_->buf4_sort(&Gaa, PSIF_DCT_DENSITY, rspq, ID("[V,V]"), ID("[O,O]"), varname("<VV|OO>"));
    global_dpd_->buf4_close(&Gaa);

    psio_->close(PSIF_DCT_DENSITY, 1);
}

void DCTSolver::compute_unrelaxed_density_OVOV_RHF(bool cumulant_only) {
    psio_->open(PSIF_DCT_DENSITY, PSIO_OPEN_OLD);

    dpdbuf4 LLaa, LLab, LLbb, Laa, Lab, Lbb, Gaa, Gab, Gba, Gbb, Tab;

    const std::string density_variable = cumulant_only ? "Lambda " : "Gamma ";
    auto varname = [&density_variable](const std::string& x) { return (density_variable + x); };

    /*
     * The OVOV block
     */

    // There are three unique spin cases for closed-shell systems: Г<IAJB>, Г<IaJb>, Г<IajB>

    // Г<IAJB> spin case

    global_dpd_->buf4_init(&Gaa, PSIF_DCT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           varname("(OV|OV)"));
    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Amplitude (OV|OV)");
    global_dpd_->buf4_init(&LLaa, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Amplitude (OV|OV)");
    global_dpd_->contract444(&Laa, &LLaa, &Gaa, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&LLaa);

    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Amplitude SF (OV|OV):(OV|ov)");  // Amplitude (OV|ov)
    global_dpd_->buf4_init(&LLab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Amplitude SF (OV|OV):(OV|ov)");  // Amplitude (OV|ov)

    global_dpd_->contract444(&Lab, &LLab, &Gaa, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&LLab);
    global_dpd_->buf4_close(&Gaa);

    // Resort Г(OV|OV) to the Г<OV|OV>
    global_dpd_->buf4_init(&Gaa, PSIF_DCT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           varname("(OV|OV)"));
    global_dpd_->buf4_sort(&Gaa, PSIF_DCT_DENSITY, psrq, ID("[O,V]"), ID("[O,V]"), varname("<OV|OV>"));
    global_dpd_->buf4_close(&Gaa);

    // Г<IaJb> spin case:

    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Amplitude SF (OV|OV):(Ov|oV)");  // Amplitude (Ov|oV)
    global_dpd_->buf4_init(&LLab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Amplitude SF (OV|OV):(Ov|oV)");  // ambda (Ov|oV)

    global_dpd_->buf4_init(&Gab, PSIF_DCT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           varname("SF <OV|OV>:<Ov|Ov>"));  // Gamma <Ov|Ov>
    global_dpd_->contract444(&Lab, &LLab, &Gab, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_close(&Gab);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&LLab);

    // Г<IajB> spin case:

    global_dpd_->buf4_init(&Tab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "-SF (OV|OV)");  // Temp (OV|ov)
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Amplitude SF (OV|OV):(OV|ov)");  // Amplitude (OV|ov)
    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Amplitude (OV|OV)");
    global_dpd_->contract444(&Laa, &Lab, &Tab, 0, 1, -1.0, 0.0);
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_init(&LLaa, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Amplitude (OV|OV)");  // Amplitude (ov|ov)
    global_dpd_->contract444(&Lab, &LLaa, &Tab, 0, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&LLaa);
    global_dpd_->buf4_close(&Tab);
    global_dpd_->buf4_init(&Tab, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "-SF (OV|OV)");  // Temp (OV|ov)
    global_dpd_->buf4_sort(&Tab, PSIF_DCT_DENSITY, psrq, ID("[O,V]"), ID("[O,V]"),
                           varname("SF <OV|OV>:<Ov|oV>"));  // Gamma <Ov|oV>
    global_dpd_->buf4_close(&Tab);
    global_dpd_->buf4_close(&Lab);

    if (!cumulant_only) {
        compute_unrelaxed_separable_density_OVOV_RHF();
    };

    psio_->close(PSIF_DCT_DENSITY, 1);
}

void DCTSolver::compute_unrelaxed_separable_density_OVOV_RHF() {
    dpdbuf4 Gaa, Gab;

    global_dpd_->buf4_init(&Gaa, PSIF_DCT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Gamma <OV|OV>");
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&Gaa, h);
        global_dpd_->buf4_mat_irrep_rd(&Gaa, h);

#pragma omp parallel for
        for (long int ia = 0; ia < Gaa.params->rowtot[h]; ++ia) {
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

    global_dpd_->buf4_init(&Gab, PSIF_DCT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "Gamma SF <OV|OV>:<Ov|Ov>");  // Gamma <Ov|Ov>
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&Gab, h);
        global_dpd_->buf4_mat_irrep_rd(&Gab, h);

#pragma omp parallel for
        for (long int ia = 0; ia < Gab.params->rowtot[h]; ++ia) {
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
                        (kappa_mo_a_->get(Gi, i, j) + aocc_tau_.get(Gi, i, j)) * avir_tau_.get(Ga, a, b);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&Gab, h);
        global_dpd_->buf4_mat_irrep_close(&Gab, h);
    }

    global_dpd_->buf4_close(&Gab);
}

void DCTSolver::compute_unrelaxed_density_VVVV_RHF(bool cumulant_only) {
    psio_->open(PSIF_DCT_DENSITY, PSIO_OPEN_OLD);

    dpdbuf4 LLaa, Laa, Gaa, Gab;

    /*
     * The VVVV block
     */

    const std::string density_variable = cumulant_only ? "Lambda " : "Gamma ";
    auto varname = [&density_variable](const std::string& x) { return (density_variable + x); };

    // Gamma_abcd = 1/16 (Amplitude_ijab * Amplitude_ijcd + Amplitude_ijab * Amplitude_ijcd)
    global_dpd_->buf4_init(&Gaa, PSIF_DCT_DENSITY, 0, ID("[V,V]"), ID("[V,V]"), ID("[V,V]"), ID("[V,V]"), 0,
                           varname("SF <VV|VV>"));
    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Amplitude SF <OO|VV>");
    global_dpd_->buf4_init(&LLaa, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Amplitude SF <OO|VV>");
    global_dpd_->contract444(&Laa, &LLaa, &Gaa, 1, 1, 0.25, 0.0);
    global_dpd_->buf4_close(&LLaa);
    global_dpd_->buf4_close(&Gaa);
    global_dpd_->buf4_close(&Laa);

    if (!cumulant_only) {
        compute_unrelaxed_separable_density_VVVV_RHF();
    };

    // Gamma <AB|CD> = Gamma<Ab|Cd> - Gamma<Ab|Dc>
    global_dpd_->buf4_init(&Gab, PSIF_DCT_DENSITY, 0, ID("[V,V]"), ID("[V,V]"), ID("[V,V]"), ID("[V,V]"), 1,
                           varname("SF <VV|VV>"));
    global_dpd_->buf4_copy(&Gab, PSIF_DCT_DENSITY, varname("<VV|VV>"));
    global_dpd_->buf4_close(&Gab);

    psio_->close(PSIF_DCT_DENSITY, 1);
}

void DCTSolver::compute_unrelaxed_separable_density_VVVV_RHF() {
    dpdbuf4 Gab;

    global_dpd_->buf4_init(&Gab, PSIF_DCT_DENSITY, 0, ID("[V,V]"), ID("[V,V]"), ID("[V,V]"), ID("[V,V]"), 0,
                           "Gamma SF <VV|VV>");  // Gamma <Vv|Vv>
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&Gab, h);
        global_dpd_->buf4_mat_irrep_rd(&Gab, h);

#pragma omp parallel for
        for (long int ab = 0; ab < Gab.params->rowtot[h]; ++ab) {
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
                if (Ga == Gc && Gb == Gd) tpdm += 0.25 * avir_tau_.get(Ga, a, c) * avir_tau_.get(Gb, b, d);
                Gab.matrix[h][ab][cd] += tpdm;
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&Gab, h);
        global_dpd_->buf4_mat_irrep_close(&Gab, h);
    }

    global_dpd_->buf4_close(&Gab);
}

void DCTSolver::construct_oo_density_RHF() {
    Da_ = std::make_shared<Matrix>(std::move(construct_oo_density(aocc_tau_, avir_tau_, *kappa_mo_a_, *Ca_)));
    Db_ = Da_->clone();
}

}  // namespace dct
}  // namespace psi
