/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#include "psi4/libdpd/dpd.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libtrans/integraltransform.h"

namespace psi {
namespace dct {

/**
 * Updates the MO coefficients, transforms the integrals into both chemists'
 * and physcists' notation
 * for RHF reference.
 */
void DCTSolver::transform_integrals_RHF() {
    dct_timer_on("DCTSolver::transform_integrals()");

    if (options_.get_str("DCT_TYPE") == "DF") {
        // Transform b(Q|mn) to b(Q|pq) in MO basis
        transform_b_so2mo();

        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

        /*- Transform g(OV|OV) -*/
        form_df_g_ovov();
        /*- Transform g(OO|OO) -*/
        form_df_g_oooo();
        /*- Transform g(VV|OO) -*/
        form_df_g_vvoo();

        if (orbital_optimized_) {
            /*- Transform g(VO|OO) -*/
            form_df_g_vooo();
            /*- Transform g(OV|VV) -*/
            form_df_g_ovvv();
        }

        psio_->close(PSIF_LIBTRANS_DPD, 1);

        _ints->update_orbitals();
    } else {
        _ints->update_orbitals();

        if (print_ > 1) {
            outfile->Printf("\tTransforming integrals...\n");
        }
        _ints->set_print(print_ - 2 >= 0 ? print_ - 2 : 0);

        // Generate the integrals in various spaces in chemists' notation

        dct_timer_on("DCTSolver::transform_OVOV");
        _ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);
        dct_timer_off("DCTSolver::transform_OVOV");

        dct_timer_on("DCTSolver::transform_OOOO");
        _ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::occ, MOSpace::occ);
        dct_timer_off("DCTSolver::transform_OOOO");

        dct_timer_on("DCTSolver::transform_OOVV");
        _ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::occ, MOSpace::occ);
        dct_timer_off("DCTSolver::transform_OOVV");

        if ((options_.get_str("AO_BASIS") == "NONE")) {
            _ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::vir, MOSpace::vir);
        }
        if ((options_.get_str("ALGORITHM") == "QC" && options_.get_bool("QC_COUPLING") &&
             options_.get_str("QC_TYPE") == "SIMULTANEOUS") ||
            orbital_optimized_) {
            // Compute the integrals needed for the MO Hessian
            dct_timer_on("DCTSolver::transform_VOOO");
            _ints->transform_tei(MOSpace::vir, MOSpace::occ, MOSpace::occ, MOSpace::occ);
            dct_timer_off("DCTSolver::transform_VOOO");

            dct_timer_on("DCTSolver::transform_OVVV");
            _ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::vir, MOSpace::vir);
            dct_timer_off("DCTSolver::transform_OVVV");
        }
    }

    /*
     * Re-sort the chemists' notation integrals to physisicts' notation
     * (pq|rs) = <pr|qs>
     */

    // The integral object closes this file - we need to re-open it.
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    sort_OVOV_integrals_RHF();

    sort_OOOO_integrals_RHF();

    sort_OOVV_integrals_RHF();

    if (options_.get_str("AO_BASIS") == "NONE" && options_.get_str("DCT_TYPE") == "CONV") sort_VVVV_integrals_RHF();

    // VVVO and OOOV integrals are needed for the QC algorithm
    if ((options_.get_str("ALGORITHM") == "QC" && options_.get_bool("QC_COUPLING") &&
         options_.get_str("QC_TYPE") == "SIMULTANEOUS") ||
        orbital_optimized_) {
        sort_OOOV_integrals_RHF();

        sort_OVVV_integrals_RHF();
    }

    if (orbital_optimized_) transform_core_integrals_RHF();

    // The integral transformation object also provided us with the Fock matrix
    // in the current basis, from which we can get the new denominator matrices.
    // N.B. These are not neccesarily the eigenvalues, rather they are the diagonal
    // elements of F0 in the current basis; F is diagonal, not F0.
    build_denominators_RHF();

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    dct_timer_off("DCTSolver::transform_integrals()");
}

void DCTSolver::sort_OVOV_integrals_RHF() {
    dpdbuf4 I;
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "MO Ints (OV|OV)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[O,O]"), ID("[V,V]"),
                           "MO Ints <OO|VV>");  // MO Ints <Oo|Vv>
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, psrq, ID("[O,V]"), ID("[O,V]"), "MO Ints SF <OV|OV>:<Ov|oV>");
    global_dpd_->buf4_close(&I);
}

void DCTSolver::sort_OOOO_integrals_RHF() {
    dpdbuf4 I;
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"), ID("[O>=O]+"), ID("[O>=O]+"), 0,
                           "MO Ints (OO|OO)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[O,O]"), ID("[O,O]"),
                           "MO Ints <OO|OO>");  // MO Ints <Oo|Oo>
    global_dpd_->buf4_close(&I);
}

void DCTSolver::sort_OOVV_integrals_RHF() {
    dpdbuf4 I, Irs, Isr;

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,O]"), ID("[V>=V]+"), ID("[O>=O]+"), 0,
                           "MO Ints (VV|OO)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, sqrp, ID("[O,V]"), ID("[O,V]"),
                           "MO Ints <OV|OV>");  // MO Ints <oV|oV>, MO Ints <Ov|Ov>, MO Ints <ov|ov>
    global_dpd_->buf4_close(&I);

    /*
     * Antisymmetrize the <OV|OV> integrals
     */

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "MO Ints <OV|OV>");
    global_dpd_->buf4_copy(&I, PSIF_LIBTRANS_DPD, "MO Ints <OV|OV> - <OV|VO>");
    global_dpd_->buf4_close(&I);
    // Resort (OV|OV) to the physists' notation
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "MO Ints (OV|OV)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, psrq, ID("[O,V]"), ID("[O,V]"), "MO Ints <PS|RQ>");
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_init(&Irs, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "MO Ints <OV|OV> - <OV|VO>");
    global_dpd_->buf4_init(&Isr, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           "MO Ints <PS|RQ>");
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&Irs, h);
        global_dpd_->buf4_mat_irrep_init(&Isr, h);
        global_dpd_->buf4_mat_irrep_rd(&Irs, h);
        global_dpd_->buf4_mat_irrep_rd(&Isr, h);
        for (int row = 0; row < Irs.params->rowtot[h]; ++row) {
            for (int col = 0; col < Irs.params->coltot[h]; ++col) {
                Irs.matrix[h][row][col] -= Isr.matrix[h][row][col];
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&Irs, h);
        global_dpd_->buf4_mat_irrep_close(&Irs, h);
        global_dpd_->buf4_mat_irrep_close(&Isr, h);
    }
    global_dpd_->buf4_close(&Isr);
    global_dpd_->buf4_close(&Irs);
}

void DCTSolver::sort_VVVV_integrals_RHF() {
    dpdbuf4 I;
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"), ID("[V>=V]+"), ID("[V>=V]+"), 0,
                           "MO Ints (VV|VV)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[V,V]"), ID("[V,V]"), "MO Ints <VV|VV>");
    global_dpd_->buf4_close(&I);
}

void DCTSolver::sort_OOOV_integrals_RHF() {
    dpdbuf4 I;
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"), ID("[V,O]"), ID("[O>=O]+"), 0,
                           "MO Ints (VO|OO)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, rpsq, ID("[O,V]"), ID("[O,O]"),
                           "MO Ints <OV|OO>");  // MO Ints <oV|oO>, MO Ints <Ov|Oo>
    // Intermediate MO_SF <OV|OO> = MO <Ov|oO>
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, rpqs, ID("[O,V]"), ID("[O,O]"), "MO Ints SF <OV|OO>");
    global_dpd_->buf4_close(&I);
}

void DCTSolver::sort_OVVV_integrals_RHF() {
    dpdbuf4 I;
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"), ID("[O,V]"), ID("[V>=V]+"), 0,
                           "MO Ints (OV|VV)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[O,V]"), ID("[V,V]"),
                           "MO Ints <OV|VV>");  // MO Ints <Ov|Vv>
    // Intermediate MO_SF <OV|VV> = MO <oV|Vv>
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, prsq, ID("[O,V]"), ID("[V,V]"), "MO Ints SF <OV|VV>");
    global_dpd_->buf4_close(&I);
}

void DCTSolver::transform_core_integrals_RHF() {
    // Transform one-electron integrals to the MO basis and store them in the DPD file
    dpdfile2 H;
    Matrix aH(so_h_);
    aH.transform(Ca_);

    global_dpd_->file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "H <O|O>");
    const auto& O_slice = slices_.at("ACTIVE_OCC_A");
    auto temp = *aH.get_block(O_slice);
    temp.write_to_dpdfile2(&H);
    global_dpd_->file2_close(&H);

    global_dpd_->file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "H <V|V>");
    const auto &V_slice = slices_.at("ACTIVE_VIR_A");
    temp = *aH.get_block(V_slice);
    temp.write_to_dpdfile2(&H);
    global_dpd_->file2_close(&H);

    global_dpd_->file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('V'), "H <O|V>");
    temp = *aH.get_block(O_slice, V_slice);
    temp.write_to_dpdfile2(&H);
    global_dpd_->file2_close(&H);
}

/**
 * Consructs the denominators using the diagonal elements of the unperturbed Fock matrix.
 * Also builds the MO coefficient tensors in the alpha/beta, occ/vir spaces
 * for RHF reference
 */
void DCTSolver::build_denominators_RHF() {
    dct_timer_on("DCTSolver::build_denominators()");

    dpdbuf4 D;
    dpdfile2 F;

    auto aOccEvals = std::vector<double>(nalpha_);
    auto aVirEvals = std::vector<double>(navir_);
    // Pick out the diagonal elements of the Fock matrix, making sure that they are in the order
    // used by the DPD library, i.e. starting from zero for each space and ordering by irrep
    int aOccCount = 0, aVirCount = 0;

    dpdfile2 T_OO, T_VV;
    global_dpd_->file2_init(&T_OO, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    global_dpd_->file2_init(&T_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    global_dpd_->file2_mat_init(&T_OO);
    global_dpd_->file2_mat_init(&T_VV);
    global_dpd_->file2_mat_rd(&T_OO);
    global_dpd_->file2_mat_rd(&T_VV);

    // Diagonal elements of the Fock matrix
    // Alpha spin
    aocc_c_ = Ca_->get_block(slices_.at("SO"), slices_.at("ACTIVE_OCC_A"));
    avir_c_ = Ca_->get_block(slices_.at("SO"), slices_.at("ACTIVE_VIR_A"));
    for (int h = 0; h < nirrep_; ++h) {
        for (int i = 0; i < naoccpi_[h]; ++i) {
            if (!exact_tau_) {
                aOccEvals[aOccCount++] = moFa_->get(h, i, i);
            } else {
                aOccEvals[aOccCount++] = moFa_->get(h, i, i) / (1.0 + 2.0 * T_OO.matrix[h][i][i]);
            }
        }

        for (int a = 0; a < navirpi_[h]; ++a) {
            if (!exact_tau_) {
                aVirEvals[aVirCount++] = moFa_->get(h, naoccpi_[h] + a, naoccpi_[h] + a);
            } else {
                aVirEvals[aVirCount++] =
                    moFa_->get(h, a + naoccpi_[h], a + naoccpi_[h]) / (1.0 - 2.0 * T_VV.matrix[h][a][a]);
            }
        }
    }

    global_dpd_->file2_mat_close(&T_OO);
    global_dpd_->file2_mat_close(&T_VV);
    global_dpd_->file2_close(&T_OO);
    global_dpd_->file2_close(&T_VV);

    // Elements of the Fock matrix
    if (!exact_tau_) {
        // Alpha occupied
        global_dpd_->file2_init(&F, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "F <O|O>");
        global_dpd_->file2_mat_init(&F);
        for (int h = 0; h < nirrep_; ++h) {
            for (int i = 0; i < naoccpi_[h]; ++i) {
                for (int j = 0; j < naoccpi_[h]; ++j) {
                    F.matrix[h][i][j] = moFa_->get(h, i, j);
                }
            }
        }
        global_dpd_->file2_mat_wrt(&F);
        global_dpd_->file2_mat_close(&F);
        global_dpd_->file2_close(&F);

        // Alpha Virtual
        global_dpd_->file2_init(&F, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "F <V|V>");
        global_dpd_->file2_mat_init(&F);
        for (int h = 0; h < nirrep_; ++h) {
            for (int i = 0; i < navirpi_[h]; ++i) {
                for (int j = 0; j < navirpi_[h]; ++j) {
                    F.matrix[h][i][j] = moFa_->get(h, i + naoccpi_[h], j + naoccpi_[h]);
                }
            }
        }
        global_dpd_->file2_mat_wrt(&F);
        global_dpd_->file2_mat_close(&F);
        global_dpd_->file2_close(&F);
    }

    global_dpd_->buf4_init(&D, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>=O]+"), ID("[V>=V]+"), 0, "D <OO|VV>");
    for (int h = 0; h < nirrep_; ++h) {
        global_dpd_->buf4_mat_irrep_init(&D, h);
        for (int row = 0; row < D.params->rowtot[h]; ++row) {
            int i = D.params->roworb[h][row][0];
            int j = D.params->roworb[h][row][1];
            for (int col = 0; col < D.params->coltot[h]; ++col) {
                int a = D.params->colorb[h][col][0];
                int b = D.params->colorb[h][col][1];
                D.matrix[h][row][col] =
                    1.0 / (regularizer_ + aOccEvals[i] + aOccEvals[j] - aVirEvals[a] - aVirEvals[b]);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&D, h);
        global_dpd_->buf4_mat_irrep_close(&D, h);
    }
    global_dpd_->buf4_close(&D);

    dct_timer_off("DCTSolver::build_denominators()");
}

}  // namespace dct
}  // namespace psi
