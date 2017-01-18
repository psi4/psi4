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
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libtrans/integraltransform.h"
#include "defines.h"



namespace psi{ namespace dcft{

/**
 * Updates the MO coefficients, transforms the integrals into both chemists'
 * and physcists' notation.
 */
void
DCFTSolver::transform_integrals()
{
    dcft_timer_on("DCFTSolver::transform_integrals()");

    if (options_.get_str("DCFT_TYPE") == "DF"){
        // Transform b(Q|mn) to b(Q|pq) in MO basis
        transform_b();
        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
        /*- Transform g(OV|OV) -*/
        form_df_g_ovov();
        /*- Transform g(OO|OO) -*/
        form_df_g_oooo();
        /*- Transform g(VV|OO) -*/
        form_df_g_vvoo();

        if(orbital_optimized_){
            /*- Transform g(VO|OO) -*/
            form_df_g_vooo();
            /*- Transform g(OV|VV) -*/
            form_df_g_ovvv();
        }

        psio_->close(PSIF_LIBTRANS_DPD, 1);

        _ints->update_orbitals();
    }
    else{
        _ints->update_orbitals();
        if(print_ > 1){
            outfile->Printf( "\tTransforming integrals...\n");

        }
        _ints->set_print(print_ - 2 >= 0 ? print_ - 2 : 0);

        // Generate the integrals in various spaces in chemists' notation
        _ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);
        _ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::occ, MOSpace::occ);
        _ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::occ, MOSpace::occ);
        _ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::vir);
        if((options_.get_str("AO_BASIS") == "NONE")){
            _ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::vir, MOSpace::vir);
        }
        if ((options_.get_str("ALGORITHM") == "QC" && options_.get_bool("QC_COUPLING")
                                                   && options_.get_str("QC_TYPE") == "SIMULTANEOUS")
                                                   || orbital_optimized_) {
            // Compute the integrals needed for the MO Hessian
            _ints->transform_tei(MOSpace::vir, MOSpace::occ, MOSpace::occ, MOSpace::occ);
            _ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::occ);
            _ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::vir, MOSpace::vir);
            _ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::occ, MOSpace::vir);

        }
    }

    /*
     * Re-sort the chemists' notation integrals to physisists' notation
     * (pq|rs) = <pr|qs>
     */

    // The integral object closes this file - we need to re-open it.
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    sort_OVOV_integrals();

    sort_OOOO_integrals();

    sort_OOVV_integrals();

    if(options_.get_str("AO_BASIS") == "NONE" && options_.get_str("DCFT_TYPE") == "CONV")
        sort_VVVV_integrals();

    // VVVO and OOOV integrals are needed for the QC algorithm
    if ((options_.get_str("ALGORITHM") == "QC" && options_.get_bool("QC_COUPLING")
                                               && options_.get_str("QC_TYPE") == "SIMULTANEOUS")
                                               || orbital_optimized_) {

        sort_OOOV_integrals();

        sort_OVVV_integrals();

    }

    if (orbital_optimized_) transform_core_integrals();

    // The integral transformation object also provided us with the Fock matrix
    // in the current basis, from which we can get the new denominator matrices.
    // N.B. These are not neccesarily the eigenvalues, rather they are the diagonal
    // elements of F0 in the current basis; F is diagonal, not F0.
    build_denominators();

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    dcft_timer_off("DCFTSolver::transform_integrals()");
}

void
DCFTSolver:: sort_OVOV_integrals() {

    dpdbuf4 I;

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"),0, "MO Ints (OV|OV)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "MO Ints <OO|VV>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"),0, "MO Ints <OO|VV>");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[V,V]"), ID("[O,O]"), "MO Ints <VV|OO>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[O,o]"), ID("[V,v]"), "MO Ints <Oo|Vv>");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, psrq, ID("[O,v]"), ID("[o,V]"), "MO Ints <Ov|oV>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov|ov)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[o,o]"), ID("[v,v]"), "MO Ints <oo|vv>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo|vv>");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[v,v]"), ID("[o,o]"), "MO Ints <vv|oo>");
    global_dpd_->buf4_close(&I);

    if ((options_.get_str("ALGORITHM") == "QC" && options_.get_bool("QC_COUPLING")
                                               && options_.get_str("QC_TYPE") == "SIMULTANEOUS") || orbital_optimized_) {

        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                      ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
        global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[o,v]"), ID("[O,V]"), "MO Ints (ov|OV)");
        global_dpd_->buf4_close(&I);

        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
        global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, rqps, ID("[V,O]"), ID("[O,V]"), "MO Ints <VO|OV>");
        global_dpd_->buf4_close(&I);

        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo|vv>");
        global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, rqps, ID("[v,o]"), ID("[o,v]"), "MO Ints <vo|ov>");
        global_dpd_->buf4_close(&I);

    }

    if (options_.get_str("DCFT_FUNCTIONAL") == "ODC-13") {

        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
        global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, qprs, ID("[V,O]"), ID("[O,V]"), "MO Ints (VO|OV)");
        global_dpd_->buf4_close(&I);

        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                      ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
        global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, qprs, ID("[V,O]"), ID("[o,v]"), "MO Ints (VO|ov)");
        global_dpd_->buf4_close(&I);

        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                      ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
        global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, srpq, ID("[v,o]"), ID("[O,V]"), "MO Ints (vo|OV)");
        global_dpd_->buf4_close(&I);

        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                      ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov|ov)");
        global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, qprs, ID("[v,o]"), ID("[o,v]"), "MO Ints (vo|ov)");
        global_dpd_->buf4_close(&I);

        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                      ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
        global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, psqr, ID("[O,v]"), ID("[V,o]"), "MO Ints <Ov|Vo>");
        global_dpd_->buf4_close(&I);

    }

}

void
DCFTSolver:: sort_OOOO_integrals() {

    dpdbuf4 I;

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O>=O]+"), ID("[O>=O]+"), 0, "MO Ints (OO|OO)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[O,O]"), ID("[O,O]"), "MO Ints <OO|OO>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[o,o]"),
                  ID("[O>=O]+"), ID("[o>=o]+"), 0, "MO Ints (OO|oo)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[O,o]"), ID("[O,o]"), "MO Ints <Oo|Oo>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[o,o]"),
                  ID("[O>=O]+"), ID("[o>=o]+"), 0, "MO Ints (OO|oo)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[o,o]"), ID("[O,O]"), "MO Ints (oo|OO)");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o>=o]+"), ID("[o>=o]+"), 0, "MO Ints (oo|oo)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[o,o]"), ID("[o,o]"), "MO Ints <oo|oo>");
    global_dpd_->buf4_close(&I);

}

void
DCFTSolver:: sort_OOVV_integrals() {

    dpdbuf4 I, Irs, Isr;

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[o,o]"),
                  ID("[V>=V]+"), ID("[o>=o]+"), 0, "MO Ints (VV|oo)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, sqrp, ID("[o,V]"), ID("[o,V]"), "MO Ints <oV|oV>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[o,o]"),
                  ID("[V>=V]+"), ID("[o>=o]+"), 0, "MO Ints (VV|oo)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[o,o]"), ID("[V,V]"), "MO Ints (oo|VV)");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[v,v]"),
                  ID("[O>=O]+"), ID("[v>=v]+"), 0, "MO Ints (OO|vv)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[v,v]"), ID("[O,O]"), "MO Ints (vv|OO)");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>=O]+"), ID("[V>=V]+"), 0, "MO Ints (OO|VV)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[O,V]"), ID("[O,V]"), "MO Ints <OV|OV>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[v,v]"),
                  ID("[O>=O]+"), ID("[v>=v]+"), 0, "MO Ints (OO|vv)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[O,v]"), ID("[O,v]"), "MO Ints <Ov|Ov>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>=o]+"), ID("[v>=v]+"), 0, "MO Ints (oo|vv)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[o,v]"), ID("[o,v]"), "MO Ints <ov|ov>");
    global_dpd_->buf4_close(&I);

    /*
     * Antisymmetrize the <OV|OV> and <ov|ov> integrals
     */
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>");
    global_dpd_->buf4_copy(&I, PSIF_LIBTRANS_DPD, "MO Ints <OV|OV> - <OV|VO>");
    global_dpd_->buf4_close(&I);
    // Resort (OV|OV) to the physists' notation
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, psrq, ID("[O,V]"), ID("[O,V]"), "MO Ints <PS|RQ>");
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_init(&Irs, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV> - <OV|VO>");
    global_dpd_->buf4_init(&Isr, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <PS|RQ>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&Irs, h);
        global_dpd_->buf4_mat_irrep_init(&Isr, h);
        global_dpd_->buf4_mat_irrep_rd(&Irs, h);
        global_dpd_->buf4_mat_irrep_rd(&Isr, h);
        for(int row = 0; row < Irs.params->rowtot[h]; ++row){
            for(int col = 0; col < Irs.params->coltot[h]; ++col){
                Irs.matrix[h][row][col] -= Isr.matrix[h][row][col];
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&Irs, h);
        global_dpd_->buf4_mat_irrep_close(&Irs, h);
        global_dpd_->buf4_mat_irrep_close(&Isr, h);
    }

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov|ov>");
    global_dpd_->buf4_copy(&I, PSIF_LIBTRANS_DPD, "MO Ints <ov|ov> - <ov|vo>");
    global_dpd_->buf4_close(&I);
    // Resort (ov|ov) to the physists' notation
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov|ov)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, psrq, ID("[o,v]"), ID("[o,v]"), "MO Ints <ps|rq>");
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_init(&Irs, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov|ov> - <ov|vo>");
    global_dpd_->buf4_init(&Isr, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ps|rq>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&Irs, h);
        global_dpd_->buf4_mat_irrep_init(&Isr, h);
        global_dpd_->buf4_mat_irrep_rd(&Irs, h);
        global_dpd_->buf4_mat_irrep_rd(&Isr, h);
        for(int row = 0; row < Irs.params->rowtot[h]; ++row){
            for(int col = 0; col < Irs.params->coltot[h]; ++col){
                Irs.matrix[h][row][col] -= Isr.matrix[h][row][col];
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&Irs, h);
        global_dpd_->buf4_mat_irrep_close(&Irs, h);
        global_dpd_->buf4_mat_irrep_close(&Isr, h);
    }

}

void
DCFTSolver:: sort_VVVV_integrals() {

    dpdbuf4 I;

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V>=V]+"), ID("[V>=V]+"), 0, "MO Ints (VV|VV)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[V,V]"), ID("[V,V]"), "MO Ints <VV|VV>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[v,v]"),
                  ID("[V>=V]+"), ID("[v>=v]+"), 0, "MO Ints (VV|vv)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[V,v]"), ID("[V,v]"), "MO Ints <Vv|Vv>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[v,v]"),
                  ID("[V>=V]+"), ID("[v>=v]+"), 0, "MO Ints (VV|vv)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[v,v]"), ID("[V,V]"), "MO Ints (vv|VV)");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v>=v]+"), ID("[v>=v]+"), 0, "MO Ints (vv|vv)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[v,v]"), ID("[v,v]"), "MO Ints <vv|vv>");
    global_dpd_->buf4_close(&I);

}

void
DCFTSolver:: sort_OOOV_integrals() {

    // <VO|OO> type

    dpdbuf4 I;

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O>=O]+"), 0, "MO Ints (VO|OO)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[V,O]"), ID("[O,O]"), "MO Ints <VO|OO>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O>=O]+"), 0, "MO Ints (VO|OO)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, rsqp, ID("[O,O]"), ID("[O,V]"), "MO Ints (OO|OV)");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 0, "MO Ints <VO|OO>");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, qpsr, ID("[O,V]"), ID("[O,O]"), "MO Ints <OV|OO>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O,O]"), 0, "MO Ints <OV|OO>");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[O,O]"), ID("[O,V]"), "MO Ints <OO|OV>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[o,o]"),
                  ID("[V,O]"), ID("[o>=o]+"), 0, "MO Ints (VO|oo)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[V,o]"), ID("[O,o]"), "MO Ints <Vo|Oo>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[o,o]"),
                  ID("[V,O]"), ID("[o>=o]+"), 0, "MO Ints (VO|oo)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[o,o]"), ID("[V,O]"), "MO Ints (oo|VO)");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[O,o]"),
                  ID("[V,o]"), ID("[O,o]"), 0, "MO Ints <Vo|Oo>");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, qpsr, ID("[o,V]"), ID("[o,O]"), "MO Ints <oV|oO>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[o,O]"),
                  ID("[o,V]"), ID("[o,O]"), 0, "MO Ints <oV|oO>");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, pqsr, ID("[o,V]"), ID("[O,o]"), "MO Ints <oV|Oo>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[v,o]"),
                  ID("[O>=O]+"), ID("[v,o]"), 0, "MO Ints (OO|vo)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, rqsp, ID("[v,O]"), ID("[o,O]"), "MO Ints <vO|oO>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,O]"), ID("[o,O]"),
                  ID("[v,O]"), ID("[o,O]"), 0, "MO Ints <vO|oO>");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, qpsr, ID("[O,v]"), ID("[O,o]"), "MO Ints <Ov|Oo>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "MO Ints <Ov|Oo>");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, pqsr, ID("[O,v]"), ID("[o,O]"), "MO Ints <Ov|oO>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O>=O]+"), ID("[v,o]"),
                  ID("[O>=O]+"), ID("[v,o]"), 0, "MO Ints (OO|vo)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[v,o]"), ID("[O>=O]+"), "MO Ints (vo|OO)");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o>=o]+"), 0, "MO Ints (vo|oo)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[v,o]"), ID("[o,o]"), "MO Ints <vo|oo>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o>=o]+"), 0, "MO Ints (vo|oo)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, rsqp, ID("[o,o]"), ID("[o,v]"), "MO Ints (oo|ov)");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o,o]"), 0, "MO Ints <vo|oo>");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, qpsr, ID("[o,v]"), ID("[o,o]"), "MO Ints <ov|oo>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,o]"),
                  ID("[o,v]"), ID("[o,o]"), 0, "MO Ints <ov|oo>");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[o,o]"), ID("[o,v]"), "MO Ints <oo|ov>");
    global_dpd_->buf4_close(&I);

}

void
DCFTSolver:: sort_OVVV_integrals() {

    // <OV|VV> type

    dpdbuf4 I;

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V>=V]+"), 0, "MO Ints (OV|VV)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[O,V]"), ID("[V,V]"), "MO Ints <OV|VV>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[V,V]"), ID("[O,V]"), "MO Ints <VV|OV>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[v,v]"),
                  ID("[O,V]"), ID("[v>=v]+"), 0, "MO Ints (OV|vv)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[O,v]"), ID("[V,v]"), "MO Ints <Ov|Vv>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[v,v]"),
                  ID("[O,V]"), ID("[v>=v]+"), 0, "MO Ints (OV|vv)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[v,v]"), ID("[O,V]"), "MO Ints (vv|OV)");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv>");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, pqsr, ID("[O,v]"), ID("[v,V]"), "MO Ints <Ov|vV>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[o,v]"),
                  ID("[V>=V]+"), ID("[o,v]"), 0, "MO Ints (VV|ov)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, rqsp, ID("[o,V]"), ID("[v,V]"), "MO Ints <oV|vV>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[v,V]"),
                  ID("[o,V]"), ID("[v,V]"), 0, "MO Ints <oV|vV>");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, pqsr, ID("[o,V]"), ID("[V,v]"), "MO Ints <oV|Vv>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[v,V]"),
                  ID("[o,V]"), ID("[v,V]"), 0, "MO Ints <oV|vV>");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, qpsr, ID("[V,o]"), ID("[V,v]"), "MO Ints <Vo|Vv>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V>=V]+"), ID("[o,v]"),
                  ID("[V>=V]+"), ID("[o,v]"), 0, "MO Ints (VV|ov)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[o,v]"), ID("[V>=V]+"), "MO Ints (ov|VV)");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v>=v]+"), 0, "MO Ints (ov|vv)");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[o,v]"), ID("[v,v]"), "MO Ints <ov|vv>");
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov|vv>");
    global_dpd_->buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[v,v]"), ID("[o,v]"), "MO Ints <vv|ov>");
    global_dpd_->buf4_close(&I);

}

void DCFTSolver::transform_core_integrals() {

    // Transform one-electron integrals to the MO basis and store them in the DPD file
    dpdfile2 H;
    Matrix aH(so_h_);
    Matrix bH(so_h_);
    aH.transform(Ca_);
    bH.transform(Cb_);

    global_dpd_->file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "H <O|O>");
    global_dpd_->file2_mat_init(&H);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < naoccpi_[h]; ++i){
            for(int j = 0 ; j < naoccpi_[h]; ++j){
                H.matrix[h][i][j] = aH.get(h, i, j);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&H);
    global_dpd_->file2_close(&H);

    global_dpd_->file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "H <V|V>");
    global_dpd_->file2_mat_init(&H);
    for(int h = 0; h < nirrep_; ++h){
        for(int a = 0 ; a < navirpi_[h]; ++a){
            for(int b = 0 ; b < navirpi_[h]; ++b){
                H.matrix[h][a][b] = aH.get(h, naoccpi_[h] + a, naoccpi_[h] + b);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&H);
    global_dpd_->file2_close(&H);

    global_dpd_->file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "H <o|o>");
    global_dpd_->file2_mat_init(&H);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < nboccpi_[h]; ++i){
            for(int j = 0 ; j < nboccpi_[h]; ++j){
                H.matrix[h][i][j] = bH.get(h, i, j);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&H);
    global_dpd_->file2_close(&H);

    global_dpd_->file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "H <v|v>");
    global_dpd_->file2_mat_init(&H);
    for(int h = 0; h < nirrep_; ++h){
        for(int a = 0 ; a < nbvirpi_[h]; ++a){
            for(int b = 0 ; b < nbvirpi_[h]; ++b){
                H.matrix[h][a][b] = bH.get(h, nboccpi_[h] + a, nboccpi_[h] + b);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&H);
    global_dpd_->file2_close(&H);

    global_dpd_->file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('V'), "H <O|V>");
    global_dpd_->file2_mat_init(&H);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < naoccpi_[h]; ++i){
            for(int j = 0 ; j < navirpi_[h]; ++j){
                H.matrix[h][i][j] = aH.get(h, i, naoccpi_[h] + j);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&H);
    global_dpd_->file2_close(&H);

    global_dpd_->file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('v'), "H <o|v>");
    global_dpd_->file2_mat_init(&H);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < nboccpi_[h]; ++i){
            for(int j = 0 ; j < nbvirpi_[h]; ++j){
                H.matrix[h][i][j] = bH.get(h, i, nboccpi_[h] + j);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&H);
    global_dpd_->file2_close(&H);

}

/**
 * Consructs the denominators using the diagonal elements of the unperturbed Fock matrix.
 * Also builds the MO coefficient tensors in the alpha/beta, occ/vir spaces.
 */
void
DCFTSolver::build_denominators()
{
    dcft_timer_on("DCFTSolver::build_denominators()");

    dpdbuf4 D;
    dpdfile2 F;

    double *aOccEvals = new double [nalpha_];
    double *bOccEvals = new double [nbeta_];
    double *aVirEvals = new double [navir_];
    double *bVirEvals = new double [nbvir_];
    // Pick out the diagonal elements of the Fock matrix, making sure that they are in the order
    // used by the DPD library, i.e. starting from zero for each space and ordering by irrep
    int aOccCount = 0, bOccCount = 0, aVirCount = 0, bVirCount = 0;

    dpdfile2 T_OO, T_oo, T_VV, T_vv;
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

    //Diagonal elements of the Fock matrix
    //Alpha spin
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < naoccpi_[h]; ++i){
            if (!exact_tau_) {
                aOccEvals[aOccCount++] = moFa_->get(h, i, i);
            }
            else {
                aOccEvals[aOccCount++] = moFa_->get(h, i, i) / (1.0 + 2.0 * T_OO.matrix[h][i][i]);
            }
            for(int mu = 0; mu < nsopi_[h]; ++mu)
                aocc_c_->set(h, mu, i, Ca_->get(h, mu, i));
        }

        for(int a = 0; a < navirpi_[h]; ++a){
            if (!exact_tau_) {
                aVirEvals[aVirCount++] = moFa_->get(h, naoccpi_[h] + a, naoccpi_[h] + a);
            }
            else {
                aVirEvals[aVirCount++] = moFa_->get(h, a + naoccpi_[h], a + naoccpi_[h]) / (1.0 - 2.0 * T_VV.matrix[h][a][a]);
            }
            for(int mu = 0; mu < nsopi_[h]; ++mu)
                avir_c_->set(h, mu, a, Ca_->get(h, mu, naoccpi_[h] + a));
        }
    }

    //Elements of the Fock matrix
    //Alpha occupied

    if (!exact_tau_) {
        global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
        global_dpd_->file2_mat_init(&F);
        int offset = 0;
        for(int h = 0; h < nirrep_; ++h){
            offset += frzcpi_[h];
            for(int i = 0 ; i < naoccpi_[h]; ++i){
                for(int j = 0 ; j < naoccpi_[h]; ++j){
                    F.matrix[h][i][j] = moFa_->get(h, i, j);
                }
            }
            offset += nmopi_[h];
        }
        global_dpd_->file2_mat_wrt(&F);
        global_dpd_->file2_close(&F);

        //Alpha Virtual
        global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");
        global_dpd_->file2_mat_init(&F);
        offset = 0;
        for(int h = 0; h < nirrep_; ++h){
            offset += naoccpi_[h];
            for(int i = 0 ; i < navirpi_[h]; ++i){
                for(int j = 0 ; j < navirpi_[h]; ++j){
                    F.matrix[h][i][j] = moFa_->get(h, i + naoccpi_[h], j + naoccpi_[h]);
                }
            }
            offset += nmopi_[h] - naoccpi_[h];
        }
        global_dpd_->file2_mat_wrt(&F);
        global_dpd_->file2_close(&F);
    }

    //Diagonal elements of the Fock matrix
    //Beta spin
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < nboccpi_[h]; ++i){
            if (!exact_tau_) {
                bOccEvals[bOccCount++] = moFb_->get(h, i, i);
            }
            else {
                bOccEvals[bOccCount++] = moFb_->get(h, i, i) / (1.0 + 2.0 * T_oo.matrix[h][i][i]);
            }
            for(int mu = 0; mu < nsopi_[h]; ++mu)
                bocc_c_->set(h, mu, i, Cb_->get(h, mu, i));
        }
        for(int a = 0; a < nbvirpi_[h]; ++a){
            if (!exact_tau_) {
                bVirEvals[bVirCount++] = moFb_->get(h, nboccpi_[h] + a, nboccpi_[h] + a);
            }
            else {
                bVirEvals[bVirCount++] = moFb_->get(h, a + nboccpi_[h], a + nboccpi_[h]) / (1.0 - 2.0 * T_vv.matrix[h][a][a]);
            }
            for(int mu = 0; mu < nsopi_[h]; ++mu)
                bvir_c_->set(h, mu, a, Cb_->get(h, mu, nboccpi_[h] + a));
        }
    }

    //Off-diagonal elements of the Fock matrix
    //Beta Occupied

    if (!exact_tau_) {
        global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "F <o|o>");
        global_dpd_->file2_mat_init(&F);
        int offset = 0;
        for(int h = 0; h < nirrep_; ++h){
            offset += frzcpi_[h];
            for(int i = 0 ; i < nboccpi_[h]; ++i){
                for(int j = 0 ; j < nboccpi_[h]; ++j){
                    F.matrix[h][i][j] = moFb_->get(h, i, j);
                }
            }
            offset += nmopi_[h];
        }
        global_dpd_->file2_mat_wrt(&F);
        global_dpd_->file2_close(&F);

        //Beta Virtual
        global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "F <v|v>");
        global_dpd_->file2_mat_init(&F);
        offset = 0;
        for(int h = 0; h < nirrep_; ++h){
            offset += nboccpi_[h];
            for(int i = 0 ; i < nbvirpi_[h]; ++i){
                for(int j = 0 ; j < nbvirpi_[h]; ++j){
                    F.matrix[h][i][j] = moFb_->get(h, i + nboccpi_[h], j + nboccpi_[h]);
                }
            }
            offset += nmopi_[h] - nboccpi_[h];
        }
        global_dpd_->file2_mat_wrt(&F);
        global_dpd_->file2_close(&F);
    }

    global_dpd_->file2_close(&T_OO);
    global_dpd_->file2_close(&T_oo);
    global_dpd_->file2_close(&T_VV);
    global_dpd_->file2_close(&T_vv);

    ///////////////////////////////
    // The alpha-alpha spin case //
    ///////////////////////////////
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>=O]+"), ID("[V>=V]+"), 0, "D <OO|VV>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            int i = D.params->roworb[h][row][0];
            int j = D.params->roworb[h][row][1];
            for(int col = 0; col < D.params->coltot[h]; ++col){
                int a = D.params->colorb[h][col][0];
                int b = D.params->colorb[h][col][1];
                D.matrix[h][row][col] = 1.0/
                                    (regularizer_ + aOccEvals[i] + aOccEvals[j] - aVirEvals[a] - aVirEvals[b]);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&D, h);
        global_dpd_->buf4_mat_irrep_close(&D, h);
    }
    global_dpd_->buf4_close(&D);

    //////////////////////////////
    // The alpha-beta spin case //
    //////////////////////////////
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "D <Oo|Vv>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            int i = D.params->roworb[h][row][0];
            int j = D.params->roworb[h][row][1];
            for(int col = 0; col < D.params->coltot[h]; ++col){
                int a = D.params->colorb[h][col][0];
                int b = D.params->colorb[h][col][1];
                D.matrix[h][row][col] = 1.0/
                                    (regularizer_ + aOccEvals[i] + bOccEvals[j] - aVirEvals[a] - bVirEvals[b]);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&D, h);
        global_dpd_->buf4_mat_irrep_close(&D, h);
    }
    global_dpd_->buf4_close(&D);

    /////////////////////////////
    // The beta-beta spin case //
    /////////////////////////////
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>=o]+"), ID("[v>=v]+"), 0, "D <oo|vv>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            int i = D.params->roworb[h][row][0];
            int j = D.params->roworb[h][row][1];
            for(int col = 0; col < D.params->coltot[h]; ++col){
                int a = D.params->colorb[h][col][0];
                int b = D.params->colorb[h][col][1];
                D.matrix[h][row][col] = 1.0/
                                    (regularizer_ + bOccEvals[i] + bOccEvals[j] - bVirEvals[a] - bVirEvals[b]);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&D, h);
        global_dpd_->buf4_mat_irrep_close(&D, h);
    }
    global_dpd_->buf4_close(&D);

    delete [] aOccEvals;
    delete [] bOccEvals;
    delete [] aVirEvals;
    delete [] bVirEvals;

    dcft_timer_off("DCFTSolver::build_denominators()");
}

void DCFTSolver::build_gtau()
{
    dcft_timer_on("DCFTSolver::build_gtau()");

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdfile2 GT_OO, GT_oo, GT_VV, GT_vv, T_OO, T_oo, T_VV, T_vv;
    dpdbuf4 I;

    // Compute G * Tau contribution

    global_dpd_->file2_init(&T_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    global_dpd_->file2_init(&T_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
    global_dpd_->file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    global_dpd_->file2_init(&T_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");

    global_dpd_->file2_init(&GT_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "GTau <V|V>");

    // GT_AB = (AB|CD) Tau_CD
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V>=V]+"), ID("[V>=V]+"), 0, "MO Ints (VV|VV)");
    global_dpd_->contract422(&I, &T_VV, &GT_VV, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&I);
    // GT_AB -= <AB|CD> Tau_CD
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV|VV>");
    global_dpd_->contract422(&I, &T_VV, &GT_VV, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&I);
    // GT_AB += (AB|cd) Tau_cd
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[v,v]"),
                  ID("[V>=V]+"), ID("[v>=v]+"), 0, "MO Ints (VV|vv)");
    global_dpd_->contract422(&I, &T_vv, &GT_VV, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&I);

    // GT_AB = +(AB|IJ) Tau_IJ
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,O]"),
                  ID("[V>=V]+"), ID("[O>=O]+"), 0, "MO Ints (VV|OO)");
    global_dpd_->contract422(&I, &T_OO, &GT_VV, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&I);
    // GT_AB -= <AB|IJ> Tau_IJ
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,O]"),
                  ID("[V,V]"), ID("[O,O]"), 0, "MO Ints <VV|OO>");
    global_dpd_->contract422(&I, &T_OO, &GT_VV, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&I);
    // GT_AB += (AB|ij) Tau_ij
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[o,o]"),
                  ID("[V>=V]+"), ID("[o>=o]+"), 0, "MO Ints (VV|oo)");
    global_dpd_->contract422(&I, &T_oo, &GT_VV, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&GT_VV);

    global_dpd_->file2_init(&GT_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "GTau <v|v>");

    // GT_ab = +(ab|cd) Tau_cd
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v>=v]+"), ID("[v>=v]+"), 0, "MO Ints (vv|vv)");
    global_dpd_->contract422(&I, &T_vv, &GT_vv, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&I);
    // GT_ab -= <ab|cd> Tau_cd
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "MO Ints <vv|vv>");
    global_dpd_->contract422(&I, &T_vv, &GT_vv, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&I);
    // GT_ab += (ab|CD) Tau_CD
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[V,V]"),
                  ID("[v,v]"), ID("[V,V]"), 0, "MO Ints (vv|VV)");
    global_dpd_->contract422(&I, &T_VV, &GT_vv, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&I);

    // GT_ab = +(ab|ij) Tau_ij
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[o,o]"),
                  ID("[v>=v]+"), ID("[o>=o]+"), 0, "MO Ints (vv|oo)");
    global_dpd_->contract422(&I, &T_oo, &GT_vv, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&I);
    // GT_ab -= <ab|ij> Tau_ij
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[o,o]"),
                  ID("[v,v]"), ID("[o,o]"), 0, "MO Ints <vv|oo>");
    global_dpd_->contract422(&I, &T_oo, &GT_vv, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&I);
    // GT_ab += (ab|IJ) Tau_IJ
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[O,O]"),
                  ID("[v,v]"), ID("[O,O]"), 0, "MO Ints (vv|OO)");
    global_dpd_->contract422(&I, &T_OO, &GT_vv, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&GT_vv);

    global_dpd_->file2_init(&GT_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "GTau <O|O>");
    // GT_IJ = +(IJ|AB) Tau_AB
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>=O]+"), ID("[V>=V]+"), 0, "MO Ints (OO|VV)");
    global_dpd_->contract422(&I, &T_VV, &GT_OO, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&I);

    // GT_IJ -= <IJ|AB> Tau_AB
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    global_dpd_->contract422(&I, &T_VV, &GT_OO, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&I);
    // GT_IJ += (IJ|ab) Tau_ab
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[v,v]"),
                  ID("[O>=O]+"), ID("[v>=v]+"), 0, "MO Ints (OO|vv)");
    global_dpd_->contract422(&I, &T_vv, &GT_OO, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&I);
    // GT_IJ = +(IJ|KL) Tau_KL
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O>=O]+"), ID("[O>=O]+"), 0, "MO Ints (OO|OO)");
    global_dpd_->contract422(&I, &T_OO, &GT_OO, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&I);
    // GT_IJ -= <IJ|KL> Tau_KL
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "MO Ints <OO|OO>");
    global_dpd_->contract422(&I, &T_OO, &GT_OO, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&I);
    // GT_IJ += (IJ|kl) Tau_kl
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[o,o]"),
                  ID("[O>=O]+"), ID("[o>=o]+"), 0, "MO Ints (OO|oo)");
    global_dpd_->contract422(&I, &T_oo, &GT_OO, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&GT_OO);


    global_dpd_->file2_init(&GT_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "GTau <o|o>");
    // GT_ij = +(ij|ab) Tau_ab
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>=o]+"), ID("[v>=v]+"), 0, "MO Ints (oo|vv)");
    global_dpd_->contract422(&I, &T_vv, &GT_oo, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&I);
    // GT_ij -= <ij|ab> Tau_ab
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo|vv>");
    global_dpd_->contract422(&I, &T_vv, &GT_oo, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&I);
    // GT_ij += (ij|AB) Tau_AB
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[V,V]"),
                  ID("[o,o]"), ID("[V,V]"), 0, "MO Ints (oo|VV)");
    global_dpd_->contract422(&I, &T_VV, &GT_oo, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&I);
    // GT_ij = +(ij|kl) Tau_kl
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o>=o]+"), ID("[o>=o]+"), 0, "MO Ints (oo|oo)");
    global_dpd_->contract422(&I, &T_oo, &GT_oo, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&I);
    // GT_ij -= <ij|kl> Tau_kl
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "MO Ints <oo|oo>");
    global_dpd_->contract422(&I, &T_oo, &GT_oo, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&I);
    // GT_IJ += (ij|KL) Tau_KL
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[O,O]"),
                  ID("[o,o]"), ID("[O,O]"), 0, "MO Ints (oo|OO)");
    global_dpd_->contract422(&I, &T_OO, &GT_oo, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&GT_oo);

    global_dpd_->file2_close(&T_OO);
    global_dpd_->file2_close(&T_oo);
    global_dpd_->file2_close(&T_VV);
    global_dpd_->file2_close(&T_vv);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    dcft_timer_off("DCFTSolver::build_gtau()");
}

}} // Namespaces


