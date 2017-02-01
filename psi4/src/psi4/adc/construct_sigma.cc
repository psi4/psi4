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

#include "psi4/psi4-dec.h"
#include "psi4/libtrans/integraltransform.h"
#include "adc.h"

namespace psi{ namespace adc{

//
//  DOV   : Tensors that contain the contribution from a 3h-3p diagram that cannot be
//          be evaluated independently from the trial vector, B.
//  EOV   : An intermediate introduced in order to retain the symmetricity of the response matrix.
//  XOOVV : A 2h-2p intermediate.
//  YOOVV : A 2h-2p intermediate.
//  ZOOVV : A 2h-2p intermediate.
//  BOOVV : A 2h-2p intermediate.
//

void
ADCWfn::rhf_construct_sigma(int irrep, int root)
{
    bool do_pr = options_.get_bool("PR");
    char lbl[32], ampname[32];
    dpdfile2 B, S, D, E, Bt, C;
    dpdbuf4 A, V, K, Z, BT, XT;

    sprintf(lbl, "S^(%d)_[%d]12", root, irrep);
    global_dpd_->file2_init(&S, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
    sprintf(lbl, "B^(%d)_[%d]12", root, irrep);
    global_dpd_->file2_init(&B, PSIF_ADC,     irrep, ID('O'), ID('V'), lbl);

    // CIS term and the two 3h-3p diagrams are summed into the sigma vector.
    global_dpd_->buf4_init(&A, PSIF_ADC_SEM, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "A3h3p1234");
    global_dpd_->contract422(&A, &B, &S, 0, 0, 1, 0);
    global_dpd_->buf4_close(&A);

    // Evaluation of the remaining one 3h-3p diagram
    if(do_pr) strcpy(ampname, "tilde 2 K1234 - K1243");
    else      strcpy(ampname, "2 K1234 - K1243");
    global_dpd_->buf4_init(&K, PSIF_ADC,          0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, ampname);
    global_dpd_->buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "MO Ints 2 V1234 - V1243");

    sprintf(lbl, "DOV_[%d]12", irrep);
    global_dpd_->file2_init(&D, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
    // D_{ia} <-- \sum_{jb} (2 <ij|ab> - <ij|ba>) b_{jb}
    global_dpd_->dot24(&B, &V, &D, 0, 0, 1, 0);
    // \sigma_{ia} <-- 0.5 \sum_{jb} (2 K_{ijab} - K_{ijba}) D_{jb}
    global_dpd_->dot24(&D, &K, &S, 0, 0, 0.5, 1);
    global_dpd_->file2_close(&D);

    sprintf(lbl, "EOV_[%d]12", irrep);
    global_dpd_->file2_init(&E, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
    // E_{ia} <-- \sum_{jb} (2 K_{ijab} - K_{ijba}) b_{jb}
    global_dpd_->dot24(&B, &K, &E, 0, 0, 1, 0);
    // \sigma_{ia} <-- \sum_{jb} (2 <ij|ab> - <ij|ba>) E_{jb}
    global_dpd_->dot24(&E, &V, &S, 0, 0, 0.5, 1);
    global_dpd_->file2_close(&E);

#if DEBUG_
    outfile->Printf( ">> In construction of sigma <<\n");
    global_dpd_->buf4_print(&K, outfile, 1);
    //abort();
#endif

    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"), ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
    sprintf(lbl, "ZOOVV_[%d]1234", irrep);
    global_dpd_->buf4_init(&Z, PSIF_ADC_SEM, irrep, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, lbl);
    // ZOVOV_{jiab} <--  \sum_{c} <jc|ab> b_{ic}
    global_dpd_->contract424(&V, &B, &Z, 1, 1, 1,  1, 0);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,O]"), ID("[O,O]"), ID("[V,O]"), 0, "MO Ints <OO|VO>");
    // ZOVOV_{ijab} <-- - \sum_{k} <ij|ak> b_{kb}
    global_dpd_->contract424(&V, &B, &Z, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&V);

    // B_{iajb} <-- (2Z_{ijab}-Z_{ijba}+2Z_{jiab}-Z_{jiba}) / (\omega+e_i-e_a+e_j-e_b)
    sprintf(lbl, "BOOVV_[%d]1234", irrep);
    global_dpd_->buf4_scmcopy(&Z, PSIF_ADC_SEM, lbl, 2.0);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_ADC_SEM, pqsr, ID("[O,O]"), ID("[V,V]"), lbl, -1.0);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_ADC_SEM, qprs, ID("[O,O]"), ID("[V,V]"), lbl, -1.0);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_ADC_SEM, qpsr, ID("[O,O]"), ID("[V,V]"), lbl,  2.0);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z, PSIF_ADC_SEM, irrep, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, lbl);
    sprintf(lbl, "D_[%d]1234", irrep);
    global_dpd_->buf4_init(&A, PSIF_ADC_SEM, irrep, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, lbl);
    global_dpd_->buf4_dirprd(&A, &Z);
    global_dpd_->buf4_close(&A);

    global_dpd_->buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"), ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
    // \sigma_{ia} <-- \sum_{jbc} B_{jicb} <ja|cb>
    global_dpd_->contract442(&Z, &V, &S, 1, 1, 1, 1);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,O]"), ID("[O,O]"), ID("[V,O]"), 0, "MO Ints <OO|VO>");
    // \sigma_{ia} <-- - \sum_{jkb} <kj|bi> B_{jkab}
    global_dpd_->contract442(&V, &Z, &S, 3, 3, -1, 1); //This is genuine
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&Z);

    global_dpd_->file2_close(&S);
    global_dpd_->file2_close(&B);
}

}} // End Namespaces
