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
//  V     : The converged eigenvector in the S manifold.
//  D     : This gives contribution such that \frac{\partial A(\omega)}{\partial \omega}V(\omega).
//  XOVOV : A 2h-2p intermediate.
//  YOVOV : A 2h-2p intermediate.
//  ZOVOV : A 2h-2p intermediate.
//  BOVOV : A 2h-2p intermediate.
//

double
ADCWfn::rhf_differentiate_omega(int irrep, int root)
{
    char lbl[32];
    double dot;
    dpdfile2 S, D;
    dpdbuf4 A, V, K, Z;

    sprintf(lbl, "V^(%d)_[%d]12", root, irrep);
    global_dpd_->file2_init(&S, PSIF_ADC, irrep, ID('O'), ID('V'), lbl);
    sprintf(lbl, "D^(%d)_[%d]12", root, irrep);
    global_dpd_->file2_init(&D, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);

    global_dpd_->buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"), ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
    sprintf(lbl, "ZOOVV_[%d]1234", irrep);
    global_dpd_->buf4_init(&Z, PSIF_ADC_SEM, irrep, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, lbl);
    // ZOVOV_{jiab} <--   \sum_{c} <jc|ab> v_{ic}
    global_dpd_->contract424(&V, &S, &Z, 1, 1, 1,  1, 0);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,O]"), ID("[O,O]"), ID("[V,O]"), 0, "MO Ints <OO|VO>");
    // ZOVOV_{ijab} <-- - \sum_{k} <ij|ak> v_{kb}
    global_dpd_->contract424(&V, &S, &Z, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&V);

    // dB_{iajb} <-- - (2Z_{ijab}-Z_{ijba}+2Z_{jiab}-Z_{jiba}) / (\omega+e_i-e_a+e_j-e_b)^2
    sprintf(lbl, "BOOVV_[%d]1234", irrep);
    global_dpd_->buf4_scmcopy(&Z, PSIF_ADC_SEM, lbl, -2.0);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_ADC_SEM, pqsr, ID("[O,O]"), ID("[V,V]"), lbl,  1.0);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_ADC_SEM, qprs, ID("[O,O]"), ID("[V,V]"), lbl,  1.0);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_ADC_SEM, qpsr, ID("[O,O]"), ID("[V,V]"), lbl, -2.0);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z, PSIF_ADC_SEM, irrep, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, lbl);
    sprintf(lbl, "D_[%d]1234", irrep);
    global_dpd_->buf4_init(&A, PSIF_ADC_SEM, irrep, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, lbl);
    global_dpd_->buf4_dirprd(&A, &Z);
    global_dpd_->buf4_dirprd(&A, &Z);
    global_dpd_->buf4_close(&A);

    global_dpd_->buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"), ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
    // \dV_{ia} <-- \sum_{jbc} B_{jicb} <ja|cb>
    global_dpd_->contract442(&Z, &V, &D, 1, 1, 1, 0);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,O]"), ID("[O,O]"), ID("[V,O]"), 0, "MO Ints <OO|VO>");
    // \dV_{ia} <-- - \sum_{jkb} <kj|bi> B_{jkab}
    global_dpd_->contract442(&V, &Z, &D, 3, 3, -1, 0); //This is genuine
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&Z);

    // \frac{\partial \omega^{eigen}}{\partial omega} = V^t\frac{\partial A(\omega)}{\partial \omega}V
    dot = global_dpd_->file2_dot(&S, &D);

    global_dpd_->file2_close(&S);
    global_dpd_->file2_close(&D);

    //std::cout << "DIFF:: " << dot << std::endl;

    return dot;
}

}} // End Namespaces
