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

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {
/* replace R2L2_oOoO and tau2l2_oOoO with sort later */
/* just removed cleaning  - shouldn't need this */

extern void c_clean_CIJAB(dpdbuf4 *CMNEF);

void x_te_intermediates(void)
{
  dpdfile2 R1, L1;
  dpdbuf4 V, L, R, T2;
  int G_irr, L_irr, R_irr;
  G_irr = params.G_irr; L_irr = params.L_irr; R_irr = params.R_irr;

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    /* R2L2_OOOO = 0.5 Rijef Lklef */
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 2, 2, 2, 2, 0, "R2L2_OOOO");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 2, 7, 2, 7, 0, "RIJAB");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    global_dpd_->contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 2, 2, 2, 2, 0, "R2L2_oooo");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 2, 7, 2, 7, 0, "Rijab");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 2, 7, 2, 7, 0, "Lijab");
    global_dpd_->contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, 0, 0, 0, 0, 0, 0, "R2L2_OoOo");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, qpsr, 0, 0, "R2L2_oOoO");
    global_dpd_->buf4_close(&V);

    /* Tau2L2_OOOO = 0.5 Tau_ijef Lklef */
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 2, 2, 2, 2, 0, "Tau2L2_OOOO");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    global_dpd_->contract444(&T2, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 2, 2, 2, 2, 0, "Tau2L2_oooo");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 2, 7, 2, 7, 0, "Lijab");
    global_dpd_->contract444(&T2, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 0, 0, 0, 0, 0, "Tau2L2_OoOo");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->contract444(&T2, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, qpsr, 0, 0, "Tau2L2_oOoO");
    global_dpd_->buf4_close(&V);

    /* R2L2_OVOV = Rimae Ljmbe */
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVOV");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 10, 10, 10, 10, 0, "RIAJB");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAJB");
    global_dpd_->contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 10, 10, 10, 10, 0, "RIAjb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    global_dpd_->contract444(&R, &L, &V, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovov");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 10, 10, 10, 10, 0, "Riajb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "Liajb");
    global_dpd_->contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 10, 10, 10, 10, 0, "RIAjb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    global_dpd_->contract444(&R, &L, &V, 1, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVov");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 10, 10, 10, 10, 0, "RIAjb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "Liajb");
    global_dpd_->contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 10, 10, 10, 10, 0, "RIAJB");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    global_dpd_->contract444(&R, &L, &V, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovOV");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 10, 10, 10, 10, 0, "RiaJB");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAJB");
    global_dpd_->contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 10, 10, 10, 10, 0, "Riajb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    global_dpd_->contract444(&R, &L, &V, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_oVoV");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 10, 10, 10, 10, 0, "RjAIb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIbjA");
    global_dpd_->contract444(&R, &L, &V, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OvOv");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 10, 10, 10, 10, 0, "RIbjA");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LjAIb");
    global_dpd_->contract444(&R, &L, &V, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&V);

    /* L2R1_OOVO = Lijae Rke */
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 2, 11, 2, 11, 0, "L2R1_OOVO");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 2, 5, 2, 7, 0, "LIJAB");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    global_dpd_->contract424(&L, &R1, &V, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 2, 11, 2, 11, 0, "L2R1_oovo");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 2, 5, 2, 7, 0, "Lijab");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "Ria");
    global_dpd_->contract424(&L, &R1, &V, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L2R1_OoVo");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "Ria");
    global_dpd_->contract424(&L, &R1, &V, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L2R1_OovO");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjaB");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    global_dpd_->contract424(&L, &R1, &V, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 2, 11, 2, 11, 0, "L2R1_OOVO");
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 2, 10, "L2R1_OOVO(pqsr)");
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 2, 11, 2, 11, 0, "L2R1_oovo");
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 2, 10, "L2R1_oovo(pqsr)");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L2R1_OoVo");
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 0, 10, "L2R1_OoVo(pqsr)");
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, qprs, 0, 11, "L2R1_OoVo(qprs)");
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 0, 10, 0, 10, 0, "L2R1_OoVo(pqsr)");
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, qprs, 0, 10, "L2R1_OoVo(qpsr)");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L2R1_OovO");
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 0, 10, "L2R1_OovO(pqsr)");
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, qprs, 0, 11, "L2R1_OovO(qprs)");
    global_dpd_->buf4_close(&V);

    /* L1R2_OOVO = Rijae Lke */
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 2, 11, 2, 11, 0, "L1R2_OOVO");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 2, 5, 2, 7, 0, "RIJAB");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "LIA");
    global_dpd_->contract424(&R, &L1, &V, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 2, 11, 2, 11, 0, "L1R2_oovo");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 2, 5, 2, 7, 0, "Rijab");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "Lia");
    global_dpd_->contract424(&R, &L1, &V, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L1R2_OoVo");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "Lia");
    global_dpd_->contract424(&R, &L1, &V, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L1R2_OovO");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjaB");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "LIA");
    global_dpd_->contract424(&R, &L1, &V, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 2, 11, 2, 11, 0, "L1R2_OOVO");
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 2, 10, "L1R2_OOVO(pqsr)");
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 2, 11, 2, 11, 0, "L1R2_oovo");
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 2, 10, "L1R2_oovo(pqsr)");
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L1R2_OoVo");
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 0, 10, "L1R2_OoVo(pqsr)");
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 0, 10, 0, 10, 0, "L1R2_OoVo(pqsr)");
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, qprs, 0, 10, "L1R2_OoVo(qpsr)");
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L1R2_OovO");
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 0, 10, "L1R2_OovO(pqsr)");
    global_dpd_->buf4_close(&V);

    /* L2R1_VVOV = Limab Rmc */
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 7, 10, 7, 10, 0, "L2R1_VVOV");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 0, 7, 2, 7, 0, "LIJAB");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    global_dpd_->contract424(&L, &R1, &V, 1, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 7, 11, "L2R1_VVOV(pqsr)");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 7, 10, 7, 10, 0, "L2R1_vvov");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 0, 7, 2, 7, 0, "Lijab");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "Ria");
    global_dpd_->contract424(&L, &R1, &V, 1, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 7, 11, "L2R1_vvov(pqsr)");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 5, 10, 5, 10, 0, "L2R1_VvOv");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "Ria");
    global_dpd_->contract424(&L, &R1, &V, 1, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 5, 11, "L2R1_VvOv(pqsr)");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 5, 10, 5, 10, 0, "L2R1_VvoV");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJAb");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    global_dpd_->contract424(&L, &R1, &V, 1, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 5, 11, "L2R1_VvoV(pqsr)");
    global_dpd_->buf4_close(&V);

    /* R2L1_VVOV = Rimab Lmc */
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 7, 10, 7, 10, 0, "R2L1_VVOV");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 0, 7, 2, 7, 0, "RIJAB");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "LIA");
    global_dpd_->contract424(&R, &L1, &V, 1, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 7, 11, "R2L1_VVOV(pqsr)");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 7, 10, 7, 10, 0, "R2L1_vvov");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 0, 7, 2, 7, 0, "Rijab");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "Lia");
    global_dpd_->contract424(&R, &L1, &V, 1, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 7, 11, "R2L1_vvov(pqsr)");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 5, 10, 5, 10, 0, "R2L1_VvOv");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "Lia");
    global_dpd_->contract424(&R, &L1, &V, 1, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 5, 11, "R2L1_VvOv(pqsr)");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 5, 10, 5, 10, 0, "R2L1_VvoV");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 0, 5, 0, 5, 0, "RiJAb");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "LIA");
    global_dpd_->contract424(&R, &L1, &V, 1, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 5, 11, "R2L1_VvoV(pqsr)");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->buf4_sort(&T2, PSIF_CC_TAMPS, qprs, 0, 5, "taujIAb");
    global_dpd_->buf4_close(&T2);

  }
  else if(params.ref == 2) { /** UHF **/
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 2, 2, 2, 2, 0, "R2L2_OOOO");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 2, 7, 2, 7, 0, "RIJAB");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    global_dpd_->contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 12, 12, 12, 12, 0, "R2L2_oooo");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 12, 17, 12, 17, 0, "Rijab");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 12, 17, 12, 17, 0, "Lijab");
    global_dpd_->contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 22, 22, 22, 22, 0, "R2L2_OoOo");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 22, 28, 22, 28, 0, "RIjAb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, qpsr, 23, 23, "R2L2_oOoO");
    global_dpd_->buf4_close(&V);

    /* Tau2L2_OOOO = 0.5 Tau_ijef Lklef */
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 2, 2, 2, 2, 0, "Tau2L2_OOOO");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    global_dpd_->contract444(&T2, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 12, 12, 12, 12, 0, "Tau2L2_oooo");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 12, 17, 12, 17, 0, "Lijab");
    global_dpd_->contract444(&T2, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 22, 22, 22, 22, 0, "Tau2L2_OoOo");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->contract444(&T2, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, qpsr, 23, 23, "Tau2L2_oOoO");
    global_dpd_->buf4_close(&V);

    /* R2L2_OVOV = Rimae Ljmbe */
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 20, 20, 20, 20, 0, "R2L2_OVOV");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 20, 20, 20, 20, 0, "RIAJB");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 20, 20, 20, 20, 0, "LIAJB");
    global_dpd_->contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 20, 30, 20, 30, 0, "RIAjb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    global_dpd_->contract444(&R, &L, &V, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 30, 30, 30, 30, 0, "R2L2_ovov");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 30, 30, 30, 30, 0, "Riajb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 30, 30, 30, 30, 0, "Liajb");
    global_dpd_->contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 20, 30, 20, 30, 0, "RIAjb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    global_dpd_->contract444(&R, &L, &V, 1, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 20, 30, 20, 30, 0, "R2L2_OVov");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 20, 30, 20, 30, 0, "RIAjb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 30, 30, 30, 30, 0, "Liajb");
    global_dpd_->contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 20, 20, 20, 20, 0, "RIAJB");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    global_dpd_->contract444(&R, &L, &V, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 30, 20, 30, 20, 0, "R2L2_ovOV");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 30, 20, 30, 20, 0, "RiaJB");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 20, 20, 20, 20, 0, "LIAJB");
    global_dpd_->contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 30, 30, 30, 30, 0, "Riajb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    global_dpd_->contract444(&R, &L, &V, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 27, 27, 27, 27, 0, "R2L2_oVoV");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 24, 27, 24, 27, 0, "LIbjA");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 27, 24, 27, 24, 0, "RjAIb");
    global_dpd_->contract444(&R, &L, &V, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 24, 24, 24, 24, 0, "R2L2_OvOv");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 27, 24, 27, 24, 0, "LjAIb");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 24, 27, 24, 27, 0, "RIbjA");
    global_dpd_->contract444(&R, &L, &V, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    /* L2R1_OOVO = Lijae Rke */
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 2, 21, 2, 21, 0, "L2R1_OOVO");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 2, 5, 2, 7, 0, "LIJAB");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    global_dpd_->contract424(&L, &R1, &V, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 12, 31, 12, 31, 0, "L2R1_oovo");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 12, 15, 12, 17, 0, "Lijab");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
    global_dpd_->contract424(&L, &R1, &V, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 22, 26, 22, 26, 0, "L2R1_OoVo");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
    global_dpd_->contract424(&L, &R1, &V, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 22, 25, 22, 25, 0, "L2R1_OovO");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 22, 29, 22, 29, 0, "LIjaB");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    global_dpd_->contract424(&L, &R1, &V, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 2, 21, 2, 21, 0, "L2R1_OOVO");
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 2, 20, "L2R1_OOVO(pqsr)");
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 12, 31, 12, 31, 0, "L2R1_oovo");
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 12, 30, "L2R1_oovo(pqsr)");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 22, 26, 22, 26, 0, "L2R1_OoVo");
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 22, 27, "L2R1_OoVo(pqsr)");
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, qprs, 23, 26, "L2R1_OoVo(qprs)");
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 22, 27, 22, 27, 0, "L2R1_OoVo(pqsr)");
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, qprs, 23, 27, "L2R1_OoVo(qpsr)");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 22, 25, 22, 25, 0, "L2R1_OovO");
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 22, 24, "L2R1_OovO(pqsr)");
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, qprs, 23, 25, "L2R1_OovO(qprs)");
    global_dpd_->buf4_close(&V);

    /* L1R2_OOVO = Rijae Lke */
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 2, 21, 2, 21, 0, "L1R2_OOVO");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 2, 5, 2, 7, 0, "RIJAB");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "LIA");
    global_dpd_->contract424(&R, &L1, &V, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 12, 31, 12, 31, 0, "L1R2_oovo");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 12, 15, 12, 17, 0, "Rijab");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 2, 3, "Lia");
    global_dpd_->contract424(&R, &L1, &V, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 22, 26, 22, 26, 0, "L1R2_OoVo");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 22, 28, 22, 28, 0, "RIjAb");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 2, 3, "Lia");
    global_dpd_->contract424(&R, &L1, &V, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 22, 25, 22, 25, 0, "L1R2_OovO");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 22, 29, 22, 29, 0, "RIjaB");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "LIA");
    global_dpd_->contract424(&R, &L1, &V, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 2, 21, 2, 21, 0, "L1R2_OOVO");
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 2, 20, "L1R2_OOVO(pqsr)");
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 12, 31, 12, 31, 0, "L1R2_oovo");
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 12, 30, "L1R2_oovo(pqsr)");
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 22, 26, 22, 26, 0, "L1R2_OoVo");
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 22, 27, "L1R2_OoVo(pqsr)");
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 22, 27, 22, 27, 0, "L1R2_OoVo(pqsr)");
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, qprs, 23, 27, "L1R2_OoVo(qpsr)");
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 22, 25, 22, 25, 0, "L1R2_OovO");
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 22, 24, "L1R2_OovO(pqsr)");
    global_dpd_->buf4_close(&V);

    /* L2R1_VVOV = Limab Rmc */
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 7, 20, 7, 20, 0, "L2R1_VVOV");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 0, 7, 2, 7, 0, "LIJAB");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    global_dpd_->contract424(&L, &R1, &V, 1, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 7, 21, "L2R1_VVOV(pqsr)");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 17, 30, 17, 30, 0, "L2R1_vvov");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 17, 12, 17, 0, "Lijab");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
    global_dpd_->contract424(&L, &R1, &V, 1, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 17, 31, "L2R1_vvov(pqsr)");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 28, 24, 28, 24, 0, "L2R1_VvOv");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
    global_dpd_->contract424(&L, &R1, &V, 1, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 28, 25, "L2R1_VvOv(pqsr)");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 28, 27, 28, 27, 0, "L2R1_VvoV");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 23, 28, 23, 28, 0, "LiJAb");
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    global_dpd_->contract424(&L, &R1, &V, 1, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 28, 26, "L2R1_VvoV(pqsr)");
    global_dpd_->buf4_close(&V);

    /* R2L1_VVOV = Rimab Lmc */
    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 7, 20, 7, 20, 0, "R2L1_VVOV");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 0, 7, 2, 7, 0, "RIJAB");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "LIA");
    global_dpd_->contract424(&R, &L1, &V, 1, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 7, 21, "R2L1_VVOV(pqsr)");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 17, 30, 17, 30, 0, "R2L1_vvov");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 10, 17, 12, 17, 0, "Rijab");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 2, 3, "Lia");
    global_dpd_->contract424(&R, &L1, &V, 1, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 17, 31, "R2L1_vvov(pqsr)");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 28, 24, 28, 24, 0, "R2L1_VvOv");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 22, 28, 22, 28, 0, "RIjAb");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 2, 3, "Lia");
    global_dpd_->contract424(&R, &L1, &V, 1, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 28, 25, "R2L1_VvOv(pqsr)");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 28, 27, 28, 27, 0, "R2L1_VvoV");
    global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 23, 28, 23, 28, 0, "RiJAb");
    global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "LIA");
    global_dpd_->contract424(&R, &L1, &V, 1, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&R);
    global_dpd_->buf4_sort(&V, PSIF_EOM_TMP, pqsr, 28, 26, "R2L1_VvoV(pqsr)");
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    global_dpd_->buf4_sort(&T2, PSIF_CC_TAMPS, qprs, 23, 28, "taujIAb");
    global_dpd_->buf4_close(&T2);
  }
}

}} // namespace psi::ccdensity
