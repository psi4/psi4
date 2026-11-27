/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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

/*! \file
    \ingroup CCRESPONSE
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccresponse {

void G_build(const char *pert, int irrep, double omega) {
    dpdbuf4 tIjAb, Y2;
    dpdfile2 GAE, GMI;
    char lbl[32];
    double Y1_norm;

        sprintf(lbl, "G_%s_IA (%5.3f)", pert, omega);
        global_dpd_->file2_init(&GMI, PSIF_CC_OEI, irrep, 0, 0, lbl);

        /* Y(Mj,Ab) * [ 2 Y(Ij,Ab) - Y(Ij,Ba) ] --> G(M,I) */
        global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
        //sprintf(lbl, "Y_%s_(2IjAb-IjbA) (%5.3f)", pert, omega);
        sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert, omega);
        global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
        global_dpd_->contract442(&tIjAb, &Y2, &GMI, 0, 0, 1, 0);
        global_dpd_->buf4_close(&tIjAb);
        global_dpd_->buf4_close(&Y2);
        global_dpd_->file2_close(&GMI);

        sprintf(lbl, "G_%s_AE (%5.3f)", pert, omega);
        global_dpd_->file2_init(&GAE, PSIF_CC_OEI, irrep, 1, 1, lbl);
        // Y(Ij,Eb) * [ 2 Y(Ij,Ab) - Y(Ij,Ba) ] --> G(A,E) 
        global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
        //sprintf(lbl, "Y_%s_(2IjAb-IjbA) (%5.3f)", pert, omega);
	sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert, omega);
        global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
        global_dpd_->contract442(&Y2, &tIjAb, &GAE, 2, 2, -1, 0);
        global_dpd_->buf4_close(&tIjAb);
        global_dpd_->buf4_close(&Y2);
        global_dpd_->file2_close(&GAE);

    return;
}

}  // namespace ccresponse
}  // namespace psi
