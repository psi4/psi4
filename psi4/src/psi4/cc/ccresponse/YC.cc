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
    \ingroup ccresponse
      
    Computes YCX contributions to Quadratic Response Functions.
    
    Author: Monika Kodrycka
*/

#include <cstdio>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccresponse {

double YC(const char *pert_y, int irrep_y, double omega_y, const char *pert_c, int irrep_c, double omega_c) {
    double polar = 0.0;
    dpdfile2 Y1, mu1;
    dpdbuf4 Y2, mu2;
    char lbl[32];
    double Y1_norm;

    /*** <0|Y1(B) * A_bar|0> ***/

    sprintf(lbl, "%sBAR_IA", pert_c);
    global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep_c, 0, 1, lbl); 
    sprintf(lbl, "Y_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&Y1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    
    polar += global_dpd_->file2_dot(&mu1, &Y1);

    global_dpd_->file2_close(&mu1);
    global_dpd_->file2_close(&Y1);

    /*** <0|Y2(B) * A_bar|0> ***/

    sprintf(lbl, "%sBAR_IjAb", pert_c);
    global_dpd_->buf4_init(&mu2, PSIF_CC_LR, irrep_c, 0, 5, 0, 5, 0, lbl);  
    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);

    polar +=  0.5 * global_dpd_->buf4_dot(&mu2, &Y2);
   
    global_dpd_->buf4_close(&mu2);
    global_dpd_->buf4_close(&Y2); 

    return polar;
}

}  // namespace ccresponse
}  // namespace psi
