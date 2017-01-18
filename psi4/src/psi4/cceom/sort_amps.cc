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
    \ingroup CCEOM
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

/*
 copied from cclamba to make consistent copies of R and L
 changes should be made to sort_amps in both programs as
 cceom_density is written
*/

void sort_amps(void)
{
  dpdbuf4 R2;
  int R_irr;

  /* calculate irrep of R, the irrep for root of interest. */
  R_irr = eom_params.prop_sym^moinfo.sym;

  if(params.ref == 0 || params.ref == 1) { /* RHF/ROHF */
    /* Build R2iJaB list */
    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, R_irr, 0, 5, 0, 5, 0, "RIjAb");
    global_dpd_->buf4_sort(&R2, PSIF_CC_RAMPS, qpsr, 0, 5, "RiJaB");
    global_dpd_->buf4_close(&R2);

    /* Build R2IAJB list */
    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, R_irr, 0, 5, 2, 7, 0, "RIJAB");
    global_dpd_->buf4_sort(&R2, PSIF_CC_RAMPS, prqs, 10, 10, "RIAJB");
    global_dpd_->buf4_close(&R2);

    /* Build R2iajb list */
    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, R_irr, 0, 5, 2, 7, 0, "Rijab");
    global_dpd_->buf4_sort(&R2, PSIF_CC_RAMPS, prqs, 10, 10, "Riajb");
    global_dpd_->buf4_close(&R2);

    /* Build R2IAjb list */
    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, R_irr, 0, 5, 0, 5, 0, "RIjAb");
    global_dpd_->buf4_sort(&R2, PSIF_CC_RAMPS, prqs, 10, 10, "RIAjb");
    global_dpd_->buf4_close(&R2);

    /* Build R2iaJB list */
    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, R_irr, 0, 5, 0, 5, 0, "RiJaB");
    global_dpd_->buf4_sort(&R2, PSIF_CC_RAMPS, prqs, 10, 10, "RiaJB");
    global_dpd_->buf4_close(&R2);

    /* Build R2IbjA and R2 jAIb list */
    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, R_irr, 10, 10, 10, 10, 0, "RIAjb");
    global_dpd_->buf4_sort(&R2, PSIF_CC_RAMPS, psrq, 10, 10, "RIbjA");
    global_dpd_->buf4_sort(&R2, PSIF_CC_RAMPS, rqps, 10, 10, "RjAIb");
    global_dpd_->buf4_close(&R2);
  }
  else if(params.ref == 2) { /* UHF */

    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, R_irr, 22, 28, 22, 28, 0, "RIjAb");
    global_dpd_->buf4_sort(&R2, PSIF_CC_RAMPS, qpsr, 23, 29, "RiJaB");
    global_dpd_->buf4_close(&R2);

    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, R_irr, 0, 5, 2, 7, 0, "RIJAB");
    global_dpd_->buf4_sort(&R2, PSIF_CC_RAMPS, prqs, 20, 20, "RIAJB");
    global_dpd_->buf4_close(&R2);

    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, R_irr, 10, 15, 12, 17, 0, "Rijab");
    global_dpd_->buf4_sort(&R2, PSIF_CC_RAMPS, prqs, 30, 30, "Riajb");
    global_dpd_->buf4_close(&R2);

    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, R_irr, 22, 28, 22, 28, 0, "RIjAb");
    global_dpd_->buf4_sort(&R2, PSIF_CC_RAMPS, prqs, 20, 30, "RIAjb");
    global_dpd_->buf4_close(&R2);

    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, R_irr, 23, 29, 23, 29, 0, "RiJaB");
    global_dpd_->buf4_sort(&R2, PSIF_CC_RAMPS, prqs, 30, 20, "RiaJB");
    global_dpd_->buf4_close(&R2);

    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, R_irr, 20, 30, 20, 30, 0, "RIAjb");
    global_dpd_->buf4_sort(&R2, PSIF_CC_RAMPS, psrq, 24, 27, "RIbjA");
    global_dpd_->buf4_sort(&R2, PSIF_CC_RAMPS, rqps, 27, 24, "RjAIb");
    global_dpd_->buf4_close(&R2);
  }

}


}} // namespace psi::cceom
