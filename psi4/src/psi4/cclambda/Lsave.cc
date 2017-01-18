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
    \ingroup CCLAMBDA
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

void Lsave(int L_irr)
{
  dpdfile2 L1;
  dpdbuf4 L2;

  if(params.ref == 0 || params.ref == 1) { /** ROHF **/

    global_dpd_->file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    global_dpd_->file2_copy(&L1, PSIF_CC_LAMBDA, "LIA");
    global_dpd_->file2_close(&L1);

    global_dpd_->file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 0, 1, "New Lia");
    global_dpd_->file2_copy(&L1, PSIF_CC_LAMBDA, "Lia");
    global_dpd_->file2_close(&L1);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    global_dpd_->buf4_copy(&L2, PSIF_CC_LAMBDA, "LIJAB");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
    global_dpd_->buf4_copy(&L2, PSIF_CC_LAMBDA, "Lijab");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    global_dpd_->buf4_copy(&L2, PSIF_CC_LAMBDA, "LIjAb");
    global_dpd_->buf4_close(&L2);
  }
  else if(params.ref == 2) { /** UHF **/
    global_dpd_->file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    global_dpd_->file2_copy(&L1, PSIF_CC_LAMBDA, "LIA");
    global_dpd_->file2_close(&L1);

    global_dpd_->file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 2, 3, "New Lia");
    global_dpd_->file2_copy(&L1, PSIF_CC_LAMBDA, "Lia");
    global_dpd_->file2_close(&L1);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    global_dpd_->buf4_copy(&L2, PSIF_CC_LAMBDA, "LIJAB");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "New Lijab");
    global_dpd_->buf4_copy(&L2, PSIF_CC_LAMBDA, "Lijab");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    global_dpd_->buf4_copy(&L2, PSIF_CC_LAMBDA, "LIjAb");
    global_dpd_->buf4_close(&L2);

  }
}


}} // namespace psi::cclambda
