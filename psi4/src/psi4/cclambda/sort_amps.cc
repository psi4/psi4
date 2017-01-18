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

void sort_amps(int L_irr)
{
  dpdbuf4 L2;

  if(params.ref == 0) {
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_scmcopy(&L2, PSIF_CC_LAMBDA, "2 LIjAb - LIjBa", 2);
    global_dpd_->buf4_sort_axpy(&L2, PSIF_CC_LAMBDA, pqsr, 0, 5, "2 LIjAb - LIjBa", -1);
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMBDA, prqs, 10, 10, "LIAjb");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMBDA, psqr, 10, 10, "LIbjA");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    global_dpd_->buf4_scmcopy(&L2, PSIF_CC_LAMBDA, "2 LIAjb - LIbjA", 2);
    global_dpd_->buf4_sort_axpy(&L2, PSIF_CC_LAMBDA, psrq, 10, 10, "2 LIAjb - LIbjA", -1);
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMBDA, qpsr, 0, 5, "LiJaB");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMBDA, prqs, 10, 10, "LiaJB");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMBDA, rqps, 10, 10, "LjAIb");
    global_dpd_->buf4_close(&L2);
  }

  if(params.ref == 1) { /** ROHF **/
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMBDA, prqs, 10, 10, "LIAjb");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMBDA, psqr, 10, 10, "LIbjA");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMBDA, qpsr, 0, 5, "LiJaB");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMBDA, prqs, 10, 10, "LiaJB");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMBDA, rqps, 10, 10, "LjAIb");
    global_dpd_->buf4_close(&L2);

    /* Build L2IAJB List */
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMBDA, prqs, 10, 10, "LIAJB");
    global_dpd_->buf4_close(&L2);
    /* Build L2iajb List */
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 2, 7, 0, "Lijab");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMBDA, prqs, 10, 10, "Liajb");
    global_dpd_->buf4_close(&L2);
  }
  else if(params.ref == 2) { /** UHF **/

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMBDA, qpsr, 23, 29, "LiJaB");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMBDA, prqs, 20, 20, "LIAJB");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 10, 15, 12, 17, 0, "Lijab");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMBDA, prqs, 30, 30, "Liajb");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMBDA, prqs, 20, 30, "LIAjb");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMBDA, prqs, 30, 20, "LiaJB");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMBDA, psrq, 24, 27, "LIbjA");
    global_dpd_->buf4_sort(&L2, PSIF_CC_LAMBDA, rqps, 27, 24, "LjAIb");
    global_dpd_->buf4_close(&L2);
  }

}


}} // namespace psi::cclambda
