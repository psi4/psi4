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

void local_filter_T1(dpdfile2 *T1);
void local_filter_T2(dpdbuf4 *T2);

void dijabL2(int L_irr)
{
  dpdbuf4 L2, newLIJAB, newLijab, newLIjAb;
  dpdbuf4 d2, dIJAB, dijab, dIjAb;

  if(params.ref == 0) { /** RHF **/
    global_dpd_->buf4_init(&newLIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    global_dpd_->buf4_copy(&newLIjAb, PSIF_CC_LAMBDA, "New LIjAb Increment");
    global_dpd_->buf4_close(&newLIjAb);

    global_dpd_->buf4_init(&newLIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb Increment");
    if(params.local) local_filter_T2(&newLIjAb);
    else {
      global_dpd_->buf4_init(&dIjAb, PSIF_CC_DENOM, L_irr, 0, 5, 0, 5, 0, "dIjAb");
      global_dpd_->buf4_dirprd(&dIjAb, &newLIjAb);
      global_dpd_->buf4_close(&dIjAb);
    }
    global_dpd_->buf4_close(&newLIjAb);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_copy(&L2, PSIF_CC_LAMBDA, "New LIjAb");
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&newLIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb Increment");
    global_dpd_->buf4_axpy(&L2, &newLIjAb, 1);
    global_dpd_->buf4_close(&L2);
    /*dpd_buf4_print(&newLIjAb,outfile,1);*/
    global_dpd_->buf4_close(&newLIjAb);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 0, 5, 1, "New LIjAb");
    global_dpd_->buf4_copy(&L2, PSIF_CC_LAMBDA, "New LIJAB");
    global_dpd_->buf4_copy(&L2, PSIF_CC_LAMBDA, "New Lijab");
    global_dpd_->buf4_close(&L2);
  }
  else if(params.ref == 1) { /** ROHF **/

    global_dpd_->buf4_init(&newLIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    global_dpd_->buf4_copy(&newLIJAB, PSIF_CC_LAMBDA, "New LIJAB Increment");
    global_dpd_->buf4_close(&newLIJAB);

    global_dpd_->buf4_init(&newLIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB Increment");
    global_dpd_->buf4_init(&dIJAB, PSIF_CC_DENOM, L_irr, 1, 6, 1, 6, 0, "dIJAB");
    global_dpd_->buf4_dirprd(&dIJAB, &newLIJAB);
    global_dpd_->buf4_close(&dIJAB);
    global_dpd_->buf4_close(&newLIJAB);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_copy(&L2, PSIF_CC_LAMBDA, "New LIJAB");
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&newLIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB Increment");
    global_dpd_->buf4_axpy(&L2, &newLIJAB, 1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&newLIJAB);

    global_dpd_->buf4_init(&newLijab, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
    global_dpd_->buf4_copy(&newLijab, PSIF_CC_LAMBDA, "New Lijab Increment");
    global_dpd_->buf4_close(&newLijab);

    global_dpd_->buf4_init(&newLijab, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab Increment");
    global_dpd_->buf4_init(&dijab, PSIF_CC_DENOM, L_irr, 1, 6, 1, 6, 0, "dijab");
    global_dpd_->buf4_dirprd(&dijab, &newLijab);
    global_dpd_->buf4_close(&dijab);
    global_dpd_->buf4_close(&newLijab);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
    global_dpd_->buf4_copy(&L2, PSIF_CC_LAMBDA, "New Lijab");
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&newLijab, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab Increment");
    global_dpd_->buf4_axpy(&L2, &newLijab, 1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&newLijab);

    global_dpd_->buf4_init(&newLIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    global_dpd_->buf4_copy(&newLIjAb, PSIF_CC_LAMBDA, "New LIjAb Increment");
    global_dpd_->buf4_close(&newLIjAb);

    global_dpd_->buf4_init(&newLIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb Increment");
    global_dpd_->buf4_init(&dIjAb, PSIF_CC_DENOM, L_irr, 0, 5, 0, 5, 0, "dIjAb");
    global_dpd_->buf4_dirprd(&dIjAb, &newLIjAb);
    global_dpd_->buf4_close(&dIjAb);
    global_dpd_->buf4_close(&newLIjAb);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_copy(&L2, PSIF_CC_LAMBDA, "New LIjAb");
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&newLIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb Increment");
    global_dpd_->buf4_axpy(&L2, &newLIjAb, 1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&newLIjAb);
  }
  else if(params.ref == 2) { /** UHF **/

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    global_dpd_->buf4_init(&d2, PSIF_CC_DENOM, L_irr, 1, 6, 1, 6, 0, "dIJAB");
    global_dpd_->buf4_dirprd(&d2, &L2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&d2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "New Lijab");
    global_dpd_->buf4_init(&d2, PSIF_CC_DENOM, L_irr, 11, 16, 11, 16, 0, "dijab");
    global_dpd_->buf4_dirprd(&d2, &L2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&d2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    global_dpd_->buf4_init(&d2, PSIF_CC_DENOM, L_irr, 22, 28, 22, 28, 0, "dIjAb");
    global_dpd_->buf4_dirprd(&d2, &L2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&d2);

  }
}


}} // namespace psi::cclambda
