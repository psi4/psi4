/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

/*! \file
    \ingroup CCLAMBDA
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
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
    dpd_->buf4_init(&newLIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_->buf4_copy(&newLIjAb, PSIF_CC_LAMBDA, "New LIjAb Increment");
    dpd_->buf4_close(&newLIjAb);

    dpd_->buf4_init(&newLIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb Increment");
    if(params.local) local_filter_T2(&newLIjAb);
    else {
      dpd_->buf4_init(&dIjAb, PSIF_CC_DENOM, L_irr, 0, 5, 0, 5, 0, "dIjAb");
      dpd_->buf4_dirprd(&dIjAb, &newLIjAb);
      dpd_->buf4_close(&dIjAb);
    }
    dpd_->buf4_close(&newLIjAb);

    dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_->buf4_copy(&L2, PSIF_CC_LAMBDA, "New LIjAb");
    dpd_->buf4_close(&L2);
    dpd_->buf4_init(&newLIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb Increment");
    dpd_->buf4_axpy(&L2, &newLIjAb, 1);
    dpd_->buf4_close(&L2);
    /*dpd_buf4_print(&newLIjAb,outfile,1);*/
    dpd_->buf4_close(&newLIjAb);

    dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 0, 5, 1, "New LIjAb");
    dpd_->buf4_copy(&L2, PSIF_CC_LAMBDA, "New LIJAB");
    dpd_->buf4_copy(&L2, PSIF_CC_LAMBDA, "New Lijab");
    dpd_->buf4_close(&L2);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_->buf4_init(&newLIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_->buf4_copy(&newLIJAB, PSIF_CC_LAMBDA, "New LIJAB Increment");
    dpd_->buf4_close(&newLIJAB);

    dpd_->buf4_init(&newLIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB Increment");
    dpd_->buf4_init(&dIJAB, PSIF_CC_DENOM, L_irr, 1, 6, 1, 6, 0, "dIJAB");
    dpd_->buf4_dirprd(&dIJAB, &newLIJAB);
    dpd_->buf4_close(&dIJAB);
    dpd_->buf4_close(&newLIJAB);

    dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_->buf4_copy(&L2, PSIF_CC_LAMBDA, "New LIJAB");
    dpd_->buf4_close(&L2);
    dpd_->buf4_init(&newLIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB Increment");
    dpd_->buf4_axpy(&L2, &newLIJAB, 1);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&newLIJAB);

    dpd_->buf4_init(&newLijab, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
    dpd_->buf4_copy(&newLijab, PSIF_CC_LAMBDA, "New Lijab Increment");
    dpd_->buf4_close(&newLijab);

    dpd_->buf4_init(&newLijab, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab Increment");
    dpd_->buf4_init(&dijab, PSIF_CC_DENOM, L_irr, 1, 6, 1, 6, 0, "dijab");
    dpd_->buf4_dirprd(&dijab, &newLijab);
    dpd_->buf4_close(&dijab);
    dpd_->buf4_close(&newLijab);

    dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
    dpd_->buf4_copy(&L2, PSIF_CC_LAMBDA, "New Lijab");
    dpd_->buf4_close(&L2);
    dpd_->buf4_init(&newLijab, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
    dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab Increment");
    dpd_->buf4_axpy(&L2, &newLijab, 1);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&newLijab);

    dpd_->buf4_init(&newLIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_->buf4_copy(&newLIjAb, PSIF_CC_LAMBDA, "New LIjAb Increment");
    dpd_->buf4_close(&newLIjAb);

    dpd_->buf4_init(&newLIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb Increment");
    dpd_->buf4_init(&dIjAb, PSIF_CC_DENOM, L_irr, 0, 5, 0, 5, 0, "dIjAb");
    dpd_->buf4_dirprd(&dIjAb, &newLIjAb);
    dpd_->buf4_close(&dIjAb);
    dpd_->buf4_close(&newLIjAb);

    dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_->buf4_copy(&L2, PSIF_CC_LAMBDA, "New LIjAb");
    dpd_->buf4_close(&L2);
    dpd_->buf4_init(&newLIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb Increment");
    dpd_->buf4_axpy(&L2, &newLIjAb, 1);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&newLIjAb);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_->buf4_init(&d2, PSIF_CC_DENOM, L_irr, 1, 6, 1, 6, 0, "dIJAB");
    dpd_->buf4_dirprd(&d2, &L2);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&d2);

    dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "New Lijab");
    dpd_->buf4_init(&d2, PSIF_CC_DENOM, L_irr, 11, 16, 11, 16, 0, "dijab");
    dpd_->buf4_dirprd(&d2, &L2);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&d2);

    dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    dpd_->buf4_init(&d2, PSIF_CC_DENOM, L_irr, 22, 28, 22, 28, 0, "dIjAb");
    dpd_->buf4_dirprd(&d2, &L2);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&d2);

  }
}


}} // namespace psi::cclambda
