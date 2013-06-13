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

void Lsave(int L_irr)
{
  dpdfile2 L1;
  dpdbuf4 L2;

  if(params.ref == 0 || params.ref == 1) { /** ROHF **/

    dpd_->file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    dpd_->file2_copy(&L1, PSIF_CC_LAMBDA, "LIA");
    dpd_->file2_close(&L1);

    dpd_->file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 0, 1, "New Lia");
    dpd_->file2_copy(&L1, PSIF_CC_LAMBDA, "Lia");
    dpd_->file2_close(&L1);

    dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_->buf4_copy(&L2, PSIF_CC_LAMBDA, "LIJAB");
    dpd_->buf4_close(&L2);

    dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
    dpd_->buf4_copy(&L2, PSIF_CC_LAMBDA, "Lijab");
    dpd_->buf4_close(&L2);

    dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_->buf4_copy(&L2, PSIF_CC_LAMBDA, "LIjAb");
    dpd_->buf4_close(&L2);
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_->file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
    dpd_->file2_copy(&L1, PSIF_CC_LAMBDA, "LIA");
    dpd_->file2_close(&L1);

    dpd_->file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 2, 3, "New Lia");
    dpd_->file2_copy(&L1, PSIF_CC_LAMBDA, "Lia");
    dpd_->file2_close(&L1);

    dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_->buf4_copy(&L2, PSIF_CC_LAMBDA, "LIJAB");
    dpd_->buf4_close(&L2);

    dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "New Lijab");
    dpd_->buf4_copy(&L2, PSIF_CC_LAMBDA, "Lijab");
    dpd_->buf4_close(&L2);

    dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    dpd_->buf4_copy(&L2, PSIF_CC_LAMBDA, "LIjAb");
    dpd_->buf4_close(&L2);

  }
}


}} // namespace psi::cclambda
