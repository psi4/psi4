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
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libdpd/dpd.h>
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void tsave(void)
{
  dpdfile2 t1;
  dpdbuf4 t2;

  if(params.ref == 0) { /** RHF **/
    dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    dpd_->file2_copy(&t1, PSIF_CC_OEI, "tIA");
    dpd_->file2_close(&t1);

    dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_->buf4_copy(&t2, PSIF_CC_TAMPS, "tIjAb");
    dpd_->buf4_close(&t2);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    dpd_->file2_copy(&t1, PSIF_CC_OEI, "tIA");
    dpd_->file2_close(&t1);

    dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "New tia");
    dpd_->file2_copy(&t1, PSIF_CC_OEI, "tia");
    dpd_->file2_close(&t1);

    dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    dpd_->buf4_copy(&t2, PSIF_CC_TAMPS, "tIJAB");
    dpd_->buf4_close(&t2);

    dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tijab");
    dpd_->buf4_copy(&t2, PSIF_CC_TAMPS, "tijab");
    dpd_->buf4_close(&t2);

    dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_->buf4_copy(&t2, PSIF_CC_TAMPS, "tIjAb");
    dpd_->buf4_close(&t2);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    dpd_->file2_copy(&t1, PSIF_CC_OEI, "tIA");
    dpd_->file2_close(&t1);

    dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 2, 3, "New tia");
    dpd_->file2_copy(&t1, PSIF_CC_OEI, "tia");
    dpd_->file2_close(&t1);

    dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    dpd_->buf4_copy(&t2, PSIF_CC_TAMPS, "tIJAB");
    dpd_->buf4_close(&t2);

    dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "New tijab");
    dpd_->buf4_copy(&t2, PSIF_CC_TAMPS, "tijab");
    dpd_->buf4_close(&t2);

    dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
    dpd_->buf4_copy(&t2, PSIF_CC_TAMPS, "tIjAb");
    dpd_->buf4_close(&t2);

  }
}
}} // namespace psi::ccenergy
