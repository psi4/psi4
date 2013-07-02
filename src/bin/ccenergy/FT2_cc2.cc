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
#include <libdpd/dpd.h>
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void FT2_CC2(void)
{
  dpdbuf4 newT2, T2, Z;
  dpdfile2 F;

  if(params.ref == 0) { /* RHF */
    global_dpd_->buf4_init(&newT2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->contract424(&T2, &F, &newT2, 3, 1, 0, 1, 1);
    global_dpd_->contract244(&F, &T2, &newT2, 1, 2, 1, 1, 1);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    global_dpd_->contract424(&T2, &F, &newT2, 1, 0, 1, -1, 1);
    global_dpd_->contract244(&F, &T2, &newT2, 0, 0, 0, -1, 1);
    global_dpd_->file2_close(&F);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&newT2);
  }
  else if(params.ref == 1) { /* ROHF */
    global_dpd_->buf4_init(&newT2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "Z(I>J,AB)");
    global_dpd_->contract424(&T2, &F, &Z, 3, 1, 0, 1, 0);
    global_dpd_->contract244(&F, &T2, &Z, 1, 2, 1, 1, 1);
    global_dpd_->buf4_axpy(&Z, &newT2, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->file2_close(&F);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&newT2);

    global_dpd_->buf4_init(&newT2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tijab");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "fab");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "Z(I>J,AB)");
    global_dpd_->contract424(&T2, &F, &Z, 3, 1, 0, 1, 0);
    global_dpd_->contract244(&F, &T2, &Z, 1, 2, 1, 1, 1);
    global_dpd_->buf4_axpy(&Z, &newT2, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->file2_close(&F);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&newT2);

    global_dpd_->buf4_init(&newT2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "fab");
    global_dpd_->contract424(&T2, &F, &newT2, 3, 1, 0, 1, 1);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->contract244(&F, &T2, &newT2, 1, 2, 1, 1, 1);
    global_dpd_->file2_close(&F);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&newT2);
  }
  else if(params.ref == 2) { /* UHF */

    global_dpd_->buf4_init(&newT2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "Z(I>J,AB)");
    global_dpd_->contract424(&T2, &F, &Z, 3, 1, 0, 1, 0);
    global_dpd_->contract244(&F, &T2, &Z, 1, 2, 1, 1, 1);
    global_dpd_->buf4_axpy(&Z, &newT2, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->file2_close(&F);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&newT2);

    global_dpd_->buf4_init(&newT2, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "New tijab");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 3, 3, "fab");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 12, 15, 12, 15, 0, "Z(i>j,ab)");
    global_dpd_->contract424(&T2, &F, &Z, 3, 1, 0, 1, 0);
    global_dpd_->contract244(&F, &T2, &Z, 1, 2, 1, 1, 1);
    global_dpd_->buf4_axpy(&Z, &newT2, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->file2_close(&F);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&newT2);

    global_dpd_->buf4_init(&newT2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 3, 3, "fab");
    global_dpd_->contract424(&T2, &F, &newT2, 3, 1, 0, 1, 1);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->contract244(&F, &T2, &newT2, 1, 2, 1, 1, 1);
    global_dpd_->file2_close(&F);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&newT2);

  }
}
}} // namespace psi::ccenergy
